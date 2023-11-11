import argparse
import logging
import aminos
from collections import Counter
import tqdm

import cProfile
import pstats

def main():
    parser = argparse.ArgumentParser(description="aminos")
    parser.add_argument('--vcf', help='Path to VCF file', required=True)
    parser.add_argument('--gff', help='Path to GFF file', required=True)
    parser.add_argument('--fasta', help='Path to FASTA file', required=True)
    parser.add_argument('--output', help='Output directory', required=True)
    parser.add_argument('--debug', help='Enable debug mode', action='store_true', default=False)

    args = parser.parse_args()

    # Configure logging
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    vcf = aminos.io.VCF(args.vcf)
    gff = aminos.io.GFF(args.gff)
    transcript_references = aminos.io.fasta.read_transcript_references(args.fasta)

    gff_transcripts = gff.get_unique_transcripts()

    # Counters for stats
    total_transcripts_seen = 0
    total_mutations_seen = 0
    accepted_mutations = 0

    for transcript_id in tqdm.tqdm(gff_transcripts, desc="Iterating over transcripts"):
        
        if transcript_id not in transcript_references:
            continue
        else:
            transcript_reference = transcript_references[transcript_id]
            total_transcripts_seen += 1

        mutations = aminos.processing.mutations.Mutations(transcript_reference)
        chr, start, end = gff.get_transcript_range(transcript_id)

        logging.debug(f"Processing transcript: {transcript_id} ({chr}:{start}-{end})")

        for record in vcf(f'{chr}:{start}-{end}'):
            bcsq = record.INFO.get('BCSQ')
            if bcsq is None:
                continue

            for csq in bcsq.split(','):
                mutation = aminos.processing.csq.process_csq(transcript_id, csq)

                if not mutation:
                    continue

                for individual_call, individual in zip(record.genotypes, vcf.samples):
                    for haplotype in [0, 1]:
                        if individual_call[haplotype] == 1:
                            mutations.add_sample_mutation(f'{individual}_{haplotype}', mutation)
                            total_mutations_seen += 1

        if mutations.accepted_mutations == 0:
            continue

        for individual in vcf.samples:
            for haplotype in [0, 1]:
                sample_name = f'{individual}_{haplotype}'
                mutations.concat_mutations(sample_name)

        file = aminos.io.fasta.Writer(args.output, transcript_id)

        file.write(mutations)
        file.close()

        accepted_mutations += mutations.accepted_mutations

    logging.info(f"Accepted mutations total: {accepted_mutations} ({100 * accepted_mutations/total_mutations_seen:.2f}% accepted)")
    logging.info(f"Total transcripts seen: {total_transcripts_seen}")

# example:
# python3 aminos.py --vcf data/ALL_GGVP.chr21.vcf.gz --gff data/Homo_sapiens.GRCh38.110.chromosome.21.gff3.gz --fasta data/reference_sequences.fasta.gz --output data/test

if __name__ == "__main__":
    with cProfile.Profile() as pr:
        main()

    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.CUMULATIVE)  # Sorting by time spent
    stats.print_stats(50)
