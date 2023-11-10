import argparse
import logging
import aminos
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description="Aminos Genetic Data Processing")
    parser.add_argument('--vcf', help='Path to VCF file', required=True)
    parser.add_argument('--gff', help='Path to GFF file', required=True)
    parser.add_argument('--fasta', help='Path to FASTA file', required=True)
    parser.add_argument('--output', help='Output directory', required=True)

    args = parser.parse_args()

    total_mutations = 0
    accepted_mutations = Counter()

    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    vcf = aminos.io.VCF(args.vcf)
    gff = aminos.io.GFF(args.gff)
    transcript_references = aminos.io.fasta.read_transcript_references(args.fasta)

    gff_transcripts = gff.get_unique_transcripts()

    for transcript_id in gff_transcripts:
        
        if transcript_id not in transcript_references:
            continue
        else:
            transcript_reference = transcript_references[transcript_id]

        mutations = aminos.processing.mutations.Mutations(transcript_reference)
        chr, start, end = gff.get_transcript_range(transcript_id)

        logging.info(f"Processing transcript: {transcript_id} ({chr}:{start}-{end})")

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

        file = aminos.io.fasta.Writer(args.output, transcript_id)

        for individual in vcf.samples:
            for haplotype in [0, 1]:
                sample_name = f'{individual}_{haplotype}'
                seq = mutations.concat_mutations(sample_name)

                file.write_header(sample_name)
                file.write_sequence(seq)

        file.close()

        total_mutations += mutations.total_mutations
        accepted_mutations += mutations.accepted_mutations

    logging.info(f"Accepted mutations total: {sum(accepted_mutations.values())} ({sum(accepted_mutations.values())/total_mutations:.2f}% accepted)")
    logging.info(f"Accepted mutations by type: {accepted_mutations}")

# example:
# python3 aminos.py --vcf data/ALL_GGVP.chr21.vcf.gz --gff data/Homo_sapiens.GRCh38.110.chromosome.21.gff3.gz --fasta data/reference_sequences.fasta.gz --output data/test

if __name__ == "__main__":
    main()