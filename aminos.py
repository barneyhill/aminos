import argparse
import logging
from tqdm_loggable.auto import tqdm

import os
import numpy as np
import sys

import aminos

import cProfile
import pstats

def run(args):

    # Configure logging
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    vcf = aminos.io.VCF(args.vcf, args.threads)
    gff = aminos.io.GFF(args.gff, args.set_chr)
    transcript_references = aminos.io.fasta.read_transcript_references(args.fasta, args.threads)

    gff_transcripts = gff.get_unique_transcripts()

    # Counters for stats
    total_transcripts_seen = 0
    total_mutations_seen = 0
    accepted_mutations = 0

    samples = np.array([f"{individual}_{haplotype}" for individual in vcf.samples for haplotype in [0, 1]])

    for transcript_id in tqdm(gff_transcripts, desc="Iterating over transcripts"):
        
        if transcript_id not in transcript_references:
            continue
        else:
            transcript_reference = transcript_references[transcript_id]
            total_transcripts_seen += 1

        mutations = aminos.processing.mutations.Mutations(transcript_reference)
        chr, start, end = gff.get_transcript_range(transcript_id)

        if args.set_chr:
            chr = args.set_chr

        logging.debug(f"Processing transcript: {transcript_id} ({chr}:{start}-{end})")

        for record in vcf(f'{chr}:{start}-{end}'):
            bcsq = record.INFO.get('BCSQ')
            if bcsq is None:
                continue

            for csq in bcsq.split(','):
                mutation = aminos.processing.csq.process_csq(transcript_id, csq)

                if not mutation:
                    continue

                mutation_id = mutations.get_mutation_id(mutation)
                for mutated_sample in samples[(record.genotype.array()[:,0:2] == 1).flatten()]:
                    mutations.add_sample_mutation(mutated_sample, mutation_id)
                    total_mutations_seen += 1

        for sample in samples:
            mutations.concat_mutations(sample)

        if mutations.accepted_mutations == 0:
            continue

        file = aminos.io.fasta.Writer(args.output, transcript_id, args.threads)

        file.write(mutations)
        file.close()

        accepted_mutations += mutations.accepted_mutations

    logging.info(f"Accepted mutations total: {accepted_mutations} ({100 * accepted_mutations/total_mutations_seen:.2f}% accepted)")
    logging.info(f"Total transcripts seen: {total_transcripts_seen}")

# example:
# python3 aminos.py --vcf data/ALL_GGVP.chr21.vcf.gz --gff data/Homo_sapiens.GRCh38.110.chromosome.21.gff3.gz --fasta data/reference_sequences.fasta.gz --output data/test

def main():

    parser = argparse.ArgumentParser(description="aminos")
    parser.add_argument('--vcf', help='Path to VCF file', required=True)
    parser.add_argument('--gff', help='Path to GFF file', required=True)
    parser.add_argument('--fasta', help='Path to FASTA file', required=True)
    parser.add_argument('--output', help='Output directory', required=True)
    parser.add_argument('--debug', help='Enable debug mode', action='store_true', default=False)
    parser.add_argument('--cprofile', help='profile code and print summary', action='store_true')
    parser.add_argument('--set-chr', help='force a chr value', required=False)
    parser.add_argument('--threads', help='Number of threads to use for VCF reader', default=os.cpu_count(), type=int)

    args = parser.parse_args()

    if args.cprofile:
        with cProfile.Profile() as pr:
            run(args)

        stats = pstats.Stats(pr)
        stats.sort_stats(pstats.SortKey.CUMULATIVE)  # Sorting by time spent
        stats.print_stats(50)
    else:
        run(args)


if __name__ == "__main__":
    main()