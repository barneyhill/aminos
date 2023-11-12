import cyvcf2
import logging

def VCF(vcf_path, n_threads=1):
    vcf = cyvcf2.VCF(vcf_path, lazy=True, threads=n_threads)
    logging.info(f"Successfully read and processed {vcf_path} with {n_threads} threads")
    return vcf