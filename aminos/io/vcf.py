import cyvcf2
import logging

def VCF(vcf_path):
    vcf = cyvcf2.VCF(vcf_path)
    logging.info(f"Successfully read and processed {vcf_path}")
    return vcf