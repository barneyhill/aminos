import cyvcf2

def read_vcf(vcf_path):
    return cyvcf2.VCF(vcf_path)