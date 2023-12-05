import cyvcf2
import logging
import numpy as np
from typing import Tuple

def interleave_arrays(haplo1_set: np.ndarray, haplo2_set: np.ndarray) -> np.ndarray:
    # Stack arrays depth-wise and then flatten
    interleaved = np.dstack((haplo1_set, haplo2_set)).flatten()
    return interleaved

def get_haplotypes_at_index(bitmask_array: np.ndarray, index: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns haplotype presence information for a given index across a NumPy array of bitmasks.

    Parameters:
    bitmask_array (np.ndarray): Array of bitmasks.
    index (int): The index position to check in each bitmask.

    Returns:
    Tuple[np.ndarray, np.ndarray]: Two boolean arrays indicating the presence of haplotype 1 and 2 at the given index.
    """

    # Masks for checking the specific bit position for each haplotype
    haplo1_mask = 1 << (2 * index)
    haplo2_mask = 2 << (2 * index)

    # Check if the haplotype bits are set across the array
    haplo1_set = (bitmask_array & haplo1_mask) != 0
    haplo2_set = (bitmask_array & haplo2_mask) != 0

    return haplo1_set, haplo2_set

def VCF(vcf_path: str, n_threads=1) -> cyvcf2.VCF:
    vcf = cyvcf2.VCF(vcf_path, lazy=True, threads=n_threads)
    logging.info(f"Successfully read and processed {vcf_path} with {n_threads} threads")
    return vcf