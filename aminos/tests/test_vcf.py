import unittest
import numpy as np
from aminos.io.vcf import get_haplotypes_at_index

class TestGetHaplotypesAtIndex(unittest.TestCase):

    def test_single_value(self):
        bitmask_array = np.array([3], dtype=np.int32)  # 3 = 11 in binary
        index = 0
        expected_haplo1 = np.array([True])
        expected_haplo2 = np.array([True])
        haplo1_set, haplo2_set = get_haplotypes_at_index(bitmask_array, index)
        np.testing.assert_array_equal(haplo1_set, expected_haplo1)
        np.testing.assert_array_equal(haplo2_set, expected_haplo2)

    def test_multiple_values(self):
        bitmask_array = np.array([3, 1, 2, 0], dtype=np.int32)  # 1 = 01, 2 = 10 in binary
        index = 0
        expected_haplo1 = np.array([True, True, False, False])
        expected_haplo2 = np.array([True, False, True, False])
        haplo1_set, haplo2_set = get_haplotypes_at_index(bitmask_array, index)
        np.testing.assert_array_equal(haplo1_set, expected_haplo1)
        np.testing.assert_array_equal(haplo2_set, expected_haplo2)

    def test_no_haplotypes_set(self):
        bitmask_array = np.array([0, 0, 0], dtype=np.int32)
        index = 1
        expected_haplo1 = np.array([False, False, False])
        expected_haplo2 = np.array([False, False, False])
        haplo1_set, haplo2_set = get_haplotypes_at_index(bitmask_array, index)
        np.testing.assert_array_equal(haplo1_set, expected_haplo1)
        np.testing.assert_array_equal(haplo2_set, expected_haplo2)