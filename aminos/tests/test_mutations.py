import unittest
from aminos.processing.mutations import Mutations, Mutation  # replace 'your_module' with the actual module name

class TestConcatMutations(unittest.TestCase):

    def setUp(self):
        # Mock transcript reference
        self.transcript_reference = "ACGTACGTACGT"
        self.mutations = Mutations(self.transcript_reference)

    def test_simple_mutation(self):
        self.mutations.add_sample_mutation("sample1", Mutation("M", 3, "G", 3, "T"))
        self.assertEqual(self.mutations.concat_mutations("sample1"), "ACTTACGTACGT")

    def test_multiple_non_overlapping_mutations(self):
        self.mutations.add_sample_mutation("sample2", Mutation("M", 2, "C", 2, "G"))
        self.mutations.add_sample_mutation("sample2", Mutation("M", 5, "A", 5, "C"))
        self.assertEqual(self.mutations.concat_mutations("sample2"), "AGGTCCGTACGT")

    def test_mismatched_reference_sequence(self):
        self.mutations.add_sample_mutation("sample4", Mutation("M", 6, "T", 6, "G"))
        # Expecting warning and unchanged sequence due to mismatch
        self.assertEqual(self.mutations.concat_mutations("sample4"), self.transcript_reference)
