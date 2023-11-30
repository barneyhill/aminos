import unittest
from aminos.processing.mutations import Mutations, Mutation  # replace 'your_module' with the actual module name

class TestConcatMutations(unittest.TestCase):

    def setUp(self):
        # Mock transcript reference
        self.transcript_reference = "ACGTACGTACGT"
        self.mutations = Mutations(self.transcript_reference)

    def test_simple_mutation(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 3, "G", 3, "T"))
        self.mutations.add_sample_mutation("sample1", mut_id)
        self.assertEqual(self.mutations.concat_mutations("sample1"), "ACTTACGTACGT")

    def test_multiple_non_overlapping_mutations(self):
        mut_id1 = self.mutations.get_mutation_id(Mutation("M", 2, "C", 2, "G"))
        mut_id2 = self.mutations.get_mutation_id(Mutation("M", 5, "A", 5, "C"))

        self.mutations.add_sample_mutation("sample2", mut_id1)
        self.mutations.add_sample_mutation("sample2", mut_id2)
        self.assertEqual(self.mutations.concat_mutations("sample2"), "AGGTCCGTACGT")

    def test_mismatched_reference_sequence(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 6, "T", 6, "G"))

        self.mutations.add_sample_mutation("sample4", mut_id)
        # Expecting warning and unchanged sequence due to mismatch
        self.assertEqual(self.mutations.concat_mutations("sample4"), self.transcript_reference)

    def test_out_of_bounds(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 13, "T", 13, "GG"))

        self.mutations.add_sample_mutation("sample5", mut_id)
        # Expecting warning and unchanged sequence due to out-of-bounds
        self.assertEqual(self.mutations.concat_mutations("sample5"), self.transcript_reference)

    def test_out_of_bounds(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 13, "T", 13, "GG"))

        self.mutations.add_sample_mutation("sample5", mut_id)
        # Expecting warning and unchanged sequence due to out-of-bounds
        self.assertEqual(self.mutations.concat_mutations("sample5"), self.transcript_reference)