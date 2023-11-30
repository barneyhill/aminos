import unittest
from aminos.processing.mutations import Mutations, Mutation  # replace 'your_module' with the actual module name

class TestConcatMutations(unittest.TestCase):

    def setUp(self):
        # Mock transcript reference
        self.transcript_reference = "ACGTACGTACGT"
        self.mutations = Mutations(self.transcript_reference)

    def test_missense(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 3, "G", 3, "T"))
        self.mutations.add_sample_mutation("sample1", mut_id)
        self.assertEqual(self.mutations.concat_mutations("sample1"), "A,C,T,T,A,C,G,T,A,C,G,T")

    def test_two_nonoverlapping_missense(self):
        mut_id1 = self.mutations.get_mutation_id(Mutation("M", 2, "C", 2, "G"))
        mut_id2 = self.mutations.get_mutation_id(Mutation("M", 5, "A", 5, "C"))

        self.mutations.add_sample_mutation("sample2", mut_id1)
        self.mutations.add_sample_mutation("sample2", mut_id2)
        self.assertEqual(self.mutations.concat_mutations("sample2"), "A,G,G,T,C,C,G,T,A,C,G,T")

    def test_mismatched_reference_sequence(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 6, "T", 6, "G"))

        self.mutations.add_sample_mutation("sample4", mut_id)
        # Expecting warning and unchanged sequence due to mismatch
        self.assertEqual(self.mutations.concat_mutations("sample4"), ','.join(list(self.transcript_reference)))

    def test_out_of_bounds(self):
        mut_id = self.mutations.get_mutation_id(Mutation("M", 13, "T", 13, "GG"))

        self.mutations.add_sample_mutation("sample5", mut_id)
        # Expecting warning and unchanged sequence due to out-of-bounds
        self.assertEqual(self.mutations.concat_mutations("sample5"), ','.join(list(self.transcript_reference)))

    def test_inframe_insertion(self):
        mut_id = self.mutations.get_mutation_id(Mutation("I", 3, "G", 3, "GT"))

        self.mutations.add_sample_mutation("sample6", mut_id)
        self.assertEqual(self.mutations.concat_mutations("sample6"), "A,C,GT,T,A,C,G,T,A,C,G,T")

    def test_inframe_deletion(self):
        mut_id = self.mutations.get_mutation_id(Mutation("I", 3, "GT", 3, "G"))

        self.mutations.add_sample_mutation("sample7", mut_id)
        self.assertEqual(self.mutations.concat_mutations("sample7"), "A,C,G,,A,C,G,T,A,C,G,T")