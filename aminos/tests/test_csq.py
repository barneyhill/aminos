import unittest
import aminos

csq_examples = {
    "inframe_insertion": {
        "bsq": "inframe_insertion|COPZ2|ENST00000006101|protein_coding|-|18AGRGP>18AQAGGP|46115072C>CG+46115084T>TGG",
        "expected_transcript": "ENST00000006101",
        "expected_mutation": aminos.processing.mutations.Mutation(mut_code='I', ref_pos=18, ref_seq='AGRGP', alt_pos=18, alt_seq='AQAGGP')
    },
    "inframe_deletion": {
        "bsq": "inframe_deletion|TPTE|ENST00000622113|protein_coding|+|201IR>201I|10569526TAAG>T",
        "expected_transcript": "ENST00000622113",
        "expected_mutation": aminos.processing.mutations.Mutation(mut_code='D', ref_pos=201, alt_pos=201, ref_seq='IR', alt_seq='I')
    },
    "missense": {
        "bsq": "missense|TPTE|ENST00000622113|protein_coding|+|368K>368E|10592359A>G",
        "expected_transcript": "ENST00000622113",
        "expected_mutation": aminos.processing.mutations.Mutation(mut_code='M', ref_pos=368, alt_pos=368, ref_seq='K', alt_seq='E')
    },
}


class TestCsqProcessing(unittest.TestCase):

    def test_mutations(self):
        for mutation_type, example in csq_examples.items():
            with self.subTest(mutation_type=mutation_type):
                actual_mutation = aminos.processing.csq.process_csq(example['expected_transcript'], example['bsq'])
                self.assertEqual(actual_mutation, example['expected_mutation'])