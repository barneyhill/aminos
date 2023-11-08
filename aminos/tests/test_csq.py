import unittest
import aminos

csq_examples = {
    "inframe_deletion": (
        "inframe_deletion|TPTE|ENST00000622113|protein_coding|+|201IR>201I|10569526TAAG>T",
        ("ENST00000622113", ('D', 201, 201, 'I', 2))
    ),
    "missense_1": (
        "*missense|TPTE|ENST00000622113|protein_coding|+|368K>367E|10592359A>G",
        ("ENST00000622113", ('N', 368, 367, 'E', 1))
    ),
    "missense_2": (
        "missense|TPTE|ENST00000622113|protein_coding|+|368K>368E|10592359A>G",
        ("ENST00000622113", ('M', 368, 368, 'E', 1))
    ),
    "frameshift": (
        "frameshift|DNAJC28|ENST00000381947|protein_coding|-|316LIVPILTRQKVHFDAQKEIVRAQKIYETLIKTKEVTDRNPNNLDQGEGEKTPEIKKGFLNWMNLWKFIKIRSF*>316CSHPDQAKSPF*|33488442CAATTA>C",
        ("ENST00000381947", ('F', 316, 316, 'CSHPDQAKSPF*', 74))
    )
}

class TestCsqProcessing(unittest.TestCase):

    def test_inframe_deletion(self):
        csq, expected = csq_examples["inframe_deletion"]
        transcript, instruction = aminos.processing.csq.process_csq(csq)
        self.assertEqual(transcript, expected[0])
        self.assertEqual(instruction, expected[1])

    def test_missense_1(self):
        csq, expected = csq_examples["missense_1"]
        transcript, instruction = aminos.processing.csq.process_csq(csq)
        self.assertEqual(transcript, expected[0])
        self.assertEqual(instruction, expected[1])

    def test_missense_2(self):
        csq, expected = csq_examples["missense_2"]
        transcript, instruction = aminos.processing.csq.process_csq(csq)
        self.assertEqual(transcript, expected[0])
        self.assertEqual(instruction, expected[1])

    def test_frameshift(self):
        csq, expected = csq_examples["frameshift"]
        transcript, instruction = aminos.processing.csq.process_csq(csq)
        self.assertEqual(transcript, expected[0])
        self.assertEqual(instruction, expected[1])