from collections import defaultdict

class Instructions:
    '''
    An instruction is defined as:
    (code, ref_pos, mut_pos, mut_seq, len(ref_seq))
    '''

    def __init__(self):
        self._instruction_id_counter = 0
        self._instruction_to_id = {}
        self._id_to_instruction = {}
        self.transcript_instructions = defaultdict(set)

    def _get_instruction_id(self, instruction):
        instruction_key = str(instruction)
        if instruction_key not in self._instruction_to_id:
            self._instruction_to_id[instruction_key] = self._instruction_id_counter
            self._id_to_instruction[self._instruction_id_counter] = instruction
            self._instruction_id_counter += 1
        return self._instruction_to_id[instruction_key]
    
    def _get_id_to_instruction(self, instruction_id):
        return self._id_to_instruction.get(instruction_id, None)

    def _get_instruction_to_id(self, instruction):
        return self._instruction_to_id.get(str(instruction), None)
    
    def add_sample_transcript_instruction(self, sample, transcript, instruction):
        instruction_id = self._get_instruction_id(instruction)
        self.transcript_instructions[sample].add(instruction_id)

    def validate_all_instructions(self):
        '''
        From VCF2Prot paper:

        An instruction is defined as:
        (code, ref_pos, mut_pos, mut_seq, len(ref_seq))


        As stated above, Instructions are the first step toward generating an IR and hence the validity
        of translations is of paramount importance to ensure the correct generation of IR. Hence, three
        tests were implemented to validate the correctness of input mutations and of translations, first,
        position uniqueness, second isolated boundaries, and third prohibited sequences. Positionuniqueness refers to the requirement that each instruction must have a unique starting position
        in the reference and the altered sequence and hence having two mutations at the same position
        is prohibited. “Isolated boundaries” refers to the requirement that mutations should not overlap.
        Finally, “prohibited sequences” refers to an ‘illegal’ sequence of mutations in a transcript, for
        example, after a frameshift, stop-gain, or a stop-loss no independent mutation is allowed. If
        any of the three tests fail, the input collection of mutation and the translation is considered
        invalid and hence, it is by default ignored, i.e., skipped, and an error message is printed to the
        user.
        '''

        for sample in self.transcript_instructions:

            sample_instructions = self.get_sample_instructions(sample)

            # test posititon uniqueness:

            if not self.test_position_uniqueness(sample_instructions):
                self.transcript_instructions.pop(sample)
                continue

            # test isolated boundaries:

            if not self.test_isolated_boundaries(sample_instructions):
                self.transcript_instructions.pop(sample)
                continue

            # test prohibited sequences:

            if not self.test_prohibited_sequences(sample_instructions):
                self.transcript_instructions.pop(sample)
                continue

    def test_position_uniqueness(self, instructions):

        positions = [instruction[1] for instruction in instructions]

        if len(positions) > len(set(positions)):
            return False

        return True

    def test_isolated_boundaries(self, instructions):
        return True

    def test_prohibited_sequences(self, instructions):
        return True

    def get_sample_instructions(self, sample):
        return [ self._get_id_to_instruction(id) for id in self.transcript_instructions[sample] ] 
