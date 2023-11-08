from collections import defaultdict

class Instructions:
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
        self.transcript_instructions[(sample, transcript)].add(instruction_id)

    def get_sample_transcript_instructions(self, sample, transcript):
        return [ self._get_id_to_instruction(id) for id in self.transcript_instructions[(sample, transcript)] ] 
