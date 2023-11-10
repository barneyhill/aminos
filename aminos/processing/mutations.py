from collections import defaultdict, Counter
from dataclasses import dataclass
import logging

@dataclass
class Mutation:
    mut_code: str
    ref_pos: int
    ref_seq: str
    alt_pos: int
    alt_seq: str
    
class Mutations:
    '''
    This class stores all the mutations for a given transcript.
    It can then process the transformed transcript sequences across samples.
    '''

    def __init__(self, transcript_reference):
        self._mutation_id_counter = 0
        self._mutation_to_id = {}
        self._id_to_mutation = {}
        self.mutation_ids = defaultdict(set)
        self.transcript_reference = transcript_reference
        self.total_mutations = 0
        self.accepted_mutations = Counter()

    def _get_mutation_id(self, mutation):
        mutation_key = str(mutation)
        if mutation_key not in self._mutation_to_id:
            self._mutation_to_id[mutation_key] = self._mutation_id_counter
            self._id_to_mutation[self._mutation_id_counter] = mutation
            self._mutation_id_counter += 1
        return self._mutation_to_id[mutation_key]
    
    def _get_id_to_mutation(self, mutation_id):
        return self._id_to_mutation.get(mutation_id, None)

    def _get_mutation_to_id(self, mutation):
        return self._mutation_to_id.get(mutation, None)
    
    def add_sample_mutation(self, sample, mutation):
        mutation_id = self._get_mutation_id(mutation)
        self.mutation_ids[sample].add(mutation_id)

    def get_sample_mutations(self, sample):
        return [ self._get_id_to_mutation(id) for id in self.mutation_ids[sample] ] 
        
    def concat_mutations(self, sample):

        transcript_seq = list(self.transcript_reference)

        sample_mutations = self.get_sample_mutations(sample)

        self.total_mutations += len(sample_mutations)

        for mutation in sorted(sample_mutations, key=lambda mut: mut.ref_pos):
            logging.debug(f"Applying mutation: {mutation}")

            if not self._is_valid_mutation(mutation, transcript_seq):
                return ''.join(transcript_seq)  # Early return with the original sequence

            # Apply the mutation
            for i, ref_char in enumerate(mutation.ref_seq):
                transcript_seq[mutation.ref_pos - 1 + i] = ''  # Remove ref_seq characters
            transcript_seq[mutation.ref_pos - 1] = mutation.alt_seq  # Insert alt_seq

            self.accepted_mutations[mutation.mut_code] += 1

        return ''.join(transcript_seq)

    def _is_valid_mutation(self, mutation, transcript_seq):
        for i, ref_char in enumerate(mutation.ref_seq):
            pos = mutation.ref_pos - 1 + i

            # Out-of-bounds check:
            if pos >= len(transcript_seq):
                logging.debug(f"Mutation at position {pos + 1} is out of bounds!")
                return False

            # Mismatch check:
            if transcript_seq[pos] != ref_char:
                logging.debug(f"Reference {mutation.ref_pos + i}{ref_char} does not match transcript reference {pos + 1}{transcript_seq[pos]}")
                return False
            
            # Overlap check
            if transcript_seq[pos] == '':
                logging.debug(f"Conflict! Mutation at position {pos + 1} has already been touched!")
                return False
        return True