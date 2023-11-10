from aminos.processing.mutations import Mutation

mutation_dict = {
    "inframe_deletion": "D",
    "inframe_insertion": "I",
    "missense": "M",
}

def parse_amino_acid_seq_position(input_seq):

    position = ''
    sequence = ''
    for char in input_seq:
        if char.isdigit():
            position += char
        else:
            sequence += char

    if not sequence:
        sequence = '*'

    return int(position), sequence

def parse_amino_acid_field(input_string):
    parsed_strings = input_string.split('>')
    if len(parsed_strings) != 2:
        raise ValueError(f"The parsed string has a length of: {len(parsed_strings)}, expected only two")

    try:
        ref_pos, ref_seq = parse_amino_acid_seq_position(parsed_strings[0])
    except ValueError as e:
        raise Exception(f"While extracting the sequence and the position of the reference the following error was encountered: {e}")

    try:
        mut_pos, mut_seq = parse_amino_acid_seq_position(parsed_strings[1])
    except ValueError as e:
        raise Exception(f"While extracting the sequence and the position of the reference the following error was encountered: {e}")

    return ref_pos, ref_seq, mut_pos, mut_seq

# take a bcftools csq string and return the transcript and instruction
# an instruction is defined as:
# (mutation_code, ref_pos, ref_seq, mut_pos, mut_seq)

def process_csq(transcript_id, csq):
    csq = csq.split('|')

    if len(csq) != 7:
        return None

    if transcript_id != csq[2]:
        return None

    mutation_type = csq[0]

    if mutation_type not in mutation_dict:
        return None

    mutation_code = mutation_dict[mutation_type]
    
    ref_pos, ref_seq, mut_pos, mut_seq = parse_amino_acid_field(csq[5])
    
    return Mutation(mut_code=mutation_code, ref_pos=ref_pos, ref_seq=ref_seq, alt_pos=mut_pos, alt_seq=mut_seq)