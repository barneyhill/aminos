mutation_dict = {
    "*frameshift": "R",
    "*frameshift&stop_retained": "Q",
    "*inframe_deletion": "C",
    "*inframe_insertion": "J",
    "*missense": "N",
    "*missense&inframe_altering": "K",
    "*stop_gained": "X",
    "*stop_gained&inframe_altering": "A",
    "frameshift": "F",
    "frameshift&start_lost": "V",
    "frameshift&stop_retained": "B",
    "inframe_deletion": "D",
    "inframe_deletion&stop_retained": "P",
    "inframe_insertion": "I",
    "inframe_insertion&stop_retained": "Z",
    "missense": "M",
    "missense&inframe_altering": "Y",
    "phi": "E",
    "start_lost": "0",
    "start_lost&splice_region": "U",
    "stop_gained": "G",
    "stop_gained&inframe_altering": "T",
    "stop_lost": "L",
    "stop_lost&frameshift": "W"
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
# (mutation_code, ref_pos, alt_pos, alt_seq, ref_seq_len)

def process_csq(csq):
    csq = csq.split('|')

    mutation_type = csq[0]
    if mutation_type not in mutation_dict:
        return None, None
    
    if len(csq) != 7:
        return None, None

    code = mutation_dict[mutation_type]

    transcript = csq[2]

    ref_pos, ref_seq, mut_pos, mut_seq = parse_amino_acid_field(csq[5])
    
    return transcript, (code, ref_pos, mut_pos, mut_seq, len(ref_seq))