import gzip

def init_transcript_reference(file_path):
    transcript_reference = {}
    with gzip.open(file_path, 'rt') as file:  # 'rt' mode to read as text
        name = None
        seq_list = []
        for line in file:
            line = line.strip()  # Remove trailing newline characters
            if line.startswith(">"):
                if name:  # If not the first sequence, save the previous sequence
                    transcript_reference[name] = ''.join(seq_list)
                    seq_list = []
                name = line[1:].split()[0]  # Get the name, discard the ">"
            else:
                seq_list.append(line)
        if name:  # Save the last sequence
            transcript_reference[name] = ''.join(seq_list)
    return transcript_reference