import gzip
import os

def read_transcript_references(file_path):
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

class Writer:

    def __init__(self, file_dir, transcript):
        self.file_path = os.path.join(file_dir, f'{transcript}.fa.gz')
        self.file = gzip.open(self.file_path, 'w')

    def write_header(self, individuals):
        if (type(individuals) == str) or len(individuals) == 1:
            self.file.write(f'>{individuals}\n'.encode())
        elif type(individuals) == list and len(individuals) > 1:
            self.file.write(f'>{",".join(individuals)}\n'.encode())
    
    def write_sequence(self, sequence):
        self.file.write(f'{sequence}\n'.encode())

    def close(self):
        self.file.close()