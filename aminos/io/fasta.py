import gzip
import os
import logging

def read_transcript_references(file_path):
    logging.info(f"Reading transcript references from: {file_path}")
    transcript_reference = {}
    try:
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
    except Exception as e:
        logging.error(f"Error in reading transcript references: {e}")
        raise
    return transcript_reference

class Writer:

    def __init__(self, file_dir, transcript):
        self.file_path = os.path.join(file_dir, f'{transcript}.fa.gz')
        try:
            self.file = gzip.open(self.file_path, 'w')
            logging.info(f"Initialized Writer for file: {self.file_path}")
        except Exception as e:
            logging.error(f"Error opening file for writing: {e}")
            raise

    def write_header(self, individuals):
        try:
            if (type(individuals) == str) or len(individuals) == 1:
                self.file.write(f'>{individuals}\n'.encode())
            elif type(individuals) == list and len(individuals) > 1:
                self.file.write(f'>{",".join(individuals)}\n'.encode())
        except Exception as e:
            logging.error(f"Error writing header: {e}")
            raise
    
    def write_sequence(self, sequence):
        try:
            self.file.write(f'{sequence}\n'.encode())
        except Exception as e:
            logging.error(f"Error writing sequence: {e}")
            raise

    def close(self):
        try:
            self.file.close()
        except Exception as e:
            logging.error(f"Error closing file: {e}")
            raise
