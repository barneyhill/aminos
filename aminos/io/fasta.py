import pgzip
import os
import logging

def read_transcript_references(file_path, threads):
    logging.debug(f"Reading transcript references from: {file_path}")
    transcript_reference = {}
    try:
        with pgzip.open(file_path, 'rt', thread=threads) as file:  # 'rt' mode to read as text
            name = None
            seq_list = []

            for line in file:
                line = line.strip()  # Remove trailing newline characters

                if line.startswith(">"):
                    if name:
                        transcript_reference[name] = ''.join(seq_list)
                        seq_list = []

                    info = line[1:].split('|')  # Get the name, discard the ">"
                    name = info[0]

                    #example: transcript:ENST00000621433.1
                    name = name.split(' ')[4]
                    name = name.split('.')[0]
                    # remove transcript:
                    name = name.split(':')[1]

                else:
                    seq_list.append(line)
            if name:  # Save the last sequence
                transcript_reference[name] = ''.join(seq_list)

        logging.info(f"Successfully read {len(transcript_reference)} transcript references")
    except Exception as e:
        logging.error(f"Error in reading transcript references: {e}")
        raise
    return transcript_reference

class Writer:

    def __init__(self, file_dir, transcript, threads=1):
        self.file_path = os.path.join(file_dir, f'{transcript}.fa.gz')
        try:
            self.file = pgzip.open(self.file_path, 'w', compresslevel=5, thread=threads)
            logging.debug(f"Initialized Writer for file: {self.file_path}")
        except Exception as e:
            logging.error(f"Error opening file for writing: {e}")
            raise
        self.buffer = bytearray()

    def write_header(self, individual):
        try:
            self.buffer.extend(f'>{individual}\n'.encode())
        except Exception as e:
            logging.error(f"Error writing header: {e}")
            raise
    
    def write_sequence(self, sequence):
        try:
            self.buffer.extend(f'{sequence}\n'.encode())
        except Exception as e:
            logging.error(f"Error writing sequence: {e}")
            raise

    def write(self, mutations_store):
        ids_to_sequence = mutations_store.get_mutation_ids_to_sequence()

        for mutations, samples in mutations_store.get_mutation_ids_to_samples().items():
            self.write_header(','.join(samples))
            self.write_sequence(ids_to_sequence[mutations])
        self.flush_buffer()

    def flush_buffer(self):
        self.file.write(self.buffer)
        self.buffer = bytearray()

    def close(self):
        self.flush_buffer()
        self.file.close()