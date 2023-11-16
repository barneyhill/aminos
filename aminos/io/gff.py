import pandas as pd
import logging

class GFF:
    def __init__(self, file, chrom):
        logging.debug(f"Initializing GFF with file: {file}")
        self.gff = self.read_gff(file, chrom)

    @staticmethod
    def read_gff(file, chrom):
        logging.debug(f"Reading GFF file: {file}")
        field_names = [
            'seqid', 'source', 'type', 'start', 'end', 'score', 
            'strand', 'phase', 'attributes'
        ]
        try:
            gff = pd.read_csv(file, sep='\t', comment='#', names=field_names)

            gff = gff[gff['seqid'] == chrom]

            # filter to transcripts
            gff = gff[gff['type'] == 'transcript']
            gff['ID'] = gff['attributes'].str.extract(r'ID=([^;]+)')
            gff = gff[gff['ID'].notna()]

            # throwaway version info
            gff['ID'] = gff['ID'].str.extract(r'([^\.]+)\.')

            logging.info(f"Successfully read and processed {file}")
            return gff
        except Exception as e:
            logging.error(f"Error in reading GFF file: {e}")
            raise

    def get_transcript_range(self, transcript_id):
        transcript = self.gff[self.gff['ID'] == transcript_id]

        if transcript.empty:
            logging.debug(f"No transcript found for ID: {transcript_id}")
            return None, None, None
        
        if transcript.shape[0] > 1:
            logging.error(f"More than one transcript found for ID: {transcript_id}")
            raise Exception('Unsupported more than 1 transcript')

        return (
            transcript['seqid'].values[0], 
            transcript['start'].values[0], 
            transcript['end'].values[0]
        )

    def get_unique_transcripts(self):
        unique_transcripts = self.gff['ID'].unique()
        logging.debug(f"Number of unique transcripts: {len(unique_transcripts)}")
        return unique_transcripts
