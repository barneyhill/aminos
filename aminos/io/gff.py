import pandas as pd

class GFF:
    def __init__(self, file):
        self.gff = self.read_gff(file)

    @staticmethod
    def read_gff(file):
        field_names = [
            'seqid', 'source', 'type', 'start', 'end', 'score', 
            'strand', 'phase', 'attributes'
        ]
        gff = pd.read_csv(file, sep='\t', comment='#', names=field_names)
        gff['ID'] = gff['attributes'].str.extract(r'ID=transcript:([^;]+)')
        gff = gff[gff['ID'].notna()]
        return gff

    def get_transcript_range(self, transcript_id):
        transcript = self.gff[self.gff['ID'] == transcript_id]

        if transcript.empty:
            return None, None, None
        
        if transcript.shape[0] > 1:
            raise Exception('Unsupported more than 1 transcript')

        return (
            transcript['seqid'].values[0], 
            transcript['start'].values[0], 
            transcript['end'].values[0]
        )

    def get_unique_transcripts(self):
        return self.gff['ID'].unique()