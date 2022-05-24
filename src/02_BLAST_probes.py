import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import json



def transcriptome_blast(query_fasta, organism='House mouse'):

    '''
    BLAST a multi fasta of probe candidates against a user-selected transcriptome, returns \
       a single JSON contaning the BLAST results 
    '''

    cline = NcbiblastnCommandline(query=query_fasta, task='blastn-short', db="refseq_rna", evalue=1000, \
        gapopen = 5, gapextend = 2, out="./blast/results/blast_results.json", outfmt=15, entrez_query='"{} [ORGN]"'.format(organism), \
             remote=True, reward = 1, penalty = -3, word_size = 7, strand = 'minus') 

    cline()
    return None

def parse_probes(fpath):

    '''
    A function to assemble probe candidates into a BLASTable format
    '''

    db = pd.read_csv(fpath)

    records = []

    for kmer in db.to_dict(orient="records"):

        
        record_1 = SeqRecord(seq=Seq(kmer['HCR3-P1']), id=kmer['name']+'-P1')
        record_2 = SeqRecord(seq=Seq(kmer['HCR3-P2']), id=kmer['name']+'-P2')
        records.append(record_1)
        records.append(record_2)


    SeqIO.write(records, './blast/query/blast_query.fasta','fasta')

def parse_results(fpath):

    output_1 = []
    output_2 = []

    with open(fpath) as f_in:

        results = json.load(f_in)
        
        for probe in range(len(results['BlastOutput2'])):

            kmer = results['BlastOutput2'][probe]['report']['results']['search']['query_title']
            hits = results['BlastOutput2'][probe]['report']['results']['search']['hits']
            align_lengths = [{hit['description'][0]['title']:hit['hsps'][0]['align_len']} for hit in hits if hit['hsps'][0]['qseq'] == hit['hsps'][0]['hseq']]

            if len(align_lengths) == 0:
                align_lengths = ['NA']
        
            if '-P1' in kmer:
                output_1.append(align_lengths[:20])
            else:
                output_2.append(align_lengths[:20])

    return [output_1, output_2]


if __name__ == "__main__":

    query_path = Path('./blast/query')
    results_path = Path('./blast/results')
    query_path.mkdir(parents=True, exist_ok=True)
    results_path.mkdir(parents=True, exist_ok=True)
    parse_probes('./probes/non_overlap_probes.csv')
    transcriptome_blast('./blast/query/blast_query.fasta', organism=)

    probes_db = pd.read_csv('./probes/non_overlap_probes.csv')
    probes_db['BLAST_1'] = parse_results('./blast/results/blast_results.json')[0]
    probes_db['BLAST_2'] = parse_results('./blast/results/blast_results.json')[1]
    probes_db.to_csv('./probes/BLAST_probes.csv')