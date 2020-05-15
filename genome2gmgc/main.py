import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import requests
from .alignment import identity_coverage
import subprocess
import pandas as pd



def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Genome2gmgc')
    parser.add_argument('-i', '--input',required=True,help = 'path to the input genome FASTA file.',dest='genome_fasta',
                        default = None)
    parser.add_argument('-o', '--output',required=True,help = 'path to the output file.',dest='output',
                        default = None)
    return parser.parse_args()

def gene_prediction(fasta_input,output):
    os.system('prodigal -i {0} -o {1}/gene.coords.gbk -a {1}/protein_gene.faa -d {1}/dna_gene.faa -p meta'
              .format(fasta_input,output))


def split_file(aa_path, dna_path, output_file, max_size=50):
    if not os.path.exists(aa_path) or not os.path.exists(dna_path):
        raise Exception("Not exist the file!")
    records = list(SeqIO.parse(aa_path, "fasta"))
    index = 0
    if len(records) > max_size:
        if not os.path.exists(output_file):
            os.makedirs(output_file)
        split_fasta_aa = []
        split_fasta_dna = []
        num_seq = 0
        for seq_record_aa , seq_record_dna in zip(SeqIO.parse(aa_path,'fasta'),
                                                       SeqIO.parse(dna_path,'fasta')):
            rec1 = SeqRecord(Seq(str(seq_record_aa.seq)), id=seq_record_aa.id, description='')
            split_fasta_aa.append(rec1)
            rec1 = SeqRecord(Seq(str(seq_record_dna.seq)), id=seq_record_dna.id, description='')
            split_fasta_dna.append(rec1)
            num_seq += 1
            if num_seq == max_size:
                num_seq = 0
                index += 1
                SeqIO.write(split_fasta_aa,output_file + '/protein_split_{}.fna'.format(index),'fasta')
                SeqIO.write(split_fasta_dna, output_file + '/dna_split_{}.fna'.format(index), 'fasta')
                split_fasta_aa = []
                split_fasta_dna = []
        if split_fasta_aa != []:
            index += 1
            SeqIO.write(split_fasta_aa,output_file + '/protein_split_{}.fna'.format(index), 'fasta')

        if split_fasta_dna != []:
            SeqIO.write(split_fasta_dna,output_file + '/dna_split_{}.fna'.format(index), 'fasta')

    return index




def query_gmgc(fasta_file):
    if not os.path.exists(fasta_file):
        raise Exception("Missing file '{}'".format(fasta_file))

    if len(list(SeqIO.parse(fasta_file, "fasta"))) == 0:
        raise Exception("Input FASTA file '{}' file is empty!".format(fasta_file))

    besthit = subprocess.getstatusoutput('curl  -s -F \'mode=besthit\' -F \'return_seqs=0\' -F \'fasta=@{}\'  '
              '\'http://gmgc.embl.de/api/v1.0/query/sequence\''.format(fasta_file))[1]
    besthit = eval(besthit)['results']
    return besthit


def realignment(dna_path,aa_path,besthit):
    results = []
    for record_dna, record_aa, hit_index in zip(SeqIO.parse(dna_path,'fasta'),
                                                SeqIO.parse(aa_path,'fasta'),
                                                besthit):
        hit_result = []
        query_dna = str(record_dna.seq)
        query_aa = str(record_aa.seq)
        query_name = hit_index['query_name']
        if hit_index['hits'] == []:
            continue
        hit_gene_id = hit_index['hits'][0]['gene_id']
        try:
            target_dna_request = requests.get('http://gmgc.embl.de/api/v1.0/unigene/{}/dna_sequence'.
                                             format(hit_gene_id))
        except Exception:
            raise Exception('Http request error')
        target_dna = eval(target_dna_request.content)['dna_sequence']

        try:
            target_aa_request = requests.get('http://gmgc.embl.de/api/v1.0/unigene/{}/protein_sequence'.
                                             format(hit_gene_id))
        except Exception:
            raise Exception('Http request error')
        target_aa = eval(target_aa_request.content)['protein_sequence']
        if target_aa and target_dna is not None:
            realign_result = identity_coverage(query_dna,query_aa,target_dna,target_aa)
            hit_result.extend([query_name,hit_gene_id,realign_result,target_dna,target_aa])
            results.append(hit_result)
    return results


def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    out = args.output
    if not os.path.exists(out):
        os.makedirs(out)
    gene_prediction(args.genome_fasta,out)
    num_split = split_file(out+'/protein_gene.faa',out+'/dna_gene.faa',
                           output_file=out+'/split_file')
    hit_table = []
    for index in range(num_split):
        print(index)
        besthit = query_gmgc(out+'/split_file/protein_split_{}.fna'.format(index+1))
        hit_table_index = realignment(out+'/split_file/dna_split_{}.fna'.format(index+1),
                                      out+'/split_file/protein_split_{}.fna'.format(index+1),besthit)
        hit_table.extend(hit_table_index)
    hit_table = pd.DataFrame(hit_table)
    hit_table.columns = ['query_name','gene_id','align_category','gene_dna','gene_protein']
    hit_table.to_csv(out+'/hit_table.csv')


if __name__ == '__main__':
    main(sys.argv)
