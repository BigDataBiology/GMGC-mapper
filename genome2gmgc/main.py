import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import requests
import pandas as pd
import json
import time
from safeout import safeout
from tqdm import tqdm

from .alignment import identity_coverage
from .genome2gmgc_version import __version__

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Genome2gmgc')
    parser.add_argument('-i', '--input',required=True,help = 'path to the input genome FASTA file.',dest='genome_fasta',
                        default = None)
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default = None)
    return parser.parse_args()

def gene_prediction(fasta_input,output):
    import subprocess
    from os import path
    print('Start gene prediction...')
    subprocess.check_call(
            ['prodigal',
                '-i', fasta_input,
                '-o', path.join(output, 'gene.coords.gbk'),
                '-a', path.join(output, 'prodigal_out.faa'),
                '-d', path.join(output, 'prodigal_out.fna')])

            # For short inputs, prodigal will output the warning
            #
            # ```Warning:  ideally Prodigal should be given at least 100000 bases for training.
            # You may get better results with the -p meta option.```
            #
            # Arguably, we could check this, but we default to the non `-p meta` call

    print('\nGene prediction done.\n')


def split_file(aa_path, dna_path, output_dir, max_size=50):
    if not os.path.exists(aa_path) or not os.path.exists(dna_path):
        raise Exception("Not exist the file!")
    records = list(SeqIO.parse(aa_path, "fasta"))
    index = 0
    if len(records) > max_size:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
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
                SeqIO.write(split_fasta_aa,  output_dir + '/split_{}.faa'.format(index), 'fasta')
                SeqIO.write(split_fasta_dna, output_dir + '/split_{}.fna'.format(index), 'fasta')
                split_fasta_aa = []
                split_fasta_dna = []
        if split_fasta_aa != []:
            index += 1
            SeqIO.write(split_fasta_aa, output_dir + '/split_{}.faa'.format(index), 'fasta')

        if split_fasta_dna != []:
            SeqIO.write(split_fasta_dna, output_dir + '/split_{}.fna'.format(index), 'fasta')

    return index




def query_gmgc(fasta_file,max_try = 10):
    if not os.path.exists(fasta_file):
        raise Exception("Missing file '{}'".format(fasta_file))

    if len(list(SeqIO.parse(fasta_file, "fasta"))) == 0:
        raise Exception("Input FASTA file '{}' file is empty!".format(fasta_file))

    for try_index in range(max_try):
        try:
            parameter = {'mode': 'besthit', 'return_seqs': True}
            fasta = {'fasta': open('{}'.format(fasta_file), 'rb')}
            besthit = requests.post('http://gmgc.embl.de/api/v1.0/query/sequence', data = parameter, files = fasta)
        except Exception:
            print('genome2gmgc query failed {} times!'.format(try_index+1))
            time.sleep(60)
        else:
            return besthit
    return None




def realignment(dna_path,aa_path,besthit):
    results = []
    for record_dna, record_aa, hit_index in zip(SeqIO.parse(dna_path,'fasta'),
                                                SeqIO.parse(aa_path,'fasta'),
                                                besthit):
        hit_result = []
        query_dna = str(record_dna.seq)
        query_aa = str(record_aa.seq)
        query_name = hit_index['query_name']
        if hit_index['hits'] != []:
            hit_gene_id = hit_index['hits'][0]['gene_id']
            target_dna = hit_index['hits'][0]['dna_sequence']
            target_aa =  hit_index['hits'][0]['protein_sequence']
            if target_aa and target_dna is not None:
                realign_result = identity_coverage(query_dna,query_aa,target_dna,target_aa)
                hit_result.extend([query_name,hit_gene_id,realign_result,target_dna,target_aa])
                results.append(hit_result)
        else:
            hit_result.extend([query_name, None, 'NO HIT', None, None])
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
    num_split = split_file(out+'/prodigal_out.faa',
                           out+'/prodigal_out.fna',
                           output_dir=out+'/split_file')
    hit_table = []
    print('Starting GMGC queries (total: {} batches to process)'.format(num_split))
    for index in tqdm(range(num_split)):
        besthit = query_gmgc(out+'/split_file/split_{}.faa'.format(index+1))
        if besthit is not None:
            besthit = json.loads(bytes.decode(besthit.content))['results']

            hit_table_index = realignment(out+'/split_file/split_{}.fna'.format(index+1),
                                          out+'/split_file/split_{}.faa'.format(index+1),besthit)
            hit_table.extend(hit_table_index)
    hit_table = pd.DataFrame(hit_table)
    hit_table.columns = ['query_name','gene_id','align_category','gene_dna','gene_protein']
    num_gene = hit_table.shape[0]


    summary = []
    summary.append('*'*30+'Genome2gmgc results summary table'+'*'*30)
    summary.append('- Processed {} genes'.format(num_gene))
    match_result = hit_table['align_category'].value_counts().to_dict()
    if 'EXACT' in match_result:
        summary.append(' -{0} ({1:.1%}) were found in the GMGC at above 95% nucleotide identity with 95% coverage'
                .format(match_result['EXACT'], match_result['EXACT']/num_gene))
    else:
        summary.append(' -No genes were found in the GMGC at above 95% nucleotide identity with 95% coverage')

    if 'SIMILAR' in match_result:
        summary.append(' -{0} ({1:.1%}) were found in the GMGC at above 80% nucleotide identity with 80% coverage'
                .format(match_result['SIMILAR'], match_result['SIMILAR']/num_gene))
    else:
        summary.append(' -No genes were found in the GMGC at above 80% nucleotide identity with 80% coverage')

    if 'MATCH' in match_result:
        summary.append(' -{0} ({1:.1%}) were found in the GMGC at above 50% nucleotide identity with 50% coverage'
                .format(match_result['MATCH'], match_result['MATCH']/num_gene))
    else:
        summary.append(' -No genes were found in the GMGC at above 50% nucleotide identity with 50% coverage')

    no_match = match_result.get('NO MATCH', 0.0) + match_result.get('NO HIT', 0.0)
    if no_match:
        summary.append(' -{0} ({1:.1%}) had no match in the GMGC'
                    .format(no_match, no_match/num_gene))

    with safeout(out+'/hit_table.tsv', 'wt') as ofile:
        ofile.write('# Results from Genome2gmgc v{}\n'.format(__version__))
        hit_table.to_csv(ofile, sep='\t', index=False)
    with safeout(out+'/summary.txt', 'wt') as ofile:
        for s in summary:
            print(s)
            ofile.write(s+'\n')




if __name__ == '__main__':
    main(sys.argv)
