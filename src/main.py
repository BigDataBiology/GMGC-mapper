import argparse
import sys

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Genome2gmgc')
    parser.add_argument('-i', '--input',required=True,help = 'path to the input genome FASTA file.',dest='genome_fasta',
                        default = None)

    return parser.parse_args()

def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    print(args.genome_fasta)


if __name__ == '__main__':

    main(sys.argv)