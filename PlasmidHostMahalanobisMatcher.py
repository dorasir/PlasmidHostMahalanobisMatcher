import argparse
import os
import sys

parser = argparse.ArgumentParser(description='PlasmidHostMatcher')
parser.add_argument('-q', dest='query_virus_dir', nargs=1, required=True,
                    help='Directory containing query virus genomes with .fasta or .fa suffix')
parser.add_argument('-o', dest='output_dir', nargs=1, required=True, help='Output directory')
parser.add_argument('-t', dest='num_Threads', nargs=1, type=int, default=[1],
                    help='Number of threads to use. Default = 1')
parser.add_argument('-n', dest='topN', metavar='topN', nargs=1, type=int, default=[1],
                    help='Number of top predictions written to the output files. All predictions will be output if '
                         'there is a tie in score. Default = 1')
parser.add_argument('-i', dest='intermediate_dir', nargs=1, default=['./intermediate_res'],
                    help='Directory storing intermediate result. Default = ./intermediate_res')
parser.add_argument('-l', dest='genome_list', nargs=1, default=[None],
                    help='Location of the file containing host NCBI genome names of interest')

args = parser.parse_args()
