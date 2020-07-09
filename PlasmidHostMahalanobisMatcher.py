import argparse
import os
import sys
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance

parser = argparse.ArgumentParser(description='PlasmidHostMatcher')
parser.add_argument('-q', dest='query_plasmid_dir', nargs=1, required=True,
                    help='Directory containing query plasmid genomes with .fasta or .fa suffix')
parser.add_argument('-s', dest='subject_host_dir', nargs=1, required=True,
                    help='Directory containing subject host genomes with .fasta or .fa suffix')
parser.add_argument('-o', dest='output_dir', nargs=1,
                    required=True, help='Output directory')
parser.add_argument('-t', dest='num_threads', nargs=1, type=int, default=[1],
                    help='Number of threads to use. Default = 1')
parser.add_argument('-n', dest='topN', metavar='topN', nargs=1, type=int, default=[1],
                    help='Number of top predictions written to the output files. All predictions will be output if '
                         'there is a tie in score. Default = 1')
parser.add_argument('-i', dest='intermediate_dir', nargs=1, default=['./intermediate_res'],
                    help='Directory storing intermediate result. Default = ./intermediate_res')

args = parser.parse_args()

if __name__ == "__main__":
    plasmid_dir = args.query_plasmid_dir
    host_dir = args.subject_host_dir
    temp_dir = args.intermediate_dir

    # Calculate plasmid host Mahalanobis distance
    t = MahalanobisRelativeAbundance(args.subject_host_dir, args.query_plasmid_dir,
                                     temp_directory_path=args.intermediate_dir, thread=args.num_threads)

    # Calculate blast distance

    # Calculate test-training plasmid distance and svpos
