import argparse
import os
import sys
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import util
import numpy as np
import pickle

parser = argparse.ArgumentParser(description="PlasmidHostMatcher")
parser.add_argument(
    "-q",
    dest="query_plasmid_dir",
    nargs=1,
    required=True,
    help="Directory containing query plasmid genomes with .fasta or .fa suffix",
)
parser.add_argument(
    "-s",
    dest="subject_host",
    nargs=1,
    required=True,
    help="Path to the downloaded host blast database",
)
parser.add_argument("-o", dest="output_dir", nargs=1, required=True, help="Output directory")
parser.add_argument(
    "-t", dest="num_threads", nargs=1, type=int, default=[2], help="Number of threads to use. Default = 1"
)
parser.add_argument(
    "-n",
    dest="topN",
    metavar="topN",
    nargs=1,
    type=int,
    default=[1],
    help="Number of top predictions written to the output files. All predictions will be output if "
    "there is a tie in score. Default = 1",
)
parser.add_argument(
    "-i",
    dest="intermediate_dir",
    nargs=1,
    default=["./intermediate_res"],
    help="Directory storing intermediate result. Default = ./intermediate_res",
)

args = parser.parse_args()


def construct_data(index, *features):
    data = []
    for feature in features:
        data.append(feature[np.arange(feature.shape[0]), index, np.newaxis])
    return np.concatenate(data, axis=1)


def calc_svpos(plasmid_plasmid_dist: np.ndarray, plasmid_host_interaction_indicator):
    svpos = plasmid_plasmid_dist.dot(plasmid_host_interaction_indicator) / plasmid_host_interaction_indicator.sum(axis=0)
    return np.nan_to_num(svpos)


if __name__ == "__main__":
    if os.path.isdir(args.subject_host):
        if util.check_directory_content:
            util.fasta_to_database(
                args.subject_host,
            )

    """Asset will be downloaded and stored at designated path"""
    TRAINING_ASSETS_PATH = ""
    training_assets = util.load_obj(TRAINING_ASSETS_PATH)

    training_plasmid_host = training_assets.plasmid_host
    training_metadata = training_assets.metadata
    training_interaction_indicator = np.load(TRAINING_INDICATOR_PATH)

    # Calculate plasmid host Mahalanobis distance
    t = MahalanobisRelativeAbundance(
        args.subject_host, args.query_plasmid_dir, temp_directory_path=args.intermediate_dir, thread=args.num_threads
    )
    plasmid_host = t.calc_distance(args.num_threads)

    # Calculate blast distance
    blast_results_dict = util.blast_dir_to_db(args.query_plasmid_dir, args.subject_host)
    blast_results = {}
    for key in blast_results_dict:
        blast_results[key.split(".")[0]] = blast_results_dict[key]
    blast_results_mat = np.zeros((len(metadata), len(set(host_list))))
    for i in range(len(metadata)):
        success, series = blast_results[metadata.Locus_ID[i]]
        if not success:
            continue
        else:
            for key in series.keys():
                idx = host_to_idx_dict[int(key[4:])]
                blast_results_mat[i, idx] = series[key]

    # Calculate test-training plasmid distance and svpos
    test_plasmid_to_train_plasmid = util.cosine_similarity(plasmid_host, training_plasmid_host)
    svpos = calc_svpos(test_plasmid_to_train_plasmid, training_interaction_indicator)

    model_path = ""
    model = pickle.load(model_path)

    idx = np.arange(plasmid_host.shape[0])
    data = construct_data(idx, plasmid_host, blast_results_mat, svpos)

    prediction = model.predict_proba(data)

    print("Complete!")
