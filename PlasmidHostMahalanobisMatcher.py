import argparse
import os
import sys
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import util
import numpy as np
import pickle
import copy

parser = argparse.ArgumentParser(description="PlasmidHostMatcher")
parser.add_argument(
    "-q",
    dest="query_plasmid_dir",
    nargs=1,
    required=True,
    help="Directory containing query plasmid genomes with .fasta or .fa suffix",
)
# parser.add_argument(
#     "-s",
#     dest="subject_host",
#     nargs=1,
#     required=True,
#     help="Path to the downloaded host blast database",
# )
parser.add_argument("-o", dest="output", nargs=1, required=True, default=['output.txt'], help="Output file")
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


def calc_svpos(plasmid_plasmid_dist: np.ndarray, plasmid_host_interaction_indicator):
    svpos = plasmid_plasmid_dist.dot(plasmid_host_interaction_indicator) / plasmid_host_interaction_indicator.sum(
        axis=0
    )
    return np.nan_to_num(svpos)


if __name__ == "__main__":

    query_plasmid_dir = args.query_plasmid_dir[0]
    intermediate_dir = args.intermediate_dir[0]
    num_threads = args.num_threads[0]
    topN = args.topN[0]
    output = args.output[0]

    """Asset will be downloaded and stored at designated path"""
    TRAINING_HOST_FASTA_DIR_PATH = "data_exp/hosts_no_plasmid"
    TRAINING_HOST_DB_PATH = "data/hosts_no_plasmid_complete.fna"

    TRAINING_ASSETS_PATH = "data/training_assets.pkl"
    training_assets = util.load_obj(TRAINING_ASSETS_PATH)

    training_plasmid_host = copy.copy(training_assets["plasmid_host"])
    training_metadata = copy.copy(training_assets["metadata"])
    training_interaction_indicator = copy.copy(training_assets["interaction_indicator"])

    host_list = list(set(training_metadata.Assembly_chainid))
    host_list.sort()
    host_to_idx_dict = {host: i for i, host in enumerate(host_list)}
    idx_to_host_dict = {i:f'GCF_{int(host):09}' for i, host in enumerate(host_list)}

    query_list = os.listdir(query_plasmid_dir)
    query_list.sort()
    query_list = [q[:-5] for q in query_list]

    # Calculate plasmid host Mahalanobis distance
    t = MahalanobisRelativeAbundance(
        TRAINING_HOST_FASTA_DIR_PATH,
        query_plasmid_dir,
        temp_directory_path=intermediate_dir,
        thread=num_threads,
    )
    plasmid_host = t.calc_distance(num_threads)
    plasmid_host[plasmid_host > 1000] = 1000
    plasmid_host_normalized = (plasmid_host - plasmid_host.min(axis=0)) / (
        plasmid_host.max(axis=0) - plasmid_host.min(axis=0)
    )

    # Calculate blast distance
    blast_results_dict = util.blast_dir_to_db(query_plasmid_dir, TRAINING_HOST_DB_PATH)
    blast_results = {}
    for key in blast_results_dict:
        blast_results[key.split(".")[0]] = blast_results_dict[key]
    blast_results_mat = np.zeros((len(query_list), len(set(host_list))))
    for i in range(blast_results_mat.shape[0]):
        success, series = blast_results[query_list[i]]
        if not success:
            continue
        else:
            for key in series.keys():
                idx = host_to_idx_dict[int(key[4:])]
                blast_results_mat[i, idx] = series[key]

    # Calculate test-training plasmid distance and svpos
    # test_plasmid_to_train_plasmid = util.cosine_similarity(plasmid_host_normalized, training_plasmid_host[:, :6])
    test_plasmid_to_train_plasmid = util.cosine_similarity(plasmid_host_normalized, training_plasmid_host)
    svpos = calc_svpos(test_plasmid_to_train_plasmid, training_interaction_indicator)

    model_path = "data/model.pkl"
    model = util.load_obj(model_path)

    idx = np.arange(plasmid_host.shape[0])
    # features = [plasmid_host_normalized, blast_results_mat[:, :6], svpos[:, :6]]
    features = [plasmid_host_normalized, blast_results_mat, svpos]
    combined_features = [feature.flatten()[:, None] for feature in features]
    combined_features = np.hstack(combined_features)

    prediction = model.predict_proba(combined_features)
    prediction = prediction[:, 1].reshape((-1, features[0].shape[1]))

    prediction_max = (-prediction).argsort(axis=1)
    prediction_name = [[idx_to_host_dict[prediction_max[i, j]] for j in range(topN)] for i in range(prediction.shape[0])]

    with open(output, 'w') as f:
        for i, name in enumerate(query_list):
            f.write(f'Prediction for {name}:\n')
            for j in range(topN):
                f.write(f'{prediction_name[i][j]}\t')
            f.write('\n')

    print("Complete!")
