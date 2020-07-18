# %% Initial import
import matplotlib.pyplot as plt
from numpy.lib.function_base import average
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from tools import kmer_count
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np
import pandas as pd
import os
import shutil
import util
from ete3 import NCBITaxa

import warnings
warnings.filterwarnings("ignore")


# %% Defining taxonomic accuracy function
desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']


def taxid_to_lineage(taxid):
    """
    Function for retrieving the taxonomic rank of given taxid
    :param taxid:
    :return:
    """
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    rank_to_id = {rank: id for (id, rank) in ncbi.get_rank(lineage).items()}
    rank_to_id = {desired_rank: (rank_to_id[desired_rank] if desired_rank in rank_to_id.keys() else None) for
                  desired_rank in desired_ranks}
    return rank_to_id


def taxonomic_accuracy(prediction, target, use_max=False):
    cnt = np.zeros(7)
    for i in range(prediction.shape[0]):
        if use_max:
            rank_to_id_1 = taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[prediction[i, :].argmax()]])
        else:
            rank_to_id_1 = taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[prediction[i, :].argmin()]])
        rank_to_id_2 = taxid_to_lineage(
            host_to_speciesid[idx_to_host_dict[target[i]]])
        for j, rank in enumerate(desired_ranks):
            if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
                cnt[j] += 1
    """Strain accuracy is handled separately"""
    if use_max:
        cnt[6] = np.sum(prediction.argmax(axis=1) == target)
    else:
        cnt[6] = np.sum(prediction.argmin(axis=1) == target)
    acc = cnt / prediction.shape[0]
    return acc

# %%


def taxonomic_accuracy_threshold(prediction, target, use_max=True):
    order = np.argsort(prediction.max(axis=1))
    prediction = prediction[order, :]
    target = target[order.astype(int)]

    n = prediction.shape[0]

    prediction_correct_indicator = np.zeros((n, 7))
    for i in range(n):
        if use_max:
            rank_to_id_1 = taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[prediction[i, :].argmax()]])
        else:
            rank_to_id_1 = taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[prediction[i, :].argmin()]])
        rank_to_id_2 = taxid_to_lineage(
            host_to_speciesid[idx_to_host_dict[target[i]]])
        for j, rank in enumerate(desired_ranks):
            if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
                prediction_correct_indicator[i, j] = 1

        if use_max:
            if prediction[i, :].argmax() == target[i]:
                prediction_correct_indicator[i, 6] = 1
        else:
            if prediction[i, :].argmin() == target[i]:
                prediction_correct_indicator[i, 6] = 1

    thresholded_accuracy = np.zeros(prediction_correct_indicator.shape)
    for i in range(n):
        thresholded_accuracy[i, :] = prediction_correct_indicator[i:, :].sum(
            axis=0) / (n-i)
    return thresholded_accuracy


# %% Loading metadata


def load_metadata():
    metadata = pd.read_csv('metadata.csv')
    metadata = metadata.dropna()
    metadata = metadata.reset_index(drop=True)

    for i in range(len(metadata.Locus_ID)):
        metadata.at[i, 'Locus_ID'] = metadata.at[i, 'Locus_ID'].split('.')[0]

    metadata.sort_values(by=['Locus_ID'], inplace=True)
    metadata.reset_index(drop=True, inplace=True)
    return metadata


metadata = load_metadata()

# %% Moving plasmid


def move_plasmid(metadata, plasmid_path='data/plasmids', used_plasmid_path='data/plasmids_used'):
    """Move plasmids with host information to plasmids_used directory"""
    os.makedirs(used_plasmid_path)
    plasmid_list = list(metadata.Locus_ID)
    for fn in os.listdir(plasmid_path):
        if fn.split('.')[0] in plasmid_list:
            shutil.copyfile(os.path.join('data/plasmids', fn),
                            os.path.join('data/plasmids_used', fn))


move_plasmid(metadata)

util.remove_plasmid_seq('data/hosts', 'data/hosts_no_plasmid')


# %%
def generate_data(subject_path, query_path, save_path, recalculate=False, save=False, k=3):

    if recalculate:
        t = MahalanobisRelativeAbundance(
            subject_path, query_path, temp_directory_path='temp_dir/plasmid_host', k=k)
        plasmid_host = t.calc_distance(8)
        if save:
            np.save(save_path, plasmid_host)
        plasmid_kmer_freq = t.query_kmer_freq
        plasmid_plasmid = util.cosine_similarity(
            plasmid_kmer_freq, plasmid_kmer_freq)
    else:
        plasmid_host = np.load(save_path)
        plasmid_kmer_count = 0  # TODO


# %% Defining training parameters
k = 3
thread_num = 8

plasmid_host_class_path = 'plasmid_host.pkl'
plasmid_host_dist_path = 'plasmid_host.npy'

# %% Calculate plasmid-host distance and save result
t = MahalanobisRelativeAbundance(
    'data/hosts', 'data/plasmids_used', temp_directory_path='temp_dir/plasmid_host', k=k)
plasmid_host = t.calc_distance(thread_num)
np.save(plasmid_host_dist_path, plasmid_host)
util.save_obj(t, plasmid_host_class_path)

# %% Load calculated plasmid-host distance
plasmid_host = np.load(plasmid_host_dist_path)

# %% Normalize plasmid-host distance
plasmid_host_normalized = (plasmid_host - plasmid_host.min(axis=0)) / \
    (plasmid_host.max(axis=0) - plasmid_host.min(axis=0))

# %% Calculate plasmid-wise distance
plasmid_plasmid = util.cosine_similarity(plasmid_host, plasmid_host)


# %% Construct plasmid interaction table
host_list = list(set(metadata.Assembly_chainid))
host_list.sort()
host_to_idx_dict = {host: i for i, host in enumerate(host_list)}
idx_to_host_dict = {i: host for i, host in enumerate(host_list)}

# plasmid-strain indicator
interaction_table = np.zeros((len(metadata), len(set(host_list))))
for i in range(len(metadata)):
    interaction_table[i, host_to_idx_dict[metadata.Assembly_chainid[i]]] = 1

# %% Construct plasmid interaction table based on species
host_to_speciesid = {}
speciesid_to_host = {}
for i in range(len(metadata)):
    host_to_speciesid[metadata.Assembly_chainid[i]
                      ] = metadata.Assembly_speciestaxid[i]
    if metadata.Assembly_speciestaxid[i] not in speciesid_to_host.keys():
        speciesid_to_host[metadata.Assembly_speciestaxid[i]] = {
            metadata.Assembly_chainid[i]}
    else:
        speciesid_to_host[metadata.Assembly_speciestaxid[i]].add(
            metadata.Assembly_chainid[i])

interaction_table_species = np.zeros(
    (len(metadata), len(set(host_list))))  # plasmid-strain indicator
for i in range(len(metadata)):
    for strain in speciesid_to_host[metadata.Assembly_speciestaxid[i]]:
        interaction_table_species[i, host_to_idx_dict[strain]] = 1

# def seqacc_to_hostacc(host_dir):
#     seqacc_to_hostacc_dict = {}
#
#     from Bio import SeqIO
#     host_list = os.listdir(host_dir)
#     host_list = [h for h in host_list if h.endswith('.fna')]
#     for h in host_list:
#         host_path = os.path.join(host_dir, h)
#         host_name = h.split('.')[0]
#         for record in SeqIO.parse(host_path, 'fasta'):
#             seqacc = record.description.split()[0]
#             seqacc_to_hostacc_dict[seqacc] = host_name
#     import pickle
#     with open('seqacc_to_hostacc.pkl', 'wb') as f:
#         pickle.dump(seqacc_to_hostacc_dict, f)

# seqacc_to_hostacc('data/hosts_no_plasmid')

# seqacc_to_hostacc_dict = util.load_obj('seqacc_to_hostacc.pkl')

# %% Calculate blast
blast_results = util.blast_dir_to_db(
    'data/plasmids_used', 'data/hosts_no_plasmid_complete.fna')
util.save_obj(blast_results, 'blast_results.pkl')

# blast_results = util.blast_dir_to_db('plsdb', 'data/hosts_no_plasmid_complete.fna')
# util.save_obj(blast_results, 'blast_results_plsdb.pkl')

# %% Construct blast result matrix
blast_results_old = util.load_obj('blast_results.pkl')
blast_results = {}
for key in blast_results_old:
    blast_results[key.split('.')[0]] = blast_results_old[key]
blast_results_mat = np.zeros((len(metadata), len(set(host_list))))
for i in range(len(metadata)):
    success, series = blast_results[metadata.Locus_ID[i]]
    if not success:
        continue
    else:
        for key in series.keys():
            idx = host_to_idx_dict[int(key[4:])]
            blast_results_mat[i, idx] = series[key]

# %% Calculate svpos and svneg
plasmid_plasmid = np.load('plasmid_plasmid.npy')
svpos = plasmid_plasmid.dot(interaction_table) / interaction_table.sum(axis=0)

# %% Construct species taxid to index mapping
speciesid_to_idx_l = list(set(metadata.Assembly_speciestaxid))
speciesid_to_idx_l.sort()
speciesid_to_idx = {s: i for i, s in enumerate(speciesid_to_idx_l)}
idx_to_speciesid = {i: s for i, s in enumerate(speciesid_to_idx_l)}
"""List of index of species"""
species = [speciesid_to_idx[metadata.Assembly_speciestaxid[i]]
           for i in range(len(metadata))]

# %% Construct index of positive set and negative set for each plasmid


def generate_true_false_pairs(taxonomic_level=None):
    true_idx = [host_to_idx_dict[metadata.Assembly_chainid[i]]
                for i in range(len(metadata))]
    false_idx = []
    if taxonomic_level:
        for i in range(len(metadata)):
            idx = np.random.randint(plasmid_host.shape[1])
            while idx == true_idx[i]:
                idx = np.random.randint(plasmid_host.shape[1])
            false_idx.append(idx)
    else:
        true_family = [(taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])[
            taxonomic_level] if taxonomic_level in taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[idx]]) else None) for idx in true_idx]
        for i in range(len(metadata)):
            idx = np.random.randint(plasmid_host.shape[1])
            while idx == true_idx[i]:
                idx = np.random.randint(plasmid_host.shape[1])
            false_lineage = taxid_to_lineage(
                host_to_speciesid[idx_to_host_dict[idx]])
            false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
            while idx == true_idx[i] and false_family == true_family[i]:
                idx = np.random.randint(plasmid_host.shape[1])
                false_lineage = taxid_to_lineage(
                    host_to_speciesid[idx_to_host_dict[idx]])
                false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
            false_idx.append(idx)
    return true_idx, false_idx


true_idx, false_idx = generate_true_false_pairs()

# %% Construct positive and negative dataset
X_pos = np.concatenate([plasmid_host_normalized[np.arange(plasmid_host_normalized.shape[0]), true_idx, np.newaxis],
                        svpos[np.arange(svpos.shape[0]), true_idx, np.newaxis],
                        blast_results_mat[np.arange(blast_results_mat.shape[0]), true_idx, np.newaxis]], axis=1)
X_neg = np.concatenate([plasmid_host_normalized[np.arange(plasmid_host_normalized.shape[0]), false_idx, np.newaxis],
                        svpos[np.arange(svpos.shape[0]),
                              false_idx, np.newaxis],
                        blast_results_mat[np.arange(blast_results_mat.shape[0]), false_idx, np.newaxis]], axis=1)
X = np.concatenate((X_pos, X_neg), axis=0)
y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

# %% Plot ROC curve and calculate AUC

# %%
levelwise_similarity = False
host_similarity = np.zeros((len(host_list), len(host_list)))
host_lineage = {}
for i in range(len(host_list)):
    similarity = []
    lineage_1 = taxid_to_lineage(host_to_speciesid[host_list[i]])
    host_lineage[i] = lineage_1
    for j in range(len(host_list)):
        if j == i:
            host_similarity[i, j] = 1
            continue
        if j in host_lineage:
            lineage_2 = host_lineage[j]
        else:
            lineage_2 = taxid_to_lineage(host_to_speciesid[host_list[j]])
            host_lineage[j] = lineage_2
        rank = 'species'
        if rank in lineage_1 and rank in lineage_2 and lineage_1[rank] == lineage_2[rank]:
            host_similarity[i, j] = 1
        # Compare lineage of two hosts
        if levelwise_similarity:
            for k in range(len(desired_ranks)-1, -1, -1):
                if desired_ranks[k] in lineage_1 and desired_ranks[k] in lineage_2:
                    if lineage_1[desired_ranks[k]] == lineage_2[desired_ranks[k]]:
                        host_similarity[i, j] = 1 * (10 ** (k-6))
                        break

# %% Function for calculating svpos


def calc_svpos_training(plasmid_plasmid_dist: np.ndarray, plasmid_host_interaction_indicator, train_index, test_index):
    train_test_distance = plasmid_plasmid_dist[np.ix_(test_index, train_index)]
    train_host_indicator = plasmid_host_interaction_indicator[train_index]
    svpos = train_test_distance.dot(
        train_host_indicator) / train_host_indicator.sum(axis=0)
    return np.nan_to_num(svpos)


# %% Train and test the model


def evaluate_performance(model, true_host_idx, host_similarity, repeat=100, use_max=False, test_size=0.2, threshold=False, **kwargs):
    """
    """
    acc = []

    if len(kwargs) == 0:
        raise RuntimeError("No features specified.")

    if 'train_test_indices' in kwargs.keys():
        if len(kwargs['train_test_indices']) != repeat:
            raise RuntimeError(
                "Number of repeat and train test set split should be the same.")

    if type(true_host_idx) is list:
        true_host_idx = np.array(true_host_idx)

    acc_threshold_avg = 0
    pred_max = 0

    # Indication of if we need to calculate svpos
    svpos_indicator = False
    features = []
    for key in kwargs.keys():
        if key == 'plasmid_host' or key == 'blast':
            features.append(kwargs[key])
        elif key == 'plasmid_plasmid':
            svpos_indicator = True

    n = features[0].shape[0]

    for i in range(repeat):
        if 'train_test_indicies' in kwargs.keys():
            idx_train, idx_test = kwargs['train_test_indicies'][i]
        else:
            idx_train, idx_test = train_test_split(
                np.arange(n), test_size=test_size)

        target = true_host_idx[idx_test]

        if len(kwargs) == 1:
            if not svpos_indicator:
                acc.append(taxonomic_accuracy(
                    features[0][idx_test, :], target, use_max=use_max))
                continue
            else:
                svpos_testing = calc_svpos_training(
                    kwargs['plasmid_plasmid'], interaction_table, idx_train, idx_test)
                acc.append(taxonomic_accuracy(
                    svpos_testing[idx_test, :], target, use_max=use_max))

        combined_training_features = [
            feature[idx_train, :].flatten()[:, None] for feature in features]
        if svpos_indicator:
            svpos_training = calc_svpos_training(
                kwargs['plasmid_plasmid'], interaction_table, idx_train, idx_train)
            combined_training_features.append(
                svpos_training.flatten()[:, None])

        combined_training_features = np.hstack(combined_training_features)
        # y_train = [host_similarity[int(true_host_idx[i])]
        #            for i in true_host_idx[idx_train]]
        y_train = [host_similarity[int(true_host_idx[i])] for i in idx_train]
        y_train = np.vstack(y_train).flatten()
        model.fit(combined_training_features, y_train)

        combined_testing_features = [
            feature[idx_test, :].flatten()[:, None] for feature in features]
        if svpos_indicator:
            svpos_testing = calc_svpos_training(
                kwargs['plasmid_plasmid'], interaction_table, idx_train, idx_test)
            combined_testing_features.append(svpos_testing.flatten()[:, None])

        combined_testing_features = np.hstack(combined_testing_features)
        pred = model.predict_proba(combined_testing_features)[
            :, 1].reshape((-1, features[0].shape[1]))

        if threshold:
            acc_threshold = taxonomic_accuracy_threshold(
                pred, target, use_max=True)
            acc_threshold_avg = (acc_threshold_avg * i +
                                 acc_threshold) / (i + 1)
            pred_max = (pred_max * i + pred.max(axis=1)) / (i + 1)
        else:
            acc.append(taxonomic_accuracy(pred, target, use_max=use_max))
    if threshold:
        return acc_threshold_avg, pred_max
    else:
        acc = np.vstack(acc)
        return acc


# %%
model = LogisticRegression(class_weight='balanced')
acc, pred_max = evaluate_performance(model, true_idx, host_similarity, repeat=10, use_max=True,
                                     plasmid_host=plasmid_host_normalized, blast=blast_results_mat, plasmid_plasmid=plasmid_plasmid, threshold=True)

# %%
taxonomic_level = ['phylum', 'class', 'order',
                   'family', 'genus', 'species', 'strain']
for i in range(7):
    plt.plot(np.sort(pred_max), acc[:, i], label=taxonomic_level[i])
plt.plot(np.sort(pred_max), np.arange(
    len(pred_max), 0, -1) / len(pred_max), label='recall')
plt.xlabel('Prediction score threshold')
plt.ylabel('Accuracy and recall')
plt.grid()
plt.legend(loc=(1.04, 0.5))
plt.savefig('threshold.pdf', bbox_inches='tight')


# %%
model = LogisticRegression(class_weight='balanced')
acc = evaluate_performance(model, true_idx, host_similarity, repeat=2, use_max=True,
                           plasmid_host=plasmid_host_normalized, blast=blast_results_mat, plasmid_plasmid=plasmid_plasmid)

# %%


def train_model(model, true_host_idx, host_similarity, **kwargs):
    svpos_indicator = False
    features = []
    for key in kwargs.keys():
        if key == 'plasmid_host' or key == 'blast' or key == 'svpos':
            features.append(kwargs[key])
        elif key == 'plasmid_plasmid':
            svpos_indicator = True

    combined_features = [feature.flatten()[:, None] for feature in features]
    combined_features = np.hstack(combined_features)

    y = [host_similarity[int(i)] for i in true_host_idx]
    y = np.vstack(y).flatten()
    model.fit(combined_features, y)

    return model

# %%


model = LogisticRegression(class_weight='balaced', n_jobs=8)
model_trained = train_model(model, true_idx, host_similarity,
                            plasmid_host=plasmid_host_normalized, blast=blast_results_mat, svpos=svpos)
features = [plasmid_host_normalized, blast_results_mat, svpos]

prediction = model_trained.predict_proba(features)[:, 1].reshape((-1, 2705))

# %%
acc_threshold = taxonomic_accuracy_threshold(
    prediction, np.array(true_idx), use_max=True)

# %%
taxonomic_level = ['phylum', 'class', 'order',
                   'family', 'genus', 'species', 'strain']
for i in range(7):
    plt.plot(np.sort(prediction.max(axis=1)),
             acc_threshold[:, i], label=taxonomic_level[i])
plt.plot(np.sort(prediction.max(axis=1)), np.arange(
    6510, 0, -1) / 6510, label='recall')
plt.legend(loc=(1.04, 0.5))
plt.show()


# %%
# import matplotlib.pyplot as plt
# for i in range(6):
#     plt.figure(i)
#     plt.boxplot([net_blast_acc[:, i], net_no_blast_acc[:, i], mah_acc[:, i], svpos_acc[:, i], blast_acc[:, i]], labels=['Net-blast', 'No-blast', 'Mah', 'Svpos', 'Blast'], showmeans=True)
#     plt.title(desired_ranks[i])
#     plt.savefig(desired_ranks[i]+'.png')
#
# plt.figure(6)
# plt.boxplot([net_blast_acc[:, 6], net_no_blast_acc[:, 6], mah_acc[:, 6], svpos_acc[:, 6], blast_acc[:, 6]], labels=['Net-blast', 'No-blast', 'Mah', 'Svpos', 'Blast'], showmeans=True)
# plt.title('strain')
# plt.savefig('strain.png')
