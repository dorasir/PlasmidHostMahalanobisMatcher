from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np
import pandas as pd
import os
import shutil
import util
from sklearn.neural_network import MLPClassifier

import warnings

warnings.filterwarnings("ignore")

from ete3 import NCBITaxa

from tools import kmer_count


def taxid_to_lineage(taxid):
    """
    Function for retrieving the taxonomic rank of given taxid
    :param taxid:
    :return:
    """
    desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    rank_to_id = {rank: id for (id, rank) in ncbi.get_rank(lineage).items()}
    rank_to_id = {desired_rank: (rank_to_id[desired_rank] if desired_rank in rank_to_id.keys() else None) for
                  desired_rank in desired_ranks}
    return rank_to_id


def taxonomic_accuracy(prediction, target, use_max=False):
    desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    cnt = np.zeros(7)
    for i in range(prediction.shape[0]):
        if use_max:
            rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[prediction[i, :].argmax()]])
        else:
            rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[prediction[i, :].argmin()]])
        rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[target[i]]])
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


def move_plasmid(metadata, plasmid_path, used_plasmid_path):
    """Move plasmids with host information to plasmids_used directory"""
    os.makedirs('data/plasmids_used')
    plasmid_list = list(metadata.Locus_ID)
    for fn in os.listdir('data/plasmids'):
        if fn.split('.')[0] in plasmid_list:
            shutil.copyfile(os.path.join('data/plasmids', fn), os.path.join('data/plasmids_used', fn))

# util.remove_plasmid_seq('data/hosts', 'data/hosts_no_plasmid')


def generate_data(subject_path, query_path, save_path, recalculate=False, save=False, k=3):

    if recalculate:
        t = MahalanobisRelativeAbundance(subject_path, query_path, temp_directory_path='temp_dir/plasmid_host', k=k)
        plasmid_host = t.calc_distance(8)
        if save:
            np.save(save_path, plasmid_host)
        plasmid_kmer_freq = t.query_kmer_freq
        plasmid_plasmid = util.cosine_similarity(plasmid_kmer_freq, plasmid_kmer_freq)
    else:
        plasmid_host = np.load(save_path)
        plasmid_kmer_count = #TODO



def combine_data()


def evaluate_model(model, X, y):
    pass


"""Training parameters"""


"""
Calculate plasmid-host distance
"""
# t = MahalanobisRelativeAbundance('data/hosts', 'data/plasmids_used', temp_directory_path='temp_dir/plasmid_host')
# plasmid_host = t.calc_distance(8)
# np.save('plasmid_host.npy', dist)

plasmid_host = np.load('plasmid_host.npy')
# plasmid_host = np.load('d2star.npy')

"""
Calculate plasmid-wise distance
"""
plasmid_plasmid = util.cosine_similarity(plasmid_host, plasmid_host)
# np.save('plasmid_plasmid.npy', plasmid_plasmid)

"""
Construct plasmid interaction table
"""

host_list = list(set(metadata.Assembly_chainid))
host_list.sort()
host_to_idx_dict = {host: i for i, host in enumerate(host_list)}
idx_to_host_dict = {i: host for i, host in enumerate(host_list)}

interaction_table = np.zeros((len(metadata), len(set(host_list))))  # plasmid-strain indicator
for i in range(len(metadata)):
    interaction_table[i, host_to_idx_dict[metadata.Assembly_chainid[i]]] = 1

"""
Construct plasmid interaction table based on species
"""
host_to_speciesid = {}
speciesid_to_host = {}
for i in range(len(metadata)):
    host_to_speciesid[metadata.Assembly_chainid[i]] = metadata.Assembly_speciestaxid[i]
    if metadata.Assembly_speciestaxid[i] not in speciesid_to_host.keys():
        speciesid_to_host[metadata.Assembly_speciestaxid[i]] = {metadata.Assembly_chainid[i]}
    else:
        speciesid_to_host[metadata.Assembly_speciestaxid[i]].add(metadata.Assembly_chainid[i])

interaction_table_species = np.zeros((len(metadata), len(set(host_list))))  # plasmid-strain indicator
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

"""Calculate blast"""
# blast_results = util.blast_dir_to_db('data/plasmids_used', 'data/hosts_no_plasmid_complete.fna')
# util.save_obj(blast_results, 'blast_results.pkl')

# blast_results = util.blast_dir_to_db('plsdb', 'data/hosts_no_plasmid_complete.fna')
# util.save_obj(blast_results, 'blast_results_plsdb.pkl')

"""
Construct blast result matrix
"""
blast_results_old = util.load_obj('blast_results.pkl')
# blast_results_old = util.load_obj('blast_results_plsdb.pkl')
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

"""
Construct species taxid to index mapping
"""
speciesid_to_idx_l = list(set(metadata.Assembly_speciestaxid))
speciesid_to_idx_l.sort()
speciesid_to_idx = {s: i for i, s in enumerate(speciesid_to_idx_l)}
idx_to_speciesid = {i: s for i, s in enumerate(speciesid_to_idx_l)}
"""List of index of species"""
species = [speciesid_to_idx[metadata.Assembly_speciestaxid[i]] for i in range(len(metadata))]

"""Construct index of positive set and negative set for each plasmid"""
taxonomic_level = 'phylum'
true_idx = [host_to_idx_dict[metadata.Assembly_chainid[i]] for i in range(len(metadata))]
false_idx = []
true_family = [(taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])[
                    taxonomic_level] if taxonomic_level in taxid_to_lineage(
    host_to_speciesid[idx_to_host_dict[idx]]) else None) for idx in true_idx]
for i in range(len(metadata)):
    idx = np.random.randint(plasmid_host.shape[1])
    while idx == true_idx[i]:
        idx = np.random.randint(plasmid_host.shape[1])
    # false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
    # false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
    # while idx == true_idx[i] and false_family == true_family[i]:
    #     idx = np.random.randint(plasmid_host.shape[1])
    #     false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
    #     false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
    false_idx.append(idx)

# false_idx = []
# true_family = [(taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])[taxonomic_level] if taxonomic_level in taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]]) else None) for idx in true_idx]
# for i in range(len(metadata)):
#     idxs = []
#     for j in range(1000):
#         idx = np.random.randint(plasmid_host.shape[1])
#         while idx == true_idx[i]:
#             idx = np.random.randint(plasmid_host.shape[1])
#         # false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
#         # false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
#         # while idx == true_idx[i] and false_family == true_family[i]:
#         #     idx = np.random.randint(plasmid_host.shape[1])
#         #     false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
#         #     false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
#         idxs.append(idx)
#     false_idx.append(idxs)

"""calculate svpos and svneg"""
plasmid_plasmid = np.load('plasmid_plasmid.npy')
svpos = plasmid_plasmid.dot(interaction_table) / interaction_table.sum(axis=0)
svneg = plasmid_plasmid.dot((1 - interaction_table)) / (1 - interaction_table).sum(axis=0)

"""Normalize plasmid-host distance"""
plasmid_host_normalized = (plasmid_host - plasmid_host.min(axis=0)) / (plasmid_host.max(axis=0) - plasmid_host.min(axis=0))

"""Construct positive and negative dataset"""
X_pos = np.concatenate([plasmid_host_normalized[np.arange(plasmid_host_normalized.shape[0]), true_idx, np.newaxis],
                        svpos[np.arange(svpos.shape[0]), true_idx, np.newaxis],
                        svneg[np.arange(svpos.shape[0]), true_idx, np.newaxis],
                        blast_results_mat[np.arange(blast_results_mat.shape[0]), true_idx, np.newaxis]], axis=1)
X_neg = np.concatenate([plasmid_host_normalized[np.arange(plasmid_host_normalized.shape[0]), false_idx, np.newaxis],
                        svpos[np.arange(svpos.shape[0]), false_idx, np.newaxis],
                        svneg[np.arange(svpos.shape[0]), false_idx, np.newaxis],
                        blast_results_mat[np.arange(blast_results_mat.shape[0]), false_idx, np.newaxis]], axis=1)
# X = np.concatenate((X_pos, X_neg), axis=0)
#
# y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))
#
# X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

"""Split """
from sklearn.model_selection import train_test_split

# for i in range(10):
#     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)
#     glm_binom = sm.GLM(y_train, X_train, family=sm.families.Binomial())
#     res = glm_binom.fit()


# data = np.concatenate((X, y[:, None]), axis=1)
# data_train, data_test = train_test_split(data, test_size=0.1)
# glm_binom = sm.GLM(data_train[:, -1], data_train[:, :-1], family=sm.families.Binomial())
# res = glm_binom.fit(max_iter=200)
#
# data = np.concatenate((plasmid_host_normalized.flatten()[:, None], svpos.flatten()[:, None], svneg.flatten()[:, None], blast_results_mat.flatten()[:, None]), axis=1)
# data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)
#
# pred = res.predict(data)
# pred = pred.reshape(plasmid_host_normalized.shape)
#
#
#
#
#
# acc = taxonomic_accuracy(pred, true_idx)
# print(acc)

# cnt = np.zeros(7)
# for i in range(pred.shape[0]):
#     rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[pred[i, :].argmin()]])
#     rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[true_idx[i]]])
#     for j, rank in enumerate(desired_ranks):
#         if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
#             cnt[j] += 1
# print(cnt / pred.shape[0])

# from sklearn.neural_network import MLPRegressor
# from sklearn.svm import SVC, SVR, NuSVR, LinearSVR
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

#
net_acc = []
mah_acc = []
nn_acc = []
for i in range(100):
    idx_train, idx_test = train_test_split(np.arange(plasmid_host_normalized.shape[0]), test_size=0.2)

    true_idx = np.array(true_idx)
    false_idx = np.array(false_idx)
    # false_idx = np.vstack(false_idx)
    X_pos = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]
                               , svpos[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis],
                            blast_results_mat[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
    # X_pos = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis],
    #                         svpos[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
    # X_pos = np.tile(X_pos, (1000, 1))
    X_neg = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis],
                            svpos[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis],
                            blast_results_mat[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis]], axis=1)
    # X_neg = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), false_idx[idx_train, 0], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train, 0], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), false_idx[idx_train, 0], np.newaxis]], axis=1)
    # X_neg =np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), false_idx[idx_train, 0], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train, 0], np.newaxis]], axis=1)
    # for j in range(1, 10):
    #     X_temp = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), false_idx[idx_train, j], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train, j], np.newaxis]], axis=1)
    #     X_neg = np.concatenate((X_neg, X_temp), axis=0)
    X = np.concatenate((X_pos, X_neg), axis=0)

    y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

    X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

    # glm_binom = sm.GLM(y, X, family=sm.families.Binomial())
    # res = glm_binom.fit()

    # data = np.concatenate((plasmid_host_normalized[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None], blast_results_mat[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((plasmid_host_normalized[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

    # pred = res.predict(data)
    # pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))

    t = true_idx[idx_test]

    print(i)
    # acc = taxonomic_accuracy(pred, t, use_max=True)
    # net_acc.append(acc)
    # print(acc)

    model = MLPClassifier(hidden_layer_sizes=(5,), activation='logistic')
    # model = RandomForestClassifier()
    model.fit(X, y)

    pred = model.predict_proba(data)
    pred = pred[:, 1]
    pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))
    # t = true_idx[idx_test]

    acc = taxonomic_accuracy(pred, t, use_max=True)
    # acc = np.sum(pred[np.arange(pred.shape[0]), t] == 1) / pred.shape[0]
    nn_acc.append(acc)
    print(acc)

    acc = taxonomic_accuracy(plasmid_host_normalized[idx_test, :], t)
    # acc = taxonomic_accuracy(svpos[idx_test, :], t, use_max=True)
    # acc = taxonomic_accuracy(blast_results_mat[idx_test, :], t, use_max=True)
    mah_acc.append(acc)
    print(acc)
net_acc = np.vstack(net_acc)
nn_acc = np.vstack(nn_acc)
mah_acc = np.vstack(mah_acc)
#
# for i in range(10):
#     idx_train, idx_test = train_test_split(np.arange(plasmid_host_normalized.shape[0]), test_size=0.2)
#
#     true_idx = np.array(true_idx)
#     false_idx = np.array(false_idx)
#     X_pos = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
#     X_neg = np.concatenate([plasmid_host_normalized[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis]], axis=1)
#     X = np.concatenate((X_pos, X_neg), axis=0)
#
#     y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))
#
#     X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)
#
#     # model = MLPRegressor(hidden_layer_sizes=(50, 10,), activation='logistic', max_iter=1000, solver='lbfgs')
#     model = MLPRegressor(hidden_layer_sizes=(50, 10,), activation='logistic', max_iter=1000,)
#     model.fit(X, y)
#
#     data = np.concatenate((plasmid_host_normalized[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None], svneg[idx_test, :].flatten()[:, None], blast_results_mat[idx_test, :].flatten()[:, None]), axis=1)
#     data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)
#
#     pred = model.predict(data)
#     pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))
#
#     t = true_idx[idx_test]
#
#     print(i)
#     acc = taxonomic_accuracy(pred, t, use_max=True)
#     print(acc)
#     acc = taxonomic_accuracy(plasmid_host_normalized[idx_test, :], t)
#     print(acc)

host_similarity = np.zeros((len(host_list), len(host_list)))
# for i in range(len(host_list)):
#     host_similarity[i, i] = 1
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
        # for k in range(len(desired_ranks)-1, -1, -1):
        #     if desired_ranks[k] in lineage_1 and desired_ranks[k] in lineage_2:
        #         if lineage_1[desired_ranks[k]] == lineage_2[desired_ranks[k]]:
        #             # similarity.append(9 * (10 ** (k-6)))
        #             host_similarity[i, j] = 1 * (10 ** (k-6))
        #             break

from sklearn.linear_model import LogisticRegression

net_blast_acc = []
net_no_blast_acc = []
mah_acc = []
svpos_acc = []
blast_acc = []
for i in range(100):
    idx_train, idx_test = train_test_split(np.arange(plasmid_host_normalized.shape[0]), test_size=0.2)

    true_idx = np.array(true_idx)
    false_idx = np.array(false_idx)

    data_train_blast = np.concatenate((plasmid_host_normalized[idx_train, :].flatten()[:, None],
                                       svpos[idx_train, :].flatten()[:, None],
                                       blast_results_mat[idx_train, :].flatten()[:, None]), axis=1)
    data_train = np.concatenate((plasmid_host_normalized[idx_train, :].flatten()[:, None], svpos[idx_train, :].flatten()[:, None]),
                                axis=1)
    data_train = np.concatenate((data_train, np.ones((data_train.shape[0], 1))), axis=1)
    y = [host_similarity[int(true_idx[i])] for i in true_idx[idx_train]]
    y = np.vstack(y).flatten()

    # clf = LinearSVR()
    # clf.fit(data_train, y)
    # clf = LogisticRegression(class_weight='balanced', tol=1e-10, max_iter=1000, n_jobs=4)
    clf_blast = LogisticRegression(class_weight='balanced', n_jobs=4)
    clf_blast.fit(data_train_blast, y)
    clf = LogisticRegression(class_weight='balanced', n_jobs=4)
    clf.fit(data_train, y)

    data_blast = np.concatenate((plasmid_host_normalized[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None],
                                 blast_results_mat[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((plasmid_host_normalized[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

    # pred = clf.predict(data)
    # pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))
    pred_blast = clf_blast.predict_proba(data_blast)
    pred_blast = pred_blast[:, 1]
    pred_blast = pred_blast.reshape((-1, plasmid_host_normalized.shape[1]))
    pred = clf.predict_proba(data)
    pred = pred[:, 1]
    pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))
    # pred = reg.predict(data)
    # pred = pred.reshape((-1, plasmid_host_normalized.shape[1]))

    # clf_comb = LogisticRegression(class_weight='balanced', n_jobs=4)
    # pred_train_blast = clf_blast.predict_proba(data_train_blast)[:, 1].reshape((-1, plasmid_host_normalized.shape[1]))
    # pred_train = clf.predict_proba(data_train)[:, 1].reshape((-1, plasmid_host_normalized.shape[1]))
    # X_pos = np.concatenate([pred_train_blast[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis],
    #                         pred_train[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
    # X_neg = np.concatenate([pred_train_blast[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis],
    #                         pred_train[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis]], axis=1)
    # pred_train_comb = np.concatenate((X_pos, X_neg), axis=0)
    # y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))
    # # pred_train_comb = np.concatenate((pred_train_blast.flatten()[:, None], pred_train.flatten()[:, None]), axis=1)
    # clf_comb.fit(pred_train_comb, y)
    # pred_comb = clf_comb.predict_proba(np.concatenate((pred_blast.flatten()[:, None], pred.flatten()[:, None]), axis=1))
    # pred_comb = pred_comb[:, 1]
    # pred_comb = pred_comb.reshape((-1, plasmid_host_normalized.shape[1]))

    t = true_idx[idx_test]

    # print(i)
    # acc = taxonomic_accuracy(pred_comb, t, use_max=True)
    # print(acc)

    acc = taxonomic_accuracy(pred_blast, t, use_max=True)
    net_blast_acc.append(acc)
    print("Network w/ blast: ", acc)
    acc = taxonomic_accuracy(pred, t, use_max=True)
    # acc = np.sum(pred[np.arange(pred.shape[0]), t] == 1) / pred.shape[0]
    net_no_blast_acc.append(acc)
    print("Network: ", acc)
    acc = taxonomic_accuracy(plasmid_host_normalized[idx_test, :], t)
    mah_acc.append(acc)
    print("Mah: ", acc)
    acc = taxonomic_accuracy(svpos[idx_test, :], t, use_max=True)
    svpos_acc.append(acc)
    print("Svpos: ", acc)
    acc = taxonomic_accuracy(blast_results_mat[idx_test, :], t, use_max=True)
    blast_acc.append(acc)
    print("Blast: ", acc)
net_blast_acc = np.vstack(net_blast_acc)
net_no_blast_acc = np.vstack(net_no_blast_acc)
mah_acc = np.vstack(mah_acc)
svpos_acc = np.vstack(mah_acc)
blast_acc = np.vstack(blast_acc)

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
