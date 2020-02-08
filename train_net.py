from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np
import pandas as pd
import os
import shutil
import util
from Bio.Blast.Applications import NcbiblastnCommandline


from ete3 import NCBITaxa
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
    rank_to_id = {desired_rank: (rank_to_id[desired_rank] if desired_rank in rank_to_id.keys() else None) for desired_rank in desired_ranks}
    return rank_to_id


def taxonomic_accuracy(prediction, target, use_max=False):
    cnt = np.zeros(7)
    for i in range(prediction.shape[0]):
        if max:
            rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[prediction[i, :].argmax()]])
        else:
            rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[prediction[i, :].argmin()]])
        rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[target[i]]])
        for j, rank in enumerate(desired_ranks):
            if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
                cnt[j] += 1
    if max:
        cnt[6] = np.sum(prediction.argmax(axis=1) == target)
    acc = cnt / pred.shape[0]
    return acc

metadata = pd.read_csv('metadata.csv')
metadata = metadata.dropna()
metadata = metadata.reset_index(drop=True)

for i in range(len(metadata.Locus_ID)):
    metadata.at[i, 'Locus_ID'] = metadata.at[i, 'Locus_ID'].split('.')[0]

metadata.sort_values(by=['Locus_ID'], inplace=True)
metadata.reset_index(drop=True, inplace=True)

"""Move plasmids with host information to plasmids_used directory"""
# os.makedirs('data/plasmids_used')
# plasmid_l = list(metadata.Locus_ID)
# for fn in os.listdir('data/plasmids'):
#     if fn.split('.')[0] in plasmid_l:
#         shutil.copyfile(os.path.join('data/plasmids', fn), os.path.join('data/plasmids_used', fn))

# util.remove_plasmid_seq('data/hosts', 'data/hosts_no_plasmid')

"""
Calculate plasmid-host distance
"""
# t = MahalanobisRelativeAbundance('data/hosts', 'data/plasmids_used', temp_directory_path='temp_dir/plasmid_host')
# dist = t.calc_distance(8)
# np.save('plasmid_host.npy', dist)

plasmid_host = np.load('plasmid_host.npy')

"""
Calculate plasmid-wise distance
"""
dist = np.load('plasmid_host.npy')
plasmid_plasmid_distance = util.cosine_similarity(dist, dist)
np.save('plasmid_plasmid.npy', plasmid_plasmid_distance)

"""
Construct plasmid interaction table
"""

host_list = list(set(metadata.Assembly_chainid))
host_list.sort()
host_to_idx_dict = {host: i for i, host in enumerate(host_list)}
idx_to_host_dict = {i: host for i, host in enumerate(host_list)}

interaction_table = np.zeros((len(metadata), len(set(host_list)))) #plasmid-strain indicator
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
        speciesid_to_host[metadata.Assembly_speciestaxid[i]] = set([metadata.Assembly_chainid[i]])
    else:
        speciesid_to_host[metadata.Assembly_speciestaxid[i]].add(metadata.Assembly_chainid[i])

interaction_table_species = np.zeros((len(metadata), len(set(host_list)))) #plasmid-strain indicator
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

# blast_results = util.blast_dir_to_db('data/plasmids_used', 'data/hosts_no_plasmid_complete.fna')
# util.save_obj(blast_results, 'blast_results.pkl')

"""
Construct blast result matrix
"""
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
plasmid_host = np.load('plasmid_host.npy')
taxonomic_level = 'phylum'
true_idx = [host_to_idx_dict[metadata.Assembly_chainid[i]] for i in range(len(metadata))]
false_idx = []
true_family = [(taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])[taxonomic_level] if taxonomic_level in taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]]) else None) for idx in true_idx]
for i in range(len(metadata)):
    idx = np.random.randint(plasmid_host.shape[1])
    false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
    false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
    while idx == true_idx[i] and false_family == true_family[i]:
        idx = np.random.randint(plasmid_host.shape[1])
        false_lineage = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[idx]])
        false_family = false_lineage[taxonomic_level] if taxonomic_level in false_lineage else None
    false_idx.append(idx)

"""calculate svpos and svneg"""
plasmid_plasmid = np.load('plasmid_plasmid.npy')
svpos = plasmid_plasmid.dot(interaction_table) / interaction_table.sum(axis=0)
svneg = plasmid_plasmid.dot((1 - interaction_table)) / (1 - interaction_table).sum(axis=0)

"""Normalize plasmid-host distance"""
plasmid_host = (plasmid_host - plasmid_host.min(axis=0)) / (plasmid_host.max(axis=0) - plasmid_host.min(axis=0))

"""Construct positive and negative dataset"""
X_pos = np.concatenate([plasmid_host[np.arange(plasmid_host.shape[0]), true_idx, np.newaxis], svpos[np.arange(svpos.shape[0]), true_idx, np.newaxis], svneg[np.arange(svpos.shape[0]), true_idx, np.newaxis], blast_results_mat[np.arange(blast_results_mat.shape[0]), true_idx, np.newaxis]], axis=1)
X_neg = np.concatenate([plasmid_host[np.arange(plasmid_host.shape[0]), false_idx, np.newaxis], svpos[np.arange(svpos.shape[0]), false_idx, np.newaxis], svneg[np.arange(svpos.shape[0]), false_idx, np.newaxis], blast_results_mat[np.arange(blast_results_mat.shape[0]), false_idx, np.newaxis]], axis=1)
X = np.concatenate((X_pos, X_neg), axis=0)

y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

"""Train generalized linear model to predict likelihood of interaction"""
import statsmodels.api as sm
glm_binom = sm.GLM(y, X, family=sm.families.Binomial())
res = glm_binom.fit()

"""Split """
from sklearn.model_selection import train_test_split

# for i in range(10):
#     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)
#     glm_binom = sm.GLM(y_train, X_train, family=sm.families.Binomial())
#     res = glm_binom.fit()


data = np.concatenate((X, y[:, None]), axis=1)
data_train, data_test = train_test_split(data, test_size=0.1)
glm_binom = sm.GLM(data_train[:, -1], data_train[:, :-1], family=sm.families.Binomial())
res = glm_binom.fit(max_iter=200)

data = np.concatenate((plasmid_host.flatten()[:, None], svpos.flatten()[:, None], svneg.flatten()[:, None], blast_results_mat.flatten()[:, None]), axis=1)
data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

pred = res.predict(data)
pred = pred.reshape(plasmid_host.shape)





acc = taxonomic_accuracy(pred, true_idx)
print(acc)

# cnt = np.zeros(7)
# for i in range(pred.shape[0]):
#     rank_to_id_1 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[pred[i, :].argmin()]])
#     rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[true_idx[i]]])
#     for j, rank in enumerate(desired_ranks):
#         if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
#             cnt[j] += 1
# print(cnt / pred.shape[0])

for i in range(10):
    idx_train, idx_test = train_test_split(np.arange(plasmid_host.shape[0]), test_size=0.2)

    true_idx = np.array(true_idx)
    false_idx = np.array(false_idx)
    X_pos = np.concatenate([plasmid_host[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
    X_neg = np.concatenate([plasmid_host[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis]], axis=1)
    X = np.concatenate((X_pos, X_neg), axis=0)

    y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

    X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

    glm_binom = sm.GLM(y, X, family=sm.families.Binomial())
    res = glm_binom.fit()

    data = np.concatenate((plasmid_host[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None], svneg[idx_test, :].flatten()[:, None], blast_results_mat[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

    pred = res.predict(data)
    pred = pred.reshape((-1, plasmid_host.shape[1]))

    t = true_idx[idx_test]

    acc = taxonomic_accuracy(pred, t)
    print(acc)

from sklearn.neural_network import MLPRegressor
for i in range(10):
    idx_train, idx_test = train_test_split(np.arange(plasmid_host.shape[0]), test_size=0.2)

    true_idx = np.array(true_idx)
    false_idx = np.array(false_idx)
    X_pos = np.concatenate([plasmid_host[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), true_idx[idx_train], np.newaxis]], axis=1)
    X_neg = np.concatenate([plasmid_host[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svpos[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], svneg[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis], blast_results_mat[np.arange(idx_train.shape[0]), false_idx[idx_train], np.newaxis]], axis=1)
    X = np.concatenate((X_pos, X_neg), axis=0)

    y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

    X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

    model = MLPRegressor(hidden_layer_sizes=(50, 10,), activation='logistic', max_iter=1000, solver='lbfgs')
    model.fit(X, y)

    data = np.concatenate((plasmid_host[idx_test, :].flatten()[:, None], svpos[idx_test, :].flatten()[:, None], svneg[idx_test, :].flatten()[:, None], blast_results_mat[idx_test, :].flatten()[:, None]), axis=1)
    data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

    pred = model.predict(data)
    pred = pred.reshape((-1, plasmid_host.shape[1]))

    t = true_idx[idx_test]

    acc = taxonomic_accuracy(pred, t)
    print(acc)

# cnt = np.zeros(7)
# for i in range(pred.shape[0]):
#     rank_to_id_1, rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[pred[i, :].argmax()]],
#                                                   host_to_speciesid[idx_to_host_dict[t[i]]])
#     for j, rank in enumerate(desired_ranks):
#         if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
#             cnt[j] += 1
# print(cnt / pred.shape[0])

# data = np.concatenate((plasmid_host[idx_train, :].flatten()[:, None], svpos[idx_train, :].flatten()[:, None], svneg[idx_train, :].flatten()[:, None], blast_results_mat[idx_train, :].flatten()[:, None]), axis=1)
# data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)
#
# pred = res.predict(data)
# pred = pred.reshape((-1, plasmid_host.shape[1]))
#
# t = true_idx[idx_train]
#
# cnt = np.zeros(7)
# for i in range(pred.shape[0]):
#     rank_to_id_1, rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[pred[i, :].argmax()]],
#                                                   host_to_speciesid[idx_to_host_dict[t[i]]])
#     for j, rank in enumerate(desired_ranks):
#         if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
#             cnt[j] += 1
# print(cnt / pred.shape[0])


import matplotlib.pyplot as plt
plt.boxplot([X_pos[:, 0], X_neg[:, 0]])
plt.show()

plt.boxplot([X_pos[:, 1], X_neg[:, 1]])
plt.show()

plt.boxplot([X_pos[:, 2], X_neg[:, 2]])
plt.show()

plt.boxplot([X_pos[:, 3], X_neg[:, 3]])
plt.show()


plt.hist([X_pos[:, 0], X_neg[:, 0]])
plt.show()

plt.hist([X_pos[:, 1], X_neg[:, 1]])
plt.show()

plt.hist([X_pos[:, 2], X_neg[:, 2]])
plt.show()

plt.hist([X_pos[:, 3], X_neg[:, 3]])
plt.show()
