from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np
import pandas as pd
import os
import shutil
import util
from Bio.Blast.Applications import NcbiblastnCommandline

metadata = pd.read_csv('metadata.csv')
metadata = metadata.dropna()
metadata = metadata.reset_index(drop=True)

for i in range(len(metadata.Locus_ID)):
    metadata.at[i, 'Locus_ID'] = metadata.at[i, 'Locus_ID'].split('.')[0]

metadata.sort_values(by=['Locus_ID'], inplace=True)
metadata.reset_index(drop=True, inplace=True)

# os.makedirs('data/plasmids_used')
# plasmid_l = list(metadata.Locus_ID)
# for fn in os.listdir('data/plasmids'):
#     if fn.split('.')[0] in plasmid_l:
#         shutil.copyfile(os.path.join('data/plasmids', fn), os.path.join('data/plasmids_used', fn))

# t = MahalanobisRelativeAbundance('data/plasmids_used', 'data/plasmids_used', temp_directory_path='temp_dir/plasmid_plasmid')
# dist = t.calc_distance(8)
# np.save('plasmid_plasmid.npy', dist)

# util.remove_plasmid_seq('data/hosts', 'data/hosts_no_plasmid')

"""
Calculate plasmid-host distance
"""
# t = MahalanobisRelativeAbundance('data/hosts', 'data/plasmids_used', temp_directory_path='temp_dir/plasmid_host')
# dist = t.calc_distance(8)
# np.save('plasmid_host.npy', dist)

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


speciesid_to_idx = list(set(metadata.Assembly_speciestaxid))
speciesid_to_idx.sort()
speciesid_to_idx = {s: i for i, s in enumerate(speciesid_to_idx)}

plasmid_host = np.load('plasmid_host.npy')
true_idx = [host_to_idx_dict[metadata.Assembly_chainid[i]] for i in range(len(metadata))]
false_idx = []
for i in range(len(metadata)):
    idx = np.random.randint(plasmid_host.shape[1])
    while idx == true_idx[i]:
        idx = np.random.randint(plasmid_host.shape[1])
    false_idx.append(idx)

plasmid_plasmid = np.load('plasmid_plasmid.npy')
svpos = plasmid_plasmid.dot(interaction_table) / interaction_table.sum(axis=0)
svneg = plasmid_plasmid.dot((1 - interaction_table)) / (1 - interaction_table).sum(axis=0)

plasmid_host = (plasmid_host - plasmid_host.min(axis=0)) / (plasmid_host.max(axis=0) - plasmid_host.min(axis=0))

X_pos = np.concatenate([plasmid_host[np.arange(plasmid_host.shape[0]), true_idx, np.newaxis], svpos[np.arange(svpos.shape[0]), true_idx, np.newaxis], svneg[np.arange(svpos.shape[0]), true_idx, np.newaxis], blast_results_mat[np.arange(blast_results_mat.shape[0]), true_idx, np.newaxis]], axis=1)
X_neg = np.concatenate([plasmid_host[np.arange(plasmid_host.shape[0]), false_idx, np.newaxis], svpos[np.arange(svpos.shape[0]), false_idx, np.newaxis], svneg[np.arange(svpos.shape[0]), false_idx, np.newaxis], blast_results_mat[np.arange(blast_results_mat.shape[0]), false_idx, np.newaxis]], axis=1)
X = np.concatenate((X_pos, X_neg), axis=0)

y = np.concatenate((np.ones(X_pos.shape[0]), np.zeros(X_neg.shape[0])))

X = np.concatenate((X, np.ones((X.shape[0], 1))), axis=1)

import statsmodels.api as sm
glm_binom = sm.GLM(y, X, family=sm.families.Binomial())
res = glm_binom.fit()

from sklearn.model_selection import train_test_split

data = np.concatenate((X, y[:, None]), axis=1)
data_train, data_test = train_test_split(data, test_size=0.1)
glm_binom = sm.GLM(data_train[:, -1], data_train[:, :-1], family=sm.families.Binomial())
res = glm_binom.fit(max_iter=200)
a = res.predict(data_test[:, :-1])
a[a>0.5]=1
a[a<0.5]=0
np.sum(a==data_test[:, -1]) / a.shape

data = np.concatenate((plasmid_host.flatten()[:, None], svpos.flatten()[:, None], svneg.flatten()[:, None], blast_results_mat.flatten()[:, None]), axis=1)
data = np.concatenate((data, np.ones((data.shape[0], 1))), axis=1)

pred = res.predict(data)
pred = pred.reshape(plasmid_host.shape)


from ete3 import NCBITaxa

desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
def taxid_to_lineage(taxid_1, taxid_2):
    ncbi = NCBITaxa()
    lineage_1 = ncbi.get_lineage(taxid_1)
    lineage_2 = ncbi.get_lineage(taxid_2)
    rank_to_id_1 = {rank: id for (id, rank) in ncbi.get_rank(lineage_1).items()}
    rank_to_id_2 = {rank: id for (id, rank) in ncbi.get_rank(lineage_2).items()}
    for desired_rank in desired_ranks:
        if desired_rank in rank_to_id_1.keys():
            pass
    rank_to_id_1 = {desired_rank: (rank_to_id_1[desired_rank] if desired_rank in rank_to_id_1.keys() else None) for desired_rank in desired_ranks}
    rank_to_id_2 = {desired_rank: (rank_to_id_2[desired_rank] if desired_rank in rank_to_id_2.keys() else None) for desired_rank in desired_ranks}
    return rank_to_id_1, rank_to_id_2

cnt = np.zeros(7)
for i in range(pred.shape[0]):
    rank_to_id_1, rank_to_id_2 = taxid_to_lineage(host_to_speciesid[idx_to_host_dict[pred[i, :].argmin()]],
                                                  host_to_speciesid[idx_to_host_dict[true_idx[i]]])
    for j, rank in enumerate(desired_ranks):
        if rank_to_id_1[rank] == rank_to_id_2[rank] and rank_to_id_1[rank] is not None:
            cnt[j] += 1
            # if j == 1:
            #     print(idx_to_host_dict[j])
        # if host_to_taxid[idx_to_host[score[i, :].argmin()]] == host_to_taxid[idx_to_host[plasmids_class[i]]]:
        #     cnt += 1
print(cnt / pred.shape[0])