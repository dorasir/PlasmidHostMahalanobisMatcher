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





