import os
import numpy as np
import csv
import pickle
import statsmodels.api as sm
from typing import List


plasmids_list = []
# hosts_list = os.listdir('hosts_processed')
# hosts_list.sort()
# host_list_dict = dict((name, i) for i, name in enumerate(hosts_list))

cafe_path = './cafe_linux'
plasmids_path = 'plasmids'
hosts_path = 'hosts_processed'
jellyfish_path = './jellyfish-linux'
sh_path = 'distance.sh'
dataset_metadata_path = 'Metadata.csv'




def prepare_name():
    plasmids = []
    hosts = []
    assemblyname = []
    with open(dataset_metadata_path, 'r', encoding='mac_roman') as metadata:
        for i, row in enumerate(csv.reader(metadata)):
            plasmids.append(row[1])
            hosts.append(row[5])
            assemblyname.append('_'.join(row[9].split(' ')))
        plasmids = plasmids[1:]
        hosts = hosts[1:]
        assemblyname = assemblyname[1:]

    return plasmids, hosts, assemblyname


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

# host_to_idx = load_obj('host_to_idx.pkl')
# idx_to_host = load_obj('idx_to_host.pkl')


def check_plasmid_fn(plasmid_name):
    if plasmid_name[-2] == '.':
        if os.path.isfile(plasmids_path + '/{}.fa'.format(plasmid_name)):
            plasmid_path = plasmids_path + '/{}.fa'.format(plasmid_name)
        else:
            print(plasmid_name)
            return False
    else:
        if os.path.isfile(plasmids_path + '/{}.1.fa'.format(plasmid_name)):
            plasmid_path = plasmids_path + '/{}.1.fa'.format(plasmid_name)
        elif os.path.isfile(plasmids_path + '/{}.2.fa'.format(plasmid_name)):
            plasmid_path = plasmids_path + '/{}.2.fa'.format(plasmid_name)
        else:
            print(plasmid_name)
            return False
    return plasmid_path


def host_to_species_taxid(host_name):
    pass


def lowess_adjustment(y: List[np.ndarray], x, target):
    lowess = sm.nonparametric.lowess

    z = lowess(y[0], x, return_sorted=False)
    y_new = []
    for y_ in y:
        y_new.append(y_ + target - z)

    return y


def count_nucleotide(file):
    try:
        f = open(file)
    except UnicodeDecodeError:
        print(file)
        return np.ones(4)
    occ = np.zeros(4)
    try:
        for line in f.readlines():
            if line[0] == '>':
                continue
            for c in line:
                if c == 'A':
                    occ[0] += 1
                elif c == 'C':
                    occ[1] += 1
                elif c == 'G':
                    occ[2] += 1
                elif c == 'T':
                    occ[3] += 1
        return occ
    except UnicodeDecodeError:
        print(file)
        return np.ones(4)
