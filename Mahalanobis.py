import csv
import os
import numpy as np
import helper
import pickle
import subprocess
import argparse
from Bio import SeqIO
from matplotlib import pyplot as plt
from scipy.stats import pearsonr, spearmanr
from multiprocessing import Pool
from sklearn.preprocessing import normalize
from itertools import product
import glob
import argparse

np.seterr(all='raise')

parser = argparse.ArgumentParser(description='Mahalanobis Distance')
parser.add_argument('-k', type=int, default=3, required=True, help='k-mer length')
parser.add_argument('--host_dir', type=str, required=True, help='Directory containing host genome')
parser.add_argument('--plasmid_dir', type=str, required=True, help='Directory containing plasmid genome')
parser.add_argument('--temp_dir', type=str, default='temp_dir', help='Temp directory to store intermediate calculation result')
parser.add_argument('--output', type=str, required=True)

args = parser.parse_args()


k = 3



split_dir = './hosts_split'
split_dir_list = os.listdir(split_dir)
hash_dir = './hosts_hash'
matrix_dir = './hosts_matrix'


plasmid_hash_dir = 'plasmids_hash'
plasmid_matrix_dir = 'plasmids_matrix'

host_relative_abundance_dir = 'hosts_relative_abundance'
host_relative_abundance_list = os.listdir(host_relative_abundance_dir)
# host_relative_abundance_list = [l for l in host_relative_abundance_list if '_k{}'.format(k) in l]
# host_relative_abundance_list = [l for l in host_relative_abundance_list if 'markov' in l]
host_relative_abundance_list = [l for l in host_relative_abundance_list if '_z_score_k{}.npy'.format(k) in l]
plasmid_relative_abundance_dir = 'plasmids_relative_abundance'
# plasmid_relative_abundance_list = [s + '_relative_abundance_k{}.npy'.format(k) for s in plasmid_dir_list]
# plasmid_relative_abundance_list = [s + '_relative_markov_abundance_k{}.npy'.format(k) for s in plasmid_dir_list]
plasmid_relative_abundance_list = [s + '_z_score_k{}.npy'.format(k) for s in plasmid_dir_list]

NT_TYPE_NUM = 4

nucleotides = ['A', 'C', 'G', 'T']
mer_to_idx = {}
for i, mer in enumerate(product(nucleotides, repeat=k)):
    mer_to_idx[''.join(list(mer))] = i


def count_single_nucleotide(name, output_folder='temp_dir'):
    '''
    Count the number of occurrence of single nucleotides in a genome
    :param name: The name of the organism, can be a fasta file or folder containing split genome
    :param output_folder: The path to output the counting result
    :return:
    '''
    if os.path.isdir(name):
        count_single_nt_folder(name, output_folder)
    elif os.path.isfile(name):
        count_single_nt_file(name, output_folder)


def count_single_nt_folder(name, output_folder='temp_dir'):
    matrix_filename = os.path.join(output_folder, name + '_single_nt.npy')
    if os.path.exists(matrix_filename):
        return
    full_dir_path = os.path.join(split_dir, name)
    fa_list = os.listdir(full_dir_path)
    fa_list = [os.path.join(full_dir_path, fa) for fa in fa_list]
    occ_matrix = np.zeros((len(fa_list), NT_TYPE_NUM))
    for i, fa in enumerate(fa_list):
        occ = helper.count_nucleotide(fa)
        occ_matrix[i, :] = occ
    np.save(matrix_filename, occ_matrix)


def count_single_nt_file(plasmid_name, output_folder='temp_dir'):
    nt_count_filename = os.path.join(output_folder, plasmid_name + '_single_nt.npy')
    fa = os.path.join(plasmid_dir, plasmid_name)
    if os.path.exists(nt_count_filename):
        return
    occ = helper.count_nucleotide(fa)
    np.save(nt_count_filename, occ)


def calculate_frequency_product_single(single_matrix: np.ndarray, idx, k):
    if len(single_matrix.shape) == 2:
        single_frequency_product = np.ones(single_matrix.shape[0])
    else:
        single_frequency_product = 1
    for i in range(k):
        if len(single_matrix.shape) == 2:
            single_frequency_product = single_frequency_product * single_matrix[:, int(idx % k)]
        elif len(single_matrix.shape) == 1:
            single_frequency_product = single_frequency_product * single_matrix[int(idx % k)]
        idx = idx / k
    return single_frequency_product


def calculate_frequency_product_all(multiple_matrix, single_matrix, k):
    '''
    For all k-mer, calculate the product of the frequency of each single nucleotide
    :param multiple_matrix: frequency matrix for k-mer
    :param single_matrix: frequency matrix for single nucleotide
    :param k: k-mer length
    :return:
    '''
    frequency_product = np.zeros(multiple_matrix.shape)
    for i in range(multiple_matrix.shape[-1]):
        if len(multiple_matrix.shape) == 2:
            frequency_product[:, i] = calculate_frequency_product_single(single_matrix, i, k)
        else:
            frequency_product[i] = calculate_frequency_product_single(single_matrix, i, k)
    return frequency_product


def calculate_relative_abundance(multiple_matrix, single_matrix, k):
    '''
    Calculate the relative abundance defined in the Mahalanobis paper
    :param multiple_matrix: frequency matrix for k-mer
    :param single_matrix: frequency matrix for single nucleotide
    :param k: k-mer length
    :return:
    '''
    frequency_product = calculate_frequency_product_all(multiple_matrix, single_matrix, k)
    multiple_matrix[multiple_matrix == 0] = 1
    frequency_product[frequency_product == 0] = 1
    relative_abundance = multiple_matrix / frequency_product
    return relative_abundance


def relative_abundance_host(host_name):
    relative_abundance_path = os.path.join(host_relative_abundance_dir,
                                           host_name + '_relative_abundance' + '_k{}'.format(k))
    if os.path.exists(relative_abundance_path):
        return
    matrix_filename = os.path.join(matrix_dir, host_name + '_k{}'.format(k) + '.npy')
    multiple_matrix = np.load(matrix_filename)
    k_mer_freq = multiple_matrix / multiple_matrix.sum(axis=1)[:, None]
    matrix_filename = os.path.join(matrix_dir, host_name + '_single_nt' + '.npy')
    single_matrix = np.load(matrix_filename)
    single_nt_freq = single_matrix / single_matrix.sum(axis=1)[:, None]
    relative_abundance = calculate_relative_abundance(k_mer_freq, single_nt_freq, k)
    if np.isnan(relative_abundance.sum()):
        print(host_name)
    np.save(relative_abundance_path, relative_abundance)


def fetch_index(multiple_matrix, idx, k):
    pass

def calculate_markov_expectation_single(multiple_matrix, idx, k=3):
    last_two_idx = []
    first_two_idx = []
    mid_idx = []
    for i in range(4):
        last_two_idx.append(int(idx % (NT_TYPE_NUM ** (k-1)) + i * (NT_TYPE_NUM ** (k-1))))
        first_two_idx.append(int(idx - idx % 4 + i))
        for j in range(4):
            mid_idx.append(int((idx - idx % 4 + i) % (NT_TYPE_NUM ** (k-1)) + j * (NT_TYPE_NUM ** (k-1))))
    mid_idx.sort()
    if len(multiple_matrix.shape) == 2:
        last_two_sum = multiple_matrix[:, last_two_idx].sum(axis=1)
        first_two_sum = multiple_matrix[:, first_two_idx].sum(axis=1)
        mid_sum = multiple_matrix[:, mid_idx].sum(axis=1)

        last_two_sum[last_two_sum == 0] = 1
        first_two_sum[first_two_sum == 0] = 1
        mid_sum[mid_sum == 0] = 1
    else:
        last_two_sum = multiple_matrix[last_two_idx].sum()
        first_two_sum = multiple_matrix[first_two_idx].sum()
        mid_sum = multiple_matrix[mid_idx].sum()

        last_two_sum = 1 if last_two_sum == 0 else last_two_sum
        first_two_sum = 1 if first_two_sum == 0 else first_two_sum
        mid_sum = 1 if mid_sum == 0 else mid_sum

    markov_expectation = last_two_sum * first_two_sum / mid_sum

    return markov_expectation


def calculate_markov_expectation(multiple_matrix):
    markov_ratio = np.zeros(multiple_matrix.shape)
    for i in range(multiple_matrix.shape[-1]):
        if len(multiple_matrix.shape) == 2:
            markov_ratio[:, i] = calculate_markov_expectation_single(multiple_matrix, i)
        else:
            markov_ratio[i] = calculate_markov_expectation_single(multiple_matrix, i)
    return markov_ratio


def calculate_relative_markov_abundance(multiple_matrix):
    '''
    For now the Markov abundance is only defined for 3-mer
    :param multiple_matrix: Count of k-mer
    :param single_matrix: Count of single nucleotide
    :return:
    '''
    markov_ratio = calculate_markov_expectation(multiple_matrix)
    markov_ratio[markov_ratio == 0] = 1
    multiple_matrix[multiple_matrix == 0] = 1
    relative_markov_abundance = multiple_matrix / markov_ratio
    return relative_markov_abundance
    # TODO: allow k-mer with larger length


def relative_markov_abundance_host(host_name):
    relative_abundance_path = os.path.join(host_relative_abundance_dir,
                                           host_name + '_relative_markov_abundance' + '_k{}'.format(k))
    if os.path.exists(relative_abundance_path):
        return
    matrix_filename = os.path.join(matrix_dir, host_name + '_k{}'.format(k) + '.npy')
    multiple_matrix = np.load(matrix_filename)
    relative_markov_abundance = calculate_relative_markov_abundance(multiple_matrix)
    if np.isnan(relative_markov_abundance.sum()):
        print(host_name)
    np.save(relative_abundance_path, relative_markov_abundance)

def validate_matrix(multiple_matrix, idx):
    '''
    Set the 0 in a matrix to 1 to avoid divide by 0 error
    :param multiple_matrix:
    :param idx:
    :return:
    '''
    pass


def calculate_markov_variance_single(multiple_matrix, expectation, idx, k=3):
    last_idx = []
    first_idx = []
    mid_idx = []
    for i in range(4):
        last_idx.append(int(idx % (NT_TYPE_NUM ** (k - 1)) + i * (NT_TYPE_NUM ** (k - 1))))
        first_idx.append(int(idx - idx % 4 + i))
        for j in range(4):
            mid_idx.append(int((idx - idx % 4 + i) % (NT_TYPE_NUM ** (k - 1)) + j * (NT_TYPE_NUM ** (k - 1))))
    mid_idx.sort()
    if len(multiple_matrix.shape) == 2:
        last_two_sum = multiple_matrix[:, last_idx].sum(axis=1)
        first_two_sum = multiple_matrix[:, first_idx].sum(axis=1)
        mid_sum = multiple_matrix[:, mid_idx].sum(axis=1)

        last_two_sum[last_two_sum == 0] = 1
        first_two_sum[first_two_sum == 0] = 1
        mid_sum[mid_sum == 0] = 1
    else:
        last_two_sum = multiple_matrix[last_idx].sum()
        first_two_sum = multiple_matrix[first_idx].sum()
        mid_sum = multiple_matrix[mid_idx].sum()

        last_two_sum = 1 if last_two_sum == 0 else last_two_sum
        first_two_sum = 1 if first_two_sum == 0 else first_two_sum
        mid_sum = 1 if mid_sum == 0 else mid_sum

    # expectation = calculate_markov_expectation_single(multiple_matrix, idx, k)
    variance = expectation * ((mid_sum - first_two_sum) * (mid_sum - last_two_sum) / (mid_sum ** 2))

    return variance




class TeelingZScore:
    def __init__(self, hosts_kmer_count, plasmids_kmer_count, k=3):
        self.k = k
        self.hosts_kmer_count = hosts_kmer_count
        self.plasmids_kmer_count = plasmids_kmer_count

    def calculate_markov_expectation(self, multiple_matrix):
        markov_ratio = np.zeros(multiple_matrix.shape)
        for i in range(multiple_matrix.shape[-1]):
            if len(multiple_matrix.shape) == 2:
                markov_ratio[:, i] = calculate_markov_expectation_single(multiple_matrix, i)
            else:
                markov_ratio[i] = calculate_markov_expectation_single(multiple_matrix, i)
        return markov_ratio

    def calculate_markov_variance(self, multiple_matrix, expectation, idx, k=3):
        last_idx = []
        first_idx = []
        mid_idx = []
        for i in range(4):
            last_idx.append(int(idx % (NT_TYPE_NUM ** (k - 1)) + i * (NT_TYPE_NUM ** (k - 1))))
            first_idx.append(int(idx - idx % 4 + i))
            for j in range(4):
                mid_idx.append(int((idx - idx % 4 + i) % (NT_TYPE_NUM ** (k - 1)) + j * (NT_TYPE_NUM ** (k - 1))))
        mid_idx.sort()
        if len(multiple_matrix.shape) == 2:
            last_two_sum = multiple_matrix[:, last_idx].sum(axis=1)
            first_two_sum = multiple_matrix[:, first_idx].sum(axis=1)
            mid_sum = multiple_matrix[:, mid_idx].sum(axis=1)

            last_two_sum[last_two_sum == 0] = 1
            first_two_sum[first_two_sum == 0] = 1
            mid_sum[mid_sum == 0] = 1
        else:
            last_two_sum = multiple_matrix[last_idx].sum()
            first_two_sum = multiple_matrix[first_idx].sum()
            mid_sum = multiple_matrix[mid_idx].sum()

            last_two_sum = 1 if last_two_sum == 0 else last_two_sum
            first_two_sum = 1 if first_two_sum == 0 else first_two_sum
            mid_sum = 1 if mid_sum == 0 else mid_sum

        variance = expectation * ((mid_sum - first_two_sum) * (mid_sum - last_two_sum) / (mid_sum ** 2))

        return variance

    def teeling_z_score(self, multiple_matrix, idx):
        expectation = self.calculate_markov_variance(multiple_matrix, idx, k)
        variance = calculate_markov_variance_single(multiple_matrix, expectation, idx, k)
        if len(multiple_matrix.shape) == 2:
            variance[variance == 0] = 0.1
            z_score = (multiple_matrix[:, idx] - expectation) / np.sqrt(variance)
        else:
            variance = 0.1 if variance == 0 else variance
            z_score = (multiple_matrix[idx] - expectation) / np.sqrt(variance)
        return z_score

    def teeling_z_score_all(self, multiple_matrix):
        z_score = np.zeros(multiple_matrix.shape)
        for i in range(multiple_matrix.shape[-1]):
            if len(multiple_matrix.shape) == 2:
                z_score[:, i] = self.teeling_z_score(multiple_matrix, i, k)
            else:
                z_score[i] = self.teeling_z_score(multiple_matrix, i, k)
        return z_score

    def teeling_z_score_host(self, host_name):
        z_score_path = os.path.join(host_relative_abundance_dir,
                                    host_name + '_z_score' + '_k{}'.format(k))
        if os.path.exists(z_score_path):
            return
        matrix_filename = os.path.join(matrix_dir, host_name + '_k{}'.format(k) + '.npy')
        multiple_matrix = np.load(matrix_filename)
        z_score = self.teeling_z_score_all(multiple_matrix)
        if np.isnan(z_score.sum()):
            print(host_name)
        np.save(z_score_path, z_score)

    def teeling_z_score_plasmid(self, plasmid_name):
        z_score_path = os.path.join(plasmid_relative_abundance_dir, plasmid_name + '_z_score' + '_k{}'.format(k))
        if os.path.exists(z_score_path):
            return
        matrix_filename = os.path.join(plasmid_matrix_dir, plasmid_name + '_k{}'.format(k) + '.npy')
        multiple_matrix = np.load(matrix_filename)
        z_score = self.teeling_z_score_all(multiple_matrix)
        if np.isnan(z_score.sum()):
            print(plasmid_name)
        np.save(z_score_path, z_score)

    def calculate_score(self, p=None):
        if not p:
            pass
        else:
            pass








def count_plasmid_kmer(plasmid_name):
    fa_path = os.path.join(plasmid_dir, plasmid_name)
    k_mer_count = np.zeros(NT_TYPE_NUM ** k)
    hash_path = os.path.join(plasmid_hash_dir, plasmid_name) + '_k{}'.format(k)
    if not os.path.exists('{}_0'.format(hash_path)):
        command = './jellyfish-linux count -m {} -s 200M -t 8 -o {} {}'.format(k, hash_path, fa_path)
        os.system(command)
    count = os.popen('./jellyfish-linux dump -c {}_0'.format(hash_path)).read()
    counts = count.split('\n')[:-1]
    for m_c in counts:
        (m, c) = m_c.split(' ')
        k_mer_count[mer_to_idx[m]] = c
    k_mer_count_filename = os.path.join(plasmid_matrix_dir, plasmid_name) + '_k{}'.format(k)
    np.save(k_mer_count_filename, k_mer_count)


def relative_markov_abundance_plasmid(plasmid_name):
    relative_abundance_path = os.path.join(plasmid_relative_abundance_dir,
                                           plasmid_name + '_relative_markov_abundance' + '_k{}'.format(k))
    matrix_filename = os.path.join(plasmid_matrix_dir, plasmid_name + '_k{}'.format(k) + '.npy')
    multiple_matrix = np.load(matrix_filename)
    relative_markov_abundance = calculate_relative_markov_abundance(multiple_matrix)
    np.save(relative_abundance_path, relative_markov_abundance)


def relative_abundance_plasmid(plasmid_name):
    relative_abundance_path = os.path.join(plasmid_relative_abundance_dir,
                                           plasmid_name + '_relative_abundance' + '_k{}'.format(k))
    # if os.path.exists(relative_abundance_path):
    #     return
    matrix_filename = os.path.join(plasmid_matrix_dir, plasmid_name + '_k{}'.format(k) + '.npy')
    multiple_matrix = np.load(matrix_filename)
    k_mer_freq = multiple_matrix / multiple_matrix.sum()
    matrix_filename = os.path.join(plasmid_matrix_dir, plasmid_name + '_single_nt' + '.npy')
    single_matrix = np.load(matrix_filename)
    single_nt_freq = single_matrix / single_matrix.sum()
    relative_abundance = calculate_relative_abundance(k_mer_freq, single_nt_freq, k)
    # return relative_abundance
    np.save(relative_abundance_path, relative_abundance)


def mahalanobis_distance(plasmid_relative_abundance, host_relative_abundance: np.ndarray):
    '''
    Calculate the Mahalanobis distance for a plasmid and a host
    :param plasmid_relative_abundance: (n,) numpy array
    :param host_relative_abundance: (m, n) numpy matrix
    :return:
    '''
    plasmid_relative_abundance = plasmid_relative_abundance.flatten()
    host_mean = host_relative_abundance.mean(axis=0)
    diff = plasmid_relative_abundance - host_mean
    cov = np.cov(host_relative_abundance.T)
    cov_inv = np.linalg.inv(cov)
    distance = diff[None, :].dot(cov_inv).dot(diff)
    return distance


def mahalanobis_distance_plasmid(plasmid_relative_abundance_filename):
    plasmid_relative_abundance_path = os.path.join(plasmid_relative_abundance_dir, plasmid_relative_abundance_filename)
    plasmid_relative_abundance = np.load(plasmid_relative_abundance_path)
    distances = np.zeros(len(host_relative_abundance_list))
    for i, host in enumerate(host_relative_abundance_list):
        host_relative_abundance_path = os.path.join(host_relative_abundance_dir, host)
        host_relative_abundance = np.load(host_relative_abundance_path)
        distances[i] = mahalanobis_distance(plasmid_relative_abundance, host_relative_abundance)
    print(plasmid_relative_abundance_filename)
    return distances


if __name__ == '__main__':

    plasmid_dir = 'plasmids_used'
    f = open('plasmids_used.txt')
    plasmid_dir_list = [os.path.split(l.strip())[-1] for l in f.readlines()]


    with Pool(6) as p:
        p.map(count_single_nt_folder, split_dir_list)

    a = np.load('hosts_matrix/GCF_000006605.1_ASM660v1_genomic.npy')
    b = np.load('hosts_matrix/GCF_000006605.1_ASM660v1_genomic_single_nt.npy')
    c = a / a.sum(axis=1).reshape((-1, 1))
    d = b / b.sum(axis=1).reshape((-1, 1))
    e = calculate_relative_abundance(c, d, 3)
    print(np.linalg.inv(np.cov(e)))

    with Pool(6) as p:
        p.map(relative_abundance_host, split_dir_list)

    with Pool(6) as p:
        p.map(count_plasmid_kmer, plasmid_dir_list)

    with Pool(6) as p:
        p.map(count_single_nt_file, plasmid_dir_list)

    with Pool(6) as p:
        p.map(relative_abundance_plasmid, plasmid_dir_list)

    for plasmid in plasmid_dir_list:
        r_a = relative_abundance_plasmid(plasmid)

    for plasmid in plasmid_dir_list:
        relative_abundance_plasmid(plasmid)

    distances = []
    for plasmid in plasmid_relative_abundance_list:
        distances.append(mahalanobis_distance_plasmid(plasmid))
    np.save('mahalanobis_distances_k{}'.format(k), np.vstack(distances))


    with Pool(6) as p:
        p.map(relative_markov_abundance_host, split_dir_list)
    for l in split_dir_list:
        relative_markov_abundance_host(l)

    for l in split_dir_list:
        teeling_z_score_host(l)

    for l in plasmid_dir_list:
        teeling_z_score_plasmid(l)
    distances = []
    for plasmid in plasmid_relative_abundance_list:
            distances.append(mahalanobis_distance_plasmid(plasmid))
    np.save('teeling_z_score_mahalanobis_distances_k{}'.format(k), np.vstack(distances))

    with Pool(6) as p:
        p.map(plasmid_k_mer_count, plasmid_dir_list)

    with Pool(6) as p:
        p.map(count_single_nt_file, plasmid_dir_list)

    with Pool(6) as p:
        p.map(relative_markov_abundance_plasmid, plasmid_dir_list)

    for plasmid in plasmid_dir_list:
        r_a = relative_abundance_plasmid(plasmid)

    for plasmid in plasmid_dir_list:
        relative_abundance_plasmid(plasmid)
    distances = []
    for plasmid in plasmid_relative_abundance_list:
        distances.append(mahalanobis_distance_plasmid(plasmid))
    np.save('markov_mahalanobis_distances_k{}'.format(k), np.vstack(distances))

    with Pool(4) as p:
        distances = p.map(mahalanobis_distance_plasmid, plasmid_relative_abundance_list)
    np.save('mahalanobis_distances', np.vstack(distances))
