import numpy as np
import os
import helper

NT_TYPE_NUM = 4

class MahalanobisRelativeAbundance:
    def __init__(self, host_path, plasmid_path):
        self.host_path = host_path
        self.plasmid_path = plasmid_path
        self.host_single_count = 0
        self.host_multiple_count = 0
        self.plasmid_single_count = 0
        self.plasmid_single_count = 0


    def count_single_nucleotide(self, name, output_folder='temp_dir'):
        '''
        Count the number of occurrence of single nucleotides in a genome

        :param str name: The name of the organism, can be a fasta file or folder containing split genome
        :param output_folder: The path to output the counting result
        :return:
        '''
        if os.path.isdir(name):
            return self.count_single_nt_folder(name, output_folder=output_folder)
        elif os.path.isfile(name):
            return self.count_single_nt_file(name, output_folder=output_folder)

    def count_single_nt_folder(self, path, save_result=True, output_folder='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_folder, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ_matrix = np.load(nt_count_single_path)
            return
        fa_list = os.listdir(path)
        fa_list = [os.path.join(path, fa) for fa in fa_list]
        occ_matrix = np.zeros((len(fa_list), NT_TYPE_NUM))
        for i, fa in enumerate(fa_list):
            occ = helper.count_nucleotide(fa)
            occ_matrix[i, :] = occ
        self.single_count = occ_matrix
        if save_result:
            np.save(nt_count_single_path, occ_matrix)

    def count_single_nt_file(self, path, save_result=True, output_folder='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_folder, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ = np.load(nt_count_single_path)
            return occ
        occ = helper.count_nucleotide(path)
        self.single_count = occ
        if save_result:
            np.save(nt_count_single_path, occ)

    def calculate_frequency_product_single(self, single_matrix: np.ndarray, kmer_idx, k):
        '''
        For a single k-mer, calculate the product of the frequency of each single nucleotide, f_i*f_j...

        :param ndarray single_matrix: The frequency matrix for each single nucleotide, f_i. Has shape n*k or n, if for single sequence
        :param int kmer_idx: The index representing the k-mer
        :param int k: k-mer length
        :return:
        '''
        if len(single_matrix.shape) == 2:
            single_frequency_product = np.ones(single_matrix.shape[0])
        else:
            single_frequency_product = 1
        for i in range(k):
            if len(single_matrix.shape) == 2:
                single_frequency_product = single_frequency_product * single_matrix[:, int(kmer_idx % k)]
            elif len(single_matrix.shape) == 1:
                single_frequency_product = single_frequency_product * single_matrix[int(kmer_idx % k)]
            kmer_idx = kmer_idx / k
        return single_frequency_product

    def calculate_frequency_product_all(self, multiple_matrix, single_matrix, k):
        '''
        For all k-mer, calculate the product of the frequency of each single nucleotide, f_i*f_j...

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



    def calc_distance(self):
        self.host_single_count = self.count_single_nucleotide(self.host_path)

