import numpy as np
import os
import helper
from itertools import product

NT_TYPE_NUM = 4


class MahalanobisRelativeAbundance:
    def __init__(self, host_path, plasmid_path, k=3):
        self.host_path = host_path
        self.plasmid_path = plasmid_path
        self.host_single_count = 0
        self.host_multiple_count = 0
        self.plasmid_single_count = 0
        self.plasmid_multiple_count = 0

        self.jellyfish_path = ''

        self.k = k
        nucleotides = ['A', 'C', 'G', 'T']
        self.mer_to_idx = {}
        for i, mer in enumerate(product(nucleotides, repeat=k)):
            self.mer_to_idx[''.join(list(mer))] = i

    @staticmethod
    def count_single_nucleotide(path, output_folder='temp_dir'):
        """
        Count the number of occurrence of single nucleotides in a genome

        :param str path: The name of the organism, can be a fasta file or folder containing split genome
        :param output_folder: The path to output the counting result
        :return:
        """
        if os.path.isdir(path):
            return MahalanobisRelativeAbundance.count_single_nt_folder(path, output_folder=output_folder)
        elif os.path.isfile(path):
            return MahalanobisRelativeAbundance.count_single_nt_file(path, output_folder=output_folder)

    @staticmethod
    def count_single_nt_folder(path, save_result=True, output_folder='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_folder, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ_matrix = np.load(nt_count_single_path)
            return occ_matrix
        fa_list = os.listdir(path)
        fa_list = [os.path.join(path, fa) for fa in fa_list]
        occ_matrix = np.zeros((len(fa_list), NT_TYPE_NUM))
        for i, fa in enumerate(fa_list):
            occ = helper.count_nucleotide(fa)
            occ_matrix[i, :] = occ
        if save_result:
            np.save(nt_count_single_path, occ_matrix)
        return occ_matrix

    @staticmethod
    def count_single_nt_file(path, save_result=True, output_folder='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_folder, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ = np.load(nt_count_single_path)
            return occ
        occ = helper.count_nucleotide(path)
        if save_result:
            np.save(nt_count_single_path, occ)
        return occ

    @staticmethod
    def calculate_frequency_product_single(single_frequency: np.ndarray, kmer_idx, k):
        """
        For a single k-mer, calculate the product of the frequency of each single nucleotide, f_i*f_j...

        :param ndarray single_frequency: The frequency matrix for each single nucleotide, f_i. Has shape n*k or n, if for single sequence
        :param int kmer_idx: The index representing the k-mer
        :param int k: k-mer length
        :return:
        """
        if len(single_frequency.shape) == 2:
            single_frequency_product = np.ones(single_frequency.shape[0])
        else:
            single_frequency_product = 1
        for i in range(k):
            if len(single_frequency.shape) == 2:
                single_frequency_product = single_frequency_product * single_frequency[:, int(kmer_idx % k)]
            elif len(single_frequency.shape) == 1:
                single_frequency_product = single_frequency_product * single_frequency[int(kmer_idx % k)]
            kmer_idx = kmer_idx / k
        return single_frequency_product

    @staticmethod
    def calculate_frequency_product_all(multiple_frequency, single_frequency, k):
        """
        For all k-mer, calculate the product of the frequency of each single nucleotide, f_i*f_j...

        :param multiple_frequency: frequency matrix for k-mer
        :param single_frequency: frequency matrix for single nucleotide
        :param k: k-mer length
        :return:
        """
        frequency_product = np.zeros(multiple_frequency.shape)
        for i in range(multiple_frequency.shape[-1]):
            if len(multiple_frequency.shape) == 2:
                frequency_product[:, i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(single_frequency, i, k)
            else:
                frequency_product[i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(single_frequency, i, k)
        return frequency_product

    def exec_jellyfish(self, fasta_path):
        """
        Runs jellyfish for a given fasta file and returns the count of each kmer
        :param fasta_path:
        :param hash_path:
        :param jellyfish_path:
        :return:
        """

        hash_path = fasta_path + '_hash'

        command = '{} count -m {} -s 200M -t 8 -o {} {}'.format(self.jellyfish_path, self.k, hash_path, fasta_path)
        if not os.path.exists(hash_path + '_0'):
            os.system(command)
        read_hash_command = '{} dump -c {}_0'.format(self.jellyfish_path, hash_path)
        counts = os.popen(read_hash_command).read()
        counts = counts.split('\n')[:-1]
        k_mer_count = np.zeros(NT_TYPE_NUM ** self.k)
        for mer_count in counts:
            (mer, count) = mer_count.split(' ')
            k_mer_count[self.mer_to_idx[mer]] = count
        return k_mer_count


    def count_kmer_folder(self, path, save_result=True, output_folder='temp_dir'):
        genome_name = os.path.split(path)[-1]

        kmer_count_path = os.path.join(output_folder, genome_name + '_kmer.npy')
        if os.path.exists(kmer_count_path):
            kmer_count = np.load(kmer_count_path)
            return kmer_count

        fasta_list = os.listdir(path)

        kmer_count = np.zeros((len(fasta_list), NT_TYPE_NUM ** self.k))

        for fasta_filename in fasta_list:
            fasta_path = os.path.join(path, fasta_filename)
            fasta_kmer_count = self.exec_jellyfish(fasta_path)


    def count_kmer_file(self, plasmid_name):
        plasmid_name =
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
            k_mer_count[self.mer_to_idx[m]] = c
        k_mer_count_filename = os.path.join(plasmid_matrix_dir, plasmid_name) + '_k{}'.format(k)
        np.save(k_mer_count_filename, k_mer_count)


    def calculate_relative_abundance(self, multiple_matrix, single_matrix, k):
        """
        Calculate the relative abundance, f_ij.../f_i*f_j...

        :param ndarray multiple_matrix: Frequency matrix for k-mer
        :param ndarray single_matrix: Frequency matrix for single nucleotide
        :param int k: k-mer length
        :return:
        """
        frequency_product = self.calculate_frequency_product_all(multiple_matrix, single_matrix, k)
        multiple_matrix[multiple_matrix == 0] = 1
        frequency_product[frequency_product == 0] = 1
        relative_abundance = multiple_matrix / frequency_product
        return relative_abundance



    def calc_distance(self):
        self.host_single_count = self.count_single_nucleotide(self.host_path)

