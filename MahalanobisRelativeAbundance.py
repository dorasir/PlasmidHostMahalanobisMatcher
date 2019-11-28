import numpy as np
import os
import helper
from itertools import product
from multiprocessing import Pool

NT_TYPE_NUM = 4


class MahalanobisRelativeAbundance:
    def __init__(self, host_directory_path, plasmid_directory_path, k=3, temp_directory_path='temp_dir', jellyfish_path='./jellyfish-macosx'):

        self.host_directory_path = host_directory_path
        self.plasmid_directory_path = plasmid_directory_path

        self.host_single_count = 0  # list with n (k, 4) ndarray, k is number of 5k base fragments
        self.host_kmer_count = 0  # (n, 4**k) ndarray
        self.host_single_freq = 0
        self.host_kmer_freq = 0
        self.host_relative_abundance = 0  # (n, 4**k) ndarray

        self.plasmid_single_count = 0  # (n, 4)
        self.plasmid_kmer_count = 0
        self.plasmid_single_freq = 0
        self.plasmid_kmer_freq = 0
        self.plasmid_relative_abundance = 0  # (n) dim ndarray

        self.jellyfish_path = jellyfish_path

        self.k = k
        nucleotides = ['A', 'C', 'G', 'T']
        self.mer_to_idx = {}
        for i, mer in enumerate(product(nucleotides, repeat=k)):
            self.mer_to_idx[''.join(list(mer))] = i

        self.temp_directory_path = temp_directory_path
        if not os.path.exists(temp_directory_path):
            print('Temp directory not exist, creating...')
        else:
            if not os.path.isdir(temp_directory_path):
                print('Exists file with name conflict, exiting...')

    @staticmethod
    def count_single_nucleotide(path, output_directory_path='temp_dir'):
        """
        Count the number of occurrence of single nucleotides in a genome

        :param str path: The name of the organism, can be a fasta file or folder containing split genome
        :param output_directory_path: The path to output the counting result
        :return:
        """
        if os.path.isdir(path):
            return MahalanobisRelativeAbundance.count_single_nt_directory(path,
                                                                          output_directory_path=output_directory_path)
        elif os.path.isfile(path):
            return MahalanobisRelativeAbundance.count_single_nt_file(path, output_directory_path=output_directory_path)

    @staticmethod
    def count_single_nt_directory(path, save_result=True, output_directory_path='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_directory_path, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ_matrix = np.load(nt_count_single_path)
            return occ_matrix
        fa_list = os.listdir(path)
        fa_list.sort()
        fa_list = [os.path.join(path, fa) for fa in fa_list]
        occ_matrix = np.zeros((len(fa_list), NT_TYPE_NUM))
        for i, fa in enumerate(fa_list):
            occ = helper.count_nucleotide(fa)
            occ_matrix[i, :] = occ
        if save_result:
            np.save(nt_count_single_path, occ_matrix)
        return occ_matrix

    @staticmethod
    def count_single_nt_file(path, save_result=True, output_directory_path='temp_dir'):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_directory_path, genome_name + '_single_nt.npy')
        if os.path.exists(nt_count_single_path):
            occ = np.load(nt_count_single_path)
            return occ
        occ = helper.count_nucleotide(path)
        if save_result:
            np.save(nt_count_single_path, occ)
        return occ

    def count_kmer(self, path, output_directory_path='temp_dir'):
        """
        Count the number of occurrence of kmers in a genome

        :param str path: The name of the organism, can be a fasta file or folder containing split genome
        :param output_directory_path: The path to output the counting result
        :return:
        """
        if os.path.isdir(path):
            return self.count_kmer_directory(path, output_directory_path=output_directory_path)
        elif os.path.isfile(path):
            return self.count_kmer_file(path, output_directory_path=output_directory_path)

    def exec_jellyfish(self, fasta_file_path, temp_directory_path='temp_dir'):
        """
        Runs jellyfish for a given fasta file and returns the count of each kmer
        :param fasta_file_path:
        :param hash_path:
        :param jellyfish_path:
        :return k_mer_count: a length n numpy array, where n = NT_TYPE_NUM ** k, containing the kmer count for given fasta file
        """
        genome_name = os.path.split(os.path.split(fasta_file_path)[0])
        hash_path = genome_name + '_' + os.path.split(fasta_file_path)[-1] + '_hash'
        hash_path = os.path.join(temp_directory_path, hash_path)

        command = '{} count -m {} -s 200M -t 8 -o {} {}'.format(self.jellyfish_path, self.k, hash_path, fasta_file_path)
        if not os.path.exists(hash_path + '_0'):
            os.system(command)
        read_hash_command = '{} dump -c {}_0'.format(self.jellyfish_path, hash_path)
        counts = os.popen(read_hash_command).read()
        counts = counts.split('\n')[:-1]
        k_mer_count = np.zeros(NT_TYPE_NUM ** self.k)
        for mer_count in counts:
            (mer, count) = mer_count.split(' ')
            k_mer_count[self.mer_to_idx[mer]] = count

        # os.remove(hash_path + '_0')

        return k_mer_count

    def count_kmer_directory(self, directory_path, save_result=True, output_directory_path='temp_dir',
                             overwrite_existed_result=False):

        genome_name = os.path.split(directory_path)[-1]
        kmer_count_path = os.path.join(output_directory_path, genome_name + '_kmer.npy')

        if not overwrite_existed_result and os.path.exists(kmer_count_path):
            kmer_count = np.load(kmer_count_path)
            return kmer_count

        fasta_list = os.listdir(directory_path)

        kmer_count = np.zeros((len(fasta_list), NT_TYPE_NUM ** self.k))

        for i, fasta_filename in enumerate(fasta_list):
            fasta_file_path = os.path.join(directory_path, fasta_filename)
            fasta_kmer_count = self.exec_jellyfish(fasta_file_path)
            kmer_count[i, :] = fasta_kmer_count

        if save_result:
            np.save(kmer_count_path, kmer_count)

        return kmer_count

    def count_kmer_file(self, fasta_file_path, save_result=False, output_directory_path='temp_dir',
                        overwrite_existed_result=False):

        genome_name = os.path.split(fasta_file_path)[-1]
        kmer_count_path = os.path.join(output_directory_path, genome_name + '_kmer.npy')

        if not overwrite_existed_result and os.path.exists(kmer_count_path):
            kmer_count = np.load(kmer_count_path)
            return kmer_count

        kmer_count = self.exec_jellyfish(fasta_file_path)

        if save_result:
            np.save(kmer_count_path, kmer_count)

        return kmer_count

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
                frequency_product[:, i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(
                    single_frequency, i, k)
            else:
                frequency_product[i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(single_frequency,
                                                                                                       i, k)
        return frequency_product

    def calculate_relative_abundance(self, multiple_matrix, single_matrix):
        """
        Calculate the relative abundance, f_ij.../f_i*f_j...

        :param ndarray multiple_matrix: Frequency matrix for k-mer
        :param ndarray single_matrix: Frequency matrix for single nucleotide
        :param int k: k-mer length
        :return:
        """
        frequency_product = self.calculate_frequency_product_all(multiple_matrix, single_matrix, self.k)
        multiple_matrix[multiple_matrix == 0] = 1
        frequency_product[frequency_product == 0] = 1
        relative_abundance = multiple_matrix / frequency_product
        return relative_abundance

    def calculate_mahalanobis_distance(self):
        self.mahalanobis_distance = np.zeros((len(self.plasmid_relative_abundance), len(self.host_relative_abundance)))
        for i, plasmid_ra in enumerate(self.plasmid_relative_abundance):
            for j, host_ra in enumerate(self.host_relative_abundance):
                plasmid_relative_abundance = plasmid_ra.flatten()
                host_mean = host_ra.mean(axis=0)
                diff = plasmid_relative_abundance - host_mean
                cov = np.cov(host_ra.T)
                cov_inv = np.linalg.inv(cov)
                distance = diff[None, :].dot(cov_inv).dot(diff)
                self.mahalanobis_distance[i, j] = distance

    @staticmethod
    def normalize(ndarray):
        if len(ndarray.shape) == 2:
            return ndarray / ndarray.sum(axis=1)[:, None]
        if len(ndarray.shape) == 1:
            return ndarray / ndarray.sum()

    def calc_distance(self, thread=1):
        host_directory_list = os.listdir(self.host_directory_path)
        host_directory_list = [os.path.join(self.host_directory_path, f) for f in host_directory_list if
                               not f.startswith('.')]
        host_directory_list.sort()

        plasmid_directory_list = os.listdir(self.plasmid_directory_path)
        plasmid_directory_list = [os.path.join(self.plasmid_directory_path, f) for f in plasmid_directory_list if
                                  not f.startswith('.')]
        plasmid_directory_list.sort()

        if thread == 1:
            self.host_single_count = *map(self.count_single_nucleotide, host_directory_list),
            self.host_kmer_count = *map(self.count_kmer, host_directory_list),

            self.host_single_freq = *map(self.normalize, self.host_single_count),
            self.host_kmer_freq = *map(self.normalize, self.host_kmer_count),

            self.host_relative_abundance = *map(self.calculate_relative_abundance, self.host_kmer_freq, self.host_single_freq),

            self.plasmid_single_count = *map(self.count_single_nucleotide, plasmid_directory_list),
            self.plasmid_kmer_count = *map(self.count_kmer, plasmid_directory_list),

            self.plasmid_single_freq = *map(self.normalize, self.plasmid_single_count),
            self.plasmid_kmer_freq = *map(self.normalize, self.plasmid_kmer_count),

            self.plasmid_relative_abundance = *map(self.calculate_relative_abundance, self.plasmid_kmer_freq, self.plasmid_single_freq),

            self.calculate_mahalanobis_distance()

            print('a')

        elif thread != 1:
            p = Pool(thread)

            self.host_single_count = p.map(self.count_single_nucleotide, host_directory_list)
            self.host_kmer_count = p.map(self.count_kmer, host_directory_list)

            self.host_single_freq = p.map(self.normalize, self.host_single_count)
            self.host_kmer_freq = p.map(self.normalize, self.host_kmer_count)

            self.host_relative_abundance = p.starmap(self.calculate_relative_abundance, zip(self.host_kmer_freq, self.host_single_freq))

            self.plasmid_single_count = p.map(self.count_single_nucleotide, plasmid_directory_list)
            self.plasmid_kmer_count = p.map(self.count_kmer, plasmid_directory_list)

            self.plasmid_single_freq = p.map(self.normalize, self.plasmid_single_count)
            self.plasmid_kmer_freq = p.map(self.normalize, self.plasmid_kmer_count)

            self.plasmid_relative_abundance = p.starmap(self.calculate_relative_abundance,
                                                        zip(self.plasmid_kmer_freq, self.plasmid_single_freq))

            self.calculate_mahalanobis_distance()

        return self.mahalanobis_distance

