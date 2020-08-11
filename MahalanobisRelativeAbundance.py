import numpy as np
import os
import util
from itertools import product
from multiprocessing import Pool
from tools import kmer_count
from functools import partial

# np.seterr(all='raise')

NT_TYPE_NUM = 4


class MahalanobisRelativeAbundance:
    def __init__(
        self,
        subject_directory_path,
        query_directory_path,
        k=3,
        temp_directory_path="temp_dir",
        thread=1,
        recalculate=False,
    ):
        # Initialization of variables
        self.subject_directory_path = subject_directory_path
        self.query_directory_path = query_directory_path

        self.recalculate = recalculate

        self.k = k
        nucleotides = ["A", "C", "G", "T"]
        self.mer_to_idx = {}
        for i, mer in enumerate(product(nucleotides, repeat=k)):
            self.mer_to_idx["".join(list(mer))] = i

        # Number of threads to use in counting k-mer,
        # not use in calculation of M-dist as multi threading might slow down numpy calculation
        self.thread = thread

        # Create temporary directory for storing intermediate contents
        self.temp_directory_path = temp_directory_path
        os.makedirs(temp_directory_path, exist_ok=True)

        # Create buffered result directory to avoid re-calculation on repeated experiment
        self.buffered_result_dir_path = os.path.join(temp_directory_path, f"MahalanobisRelativeAbundance_k{k}")
        os.makedirs(self.buffered_result_dir_path, exist_ok=True)

        print("Preparing splitted fasta for subject sequences")
        subject_split_directory_path = os.path.join(temp_directory_path, "subject_split")
        os.makedirs(subject_split_directory_path, exist_ok=True)

        # Split each subject into their target temp directory
        self.subject_list = os.listdir(subject_directory_path)
        self.subject_name_list = []
        for subject in self.subject_list:
            subject_path = os.path.join(self.subject_directory_path, subject)
            subject_name = os.path.splitext(subject)[0]
            self.subject_name_list.append(subject_name)
            subject_split_path = os.path.join(subject_split_directory_path, subject_name)
            if not os.path.exists(subject_split_path):
                os.makedirs(subject_split_path)
                util.split_fasta_by_size(subject_path, subject_split_path)
            else:
                continue
        # Subject directory is set to split directory, original subjects are then not used
        self.subject_directory_path = subject_split_directory_path

        # Setting variable for subject, usually be host
        # list with n (k, 4) ndarray, k is number of 5k base fragments
        self.subject_single_count = 0
        self.subject_kmer_count = 0  # (n, 4**k) ndarray
        self.subject_single_freq = 0
        self.subject_kmer_freq = 0
        self.subject_relative_abundance = 0  # (n, 4**k) ndarray

        # Setting variable for query, usually be plasmid
        self.query_single_count = 0  # (n, 4)
        self.query_kmer_count = 0
        self.query_single_freq = 0
        self.query_kmer_freq = 0
        self.query_relative_abundance = 0  # (n) dim ndarray

        self.nan_substitute = 1e-4

    def calc_distance(self, thread=1):
        subject_directory_list = os.listdir(self.subject_directory_path)
        subject_directory_list = [
            os.path.join(self.subject_directory_path, f) for f in subject_directory_list if not f.startswith(".")
        ]
        subject_directory_list.sort()

        query_directory_list = os.listdir(self.query_directory_path)
        query_directory_list = [
            os.path.join(self.query_directory_path, f) for f in query_directory_list if not f.startswith(".")
        ]
        query_directory_list.sort()

        # Prepare function with preset arguments for map and pool.map
        count_single_nucleotide_with_preset = partial(
            self.count_single_nucleotide,
            output_directory_path=self.buffered_result_dir_path,
            recalculate=self.recalculate,
        )
        count_kmer_with_preset = partial(
            self.count_kmer, output_directory_path=self.buffered_result_dir_path, recalculate=self.recalculate,
        )

        if thread == 1:
            print("Counting host nucleotide and kmer frequency...")
            self.subject_single_count = (*map(self.count_single_nucleotide, subject_directory_list),)
            self.subject_kmer_count = (*map(self.count_kmer, subject_directory_list),)

            self.subject_single_freq = (*map(self.normalize, self.subject_single_count),)
            self.subject_kmer_freq = (*map(self.normalize, self.subject_kmer_count),)

            print("Calculating host relative abundance...")
            self.subject_relative_abundance = (
                *map(self.calculate_relative_abundance, self.subject_kmer_freq, self.subject_single_freq,),
            )

            print("Counting plasmid nucleotide and kmer frequency...")
            self.query_single_count = (*map(self.count_single_nucleotide, query_directory_list),)
            self.query_kmer_count = (*map(self.count_kmer, query_directory_list),)

            self.query_single_freq = (*map(self.normalize, self.query_single_count),)
            self.query_kmer_freq = (*map(self.normalize, self.query_kmer_count),)

            print("Calculating plasmid relative abundance...")
            self.query_relative_abundance = (
                *map(self.calculate_relative_abundance, self.query_kmer_freq, self.query_single_freq,),
            )

            print("Calculating Mahalanobis distance...")
            self.calculate_mahalanobis_distance()

        elif thread != 1:
            p = Pool(thread)

            print("Counting host nucleotide and kmer frequency...")
            # self.subject_single_count = p.map(
            #     self.count_single_nucleotide, subject_directory_list)
            self.subject_single_count = p.map(count_single_nucleotide_with_preset, subject_directory_list)
            self.subject_kmer_count = p.map(count_kmer_with_preset, subject_directory_list)

            self.subject_single_freq = (*map(self.normalize, self.subject_single_count),)
            self.subject_kmer_freq = (*map(self.normalize, self.subject_kmer_count),)

            print("Calculating host relative abundance...")
            # self.subject_relative_abundance = p.starmap(self.calculate_relative_abundance, zip(self.subject_kmer_freq, self.subject_single_freq))
            self.subject_relative_abundance = (
                *map(self.calculate_relative_abundance, self.subject_kmer_freq, self.subject_single_freq,),
            )

            print("Counting plasmid nucleotide and kmer frequency...")
            self.query_single_count = p.map(count_single_nucleotide_with_preset, query_directory_list)
            self.query_kmer_count = (*map(count_kmer_with_preset, query_directory_list),)

            self.query_single_freq = (*map(self.normalize, self.query_single_count),)
            self.query_kmer_freq = (*map(self.normalize, self.query_kmer_count),)

            print("Calculating plasmid relative abundance...")
            # self.plasmid_relative_abundance = p.starmap(self.calculate_relative_abundance,
            #                                             zip(self.plasmid_kmer_freq, self.plasmid_single_freq))
            self.query_relative_abundance = (
                *map(self.calculate_relative_abundance, self.query_kmer_freq, self.query_single_freq,),
            )
            p.close()

            print("Calculating Mahalanobis distance...")
            self.calculate_mahalanobis_distance()

            p.close()

        return self.mahalanobis_distance

    @staticmethod
    def count_single_nucleotide(path, output_directory_path="temp_dir", recalculate=False):
        """
        Count the number of occurrence of single nucleotides in a genome

        :param str path: The name of the organism, can be a fasta file or folder containing split genome
        :param output_directory_path: The path to output the counting result
        :return:
        """
        if os.path.isdir(path):
            return MahalanobisRelativeAbundance.count_single_nt_directory(
                path, output_directory_path=output_directory_path, recalculate=recalculate,
            )
        elif os.path.isfile(path):
            return MahalanobisRelativeAbundance.count_single_nt_file(
                path, output_directory_path=output_directory_path, recalculate=recalculate,
            )

    @staticmethod
    def count_single_nt_directory(path, output_directory_path="temp_dir", recalculate=False):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_directory_path, f"{genome_name}_single_nt.npy")
        if os.path.exists(nt_count_single_path) and not recalculate:
            occ_matrix = np.load(nt_count_single_path)
            return occ_matrix
        fa_list = os.listdir(path)
        fa_list.sort()
        fa_list = [os.path.join(path, fa) for fa in fa_list]
        occ_matrix = np.zeros((len(fa_list), NT_TYPE_NUM))
        for i, fa in enumerate(fa_list):
            occ = util.count_nucleotide(fa)
            occ_matrix[i, :] = occ
        np.save(nt_count_single_path, occ_matrix)
        return occ_matrix

    @staticmethod
    def count_single_nt_file(path, output_directory_path="temp_dir", recalculate=False):
        genome_name = os.path.split(path)[-1]
        nt_count_single_path = os.path.join(output_directory_path, f"{genome_name}_single_nt.npy")
        if os.path.exists(nt_count_single_path) and not recalculate:
            occ = np.load(nt_count_single_path)
            return occ
        occ = util.count_nucleotide(path)
        np.save(nt_count_single_path, occ)
        return occ

    def count_kmer(self, path, output_directory_path="temp_dir", recalculate=False):
        """
        Count the number of occurrence of kmers in a genome

        :param str path: The name of the organism, can be a fasta file or folder containing split genome
        :param output_directory_path: The path to output the counting result
        :return:
        """
        if os.path.isdir(path):
            return self.count_kmer_directory(
                path, output_directory_path=output_directory_path, recalculate=recalculate,
            )
        elif os.path.isfile(path):
            return self.count_kmer_file(path, output_directory_path=output_directory_path, recalculate=recalculate,)

    def exec_kmer_count(self, fasta_file_path):
        """
        Count the number of each kmer for given fasta file
        :param fasta_file_path:
        :return k_mer_count: a length n numpy array, where n = NT_TYPE_NUM ** k, containing the kmer count for given fasta file
        """

        k_mer_count = kmer_count(fasta_file_path, self.thread, False, self.k)
        k_mer_count = np.array(k_mer_count)

        return k_mer_count

    def count_kmer_directory(self, directory_path, output_directory_path="temp_dir", recalculate=False):

        genome_name = os.path.split(directory_path)[-1]
        kmer_count_path = os.path.join(output_directory_path, f"{genome_name}_kmer.npy")

        if os.path.exists(kmer_count_path) and not recalculate:
            k_mer_count = np.load(kmer_count_path)
            return k_mer_count

        fasta_list = os.listdir(directory_path)

        k_mer_count = np.zeros((len(fasta_list), NT_TYPE_NUM ** self.k))

        for i, fasta_filename in enumerate(fasta_list):
            fasta_file_path = os.path.join(directory_path, fasta_filename)
            fasta_kmer_count = self.exec_kmer_count(fasta_file_path)
            k_mer_count[i, :] = fasta_kmer_count

        np.save(kmer_count_path, k_mer_count)

        return k_mer_count

    def count_kmer_file(self, fasta_file_path, output_directory_path="temp_dir", recalculate=False):

        genome_name = os.path.split(fasta_file_path)[-1]
        kmer_count_path = os.path.join(output_directory_path, f"{genome_name}_kmer.npy")

        if os.path.exists(kmer_count_path) and not recalculate:
            k_mer_count = np.load(kmer_count_path)
            return k_mer_count

        k_mer_count = self.exec_kmer_count(fasta_file_path)

        np.save(kmer_count_path, k_mer_count)

        return k_mer_count

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
                single_frequency_product = single_frequency_product * single_frequency[:, int(kmer_idx % 4)]
            elif len(single_frequency.shape) == 1:
                single_frequency_product = single_frequency_product * single_frequency[int(kmer_idx % 4)]
            kmer_idx = kmer_idx / 4
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
        # frequency_product = np.zeros(multiple_frequency.shape)

        frequency_product = np.ones(multiple_frequency.shape)
        kmer_index = np.vstack(list(product(range(4), repeat=k)))

        for i in range(k):
            if len(multiple_frequency.shape) == 2:
                intermediate_frequency = single_frequency[:, kmer_index[:, i]]
                frequency_product = frequency_product * intermediate_frequency
            else:
                intermediate_frequency = single_frequency[kmer_index[:, i]]
                frequency_product = frequency_product * intermediate_frequency
        return frequency_product

        for i in range(multiple_frequency.shape[-1]):
            if len(multiple_frequency.shape) == 2:
                frequency_product[:, i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(
                    single_frequency, i, k
                )
            else:
                frequency_product[i] = MahalanobisRelativeAbundance.calculate_frequency_product_single(
                    single_frequency, i, k
                )
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

    @staticmethod
    def mahalanobis(u: np.ndarray, v: np.ndarray):
        u = u.flatten()
        v_mean = v.mean(axis=0)
        diff = u - v_mean
        cov = np.cov(v.T)
        cov_inv = np.linalg.inv(cov)
        distance = diff[None, :].dot(cov_inv).dot(diff)
        return distance

    def calculate_mahalanobis_distance(self):
        self.mahalanobis_distance = np.zeros((len(self.query_relative_abundance), len(self.subject_relative_abundance)))
        for i, plasmid_ra in enumerate(self.query_relative_abundance):
            for j, host_ra in enumerate(self.subject_relative_abundance):
                plasmid_relative_abundance = plasmid_ra.flatten()
                host_mean = host_ra.mean(axis=0)
                diff = plasmid_relative_abundance - host_mean
                cov = np.cov(host_ra.T)
                try:
                    cov_inv = np.linalg.inv(cov)
                except np.linalg.LinAlgError:
                    cov_inv = np.linalg.pinv(cov)
                distance = diff[None, :].dot(cov_inv).dot(diff)
                self.mahalanobis_distance[i, j] = distance

    @staticmethod
    def normalize(ndarray: np.ndarray):
        if len(ndarray.shape) == 2:
            # return np.nan_to_num(ndarray / ndarray.sum(axis=1)[:, None])
            return np.nan_to_num(
                ndarray / ndarray.sum(axis=1, keepdims=True), nan=MahalanobisRelativeAbundance.nan_substitute
            )
        if len(ndarray.shape) == 1:
            return np.nan_to_num(ndarray / ndarray.sum(), nan=MahalanobisRelativeAbundance.nan_substitute)
