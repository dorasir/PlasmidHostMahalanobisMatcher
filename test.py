from multiprocessing import Pool
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np
import argparse
from itertools import product

# class Add:
#     def add(self, a, b):
#         return a + b
#
#     def exec_add(self):
#         a = [1, 2, 3]
#         b = [4, 5, 6]
#         p = Pool(2)
#         print(p.starmap(self.add, product(a, b)))


if __name__ == '__main__':
    # a = Add()
    # a.exec_add()

    parser = argparse.ArgumentParser()
    parser.add_argument('-q', dest='query_path', required=True)
    parser.add_argument('-s', dest='subject_path', required=True)
    parser.add_argument('-t', dest='thread_number', default=1, type=int)
    parser.add_argument('-o', dest='output_path', required=True)

    args = parser.parse_args()

    query_path = args.query_path
    subject_path = args.subject_path
    thread_number = args.thread_number
    output_path = args.output_path

    t = MahalanobisRelativeAbundance(subject_path, query_path, jellyfish_path='./jellyfish-linux', thread=thread_number)
    distance = t.calc_distance(thread=thread_number)
    np.save(output_path, distance)
    #
    # t = MahalanobisRelativeAbundance('data/hosts/splited', 'data/plasmids', jellyfish_path='./jellyfish-linux')
    # distance = t.calc_distance(thread=6)
    # np.save('result.npy', distance)