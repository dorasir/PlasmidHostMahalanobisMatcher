from multiprocessing import Pool
from MahalanobisRelativeAbundance import MahalanobisRelativeAbundance
import numpy as np

class Add:
    def add(self, a, b):
        return a + b

    def exec_add(self):
        a = [1, 2, 3]
        b = [4, 5, 6]
        p = Pool(2)
        print(p.starmap(self.add, zip(a, b)))


if __name__ == '__main__':
    t = MahalanobisRelativeAbundance('data/hosts/splited', 'data/plasmids', jellyfish_path='./jellyfish-linux')
    distance = t.calc_distance(thread=2)
    np.save('result.npy', distance)