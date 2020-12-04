import ProfilingTools as pf
from load_data import open_profile_from_tsv
import EMDUnifrac as unifrac
import numpy as np
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import itertools as it
import pandas as pd

def pairwise_unifrac(dir):
    '''
    Computes pairwise unifrac distance among profiles in a given directory
    :param dir: a directory containing profile files
    :return: a matrix of pairwise distances
    '''
    file_lst = os.listdir(dir)
    sample_lst = [os.path.splitext(profile)[0] for profile in file_lst]
    cur_dir = os.getcwd()
    os.chdir(dir)
    matrix = pd.DataFrame(columns=sample_lst, index=sample_lst)
    print(matrix)
    for pair in it.combinations(file_lst, 2):
        profile_list1 = open_profile_from_tsv(pair[0], False)
        profile_list2 = open_profile_from_tsv(pair[1], False)
        name1, metadata1, profile1 = profile_list1[0]
        name2, metadata2, profile2 = profile_list2[0]
        profile1 = pf.Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: 1)
        profile2 = pf.Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: 1)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (weighted, _) = unifrac.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        print(weighted)

    os.chdir(cur_dir)


if __name__ == '__main__':
    #Rickettsiales = pf.parse_taxid_file('data/Rickettsiales.txt')
    #pf.create_profile(Rickettsiales, 'data/testdir1', 'Rickettsiales.profile')
    pairwise_unifrac('data/testdir1')
