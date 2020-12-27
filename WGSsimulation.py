import ProfilingTools as pf
from load_data import open_profile_from_tsv
import EMDUnifrac as unifrac
import numpy as np
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import itertools as it
import pandas as pd
from skbio.diversity import beta_diversity
from skbio import TreeNode

#tree = TreeNode.read('data/gg/gg_13_5_otu/trees/99_otus.tree')
#ann_tree = TreeNode.read('data/gg/gg_13_5_otus_99_annotated.tree')

def pairwise_unifrac(dir):
    '''
    Computes pairwise unifrac distance among profiles in a given directory
    :param dir: a directory containing profile files
    :return: a matrix of pairwise distances
    '''
    cur_dir = os.getcwd()
    file_lst = os.listdir(dir) #list files in the directory
    os.chdir(dir)
    sample_lst = [os.path.splitext(profile)[0] for profile in file_lst] #remove extension
    # enumerate sample_lst, for filling matrix
    id_dict = dict()
    for i, id in enumerate(file_lst):
        id_dict[id] = i
    #initialize matrix
    dim = len(file_lst)
    dist_matrix = np.zeros(shape=(dim, dim))
    for pair in it.combinations(file_lst, 2):
        id_1,id_2 = pair[0], pair[1]
        i,j = id_dict[id_1], id_dict[id_2]
        profile_list1 = open_profile_from_tsv(id_1, False)
        profile_list2 = open_profile_from_tsv(id_2, False)
        name1, metadata1, profile1 = profile_list1[0]
        name2, metadata2, profile2 = profile_list2[0]
        profile1 = pf.Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: 1)
        profile2 = pf.Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: 1)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (weighted, _) = unifrac.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        dist_matrix[i][j] = dist_matrix[j][i] = weighted
    os.chdir(cur_dir)
    return sample_lst, dist_matrix


if __name__ == '__main__':
    #Rickettsiales = pf.parse_taxid_file('data/Rickettsiales.txt')
    #pf.create_profile(Rickettsiales, 'data/testdir1', 'Rickettsiales.profile')
    (id, dist_matrix) = pairwise_unifrac('data/gg_profile_test')
    print(id, dist_matrix)
