import ProfilingTools as pf
from load_data import open_profile_from_tsv
import EMDUnifrac as unifrac
import numpy as np
import os
import pandas as pd
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
import itertools as it
from skbio import TreeNode
import json
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix
import random

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

def parse_file_by_col(file, col):
    '''
    Parse a file and store a column into a list
    :param file:
    :param col:
    :return:
    '''
    lst = []
    with open(file, 'r') as f:
        for line in f.readlines():
            lst.append(line.split('\t')[col])
    f.close()
    return lst

def filter_otus(lst1, lst2):
    '''
    return a list of otus in lst 1 that is also in lst 2.
    to filter off otus not found on the phylogenetic tree
    :param lst1: a list of otus in concern
    :param lst2: a list of nodes present in the tree
    :return:
    '''
    return list(set(lst1) & set(lst2))

def write_list_to_file(file, lst):
    with open(file, 'w') as f:
        for x in lst:
            f.write("%s\n" %x)

def pick_otu(node, otu_ref, num):
    '''
    filter leave nodes of node n to pick only those in otu_ref, and pick num of them randomly
    :param node:
    :param otu_ref: list of otus to check against
    :param num: number of leaves to be picked, if not enough, will be replaced by the max number of nodes pickable
    :return:
    '''
    if node.is_leaf():
        print("%s is a leaf" %node)
        return
    leaves = node.get_leaves()
    leaf_names = []
    for l in leaves:
        leaf_names.append(l.name)
    filtered = filter_otus(leaf_names, otu_ref)
    if num > len(filtered):
        print("%d pickable nodes" %len(filtered))
        num = len(filtered)
    picked = []
    for i in range(num):
        r = random.randint(1, len(filtered))
        if filtered[r] not in picked:
            picked.append(filtered[r])
    return picked


#def get_taxid(otu):

if __name__ == '__main__':

    (id, dist_matrix) = pairwise_unifrac('data/gg_profile_test')
    print(id, dist_matrix)
    metadata = {
        'gg_test1':{'environment': 'A'},
        'gg_test2':{'environment': 'A'},
        'gg_test3':{'environment': 'B'},
        'gg_test4':{'environment': 'B'}
    }
    df = pd.DataFrame.from_dict(metadata, orient='index')
    dm = DistanceMatrix(dist_matrix, id)
    dist_pc = pcoa(dm)
    dist_pc.plot(df=df, column="environment", cmap='Set1')
    plt.show()

