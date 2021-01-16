import random
import pandas as pd
import os
import ProfilingTools as pf
import numpy as np
import EMDUnifrac as unifrac
from biom import Table
import ete3
from copy import deepcopy

def write_list_to_file(file, lst):
    with open(file, 'w') as f:
        for x in lst:
            f.write("%s\n" %x)

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

def pick_otu_random(node, map_dict, num):
    '''

    :param node:
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
    filtered_leaves = list(set(list(map_dict.keys())) & set(leaf_names))
    print('%d leaves are pickable' %len(filtered_leaves))
    # in case num is too big, return everything possible
    if num > len(filtered_leaves):
        print("%d pickable nodes" %len(filtered_leaves))
        return filtered_leaves
    picked = []
    picked_tax = []
    for i in range(num):
        r = random.randint(1, len(filtered_leaves))
        if filtered_leaves[r] not in picked:
            otu = filtered_leaves[r]
            picked.append(otu)
            picked_tax.append(map_dict[otu])
    return picked, picked_tax

def pick_by_dist(tree, node, map_dict, num, close_enough):
    if len(close_enough) < num:
        print("you ask for too much, I only have %d nodes satisfying your requirement" %len(close_enough))
        return
    picked=[]
    picked_tax = []
    for i in range(num):
        r = random.randint(1, len(close_enough))
        if close_enough[r] not in picked:
            otu = close_enough[r]
            picked.append(otu)
            picked_tax.append(map_dict[otu])
    return picked, picked_tax

def pick_for_two(tree, node1, node2, map_dict, num):
    picked_1 = []
    picked_tax_1 = []
    picked_2 = []
    picked_tax_2 = []
    pickable = list(map_dict.keys())
    random.shuffle(pickable)
    if len(pickable) < num*2:
        print("not enough to pick from")
        return
    print(node1.get_distance(node2))
    while len(picked_1) < num and len(picked_2) < num:
        node = pickable.pop()
        if node1.get_distance(tree&node) < node2.get_distance(tree&node):
            picked_1.append(node)
            picked_tax_1.append(map_dict[node])
        else:
            picked_2.append(node)
            picked_tax_2.append(map_dict[node])
        print("length of picked_1 is %d" %len(picked_1))
        print("length of picked_2 is %d" %len(picked_2))
    while len(picked_2) < num: #will only execute if picked_2 is not filled yet
        node = pickable.pop()
        if node2.get_distance(tree&node) < node1.get_distance(tree&node):
            picked_2.append(node)
            picked_tax_2.append(map_dict[node])
            print("length of picked_2 is %d" % len(picked_2))
    while len(picked_1) < num: #will only execute if picked_1 is not filled yet
        node = pickable.pop()
        if node1.get_distance(tree&node) < node2.get_distance(tree&node):
            picked_1.append(node)
            picked_tax_1.append(map_dict[node])
            print("length of picked_1 is %d" % len(picked_1))
    return picked_1, picked_tax_1, picked_2, picked_tax_2

def create_biom_table(sample_metadata, table_id, data, filename):
    '''
    :param sample_metadata:
    :param table_id:
    :param data: dictionary in the form of sample_id:list of otus
    :return:
    '''
    otus = []
    sample_id = [] #column index
    for key,value in list(data.items()):
        sample_id.append(key)
        otus = otus + value
    otu_ids = list(set(otus)) #row index unique otus
    df = pd.DataFrame(columns=sample_id, index=otu_ids)
    for sample in sample_id:
        sample_data = np.zeros(len(otu_ids))
        for i, otu in enumerate(otu_ids):
            if otu in list(data[sample]):
                sample_data[i] = 1
        df[sample] = sample_data #fill by column
    table = Table(data=df.to_numpy(), observation_ids=otu_ids, sample_ids=sample_id, observation_metadata=None,
                  sample_metadata=None, table_id=table_id)
    normed = table.norm(axis='sample', inplace=False)
    with open(filename, "w") as f:
        normed.to_tsv(direct_io=f)

def get_valid_otu(otu_file, taxid_file, filename):
    '''
    filter otus in otu_file against taxid_file, remove those without desired ranks and write out a file
    :param otu_file: a file with first column being otu and second column being accession
    :param taxid_file: a file with first column being accession and second column being taxid
    :return:
    '''
    # merge two dataframe by accession
    otu_acc_df = pd.read_table(otu_file, header=None, dtype=str)
    otu_acc_df.rename(columns={0:"otu", 1:"acc"}, inplace=True)
    acc_taxid_df = pd.read_table(taxid_file, header=None, usecols=[0,1], dtype=str)
    acc_taxid_df.rename(columns={0:"acc", 1:"taxid"}, inplace=True)
    merged_df = otu_acc_df.merge(acc_taxid_df, left_on='acc', right_on='acc')
    #drop rows with invalid/unsuitable taxid
    print(merged_df[0:5])
    print("size before dropping is %s" %len(merged_df))
    merged_df.dropna(inplace=True) #drop rows containing na values
    for i, row in merged_df.iterrows():
        if not pf.check_rank(row['taxid']):
            merged_df.drop([0])
    print("size after dropping is %s" %len(merged_df))
    merged_df.to_csv('otu_with_valid_taxid.txt', sep='\t', index=False)

def get_valid_otu_alternative(otu_file, taxid_file, filename):
    '''
    an alternative implementation of get_valid_otu using np vectorization
    possibly faster, to be tested and confirmed
    :param otu_file:
    :param taxid_file:
    :param filename:
    :return:
    '''
    # merge two dataframe by accession
    otu_acc_df = pd.read_table(otu_file, header=None, dtype=str)
    otu_acc_df.rename(columns={0: "otu", 1: "acc"}, inplace=True)
    acc_taxid_df = pd.read_table(taxid_file, header=None, usecols=[0, 1], dtype=str)
    acc_taxid_df.rename(columns={0: "acc", 1: "taxid"}, inplace=True)
    merged_df = otu_acc_df.merge(acc_taxid_df, left_on='acc', right_on='acc')
    # drop rows with invalid/unsuitable taxid
    print(merged_df[0:5])
    print("size before dropping is %s" % len(merged_df))
    merged_df.dropna(inplace=True)  # drop rows containing na values
    merged_df['taxid'] = _check_rank_list(merged_df['taxid'].values)
    merged_df.dropna(inplace=True)
    print("size after dropping is %s" % len(merged_df))
    merged_df.to_csv(filename, sep='\t', index=False)

def filter_against_tree(otu_file, nodes):
    '''
    remove nodes in otu_file that are not found on the tree. For the remaining nodes, parse into a dictionary
    :param otu_file:
    :param nodes: nodes on a selected tree
    :return: a valid {otu:taxid} dict
    '''
    df = pd.read_table(otu_file, dtype=str)
    #create dict
    map_dict = df.set_index('otu')['taxid'].to_dict()
    #print(df[0:5])
    print("%d pairs present" %len(map_dict))
    #filter out unwanted pairs
    filtered_dict = dict()
    all_keys = list(map_dict.keys())
    present = list(set(all_keys) & set(nodes)) #otus present on the tree
    for otu in present:
        filtered_dict[otu] = map_dict[otu]
    print('%d pairs left' %len(filtered_dict))
    return filtered_dict

def _check_rank_list(id_list):
    '''
    slight variation of filter_id_list in ProfilingTools, instead of removing invalid id, set value to np.nan
    :param id_list:
    :return:
    '''
    for i, id in enumerate(id_list):
        if not pf.check_rank(id):
            id_list[i] = np.nan
    return id_list

def search_by_distance(tree, node, map_dict, dist, num):
    '''
    search for all the nodes within the distance dist from the node, all pickable
    :param node:
    :param dist:
    :return:
    '''
    if type(node) == str:
        node = tree&node
    match = []
    match_taxid = []
    pickable_copy = deepcopy(list(map_dict.keys()))
    random.shuffle(pickable_copy)
    while len(match) < num and len(pickable_copy)>0:
        n = pickable_copy.pop()
        if node.get_distance(tree&n) <= dist
            match.append(n)
            match_taxid.append(map_dict[n])
    if len(match) < num:
        print("%d pickable nodes within distance %d" %len(match) %dist)
    return match, match_taxid

def create_data(num_sam, tree, num_org, map_dict):
    '''
    return a dictionary consisting of num_env environments, with num_sam samples each
    :param num_env:
    :param num_sam:
    :return:
    '''
    otu_dict = dict()
    taxid_dict = dict()
    node1 = random.choice(tree.get_leaves())
    node2 = random.choice(tree.get_leaves())
    while node1.get_distance(node2) < 2.5:
        node1 = random.choice(tree.get_leaves())
        node2 = random.choice(tree.get_leaves())
    print(node1, node2)
    dist = node1.get_distance(node2)/2.0
    for i in range(num_sam):
        print("creating sample round %d" %i)
        env1_otu_key = "{}{}".format('env1sam', i)
        env2_otu_key = "{}{}".format('env2sam', i)
        (env1otu, env1tax) = search_by_distance(tree, node1, map_dict, dist, num_org)
        (env2otu, env2tax) = search_by_distance(tree, node1, map_dict, dist, num_org)
        otu_dict.update({env1_otu_key:env1otu})
        otu_dict.update({env2_otu_key:env2otu})
        taxid_dict.update({env1_otu_key: env1tax})
        taxid_dict.update({env2_otu_key: env2tax})
    return otu_dict, taxid_dict

#def create_metadata(data):


def test_merge():
    '''
    test merge method of pandas. checking of get_otu_taxid_dict function
    :return:
    '''
    #differ in number of rows and order of common key
    otu_acc_df = pd.read_table('test_merge_otu.txt', header=None, dtype=str)
    acc_taxid_df = pd.read_table('test_merge_taxid.txt', header=None, usecols=[0,1], dtype=str)
    merged_df = otu_acc_df.merge(acc_taxid_df, left_on=1, right_on=0)
    print(merged_df)


if __name__ == '__main__':
    os.chdir('data/taxid_otu_conversion')
    (T, l, nodes) = unifrac.parse_tree_file('../gg/gg_13_5_otus/trees/99_otus.tree')
    tree = ete3.TreeNode('../gg/gg_13_5_otus/trees/99_otus.tree', format=1, quoted_node_names=True)
    otu_tax_dict = filter_against_tree('otu_with_valid_taxid.txt', nodes)
    os.chdir('../gg_test_2')
    #pickable_otu = list(otu_tax_dict.keys())
    (data, taxdata) = create_data(25, tree, 200, otu_tax_dict)
    create_biom_table('meta', 'gg_test2', data, 'gg_test2.biom')
    for key, value in list(taxdata.items()):
        filename = "{}{}".format(key, '.profile')
        pf.create_profile(value, 'profiles', filename)