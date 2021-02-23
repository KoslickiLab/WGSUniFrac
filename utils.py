import random
import pandas as pd
import os
import ProfilingTools as pf
import numpy as np
import EMDUnifrac as unifrac
from biom import Table
import ete3
from copy import deepcopy
from multiprocessing import Pool, Manager
import time
from itertools import repeat
import sys

class Node:
    def __init__(self, name, abundance=0, tax=''):
        self.name = name
        self.tax = tax
        self.abundance = abundance


def setup():
    os.chdir('data')
    (T, l, nodes) = unifrac.parse_tree_file('gg/gg_13_5_otus/trees/99_otus.tree')
    # (T, l, nodes) = unifrac.parse_tree_file('99_otus.tree')
    tree = ete3.TreeNode('gg/gg_13_5_otus/trees/99_otus.tree', format=1, quoted_node_names=True)
    # tree = ete3.TreeNode('99_otus.tree', format=1, quoted_node_names=True)
    otu_tax_dict = filter_against_tree('taxid_otu_conversion/otu_with_valid_taxid.txt', nodes)
    return (tree, otu_tax_dict)


def write_list_to_file(file, lst):
    with open(file, 'w') as f:
        for x in lst:
            f.write("%s\n" % x)


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
        print("%s is a leaf" % node)
        return
    leaves = node.get_leaves()
    leaf_names = []
    for l in leaves:
        leaf_names.append(l.name)
    filtered_leaves = list(set(list(map_dict.keys())) & set(leaf_names))
    print('%d leaves are pickable' % len(filtered_leaves))
    # in case num is too big, return everything possible
    if num > len(filtered_leaves):
        print("%d pickable nodes" % len(filtered_leaves))
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


def pick_by_dist(tree, node, map_dict, dist, enough):
    '''
    return all nodes within distance dist of node node.
    :param tree:
    :param node:
    :param map_dict:
    :param num:
    :param close_enough:
    :return:
    '''
    if type(node) == str:
        node = tree & node
    picked = []
    picked_tax = []
    pickable = list(map_dict.keys())
    for i in range(len(pickable)):
        if len(picked) == enough:
            print("after looking at %d nodes I finally found the miserable amount you asked" % i)
            break
        n = pickable[i]
        if node.get_distance(tree & n) <= dist:
            picked.append(n)
            picked_tax.append(map_dict[n])
            print(len(picked))
    return picked, picked_tax


def pick_for_two(tree, node1, node2, map_dict, num):
    picked_1 = []
    picked_tax_1 = []
    picked_2 = []
    picked_tax_2 = []
    pickable = list(map_dict.keys())
    random.shuffle(pickable)
    if len(pickable) < num * 2:
        print("not enough to pick from")
        return
    print(node1.get_distance(node2))
    while len(picked_1) < num and len(picked_2) < num:
        node = pickable.pop()
        if node1.get_distance(tree & node) < node2.get_distance(tree & node):
            picked_1.append(node)
            picked_tax_1.append(map_dict[node])
        else:
            picked_2.append(node)
            picked_tax_2.append(map_dict[node])
        print("length of picked_1 is %d" % len(picked_1))
        print("length of picked_2 is %d" % len(picked_2))
    while len(picked_2) < num:  # will only execute if picked_2 is not filled yet
        node = pickable.pop()
        if node2.get_distance(tree & node) < node1.get_distance(tree & node):
            picked_2.append(node)
            picked_tax_2.append(map_dict[node])
            print("length of picked_2 is %d" % len(picked_2))
    while len(picked_1) < num:  # will only execute if picked_1 is not filled yet
        node = pickable.pop()
        if node1.get_distance(tree & node) < node2.get_distance(tree & node):
            picked_1.append(node)
            picked_tax_1.append(map_dict[node])
            print("length of picked_1 is %d" % len(picked_1))
    return picked_1, picked_tax_1, picked_2, picked_tax_2

def get_valid_otu(otu_file, taxid_file, filename):
    '''
    filter otus in otu_file against taxid_file, remove those without desired ranks and write out a file
    :param otu_file: a file with first column being otu and second column being accession
    :param taxid_file: a file with first column being accession and second column being taxid
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
    for i, row in merged_df.iterrows():
        if not pf.check_rank(row['taxid']):
            merged_df.drop([0])
    print("size after dropping is %s" % len(merged_df))
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
    # merge two datafame by accessiosn
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
    # create dict
    map_dict = df.set_index('otu')['taxid'].to_dict()
    # print(df[0:5])
    print("%d pairs present" % len(map_dict))
    # filter out unwanted pairs
    filtered_dict = dict()
    all_keys = list(map_dict.keys())
    present = list(set(all_keys) & set(nodes))  # otus present on the tree
    for otu in present:
        filtered_dict[otu] = map_dict[otu]
    print('%d pairs left' % len(filtered_dict))
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
        node = tree & node
    match = []
    match_taxid = []
    pickable_copy = deepcopy(list(map_dict.keys()))
    random.shuffle(pickable_copy)
    while len(match) < num and len(pickable_copy) > 0:
        n = pickable_copy.pop()
        if node.get_distance(tree & n) <= dist:
            match.append(n)
            print(len(match))
            match_taxid.append(map_dict[n])
    if len(match) < num:
        print("%d pickable nodes within distance %d" % len(match) % dist)
    return match, match_taxid


def create_data(num_sam, tree, dist, num_org, map_dict, enough):
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
    while node1.get_distance(node2) < 2.0:
        node1 = random.choice(tree.get_leaves())
        node2 = random.choice(tree.get_leaves())
    print(node1, node2)
    print(node1.get_distance(node2))
    (otu1, tax1) = pick_by_dist(tree, node1, map_dict, dist, enough)
    (otu2, tax2) = pick_by_dist(tree, node2, map_dict, dist, enough)
    if num_org > len(otu1):
        print("environment 1 does not have enough. Only %d present" % len(otu1))
        return
    if num_org > len(otu2):
        print("environment 1 does not have enough. Only %d present" % len(otu2))
        return
    print(len(otu1))
    print(len(otu2))
    for i in range(num_sam):
        print("creating sample round %d" % i)
        env1otu = []
        env1tax = []
        env2otu = []
        env2tax = []
        while len(env1otu) < num_org:
            r = random.randint(0, len(otu1) - 1)
            if otu1[r] not in env1otu:
                env1otu.append(otu1[r])
                env1tax.append(tax1[r])
        while len(env2otu) < num_org:
            r = random.randint(0, len(otu2) - 1)
            if otu1[r] not in env1otu:
                env2otu.append(otu2[r])
                env2tax.append(tax2[r])
        env1_otu_key = "{}{}".format('env1sam', i)
        env2_otu_key = "{}{}".format('env2sam', i)
        otu_dict.update({env1_otu_key: env1otu})
        otu_dict.update({env2_otu_key: env2otu})
        taxid_dict.update({env1_otu_key: env1tax})
        taxid_dict.update({env2_otu_key: env2tax})
    return otu_dict, taxid_dict


def create_data_simple(num_org, num_sample, sample_range, distance_dict, tax_dict):
    node1 = random.choice(list(distance_dict.keys()))
    node2 = distance_dict[node1][-1]
    print(node1)
    print(node2)
    data_dict = dict()
    env1_nodes = distance_dict[node1][:sample_range - 1]
    env2_nodes = distance_dict[node2][:sample_range - 1]
    #update the 2 lists above to contain Node object instead
    for i, node in enumerate(env1_nodes):
        #create Nodes, update tax
        new_node = Node(name=node, tax=tax_dict[node])
        env1_nodes[i] = new_node
    for i, node in enumerate(env2_nodes):
        new_node = Node(name=node, tax=tax_dict[node])
        env2_nodes[i] = new_node
    #create sample
    if num_org >= sample_range:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = deepcopy(env1_nodes)
            value2 = deepcopy(env2_nodes)
            random.shuffle(value1)
            random.shuffle(value2)
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    else:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = random.sample(env1_nodes, num_org)
            value2 = random.sample(env2_nodes, num_org)
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    return data_dict

def create_biom_table(sample_metadata, table_id, data, filename, normalize=False):
    '''
    to be called after obtaining data by calling create_data
    :param sample_metadata:
    :param table_id:
    :param data: dictionary in the form of sample_id:list of Nodes
    :return:
    '''
    otus = []
    sample_id = []  # column index
    for key, value in list(data.items()):
        sample_id.append(key)
        value_name = list(map(lambda x: x.name, value))
        otus = otus + value_name
        for x, node in enumerate(value):
            ab = 1. / (2 ** x)
            ab = ab + np.random.normal(1)
            if ab < 0:
                ab = sys.float_info.epsilon
            node.abundance = ab
    otu_ids = list(set(otus))  # row index unique otus
    print('total {} otus'.format(len(otu_ids)))
    df = pd.DataFrame(columns=sample_id, index=otu_ids)
    for sample in sample_id:
        print(sample)
        sample_data = np.zeros(len(otu_ids))
        for i, otu in enumerate(otu_ids):
            otu_in_sample = list(map(lambda x: x.name, data[sample]))
            if otu in otu_in_sample:
                real_node = _get_node(otu,list(data[sample]))
                sample_data[i] = real_node.abundance
        df[sample] = sample_data  # fill by column
    table = Table(data=df.to_numpy(), observation_ids=otu_ids, sample_ids=sample_id, observation_metadata=None,
                  sample_metadata=None, table_id=table_id)
    normed = table.norm(axis='sample', inplace=False)
    with open(filename, "w") as f:
        if normalize:
            normed.to_tsv(direct_io=f)
        else:
            table.to_tsv(direct_io=f)
    return data

def _get_node(node_name, node_list):
    for node in node_list:
        if node.name == node_name:
            return node

def view_data(data):
    for k, value in data.items():
        print(k)
        for v in list(value):
            print(v.name)
            print(v.abundance)

# def create_metadata(data):

def get_sorted_distance(node, node_list, tree, file):
    '''
    get ranked pairwise distance between the nodes
    :param node_list:
    :return:
    '''
    # node_dist = dict()
    # sorted_nodes = []
    # for n in node_list:
    #     dist = node.get_distance(tree&n)
    #     node_dist[dist] = n
    # for key in sorted(node_dist):
    #     sorted_nodes.append(node_dist[key])
    # node_list = list(map(lambda n: tree&n, node_list))
    real_node = tree&node
    node_list.sort(key=lambda x: real_node.get_distance(tree & x), reverse=False)
    with open(file, 'a') as f:
        f.write("%s\t" % node)  # will be the first one anyway
        f.writelines("%s\t" % n for n in node_list)
        f.write("\n")

def filter_tree(tree, valid_list):
    leaves = tree.get_leaves()
    valid_leaves = list(filter(lambda x: (x.name in valid_list), leaves))
    return valid_leaves

def get_sorted_distance_fast(node, node_list):
    node_list.sort(key=lambda x: node.get_distance(x), reverse=False)
    node_list = list(map(lambda x:x.name, node_list))
    return node.name, node_list

def get_sorted_distance_fast_star(arg):
    return get_sorted_distance_fast(*arg)

def write_dict_to_file(result, filename):
    for pair in result:
        node, node_list = pair
        with open(filename, 'a') as f:
            f.write("%s\t" % node)  # will be the first one anyway
            f.writelines("%s\t" % n for n in node_list)
            f.write("\n")

def get_dict(file):
    '''
    read in the file line by line as a dict (node:list of nodes)
    :param file:
    :return:
    '''
    node_dict = dict()
    with open(file,'r') as f:
        for line in f.readlines():
            line = line.strip()
            nodes = line.split('\t')
            node_dict[nodes[0]] = nodes[1:]
    return node_dict

######################test##########################
def _fun(node, d, file):
    return get_sorted_distance(node, d['node_list'], d['tree'], file)


def test_parallel():
    #doesn't really run unless copy into main
    (tree, otu_tax_dict, node_list) = setup()
    node_list = node_list[:10]
    manager = Manager()
    d = manager.dict({'tree': tree, 'node_list': node_list})

    start1 = time.time()
    for each_node in node_list:
        result = get_sorted_distance_fast(each_node, node_list)
    end1 = time.time()
    print(end1 - start1)

    start2 = time.time()
    pool = Pool(10)
    result = pool.map(get_sorted_distance_fast_star, zip(valid[:5], repeat(valid)))
    pool.close()
    end2 = time.time()
    print(end2 - start2)


def test_merge():
    '''
    test merge method of pandas. checking of get_otu_taxid_dict function
    :return:
    '''
    # differ in number of rows and order of common key
    otu_acc_df = pd.read_table('test_merge_otu.txt', header=None, dtype=str)
    acc_taxid_df = pd.read_table('test_merge_taxid.txt', header=None, usecols=[0, 1], dtype=str)
    merged_df = otu_acc_df.merge(acc_taxid_df, left_on=1, right_on=0)
    print(merged_df)


if __name__ == '__main__':
    distance_dict = get_dict('data/sorted_distance_complete.txt')
    (tree, otu_tax_dict) = setup()
    # (data, taxdata) = create_data(1, tree, 0.8, 100, otu_tax_dict, 120)
    # create_biom_table('meta', 'test_1sample', data, 'test_1sample.tsv', False)
    # only if needs profiles to be created at the same time
    # for key, value in list(taxdata.items()):
    #    filename = "{}{}".format(key, '.profile')
    #    pf.create_profile(value, 'profiles', filename)
    #os.chdir('node_class_test')
    data = create_data_simple(200, 25, 500, distance_dict, otu_tax_dict)
    updated_data = create_biom_table('meta', 'test_200_out_500', data, 'test_200_out_500.tsv', False)


