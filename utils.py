import random
import pandas as pd
import os
import ProfilingTools as pf
import numpy as np
import EMDUnifrac as unifrac
from biom import Table

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

def pick_otu(node, map_dict, num):
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
    # in case num is too big, return everything possible
    if num > len(leaf_names):
        print("%d pickable nodes" %len(leaf_names))
        return leaf_names
    picked = []
    for i in range(num):
        r = random.randint(1, len(leaf_names))
        if leaf_names[r] not in picked:
            picked.append(leaf_names[r])
    return picked

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
    print(df)
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

def filter_against_tree(otu_file, tree_file):
    '''
    remove nodes in otu_file that are not found on the tree. For the remaining nodes, parse into a dictionary
    :param otu_file:
    :param tree_file:
    :return: a valid {otu:taxid} dict
    '''
    (T, l, nodes) = unifrac.parse_tree_file(tree_file)
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
    filter_against_tree('otu_with_valid_taxid.txt', '../gg/gg_13_5_otus/trees/99_otus.tree')
    #test_merge()