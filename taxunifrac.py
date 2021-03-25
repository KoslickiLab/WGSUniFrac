import os
import logging
from collections import defaultdict
from ete3 import NCBITaxa
import copy
import numpy as np
import dendropy
from biom import Table
import pandas as pd
import random
import seaborn as sns
from scipy.stats import halfnorm
from skbio import DistanceMatrix #to install: pip install scikit-bio
from sklearn.metrics import silhouette_score
import re
import itertools as it
ncbi = NCBITaxa()



# classes
class Prediction:
    def __init__(self):
        pass

    @property
    def rank(self):
        return self.__rank

    @property
    def taxid(self):
        return self.__taxid

    @property
    def percentage(self):
        return self.__percentage

    @property
    def taxpath(self):
        return self.__taxpath

    @property
    def taxpathsn(self):
        return self.__taxpathsn

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @taxid.setter
    def taxid(self, taxid):
        self.__taxid = taxid

    @percentage.setter
    def percentage(self, percentage):
        self.__percentage = percentage

    @taxpath.setter
    def taxpath(self, taxpath):
        self.__taxpath = taxpath

    @taxpathsn.setter
    def taxpathsn(self, taxpathsn):
        self.__taxpathsn = taxpathsn

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {property.split("_")[3]: value for property, value in self.__dict__.items()}

    def get_metadata(self):
        return {'rank': self.__rank, 'taxpath': self.__taxpath, 'taxpathsn': self.__taxpathsn}


#profiling
class Profile(object):
    def __init__(self, sample_metadata=None, profile=None, branch_length_fun=lambda x: 1 / x):
        self.sample_metadata = sample_metadata
        self.profile = profile
        self._data = dict()
        # Stick in the root node just to make sure everything is consistent
        self._data["-1"] = dict()
        self._data["-1"]["rank"] = None
        self._data["-1"]["tax_path"] = list()
        self._data["-1"]["tax_path_sn"] = list()
        self._data["-1"]["abundance"] = 0
        self._data["-1"]["descendants"] = list()
        self._header = list()
        self._tax_id_pos = None
        self._rank_pos = None
        self._tax_path_pos = None
        self._tax_path_sn_pos = None
        self._abundance_pos = None
        self._eps = .0000000000000001  # This is to act like zero, ignore any lines with abundance below this quantity
        self._all_keys = ["-1"]
        self._merged_flag = False
        self.root_len = 1  # the length you want between the "root" of "-1" and the superkingdom level (eg. Bacteria)
        self.branch_len_func = branch_length_fun  # Given a node n at depth d in the tree, branch_len_func(d)
        # is how long you want the branch length between n and ancestor(n) to be
        self._data["-1"]["branch_length"] = self.root_len
        self.parse_file()  # TODO: this sets all the branch lengths to 1 currentlclass Profile(object):


    def parse_file(self):
        _data = self._data
        _all_keys = self._all_keys
        _header = self._header
        for k, v in self.sample_metadata.items():
            _header.append('{}:{}'.format(k, v))

        # populate all the correct keys
        for prediction in self.profile:
            _all_keys.append(prediction.taxid.strip())

        # crawl over all profiles tax_path and create the ancestors and descendants list
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()
            tax_path = prediction.taxpath.strip().split("|")  # this will be a list, join up late
            if tax_id not in _data:
                _data[tax_id] = dict()
            else:
                raise Exception(f"Improperly formatted profile: row starting with {tax_id} shows up more than once")
            _data[tax_id]["tax_path"] = tax_path

            # populate abundance
            _data[tax_id]["abundance"] = prediction.percentage

            # populate tax path sn
            if not (prediction.taxpathsn is None):  # might not be present
                _data[tax_id]["tax_path_sn"] = prediction.taxpathsn.strip().split(
                    "|")  # this will be a list, join up later

            # populate the rank
            _data[tax_id]["rank"] = prediction.rank.strip()

            # populate the branch length
            _data[tax_id]["branch_length"] = self.tax_path_to_branch_len(tax_path, self.branch_len_func, self.root_len)

            # Find the ancestors
            if len(tax_path) <= 1:  # note, due to the format, we will never run into the case tax_path == []
                _data[tax_id]["ancestor"] = "-1"  # no ancestor, it's a root
            else:  # go from the bottom up, looking for an ancestor that is an acceptable key
                ancestor = "-1"  # this is the default
                tax_path_rev = tax_path[::-1]
                for potential_ancestor in tax_path_rev:
                    if potential_ancestor != tax_id and potential_ancestor in _all_keys:
                        ancestor = potential_ancestor
                        break  # you found the ancestor, so can quit looking
                _data[tax_id]["ancestor"] = ancestor

            # Create a placeholder descendant key initialized to [], just so each tax_id has a descendant key associated to it
            if "descendants" not in _data[tax_id]:  # if this tax_id doesn't have a descendant list,
                _data[tax_id]["descendants"] = list()  # initialize to empty list

        self._add_descendants()
        self._delete_missing()  # make sure there aren't any missing internal nodes

    def _add_descendants(self):
        """
        Idea here is to look at all the ancestors of each key, and make the key the descendant of that ancestor
        Returns
        -------
        None: modifies Profile in place
        """
        _data = self._data
        _all_keys = self._all_keys
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()  # the tax ID we are looking at
            ancestor = _data[tax_id]['ancestor']  # the tax ID's ancestor
            if tax_id not in _data[ancestor]['descendants']:
                _data[ancestor]['descendants'].append(
                    tax_id)  # so make the tax ID we're looking at the descendant of the ancestor

    def _delete_missing(self):
        """
        Deletes from the descendants all those taxids that aren't keys in the profile (i.e. there is no line that starts with that taxID)
        Returns
        -------
        none: modifies Profile in place
        """
        for key in self._data:
            clean_descendants = []
            for descendant in self._data[key]["descendants"]:
                if descendant in self._all_keys:  # if it's one of the taxids that the line starts with, add it
                    clean_descendants.append(descendant)
                else:
                    pass  # don't include the taxids that aren't actually in the final tax tree
            self._data[key]["descendants"] = clean_descendants
        return

    def write_file(self, out_file_name=None):
        if out_file_name is None:
            raise Exception
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        fid = open(out_file_name, 'w')
        # Write the header
        for head in self._header:
            fid.write("%s\n" % head)

        # Loop over length of tax_path and write data
        # always make the output tax_id, rank, tax_path, tax_path_sn, abundance in that order
        for path_length in range(1, tax_path_lengths + 1):
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length and _data[key]["abundance"] > self._eps:
                    line_data = _data[key]
                    fid.write("%s\t" % key)
                    if self._rank_pos is not None:
                        fid.write("%s\t" % line_data["rank"])
                    fid.write("%s\t" % "|".join(line_data["tax_path"]))
                    if self._tax_path_sn_pos is not None:
                        fid.write("%s\t" % "|".join(line_data["tax_path_sn"]))
                    fid.write("%f\n" % line_data["abundance"])
        fid.close()
        return

    def threshold(self, threshold=None):
        if threshold is None:
            raise Exception
        _data = self._data
        keys = _data.keys()
        for key in keys:
            if _data[key]["abundance"] < threshold:
                _data[key]["abundance"] = 0
        return

    def _subtract_down(self):
        # helper function to push all the weights up by subtracting
        # NOTE: when subtracting, need to start at root and go down
        # NOTE: when adding, need to start at leaves and go up
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        for path_length in range(1, tax_path_lengths):  # eg tax_path_lengths = 5, use 1,2,3,4 since we stop at leaves
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length:
                    descendants = _data[key]["descendants"]  # get all descendants
                    for descendant in descendants:
                        _data[key]["abundance"] -= _data[descendant]["abundance"]  # subtract the descendants abundance

    def _add_up(self):
        # helper function to push all the weights up by subtracting
        # NOTE: when subtracting, need to start at root and go down
        # NOTE: when adding, need to start at leaves and go up
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        for path_length in range(tax_path_lengths, 1,
                                 -1):  # eg tax_path_lengths = 5, use 5,4,3,2, since we stop at roots
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length:
                    ancestor = _data[key]["ancestor"]
                    if ancestor in _data:  # don't do anything if this is a/the root node
                        _data[ancestor]["abundance"] += _data[key]["abundance"]  # add the descendants abundance

    def normalize(self):
        # Need to really push it up while subtracting, then normalize, then push up wile adding
        # self._push_up(operation="subtract")
        self._subtract_down()
        _data = self._data
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] /= total_abundance
                _data[key]["abundance"] *= 100  # make back into a percentage
        # self._push_up(operation="add")
        self._add_up()
        return

    def merge(self, other):
        # Warning: not checking for taxonomic consistency
        if not isinstance(other, Profile):
            print("Only works with other Profiles")
            raise Exception
        if self._merged_flag is False:
            self._header.insert(0, "# This is a merged file, ignore files in headers below")
            self._merged_flag = True
        _data = self._data
        _other_data = other._data
        other_keys = _other_data.keys()
        for key in other_keys:
            if key in _data:
                _data[key]["abundance"] += _other_data[key]["abundance"]  # if already in there, add abundances
            else:
                _data[key] = copy.copy(_other_data[key])  # otherwise use the whole thing

    @staticmethod
    def tax_path_to_branch_len(tax_path, func, root_len=1):
        """
        This function modifies the branch lengths based on the input tax_path.
        intent is: ["2", "", "123", "456"] would result in a branch length of func(4)
        Parameters
        ----------
        tax_path : a list of strings (tax ID's)
        func : a function whose argument is the depth in the tree of a tax ID, and whose output is the branch length
               from the tax ID to its ancestor.
        root_len : how long you want the root of the tree "-1" to be to the descendants (eg. "-1" -> "Bacteria")
        Returns
        -------
        float
        """
        # eg. "-1" -> "Bacteria" should have a branch length of root_len
        if not tax_path:
            return root_len
        else:
            depth_in_tree = len(tax_path)  # this takes into account that the tax_path doesn't include the root of "-1"
            return func(depth_in_tree)

    def make_unifrac_input_and_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(
            zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, nodes_to_index, P, Q

    def make_unifrac_input_no_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(
            zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                # _data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
                pass
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        # for key in keys:
        #    if total_abundance > 0:
        #        _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                # _other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
                pass
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        # for key in keys:
        #    if total_abundance > 0:
        #        _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, nodes_to_index, P / 100., Q / 100.

class Node:
    def __init__(self, name, abundance=0, tax=''):
        self.name = name
        self.tax = tax
        self.abundance = abundance

# loading file
def open_profile_from_tsv(file_path, normalize):
    header = {}
    column_name_to_index = {}
    profile = []
    samples_list = []
    predictions_dict = {}
    reading_data = False
    got_column_indices = False

    with open(file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            line = line.rstrip('\n')

            # parse header with column indices
            if line.startswith("@@"):
                for index, column_name in enumerate(line[2:].split('\t')):
                    column_name_to_index[column_name] = index
                index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(
                    column_name_to_index)
                got_column_indices = True
                reading_data = False
                continue

            # parse header with metadata
            if line.startswith("@"):
                # if last line contained sample data and new header starts, store profile for sample
                if reading_data:
                    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
                        if len(profile) > 0:
                            samples_list.append((header['SAMPLEID'], header, profile))
                            profile = []
                            predictions_dict = {}
                    else:
                        logging.getLogger('opal').critical(
                            "Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(
                                file_path))
                        raise RuntimeError
                    header = {}
                reading_data = False
                got_column_indices = False
                key, value = line[1:].split(':', 1)
                header[key.upper()] = value.strip()
                continue

            if not got_column_indices:
                logging.getLogger('opal').critical(
                    "Header line starting with @@ in file {} is missing or at wrong position.\n".format(file_path))
                raise RuntimeError

            reading_data = True
            row_data = line.split('\t')

            taxid = row_data[index_taxid]
            # if there is already a prediction for taxon, only sum abundance
            if taxid in predictions_dict:
                prediction = predictions_dict[taxid]
                prediction.percentage += float(row_data[index_percentage])
            else:
                if float(row_data[index_percentage]) == .0:
                    continue
                prediction = Prediction()
                predictions_dict[taxid] = prediction
                prediction.taxid = row_data[index_taxid]
                prediction.rank = row_data[index_rank]
                prediction.percentage = float(row_data[index_percentage])
                prediction.taxpath = row_data[index_taxpath]
                if isinstance(index_taxpathsn, int):
                    prediction.taxpathsn = row_data[index_taxpathsn]
                else:
                    prediction.taxpathsn = None
                profile.append(prediction)

    # store profile for last sample
    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
        if reading_data and len(profile) > 0:
            samples_list.append((header['SAMPLEID'], header, profile))
    else:
        logging.getLogger('opal').critical(
            "Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(
                file_path))
        raise RuntimeError

    if normalize:
        normalize_samples(samples_list)

    return samples_list

def get_column_indices(column_name_to_index):
    if "TAXID" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("TAXID"))
        raise RuntimeError
    if "RANK" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("RANK"))
        raise RuntimeError
    if "PERCENTAGE" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("PERCENTAGE"))
        raise RuntimeError
    if "TAXPATH" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("TAXPATH"))
        raise RuntimeError
    index_taxid = column_name_to_index["TAXID"]
    index_rank = column_name_to_index["RANK"]
    index_percentage = column_name_to_index["PERCENTAGE"]
    index_taxpath = column_name_to_index["TAXPATH"]
    if "TAXPATHSN" in column_name_to_index:
        index_taxpathsn = column_name_to_index["TAXPATHSN"]
    else:
        index_taxpathsn = None
    return index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def normalize_samples(samples_list):
    for sample in samples_list:
        sample_id, sample_metadata, profile = sample
        sum_per_rank = defaultdict(float)
        for prediction in profile:
            sum_per_rank[prediction.rank] += prediction.percentage
        for prediction in profile:
            if prediction.percentage > 0:
                prediction.percentage = (prediction.percentage / sum_per_rank[prediction.rank]) * 100.0

def create_profile(node_list, outdir, filename):
    '''

    :param node_list: list of Node objects
    :param outdir:
    :param filename:
    :return:
    '''
    id_list = list(map(lambda x:x.tax, node_list))
    filtered_list = filter_id_list(id_list)
    filtered_list = list(map(int, filtered_list))
    print("length after filtering is %s " % len(filtered_list))
    if len(filtered_list) == 0:
        print("nothing is left.")
        return
    outfile = outdir + '/' + filename
    f = open(outfile, "w+")
    f.write("# Taxonomic Profiling Output\n"
            "@SampleID:SAMPLEID\n"
            "@Version:0.9.1\n"
            "@Ranks:superkingdom|phylum|class|order|family|genus|species\n"
            "@TaxonomyID:ncbi-taxonomy_DATE\n"
            "@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE\n")
    rank_list = (["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    rank_dict = dict()
    for rank in rank_list:
        rank_dict[rank] = []
    rank_list.reverse()
    #update abundance information
    print("here are the nodes")
    print(filtered_list)
    for id in filtered_list:
        lin_list = ncbi.get_lineage(id)  # get lineage
        lin_dict = ncbi.get_rank(lin_list)  # create dict id:rank
        lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
        if id != lin_list[-1]: #in case the id is obsolete and translated to another
            id = lin_list[-1]
        cur_node = _get_node_from_taxid(id, node_list)
        cur_abund = cur_node.abundance
        if lin_dict[id] != 'species':
            new_node = Node(name='species'+str(id), tax=lin_dict_reverse['species'], abundance=cur_abund)
            node_list.append(new_node)
            rank_dict["species"].append(new_node)
        else:
            rank_dict["species"].append(cur_node)
        for rank in rank_list[1:]:
            cur_taxid = lin_dict_reverse[rank]
            cur_node = _get_node_from_taxid(cur_taxid, node_list)
            if cur_node is None:
                cur_node = Node(name=rank+str(cur_taxid), tax=cur_taxid)
                node_list.append(cur_node)
            cur_node.abundance = cur_node.abundance+cur_abund
            if cur_node not in rank_dict[rank]:
                rank_dict[rank].append(cur_node)
            else:
                print(cur_node.tax)
                print(list(map(lambda x: x.tax, rank_dict[rank])))
            #cur_abund = cur_node.abundance
    rank_list.reverse()
    for k, v in rank_dict.items():
        print(k)
        print(len(list(v)))
    #print out
    for rank in rank_list:
        rank_pos = rank_list.index(rank)
        for node in rank_dict[rank]:
            lin_list = ncbi.get_lineage(node.tax)
            lin_dict = ncbi.get_rank(lin_list)
            lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
            superkingdom = lin_dict_reverse["superkingdom"]
            taxpath = str(superkingdom)
            namepath = ncbi.get_taxid_translator([superkingdom])[superkingdom]
            for r in rank_list[1:rank_pos+1]:
                taxid = lin_dict_reverse[r]
                taxpath = taxpath + "|" + str(taxid)
                name = ncbi.get_taxid_translator([taxid])[taxid]
                cur_node = _get_node_from_taxid(taxid, node_list)
                namepath = namepath + "|" + name
            f.writelines([str(node.tax), "\t", rank, "\t", taxpath, "\t", namepath, "\t", str(cur_node.abundance)])
            f.write("\n")
    f.close()
    return

def _get_node_from_taxid(taxid, node_list):
    for node in node_list:
        if int(node.tax) == int(taxid):
            return node

def filter_id_list(id_list):
    filtered_list = []
    for id in id_list:
        if check_rank(id):
            filtered_list.append(id)
    print("total %s passed"  %len(filtered_list))
    return filtered_list

def check_rank(id):
    '''
    Check if all ranks are present
    :param id: a taxid to be checked
    :return: nothing if any of the ranks is missing. otherwise return id
    '''
    rank_list = (["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    try:
        lineage = ncbi.get_lineage(id)  # get lineage dict
    except ValueError:
        return False
    ranks = ncbi.get_rank(lineage).values()
    for r in rank_list:
        if r not in ranks:
            print("rank %s not present" % r)
            return False
    return True

#unifrac
def parse_tree_file(tree_str_file, suppress_internal_node_taxa=False, suppress_leaf_node_taxa=False):
    '''
    Parse a newick tree file and return the dictionary of ancestors Tint
    :param tree_str_file:
    :param suppress_internal_node_taxa:
    :param suppress_leaf_node_taxa:
    :return: (Tint, lint, nodes_in_order)
    '''
    dtree = dendropy.Tree.get(path=tree_str_file, schema="newick", suppress_internal_node_taxa=True,
                              store_tree_weights=True, suppress_leaf_node_taxa=False)
    nodes = dtree.nodes()
    i=0
    for node in nodes:
        if node.taxon is None:
            node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
            i = i+1
    full_nodes_in_order = [item for item in dtree.levelorder_node_iter()]
    full_nodes_in_order.reverse()
    nodes_in_order = [item.taxon.label for item in full_nodes_in_order]
    Tint = dict()
    lint = dict()
    nodes_to_index = dict(zip(nodes_in_order,range(len(nodes_in_order))))
    for i in range(len(nodes_in_order)):
        node = full_nodes_in_order[i]
        parent = node.parent_node
        if parent is not None:
            Tint[i] = nodes_to_index[parent.taxon.label]
            if isinstance(node.edge.length, float):
                lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
            else:
                lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = 0.0
    return Tint, lint, nodes_in_order

def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
    '''
    (Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    '''
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return (Z, diffab)

#taxunifrac
def setup(return_dist_dict=False):
    os.chdir('data')
    (T, l, nodes) = parse_tree_file('99_otus.tree')
    otu_tax_dict = filter_against_tree('otu_with_valid_taxid.txt', nodes)
    if return_dist_dict is True:
        distance_dict = get_dist_dict('sorted_distance_complete.txt')
        return (otu_tax_dict, distance_dict)
    else:
        return (otu_tax_dict)

def get_dist_dict(file):
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

def create_data_simple(num_org, num_sample, sample_range, distance_dict, tax_dict, similarity=-1):
    '''
    :param num_org: number of organisms per sample
    :param num_sample: number of samples per environment
    :param sample_range: the spread of sample
    :param similarity: how close are the two nodes chosen
    :param distance_dict:
    :param tax_dict:
    :return:
    '''
    node1 = random.choice(list(distance_dict.keys()))
    node2 = distance_dict[node1][similarity]
    #for test only
    #node2 = random.choice(list(distance_dict.keys()))
    print(node1)
    print(node2)
    data_dict = dict()
    env1_nodes = distance_dict[node1][:sample_range - 1]
    env2_nodes = distance_dict[node2][:sample_range - 1]
    #update the 2 lists above to contain Node object instead
    for i, node in enumerate(env1_nodes):
        #create Nodes, update tax
        new_node = Node(name=node, tax=tax_dict[node])
        if int(new_node.tax) != ncbi.get_lineage(new_node.tax)[-1]:
            new_node.tax = ncbi.get_lineage(new_node.tax)[-1]
        env1_nodes[i] = new_node
    for i, node in enumerate(env2_nodes):
        new_node = Node(name=node, tax=tax_dict[node])
        if int(new_node.tax) != ncbi.get_lineage(new_node.tax)[-1]:
            new_node.tax = ncbi.get_lineage(new_node.tax)[-1]
        env2_nodes[i] = new_node
    #create sample
    if num_org >= sample_range:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(env1_nodes)
            value2 = copy.deepcopy(env2_nodes)
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

def create_biom_table(table_id, data, filename, normalize=False):
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
    otu_ids = list(set(otus))  # row index unique otus
    print('total {} otus'.format(len(otu_ids)))
    df = pd.DataFrame(columns=sample_id, index=otu_ids)
    for key, value in list(data.items()): #key = sample id, value = list of Nodes
        for x, node in enumerate(value, 1):
            ab = 100. / (1.5 ** x)
            ab = ab + halfnorm.rvs()
            df.at[node.name, key] = ab
    df = df.fillna(.0)
    print(df)
    table = Table(data=df.to_numpy(), observation_ids=otu_ids, sample_ids=sample_id, observation_metadata=None,
                  sample_metadata=None, table_id=table_id)
    normed = table.norm(axis='sample', inplace=False)
    for key, value in list(data.items()):
        for node in value:
            node.abundance = normed.get_value_by_ids(node.name, key) * 100
    with open(filename, "w") as f:
        if normalize:
            normed.to_tsv(direct_io=f)
        else:
            table.to_tsv(direct_io=f)
    return data

def run_one(dist_dict, tax_dict, num_org, num_sample, range, similarity, run):
    #create directory
    cur_dir = os.getcwd()
    dir_name = "range" + str(range) + "dist" + str(similarity) + "run" + str(run)
    if os.path.exists(dir_name):
        print("directory exists")
        return
    os.mkdir(dir_name)
    os.chdir(dir_name)
    data = create_data_simple(num_org=num_org, num_sample=num_sample, sample_range=range,
                              similarity=similarity, distance_dict=dist_dict, tax_dict=tax_dict)
    updated_data = create_biom_table(dir_name, data, 'otu_table.tsv', True)
    os.mkdir("profiles")
    for key, value in list(updated_data.items()):
        filename = "{}{}".format(key, '.profile')
        create_profile(value, 'profiles', filename)
    os.chdir(cur_dir)

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
    #create metadata
    metadata = dict()
    for name in sample_lst:
        env = name[3]
        metadata[name] = {'environment': env}
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
        profile1 = Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: 1/x)
        profile2 = Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: 1/x)
        #(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
        (weighted, _) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        dist_matrix[i][j] = dist_matrix[j][i] = weighted
    os.chdir(cur_dir)
    df = pd.DataFrame.from_dict(metadata, orient='index')
    dm = DistanceMatrix(dist_matrix, sample_lst)
    #dist_pc = pcoa(dm)
    #dist_pc.plot(df=df, column="environment", cmap="Set1", title=plot_title)
    #plt.show()
    return sample_lst, dist_matrix, metadata

def get_dataframe(dir):
    col_names = ["range", "similarity", "silhouette", "data_type", "sample_id"]
    file_lst = os.listdir(dir)[1:]
    os.chdir(dir)
    sil_score_16s = []
    sil_score_wgs = []
    Range = []
    similarity = []
    for file in file_lst:
        os.chdir(file) #individual run
        #get 16s score
        rg = int(re.findall("range(.*)dist", file)[0])
        Range.append(rg)
        sim = int(re.findall("dist(.*)run", file)[0])
        if sim == -1:
            sim = 35461
        similarity.append(sim)
        dist_matrix_16s = pd.read_table("distance-matrix.tsv", index_col=0)
        label_16s = list(map(lambda x: x[3], dist_matrix_16s))
        score_16s = silhouette_score(dist_matrix_16s, label_16s, metric="precomputed")
        sil_score_16s.append(score_16s)
        #get wgs score
        sample_lst_wgs, dist_matrix_wgs, metadata = pairwise_unifrac('profiles')
        label_wgs = list(map(lambda x: x[3], sample_lst_wgs))
        score_wgs = silhouette_score(dist_matrix_wgs, label_wgs, metric="precomputed")
        sil_score_wgs.append(score_wgs)
        os.chdir('..')
    df_16s = pd.DataFrame(columns=col_names, index=range(len(file_lst)))
    df_16s["data_type"] = "16s"
    df_16s["sample_id"] = file_lst
    df_16s["range"] = Range
    df_16s["similarity"] = similarity
    df_16s["silhouette"] = sil_score_16s
    df_wgs = pd.DataFrame(columns=col_names, index=range(len(file_lst)))
    df_wgs["datatype"] = "wgs"
    df_wgs["sample_id"] = file_lst
    df_wgs["range"] = Range
    df_wgs["similarity"] = similarity
    df_wgs["silhouette"] = sil_score_wgs
    df_combined = pd.concat([df_16s, df_wgs])
    df_combined.to_csv("combined_df.txt", sep="\t")
    return df_combined

def get_boxplot(df):
    sns.set_theme(style="ticks", palette="pastel")
    sns.boxplot(x='range', y='silhoette', hue="data_type", data=df, palette=["m", "g"])

#tests
def test_create_profile():
    nodes = ['1204725', '1036678', '183756', '867917', '66852', '633148', '638762', '49547', '483214', '186057', '242129', '1293586', '574338', '224719', '218300', '483896', '65421', '90427', '54250', '2303', '582419', '877455', '215773', '638774', '120963', '187420', '768672', '523849', '589924', '604354', '490098', '418010', '35749', '2261', '38024', '145261', '386456', '195522', '985053', '502115', '286152', '766501', '2171', '1077256', '267446', '262501', '2208', '573063', '647113', '334772', '1121009', '224720', '39152', '310064', '588319', '591019', '2309', '2269', '242697', '155321', '426368', '90909', '394295', '113653', '415426', '213231', '71280', '39441', '387631', '310083', '456320', '242129', '579137', '1184251', '267453', '190192', '487687', '183759', '39441', '387957', '498375', '419665', '281435', '183759', '183762', '536044', '49339', '176306', '399550', '694429', '430614', '39152', '523850', '122420', '54262', '190976', '269247', '129848', '218300', '868131', '328406', '660064', '710190', '2226', '227598', '187879', '242697', '647171', '190977', '1151117', '406327', '547558', '176306', '138903', '391623', '38025', '224325', '290067', '183759', '984979', '880724', '634498', '679901', '71998', '253161', '262501', '710191', '267435', '35749', '2162', '765177', '69540', '392018', '59277', '572478', '267439', '638764', '638773', '638771', '176307', '262498', '35749', '66851', '242697', '572546', '404323', '2265', '328406', '119227', '187880', '1069083', '638763', '186057', '190974', '253161', '110163', '487685', '273116', '46540', '693661', '145261', '583356', '253161', '145262', '2161', '2161', '339860', '256826', '170861', '766501', '256826', '165215', '218300', '272557', '47311', '66852', '243898', '183756', '207243', '573064', '253161', '179630', '70601', '2162', '867904', '1343739', '227597', '2226', '1094980', '56636', '269797', '52001', '259564', '2186', '392018', '565033', '523846', '110164', '59277']

    #nodes = list(set(nodes))
    #print(len(nodes))
    real_nodes = []
    for i, org in enumerate(nodes,1):
        real_node = Node(name=str(org), tax=int(org), abundance=100/(2**i))
        real_nodes.append(real_node)
    create_profile(real_nodes, "data", "test_create_profile2.profile")

def test_create_biom_table():
    dist_dict = get_dist_dict('data/sorted_distance_mini.txt')
    (T, l, nodes) = parse_tree_file('data/99_otus.tree')
    otu_tax_dict = filter_against_tree('data/otu_with_valid_taxid.txt', nodes)
    data_dict = create_data_simple(num_org=50, num_sample=5, distance_dict=dist_dict, tax_dict=otu_tax_dict, sample_range=200)

    updated_data = create_biom_table(table_id="test_create_biom_table", data=data_dict, filename="data/test_create_biom_table", normalize=True)
    for k,v in updated_data.items():
        print(k)
        print(list(map(lambda x:x.name, v)))
        print(list(map(lambda x: x.abundance, v)))

