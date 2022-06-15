import os
import logging
from collections import defaultdict
import numpy as np
import pandas as pd
import itertools as it

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
                index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
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

def just_pairwise_unifrac(dir, alpha, save_as):
    '''
    :param dir: directory containing the .profile files
    :param alpha factor for branch length function  x**alpha
    :return: a dataframe of pairwise distance matrix
    '''
    if save_as is None:
        save_as = "pairwise_WGSUniFrac_matrix.csv"
    cur_dir = os.getcwd()
    file_lst = os.listdir(dir)  # list files in the directory
    # print(file_lst)
    os.chdir(dir)
    if '.DS_Store' in file_lst:
        file_lst.remove('.DS_Store')
    sample_lst = [os.path.splitext(profile)[0].split('.')[0] for profile in file_lst] #e.g.env1sam10. i.e.filenames without extension
    #print(sample_lst)
    # enumerate sample_lst, for filling matrix
    id_dict = dict()
    for i, id in enumerate(file_lst):
        id_dict[id] = i
    # initialize matrix
    dim = len(file_lst)
    dist_matrix = np.zeros(shape=(dim, dim))
    count=0
    for pair in it.combinations(file_lst, 2): #all pairwise combinations
        #to keep the running less boring
        count+=1
        if count % 100 == 0:
            print(count, " pairs processed")
        id_1, id_2 = pair[0], pair[1]
        i, j = id_dict[id_1], id_dict[id_2]
        profile_list1 = open_profile_from_tsv(id_1, False)
        profile_list2 = open_profile_from_tsv(id_2, False)
        name1, metadata1, profile1 = profile_list1[0]
        name2, metadata2, profile2 = profile_list2[0]
        profile1 = Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: x ** alpha)
        profile2 = Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: x ** alpha)
        # (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
        (weighted, _) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        dist_matrix[i][j] = dist_matrix[j][i] = weighted
    os.chdir(cur_dir)
    pd.DataFrame(data=dist_matrix, index=sample_lst, columns=sample_lst).to_csv(save_as, sep="\t")
    return dist_matrix, sample_lst
