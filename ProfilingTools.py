import copy
import numpy as np
import sys
from ete3 import NCBITaxa
from utils import *
ncbi = NCBITaxa()


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
        self.parse_file()  # TODO: this sets all the branch lengths to 1 currently

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


def create_profile(node_list, outdir, filename):
    '''

    :param node_list: list of Node objects
    :param outdir:
    :param filename:
    :return:
    '''
    id_list = list(map(lambda x:x.tax, node_list))
    print("length before filtering is %s " %len(id_list))
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
    rank_list.reverse()
    #update abundance information
    print("all your nodes are here: ")
    print(id_list)
    print("plus your filtered nodes")
    print(filtered_list)
    for id in filtered_list:
        print("Im looking for %s" % id)
        lin_list = ncbi.get_lineage(id)  # get lineage
        print("let's see the lineage...")
        print(lin_list)
        lin_dict = ncbi.get_rank(lin_list)  # create dict id:rank
        print("now let's see the dict")
        print(lin_dict)
        lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
        print("how is the reverse dict")
        print(lin_dict_reverse)
        if id != lin_list[-1]:
            id = lin_list[-1]
        cur_node = _get_node_from_taxid(id, node_list)
        cur_abund = cur_node.abundance
        if lin_dict[id] != 'species':
            new_node = Node(name='species'+str(id), tax=lin_dict_reverse['species'], abundance=cur_abund)
            node_list.append(new_node)
        for rank in rank_list[1:]:
            cur_taxid = lin_dict_reverse[rank]
            cur_node = _get_node_from_taxid(cur_taxid, node_list)
            if cur_node is None:
                cur_node = Node(name=rank+str(cur_taxid), tax=cur_taxid)
                node_list.append(cur_node)
            cur_node.abundance = cur_node.abundance+cur_abund
            #cur_abund = cur_node.abundance
    rank_list.reverse()
    #print out
    for id in filtered_list:
        print("Im looking for %s" %id)
        lin_list = ncbi.get_lineage(id)  # get lineage
        print("let's see the lineage...")
        print(lin_list)
        lin_dict = ncbi.get_rank(lin_list)  # create dict id:rank
        print("now let's see the dict")
        print(lin_dict)
        lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
        print("how is the reverse dict")
        print(lin_dict_reverse)
        taxpath = ""
        namepath = ""
        for rank in rank_list:
            taxid = lin_dict_reverse[rank]
            name = ncbi.get_taxid_translator([taxid])[taxid]
            cur_node = _get_node_from_taxid(taxid, node_list)
            if cur_node is None:
                print("%s not found" %taxid)
            print(taxid)
            if len(taxpath) == 0:
                taxpath += str(taxid)
            else:
                taxpath = taxpath + "|" + str(taxid)
            if len(namepath) == 0:
                namepath += name
            else:
                namepath = namepath + "|" + name
            f.writelines([str(taxid), "\t", rank, "\t", taxpath, "\t", namepath, "\t", str(cur_node.abundance)])
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

def parse_taxid_file(filename):
    '''
    parses a file containing taxids into a list of taxids
    :param filename: a file with a taxid in each line
    :return: a list of taxids
    '''
    id_lst = []
    with open(filename) as file:
        for line in file:
            line = line.strip()
            id_lst.append(line)
    return id_lst

def test_normalize():
    import EMDUnifrac as EMDU
    from load_data import open_profile_from_tsv
    import os
    # test files
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/agitated_blackwell_7"
    file_path2 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/goldstandard_low_1.bin"

    # import one test profile
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile1 = profile_list[0]
    profile1 = Profile(sample_metadata=metadata1, profile=profile1)

    # import another test profile
    profile_list = open_profile_from_tsv(file_path2, False)
    name2, metadata2, profile2 = profile_list[0]
    profile2 = Profile(sample_metadata=metadata2, profile=profile2)

    print("Normalized:")
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile2)  # normalized
    print(f"P sum: {np.sum(P)}")
    print(f"Q sum: {np.sum(Q)}")
    (weighted_norm, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    (unweighted_norm, diffab) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    print(f"weighted, normalized: {weighted_norm}")
    print(f"unweighted, normalized: {unweighted_norm}")

    print("No normalized:")
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile1 = profile_list[0]
    profile1 = Profile(sample_metadata=metadata1, profile=profile1)
    # import another test profile
    profile_list = open_profile_from_tsv(file_path2, False)
    name2, metadata2, profile2 = profile_list[0]
    profile2 = Profile(sample_metadata=metadata2, profile=profile2)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_no_normalize(
        profile2)  # not normalized
    print(f"P sum: {np.sum(P)}")
    print(f"Q sum: {np.sum(Q)}")
    P_missing_mass = 1 - np.sum(P)
    (weighted_no_norm, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    (unweighted_no_norm, diffab) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    print(f"weighted, not normalized: {weighted_no_norm}")
    print(f"unweighted, not normalized: {unweighted_no_norm}")
    print(f"weighted, not normalized, but missing mass added to root node: {weighted_no_norm + P_missing_mass}")

    assert unweighted_norm == unweighted_no_norm
    assert weighted_norm != unweighted_no_norm
    return


def test_branch_lengths_all_1():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small1.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    # print(f"Tint_new: {Tint_new}")
    # print(f"lint_new: {lint_new}")
    # print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 0.20, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val

    # test with branch lengths of x
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: x)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    correct_lints = {('1', '3'): 3, ('3', '4'): 2, ('4', '-1'): 1, ('2', '4'): 3, ('0', '5'): 5, ('5', '-1'): 1}
    for key, val in correct_lints.items():
        assert lint_new[key] == correct_lints[key]
    correct_vals = {'0': 0.20, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val


def test_branch_lengths_all_2():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small2.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    # print(f"Tint_new: {Tint_new}")
    # print(f"lint_new: {lint_new}")
    # print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 1 / 9., '1': 5 / 9., '2': 2 / 9., '3': 0.0, '4': 1 / 9., '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val

    # test with branch lengths of x
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: x)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_and_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    correct_lints = {('1', '3'): 3, ('3', '4'): 2, ('4', '-1'): 1, ('2', '4'): 3, ('0', '5'): 5, ('5', '-1'): 1}
    for key, val in correct_lints.items():
        assert lint_new[key] == correct_lints[key]
    correct_vals = {'0': 1 / 9., '1': 5 / 9., '2': 2 / 9., '3': 0.0, '4': 1 / 9., '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val


def test_no_normalize():
    from load_data import open_profile_from_tsv
    import os
    # test file
    file_path1 = os.path.dirname(os.path.abspath(__file__)) + "/../../data/small2.profile"
    profile_list = open_profile_from_tsv(file_path1, False)
    name1, metadata1, profile_fernando = profile_list[0]

    # Test with branch lengths of 1
    profile1 = Profile(sample_metadata=metadata1, profile=profile_fernando, branch_length_fun=lambda x: 1)
    Tint, lint, nodes_in_order, nodes_to_index, P, Q = profile1.make_unifrac_input_no_normalize(profile1)
    index_to_nodes = dict()
    for key, val in nodes_to_index.items():
        index_to_nodes[val] = key
    Tint_new = dict()
    lint_new = dict()
    for key, val in Tint.items():
        Tint_new[index_to_nodes[key]] = index_to_nodes[val]
    for key, val in lint.items():
        lint_new[(index_to_nodes[key[0]], index_to_nodes[key[1]])] = val
    # print(f"Tint_new: {Tint_new}")
    # print(f"lint_new: {lint_new}")
    # print(f"P: {P}")
    assert Tint_new['5'] == '-1'
    assert Tint_new['4'] == '-1'
    assert Tint_new['3'] == '4'
    assert Tint_new['1'] == '3'
    assert Tint_new['2'] == '4'
    assert Tint_new['0'] == '5'
    assert set(Tint_new.keys()) == {'0', '1', '2', '3', '4', '5'}
    for val in lint.values():
        assert val == 1
    correct_vals = {'0': 0.10, '1': 0.50, '2': 0.20, '3': 0.0, '4': 0.10, '5': 0.00}
    for key, val in correct_vals.items():
        assert P[nodes_to_index[key]] == val

def test_create_profile():
    org1 = Node(name="org1", tax=33936, abundance=20)
    org2 = Node(name="org2", tax=499124, abundance=70)
    org3 = Node(name="org3", tax=187074, abundance=10)
    create_profile([org1, org2, org3], "data", "test_create_profile")

def test_create_profile2():
    nodes = ['469382', '797114', '596418', '2238', '573064', '145261', '268735', '29283', '1269828', '1269822', '1071726', '1230459', '329270', '1071739', '588898', '565033', '260463', '65421', '226001', '662480', '399555', '218300', '547559', '1089732', '1427151', '376171', '227597', '660064', '596424', '2242', '1053642', '329270', '869889', '332953', '260463', '585975', '638771', '1073081', '183759', '29282', '797209', '1126238', '489912', '768065', '399555', '755311', '2265', '1230452', '744725', '869888', '2261', '692267', '329270', '1200261', '29282', '485914', '797114', '90909', '1227487', '867904', '575194', '1214227', '797210', '453848', '1269874', '1195987', '575194', '1195984', '391623', '329270', '1269883', '260463', '399555', '183759', '404323', '267435', '425309', '114529', '717756', '1230457', '1198296', '1132509', '392421', '869886', '553467', '51589', '523846', '362891', '655464', '797303', '302484', '575195', '869887', '755311', '254206', '1017351', '268739', '2238', '61858', '523841', '43928', '1017351', '872320', '51589', '38025', '218300', '869891', '1261545', '1224621', '2238', '1269884', '1269850', '267446', '1126245', '370324', '51589', '536044', '328406', '1269835', '36738', '643748', '572546', '1126240', '504937', '1126237', '2226', '797114', '186057', '262078', '340952', '46540', '57705', '627958', '589924', '310771', '370967', '178600', '1089758', '113653', '268735', '751944', '88724', '399555', '1008196', '406551', '286152', '261291', '1053639', '195522', '329275', '88724', '399555', '231715', '160432', '699431', '267368', '328406', '335951', '558529', '1269825', '329272', '329270', '1224623', '396317', '647171', '699433', '29284', '420995', '323742', '217171', '1269863', '399555', '588319', '638774', '481736', '54262', '755307', '425309', '641561', '155321', '1109010', '1220023', '28442', '183759', '1095778', '213231', '981621', '1198324', '436949', '110163', '69493', '596417', '1158307', '406549', '229731', '122420', '260463', '1121009', '371024', '869886']

    real_nodes = []
    for org in nodes:
        real_node = Node(name=str(org), tax=int(org), abundance=0.5)
        real_nodes.append(real_node)
    create_profile(real_nodes, "data", "test_create_profile2.profile")