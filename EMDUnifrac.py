import numpy as np
import dendropy
import matplotlib
import regex as re
import warnings

def parse_tree(tree_str):
    '''
    parse a newick tree string  and return the dictionary of ancestors Tint
    Tint[i] = j => j is the ancestor of i
    lint[(i,j)] = weight of the edge (i,j)
    nodes labeled from the leaves up
    :param tree_str:
    :return: (Tint, lint, nodes_in_order)
    '''
    dtree = dendropy.Tree.get(data=tree_str, schema="newick", suppress_internal_node_taxa=False, store_tree_weights=True)
    nodes = dtree.nodes()
    i=0
    for node in nodes:
        if node.taxon == None:
            node.taxon = dendropy.datamodel.taxonmodel.Taxon(label="temp"+str(i))
            i = i+1
    full_nodes_in_order = [item for item in dtree.levelorder_node_iter()] #i in path from root to j only if i>j
    full_nodes_in_order.reverse()
    nodes_in_order = [item.taxon.label for item in full_nodes_in_order]
    Tint = dict()
    lint = dict()
    nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
    for i in range(len(nodes_in_order)):
        node = full_nodes_in_order[i]
        parent = node.parent_node
        if parent != None:
            Tint[i] = nodes_to_index[parent.taxon.label]
            lint[nodes_to_index[node.taxon.label], nodes_to_index[parent.taxon.label]] = node.edge.length
    return (Tint, lint, nodes_in_order)

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

def parse_envs(envs, nodes_in_order):
    '''
    envs_prob_dict[samples[i]] is a probability vector on the basis nodes_in_order denoting for sample i
    :param envs:
    :param nodes_in_order:
    :return: (envs_prob_dict, samples)
    '''
    nodes_in_order_dict = dict(zip(nodes_in_order,range(len(nodes_in_order))))
    for node in envs.keys():
        if node not in nodes_in_order_dict:
            print("Warning: environments contain taxa" + node + " not present in given taxonomic tree. Ignoring")
    envs_prob_dict = dict()
    for i in range(len(nodes_in_order)):
        node = nodes_in_order[i]
        if node in envs:
            samples = envs[node].keys()
            for sample in samples:
                if sample not in envs_prob_dict:
                    envs_prob_dict[sample] = np.zeros(len(nodes_in_order))
                    envs_prob_dict[sample][i] = envs[node][sample]
                else:
                    envs_prob_dict[sample][i] = envs[node][sample]
    #normaalize samples
    samples = envs_prob_dict.keys()
    for sample in samples:
        if envs_prob_dict[sample].sum() == 0:
            warnings.warn("Warning: the sample %s has non-zero counts, do not use for Unifrac calculations" %sample)
        envs_prob_dict[sample] = envs_prob_dict[sample]/envs_prob_dict[sample].sum()
    return (envs_prob_dict, samples)


def create_env(sample_file):
    '''
    :param sample_file: a file containding ids and samples
    :return: an env_dict in the form of { id: {sample:count} }
    '''
    env_dict = dict()
    with open(sample_file) as fp:
        line = fp.readline()
        while line:
            list = line.split()
            key = list.pop(0) #get key
            env_dict[key] = dict()
            for str in list:
                m = re.match(r"(\d+).(\w+)_(\d+)", str)
                sample = m.group(2)
                if sample in env_dict[key]:
                    env_dict[key][sample] += 1
                else:
                    env_dict[key][sample] = 1
            line = fp.readline()
    fp.close()
    return env_dict

def EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q):
    '''
    :param Tint:
    :param lint:
    :param nodes_in_order:
    :param P: envs_prob_dict[sample[i]]
    :param Q: envs_prob_doct[sample[j]]
    :return: (Z, F, diffab) Z = weighted Unifrac distance F = flow: a dict with keys of the form (i, j) where F[(i,j)] = num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node nodes_in_order[i] to the node nodes_in_order[j]. diffab = a dictionary with tuple keys using elements of Tint and values diffabb[(i,j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)]
    '''
    num_nodes = len(nodes_in_order)
    F = dict()
    G = dict()
    diffab = dict()
    Z = 0
    w = np.zeros(num_nodes)
    pos = dict()
    neg = dict()
    for i in range(num_nodes):
        pos[i] = set([])
        neg[i] = set([])
    for i in range(num_nodes):
        print("i=", i)
        if P[i] > 0 and Q[i] > 0:
            F[(i, i)] = np.minimum(P[i], Q[i])
        G[(i,i)] = P[i] - Q[i]
        print("G[(i,i)] = ", G[(i,i)])
        if P[i] > Q[i]:
            pos[i].add(i)
        elif P[i] < Q[i]:
            neg[i].add(i)
        posremove = set()
        negremove = set()
        for j in pos[i]:
            print("j=", j)
            print("pos=", pos)
            for k in neg[i]:
                print("neg=", neg)
                print("k=", k)
                if (j not in posremove) and (k not in negremove):
                    val = np.minimum(G[i,j], -G[(i,k)])
                    print("val = ", val)
                    if val > 0:
                        F[(j,k)] = np.minimum(G[(i,j)], -G[(i,k)])
                        G[(i,j)] = G[(i,j)] - val
                        G[(i,k)] = G[(i,k)] + val
                        Z = Z + (w[j] + w[k])*val
                        print("Z=", Z)
                    if G[(i,j)] == 0:
                        posremove.add(j)
                    if G[(i,k)] == 0:
                        negremove.add(k)
        pos[i].difference_update(posremove)
        neg[i].difference_update(negremove)
        if i < num_nodes-1:
            for j in pos[i].union(neg[i]):
                if (Tint[i],j) in G:
                    print("Tint[i],j in G. ","j=", j)
                    G[(Tint[i],j)] = G[(Tint[i],j)] + G[(i,j)]
                    diffab[(i, Tint[i])] = diffab[(i, Tint[i])] + G[(i,j)]
                    print("diffab=", diffab)
                    print("G=", G)
                else:
                    print("Tint[i],j is not in G, j=", j)
                    G[(Tint[i],j)] = G[(i,j)]
                    diffab[(i,Tint[i])] = G[(i,j)]
                    print("diffab=", diffab)
                    print("G=", G)
                print("w, before ", w)
                w[j] = w[j] + lint[i,Tint[i]]
                print("w, after ", w)
            if (i, Tint[i]) in diffab:
                print("i, Tint[i] in diffab ", "Tint[i] =", Tint[i])
                diffab[(i, Tint[i])] = lint[i, Tint[i]]*diffab[(i,Tint[i])]
                print("diffab=", diffab)
            pos[Tint[i]] |= pos[i]
            neg[Tint[i]] |= neg[i]
    return (Z, F, diffab)

def EMDUniFrac_weighted_no_flow(Tint, lint, nodes_in_order, P, Q):
    num_nodes = len(nodes_in_order)
    Z = 0 #unifrac distance
    diffab = dict()
    partial_sums = P-Q
    for i in range(num_nodes-1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]*val]
        Z += lint[i, Tint[i]*abs(val)]
    return (Z, diffab)

def get_branch_length_function(function_str):
    try:
        return eval(function_str)
    except SyntaxError as exception:
        logging.getLogger('opal').warning('Invalid function provided with -b, --branch_length_function: {}. lambda x: 1/x will be used.'.format(exception.msg))
        return eval('lambda x: 1/float(x)')

#(Tint, lint, nodes_in_order) = parse_tree('((D:0.1,C:0.1)F:0.2,(B:0.1,A:0.1)E:0.2)G;')
#P = np.array([0, 1.0/2, 1.0/2, 0, 0, 0, 0])
#Q = np.array([1.0/3, 0, 0, 1.0/3, 0, 1.0/3, 0])
#print(Tint)
#(Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q)
#print(Z)
#print