import taxunifrac as tu
import argparse
from ete3 import TreeNode


def main():
    parser = argparse.ArgumentParser(description="Get a pairwise distance matrix for a given file"
                                                 "containing otus and a tree.")
    parser.add_argument('-f', '--node_file', type=str, help='File containing otus in one of the columns.')
    parser.add_argument('-t','--tree', type=str, help='Tree file.')
    parser.add_argument('-c', '--col', type=int, help='Column in the file containing otus.')
    parser.add_argument('-o', '--output', type=str, help='Output file name.')
    args = parser.parse_args()
    file = args.node_file
    tree = TreeNode(args.tree, format=1, quoted_node_names=True)
    nodes = []
    col = args.col
    with open(file,'r') as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            if len(line) > col and line[col][0] != '>':
                nodes.append(line[col])
    real_nodes = list(map(lambda x: tree&x, nodes))
    print('nodes ready')
    with open(args.output, 'a+') as f:
        for node in real_nodes:
            (this_node, sorted_nodes) = tu._get_sorted_distance(node, real_nodes)
            f.write("%s\t" % this_node)  # will be the first one anyway
            f.writelines("%s\t" % n for n in sorted_nodes)
            f.write("\n")


if __name__ == '__main__':
    main()