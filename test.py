import EMDUnifrac as unifrac
import pickle

#setup
(Tint, lint, nodes_in_order) = unifrac.parse_tree_file('../data/tax_slv_lsu_138.1.tre')
#env_dict = unifrac.create_env('../study_232_101019-092408/processed_data/227_sortmerna_picked_otus/289_seqs_otus.txt')

print(nodes_in_order)

with open('nodes.txt', 'w') as f:
    print(nodes_in_order, file=f)

