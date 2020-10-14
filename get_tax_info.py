from ete3 import NCBITaxa
from ete3 import PhyloTree
ncbi = NCBITaxa()

lin28183 = ncbi.get_lineage(28183)

tree_28183 = ncbi.get_topology(lin28183)

#tree = ncbi.get_topology([1, 131567, 2, 203691, 203692, 1643688, 170, 171, 28183, 28078, 27339, 887, 27642, 873, 896, 815, 29195, 29432, 29426, 29502, 29355])
print(tree_28183.get_ascii(attributes=["sci_name", "rank", "taxid"]))