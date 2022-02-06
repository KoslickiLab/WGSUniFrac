#!/bash/bin
set -e

IN_DIR=$1 #input directory containing experiments 
OTU_FILE="otu_table.tsv"
TREE_PATH="reproducibles/data/gg_13_5_otus/trees"
ROOT="$(dirname $0)" #this WGSUniFrac directory

if ! test -f './reproducibles/data/99_otu_tree.qza'; then
	#get the tree. do only once.
	qiime tools import --input-path "$ROOT"/"$TREE_PATH"/99_otus.tree --output-path reproducibles/data/99_otus_tree.qza --type 'Phylogeny[Rooted]'
fi
#get distance matrix file in every experiment dir
for exp_dir in $IN_DIR/*
do
	echo $exp_dir
	biom convert -i "$exp_dir"/"$OTU_FILE" -o "$exp_dir"/"$otu_table.biom" --table-type="OTU table" --to-hdf5
	qiime tools import --input-path "$exp_dir"/"$otu_table.biom" --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path "$exp_dir"/feature-table.qza
        qiime diversity beta-phylogenetic --i-table "$exp_dir"/feature-table.qza --i-phylogeny "$ROOT"/reproducibles/data/99_otus_tree.qza --p-metric weighted_unifrac --o-distance-matrix "$exp_dir"/dmatrix.qza
        qiime tools export --input-path "$exp_dir"/dmatrix.qza --output-path "$exp_dir"
done

