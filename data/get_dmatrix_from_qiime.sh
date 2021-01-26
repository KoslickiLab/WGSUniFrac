#!/bin/bash

#input .tsv file
read FILE
echo "${FILE%.*}.biom"
TREE="../99_otu_tree.qza"

biom convert -i $FILE -o "${FILE%.*}.biom" --table-type="OTU table" --to-hdf5
qiime tools import --input-path "${FILE%.*}.biom" --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table.qza
qiime diversity beta-phylogenetic --i-table feature-table.qza --i-phylogeny $TREE --p-metric weighted_unifrac --o-distance-matrix distance_matrix_from_qiime.qza
qiime tools export --input-path distance_matrix_from_qiime.qza --output-path exported_distance_matrix

