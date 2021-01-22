#get otu tables and profiles at the same time
python utils.py

#####qiime part
#convert to biom table
biom convert -i gg_test.tsv -o gg_test.biom --table-type="OTU table" --to-hdf5
qiime tools import --input-path gg_test4.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table.qza
qiime diversity beta-phylogenetic --i-table feature-table.qza --i-phylogeny tree.qza --p-metric weighted_unifrac --o-distance-matrix distance_matrix_from_qiime.qza
qiime tools export --input-path distance_matrix_from_qiime.qza --output-path exported_distance_matrix

#get plot from exported qiime data
python WGSsimulation.py
