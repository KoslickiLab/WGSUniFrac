#!/usr/bin/bash

tmp_file="index.html"
url="https://gg-sg-web.s3-us-west-2.amazonaws.com/"

wget $url -O $tmp_file
files=$(grep -oE "<Key>[^<]*gg_13_5/[^<]*</Key>" $tmp_file | sed 's/<Key>\([^<]*\)<\/Key>/\1/' | grep ".gz$")

for f in $files; do  wget -P ./reproducibles/data "$url$f"; done
tar -xvzf ./reproducibles/data/*tar*
gunzip ./reproducibles/data/*.gz
wget -P ./reproducibles/data/ https://zenodo.org/record/5885631/files/sorted_distance_complete.txt
wget -P ./reproducibles/data/ https://zenodo.org/record/5933187/files/otu_with_valid_taxid.txt
rm $tmp_file


