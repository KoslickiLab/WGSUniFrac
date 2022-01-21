#!/usr/bin/bash

tmp_file="index.html"
url="https://gg-sg-web.s3-us-west-2.amazonaws.com/"

wget $url -O $tmp_file
files=$(grep -oE "<Key>[^<]*gg_13_5/[^<]*</Key>" $tmp_file | sed 's/<Key>\([^<]*\)<\/Key>/\1/' | grep ".gz$")

for f in $files; do  wget -P ./reproducibles/experiment1/ "$url$f"; done
tar -xvzf ./reproducibles/experiment1/*tar*
gunzip ./reproducibles/experiment1/*.gz
wget -P ./reproducibles/experiment1/ https://zenodo.org/record/5885631/files/sorted_distance_complete.txt
rm $tmp_file
