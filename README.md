[1. WGSUniFrac Reproducibles](#reproducible)

[2. WGSUniFrac User Manual](#user_manual)



### Pre-requisites

* conda or miniconda
* python 3.6 (or higher)

### Installation

```
git clone https://github.com/KoslickiLab/WGSUniFrac.git
cd WGSUniFrac
```



## 1. WGSUniFrac Reproducibles <a name="reproducible"></a>

The section headings will follow the headings in the Results section of the manuscript. In the following steps it is assumed that the current directory is *WGSUniFrac*. It is also assumed that Qiime is installed in a conda environment named *qiime2*. The installations of the following packages were not included in install_dependencies.sh as they differ depending on the operating system of the user. 

### Additional pre-requisites (install as needed)

* Grinder: for simulating amplicon and WGS reads. 
  * https://github.com/zyxue/biogrinder
* mOTUs: for profiling WGS reads into taxonomic profiles.
  * https://github.com/motu-tool/mOTUs

* wget
  * `brew install wget` with Homebrew or install via other methods
* Qiime
  * https://docs.qiime2.org/2020.8/install/native/#install-qiime-2-within-a-conda-environment

1. Create a conda environment.

```
conda create -n re-wgsunifrac python=3.6
```

2. Activate environment and install dependencies.

```
conda activate re-wgsunifrac
bash install_dependencies.sh
```

3. Install other packages listed under  **Additional pre-requisites** above as needed.
4. Create the subdirectory.

```
mkdir reproducibles
```



### 1.1 On taxonomic data converted from phylogenetic data

1. Acquire data.

```
bash get_data_exp1.sh
```

2. Generate raw data (16S biom tables and WGS profiles).

```
mkdir reproducibles/exp1
python generate_rawdata_exp1.py -o reproducibles/exp1 -dm reproducibles/data/sorted_distance_complete.txt -mf reproducibles/data/otu_with_valid_taxid.txt
```

3. Get 16s pairwise UniFrac matrices using Qiime.

```
conda activate qiime2
bash get_16s_distance_matrix_with_qiime_exp1.sh reproducibles/exp1/testRange
bash get_16s_distance_matrix_with_qiime_exp1.sh reproducibles/exp1/testDissimilarity
```

4. Get combined dataframe.

```
conda activate re-wgsunifrac
python get_combined_dataframe_exp1.py -d reproducibles/exp1/testRange -a -1 -s 'reproducibles/exp1/combined_dataframe_range.txt'
python get_combined_dataframe_exp1.py -d reproducibles/exp1/testDissimilarity -a -1 -s 'reproducibles/exp1/combined_dataframe_dissimilarity.txt'
```

5. Get boxplots from the dataframe.

```
python generate_plot_exp1.py -f reproducibles/exp1/combined_dataframe_range.txt -x range -y silhouette -s 'reproducibles/exp1/range_vs_silhouette_boxplot.png' 
python generate_plot_exp1.py -f reproducibles/exp1/combined_dataframe_dissimilarity.txt -x dissimilarity -y silhouette -s 'reproducibles/exp1/dissimilarity_vs_silhouette_boxplot.png'
```



## 2. WGSUniFrac User Manual <a name="user_manual"></a>
WGSUniFrac takes a directory containing WGS profiles as input and produces a pairwise WGSUniFrac distance matrix file as an output. Run `python get_pairwise_unifrac.py -h` to see the options.

An example is given below.

```
python get_pairwise_unifrac.py -d profile_dir -a -1 -s "my_WGSUniFrac.csv"
```







