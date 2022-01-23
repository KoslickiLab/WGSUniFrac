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

The section headings will follow the headings in the Results section of the manuscript. In the following steps it is assumed that the current directory is WGSUniFrac.

### Additional pre-requisites (install as needed)

* Grinder: for simulating amplicon and WGS reads. 
  * https://github.com/zyxue/biogrinder
* mOTUs: for profiling WGS reads into taxonomic profiles.
  * https://github.com/motu-tool/mOTUs

* wget
  * `brew install wget` with Homebrew or install via other methods
* Qiime

1. Create a conda environment.

```
conda create -n re-wgsunifrac python=3.6
```

2. Activate environment and install dependencies.

```
conda activate re-wgsunifrac
bash install_dependencies.sh
```

2. Create the subdirectory.

```
mkdir reproducibles
```



### 1.1 On taxonomic data converted from phylogenetic data

1. Acquire datasets.

```
bash get_data_exp1.sh
```

2. Generate raw data (16S biom tables and WGS profiles)

```
python generate_rawdata_exp1.py
```



## 2. WGSUniFrac User Manual <a name="user_manual"></a>





