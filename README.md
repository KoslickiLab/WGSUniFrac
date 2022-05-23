1. WGSUniFrac Reproducibles

### Pre-requisites

* conda or miniconda
* python 3.6 (or higher)

### Installation

```
git clone https://github.com/KoslickiLab/WGSUniFrac.git
cd WGSUniFrac
```

# WGSUniFrac User Manual
WGSUniFrac takes a directory containing WGS profiles as input and produces a pairwise WGSUniFrac distance matrix file as an output. Run `python get_pairwise_unifrac.py -h` to see the options.

An example is given below.

```
python get_pairwise_unifrac.py -d profile_dir -a -1 -s "my_WGSUniFrac.csv"
```







