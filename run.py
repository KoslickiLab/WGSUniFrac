import os
from WGSsimulation import get_silhouette_from_qiime
from utils import create_data_simple, create_biom_table, setup, get_dict
from ProfilingTools import create_profile

def run_one(dist_dict, tax_dict, num_org, num_sample, range, similarity, run):
    #create directory
    cur_dir = os.getcwd()
    dir_name = "org" + str(num_org) + "sample" + str(num_sample) + "range" \
               + str(range) + "dist" + str(similarity) + "run" + str(run)
    os.mkdir(dir_name)
    os.chdir(dir_name)
    data = create_data_simple(num_org=num_org, num_sample=num_sample, range=range, similarity=similarity,
                       dist_dict=dist_dict, tax_dict=tax_dict)
    updated_data = create_biom_table('meta', dir_name, data, 'otu_table.tsv', True)
    os.mkdir("profiles")
    os.chdir("profiles")
    for key, value in list(updated_data.items()):
        filename = "{}{}".format(key, '.profile')
        create_profile(value, 'profiles', filename)
    os.chdir(cur_dir)


if __name__ == '__main__':
    distance_dict = get_dict('data/sorted_distance_complete.txt')
    (tree, otu_tax_dict) = setup()
    os.chdir("data")
    ranges = [200, 500, 1000, 5000, 10000, 15000, 20000]
    similarity = [-1, 20000, 10000, 5000, 1000, 900, 800]
    for i in range(5):
        for range in ranges:
            run_one(distance_dict, otu_tax_dict, num_org=200, num_sample=25, range=range, similarity=-1, run=i)
        for sim in similarity:
            run_one(dist_dict=distance_dict, tax_dict=otu_tax_dict,
                    num_org=200, num_sample=25, range=500, similarity=sim, run=i)