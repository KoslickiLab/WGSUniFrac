import os
from WGSsimulation import get_silhouette_from_qiime
from utils import create_data_simple, create_biom_table
from ProfilingTools import create_profile

def run_one(dist_dict, tax_dict, num_org, num_sample, range, similarity):
    #create directory
    cur_dir = os.getcwd()
    dir_name = "org" + str(num_org) + "sample" + str(num_sample) + "range" + str(range) + "dist" + str(similarity)
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
    os.chdir("data")
