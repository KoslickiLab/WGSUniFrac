import taxunifrac as tu
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Generate the input data for the second part of the GTDB experiments")
    parser.add_argument('-od', '--out_dir', type=str, help="Output directory.")

    args = parser.parse_args()
    out_dir = args.out_dir
    dist_matrix_file = 'data/GTDB/bac120_distance_matrix_exp2.txt'
    distance_dict = tu.get_dist_dict(dist_matrix_file)
    tax_dict = tu.parse_taxonomy_file('data/GTDB/bac120_taxonomy_no_numeric.tsv')
    name_tax_dict = tu.get_name_taxid_dict('data/GTDB/bac120_all_taxons_taxids_valid.txt')
    GTDBid_taxid_dict = tu.get_GTDBid_taxid_dict(tax_dict, name_tax_dict)

    ranges = [200, 400, 600, 800, 1000, 1500, 2000, 2500]
    dissimilarity = [-1, 4900, 4000, 3000, 2000, 1000, 800]

    for r in ranges: #fix dissimilarity to be 4000
        for i in range(100):
            exp_id = "r" + str(r) + "d4000" + "-" + str(i)
            exp_dir = out_dir + '/' + exp_id
            if os.path.exists(exp_dir):
                print("directory exists")
                break
            os.mkdir(exp_dir)
            profile_dir1 = exp_dir + '/GTDB_profiles'
            profile_dir2 = exp_dir + '/NCBI_profiles'
            os.mkdir(profile_dir1)
            os.mkdir(profile_dir2)
            data_dict = tu.node_selection(distance_dict, dissimilarity=4000, sample_range=r, num_org=200, num_sample=25)
            for key, value in list(data_dict.items()):
                print(key)
                filename = str(key) + '.profile'
                tu.create_GTDB_profile(value, profile_dir1, filename, tax_dict, name_tax_dict, GTDBid_taxid_dict)
    for d in dissimilarity:
        for i in range(100):
            exp_id = "r600" + "d" + str(d) + "-" + str(i)
            exp_dir = out_dir + '/' + exp_id
            biom_file = exp_dir + "/distance_matrix.txt"
            if os.path.exists(exp_dir):
                print("directory exists")
                break
            os.mkdir(exp_dir)
            data_dict = tu.create_GTDB_data(distance_dict, d, 600, 200, 25)
            updated_data = tu.create_GTDB_biom_table(exp_id, data_dict, biom_file, normalize=False)
            profile_dir = exp_dir + '/profiles'
            os.mkdir(profile_dir)
            for key, value in list(updated_data.items()):
                filename = str(key) + '.profile'
                tu.create_GTDB_profile(value, profile_dir, filename, tax_dict, name_tax_dict, GTDBid_taxid_dict)


if __name__ == '__main__':
    main()