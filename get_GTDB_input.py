import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Create abundance files for grinder amplicon according to user input")
    parser.add_argument('-s', '--sample_num', type=int, help='Number of samples for each environment.')
    parser.add_argument('-o', '--organism_num', type=int, help='Number of organisms per sample.')
    parser.add_argument('-od', '--out_dir', type=str, help="Output directory.")
    parser.add_argument('-e', '--env_num', type=int, help="Number of environments.", default=2)
    parser.add_argument('-r', '--range', type=int, help="How clustered nodes are in a sample. Choose "
                                                        "a number between number of organisms to 7000",
                        default=500)
    parser.add_argument('-d', '--distance', type=int, help='A rough estimate of how distant'
                                                           'the environments are. Choose a number between '
                                                           'organism number and 6000.',
                        default=6000)
    parser.add_argument('-dmatrix', '--dist_matrix', type=str, help="File containing distance matrix")
    parser.add_argument('-tb', '--table_id', type=str, help="File containing distance matrix")

    args = parser.parse_args()
    out_dir = args.out_dir
    table_id = args.table_id
    outfile = out_dir + '/' + table_id + '.tsv'
    #env_num = args.env_num
    sample_num = args.sample_num
    org_num = args.organism_num
    rnge = args.range
    dist = args.distance
    dmatrix_file = args.dist_matrix

    distance_dict = tu.get_dist_dict(dmatrix_file)
    print(len(distance_dict))
    data_dict = tu.create_GTDB_data(distance_dict, dist, rnge, org_num, sample_num)
    tu.create_GTDB_biom_table(table_id, data_dict, outfile, normalize=False)

if __name__ == "__main__":
    main()