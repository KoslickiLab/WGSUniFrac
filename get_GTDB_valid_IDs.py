import taxunifrac as tu

def main():
    tu.get_GTDB_valid_IDs('data/GTDB/bac120_all_taxons_taxids.txt', 'data/GTDB/bac120_taxonomy_no_numeric.tsv')

if __name__ == "__main__":
        main()