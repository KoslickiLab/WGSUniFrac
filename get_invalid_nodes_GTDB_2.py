import taxunifrac as tu

def main():
    invalid_nodes = tu.get_nodes_without_complete_NCBI_lineage()
    with open('data/GTDB/invalid_nodes.txt', 'w+') as f:
        for node in invalid_nodes:
            f.write(node)
            f.write('\n')

if __name__ == '__main__':
    main()