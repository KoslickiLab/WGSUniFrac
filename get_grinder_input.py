import taxunifrac as tu

if __name__ == "__main__":
    node_list = tu.get_nodes_for_simulated_data(50, 'data/otu_taxid_with_both_data.txt')
    print(node_list)
    tu.get_abundance_file(node_list=node_list, file_name='data/abundance1.txt', abund_fun='exp', factor=2)