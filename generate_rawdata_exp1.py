import taxunifrac as tu

if __name__ == '__main__':
    (otu_tax_dict, distance_dict) = tu.setup(True)
    ranges = [200, 500, 1000, 5000, 10000, 15000, 20000]
    similarity = [-1, 30000, 20000, 10000, 5000, 1000, 900, 800]
    for r in ranges:
        for i in range(100):
            tu.run_one(distance_dict, otu_tax_dict, num_org=200, num_sample=25, range=r, similarity=-1, run=i)
    for sim in similarity:
        for i in range(100):
            tu.run_one(dist_dict=distance_dict, tax_dict=otu_tax_dict,
                    num_org=200, num_sample=25, range=500, similarity=sim, run=i)