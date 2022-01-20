import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Parse taxonomy file")
    parser.add_argument('-f', '--file', type=str, help='Taxonomy file.')
    parser.add_argument('-t', '--tax', type=int, help='Generate tax file. 1: yes. 0: no')
    parser.add_argument('-o', '--outfile', type=str, help="file name to write taxonomy file to. optional")
    args = parser.parse_args()
    file = args.file
    outfile = args.outfile
    if args.tax == 1:
        generate_tax = True
    else:
        generate_tax = False

    tu.parse_taxonomy_file(file, generate_tax, outfile)


if __name__ == "__main__":
    main()