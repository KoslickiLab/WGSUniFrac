import taxunifrac as tu
import argparse

#get a combined dataframe consisting of all experiments in a directory, with 16S UniFrac and WGS UniFrac

def main():
    parser = argparse.ArgumentParser(description="Get a combined dataframe from all experiments under a directory, saving in a file.")
    parser.add_argument('-d', '--dir', type=str, help="directory containing experiments")
    parser.add_argument('-a', '--alpha', type=float, help="The factor for branch length function. -1, -0.5, 0, 1, 5")
    parser.add_argument('-s', '--save', type=str, help="File name of the saved file.")

    args = parser.parse_args()
    dir = args.dir
    alpha=args.alpha
    save = args.save
    tu.get_dataframe(dir, alpha, save)


if __name__ == '__main__':
    main()