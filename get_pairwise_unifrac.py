import wgsunifrac as wu
import argparse

#Produces pairwise UniFrac distance matrix with user given input

def main():
    parser = argparse.ArgumentParser(description="Produces pairwise WGS UniFrac distance matrix with given user input.")
    parser.add_argument('-d', '--dir', type=str, help="A directory containing WGS profiles")
    parser.add_argument('-a', '--alpha', type=float, default=-1, help="The factor for branch length function. -1, -0.5, 0, 1, 2. Default=-1")
    parser.add_argument('-s', '--save', type=str, help="File name of the saved file. If not given, the matrix will be saved as pairwise_WGSUniFrac_matrix.csv")

    args = parser.parse_args()
    dir = args.dir
    alpha=args.alpha
    save = args.save
    wu.just_pairwise_unifrac(dir, alpha, save)


if __name__ == '__main__':
    main()


    
