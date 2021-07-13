import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get pcoa plot of taxonomic profiles in a direcory.")
    parser.add_argument('-dir', '--dir', type=str, help="The directory where profiles are located.")
    parser.add_argument('-t', '--plt_title', type=str, help="Title of the pcoa plot.")
    parser.add_argument('-a', '--alpha', type=float, help="The factor for branch length function. -1, -0.5, 0, 1, 5")

    args = parser.parse_args()
    dir = args.dir
    plt_title = args.plt_title
    alpha = args.alpha

    #sample_lst, dist_matrix, metadata = tu.pairwise_unifrac(dir, plt_title)
    tu.pairwise_unifrac(dir, plt_title, alpha)

if __name__ == "__main__":
    main()