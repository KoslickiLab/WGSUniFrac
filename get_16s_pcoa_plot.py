import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get pcoa plot given a distance matrix.")
    parser.add_argument('-f', '--file', type=str, help="Distance matrix file.")
    parser.add_argument('-t', '--plt_title', type=str, help="Title of the pcoa plot.")

    args = parser.parse_args()
    file = args.file
    plt_title = args.plt_title
    tu.get_plot_from_exported(file, plt_title)


if __name__ == "__main__":
    main()
