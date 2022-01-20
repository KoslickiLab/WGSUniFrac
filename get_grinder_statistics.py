import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get statistics from grinder dataframe")
    parser.add_argument('-e', '--env_num', type=int, help="number of environments")
    parser.add_argument('-m', '--metric', type=str, help="clustering metric. e.g. silhouette")
    parser.add_argument('-s', '--stat', type=str, help="statistics. mean or median")
    parser.add_argument('-t', '--dtype', type=str, help="16s or wgs")

    args = parser.parse_args()
    result = tu.get_grinder_statistics(args.env_num, args.metric, args.dtype, args.stat)
    print(result)


if __name__ == "__main__":
    main()
