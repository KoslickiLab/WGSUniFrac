import taxunifrac as tu
import argparse
from multiprocessing import Pool
import timeit


def main():
    parser = argparse.ArgumentParser(description="Write out pariwise matrix file given a directory containing profiles. "
                                                 "Performs in parallel")
    parser.add_argument('-dir', '--dir', type=str, help="The directory where profiles are located.")
    parser.add_argument('-t', '--plt_title', type=str, help="Title of the pcoa plot.")
    parser.add_argument('-a', '--alpha', type=float, help="The factor for branch length function. -1, -0.5, 0, 1, 5")
    parser.add_argument('-by', '--by', type=str, help="For real data. choices: bodysites, study", default="bodysites")

    args = parser.parse_args()
    dir = args.dir
    plt_title = args.plt_title
    alpha = args.alpha
    by = args.by
    p = Pool(4)

    df = tu._get_empty_df(dir) #get a dataframe with col and row being profile name
    #print(df)
    input_lst = tu._get_parallization_input(dir, alpha)
    #print(input_lst)
    start=timeit.timeit()
    result = p.starmap(tu._get_unifrac, input_lst, chunksize=100)
    #print(result)
    for tri in result: #(id1, id2, unifrac)
        id1 = tri[0]
        id2 = tri[1]
        unifrac = tri[2]
        df[id1][id2] = df[id2][id1] = unifrac
    print(df)
    df.to_csv(plt_title, sep='\t')
    end = timeit.timeit()
    print(end - start)


if __name__ == "__main__":
    main()