import taxunifrac as tu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get box plot from a dataframe file.")
    parser.add_argument('-f', '--file', type=str, help="Dataframe file.")
    parser.add_argument('-x', '--x', type=str, help="x axis.")
    parser.add_argument('-y', '--y', type=str, help="y axis.")
    parser.add_argument('-p', '--palette', type=str, nargs="+", help="Palette list. e.g. -p 'm' 'g'")
    parser.add_argument('-s', '--save', type=str, help="If wants to save df for future use, file name to save as")

    args = parser.parse_args()
    file = args.file
    x = args.x
    y = args.y
    pal = args.palette
    save = args.save
    tu.get_plot_from_file(file, x, y, pal, save)


if __name__ == '__main__':
    main()