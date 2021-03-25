import taxunifrac as tu

if __name__ == '__main__':
    df = tu.get_dataframe("data/testeverything")
    tu.get_boxplot(df)