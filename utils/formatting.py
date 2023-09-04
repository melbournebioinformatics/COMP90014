
import pandas as pd


def pretty_print_matrix(the_dict):
    
    # check that this is a dict
    try:
        assert(isinstance(the_dict, dict))
    except AssertionError:
        print('Error: input is not a dictionary.')
        return None

    names = list(the_dict.keys())
    
    # check that all keys are present and nested structure is correct 
    try:
        for name1 in names:
            for name2 in names:
                the_dict[name1][name2]
    except KeyError:
        print('Error: nested structure of dictionary is incorrect.')
        return None

    # move data to pandas dataframe for pretty printing
    df = pd.DataFrame(index=names, columns=names)
    for name1 in names:
        for name2 in names:
            df[name1][name2] = the_dict[name1][name2]

    # print it
    print(df)