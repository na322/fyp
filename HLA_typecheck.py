import pandas as pd
import numpy as np
import os
import json
mhcI, mhcII = ([] for i in range(2)) #initialisation of the HLA type selectors

def run():
    """
    Checks the excel spreadsheet containing HLA types for correct format. When an allele name cannot be found in the HLA database,
    suggestions will be shown to correct the allele name.
    """
    try:
        types = pd.read_excel("{}/spreadsheets/types.xlsx".format(os.path.dirname(__file__))).replace("n.t.", np.nan)
    except FileNotFoundError:
        print("Please ensure 'types.xlsx' is in the spreadsheets folder.")
        raise

    try:
        with open("{}/databases/HLA_groups.txt".format(os.path.dirname(__file__))) as json_file:
            groups = json.load(json_file)  # G/P group dictionary for grouped alleles
    except FileNotFoundError:
        print("Please ensure 'HLA_groups.txt' is in the databases folder.")
        raise

    try:
        with open("{}/databases/HLA_alleles.txt".format(os.path.dirname(__file__))) as json_file:
            HLA_dict = json.load(json_file)
    except FileNotFoundError:
        print("Please ensure 'HLA_alleles.txt' is in the databases folder.")
        raise

    for column in list(types):
        if column in ("A", "B", "C"):
            mhcI.append(column)
        if column in ("DRB1", "DQB1"):
            mhcII.append(column)

    grouped_members = list(groups)
    names = list(HLA_dict)

    # splits the alleles up for each HLA type in each donor
    for t in mhcI + mhcII:
        one = '{}.1'.format(t)
        two = '{}.2'.format(t)
        types[one], types[two] = types[t].str.split(', ', 1).str
        types[[one, two]] = types[[one, two]].apply(lambda x: t + x)
        types_list = []
        t = types[one].tolist()
        t2 = types[two].tolist()
        types_list.extend([t, t2])
        for n in range(2):
            for i in range(len(types_list[n])):
                test = types_list[n][i]
                if not pd.isnull(test):
                    if test not in names: #not in HLA types database
                        if test in grouped_members:
                            print("\nFor donor {}, in row {}, for the allele in position {}, {} is a grouped allele".format(i+1, i+2, n+1, types_list[n][i]))
                            print("Choice of replacement is/are {}".format(groups[test]))
                        elif test not in grouped_members: #not in grouped HLA types database
                            print("\nIn row {}, donor {}, for the allele in position {}, {} cannot be found in database of HLA alleles".format(i+1, i+2, n+1, types_list[n][i]))

if __name__ == '__main__':
    run()
    print("Can now proceed to HLA_sim_mat.py after corrections have been done, if any were required.")
    input('Press ENTER to exit')