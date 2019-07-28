import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json
from prince import CA
from scipy import sparse
cd4, cd8, mhcI, mhcI_ca, mhcII, mhcII_ca, totals, select = ([] for i in range(8)) #initialisation of the antigen selectors

bin = np.array([-1, 0, 0.001, 0.01, 1]) #bins for classifications of response strengths for specific antigens
labels = ["No response", "Weak", "Moderate", "Strong"] #labels for the binning

try:
    with open("{}/databases/MHCI_clusters.txt".format(os.path.dirname(__file__))) as json_file:
        clusters = json.load(json_file)  # cluster dictionary for MHC I
except FileNotFoundError:
    print("Please ensure 'MHCI_clusters.txt' is in the databases folder.")
    raise

try:
    with open("{}/databases/MHCII_clusters.txt".format(os.path.dirname(__file__))) as json_file:
        clusters2 = json.load(json_file) #cluster dictionary for MHC II
except FileNotFoundError:
    print("Please ensure 'MHCII_clusters.txt' is in the databases folder.")
    raise



def run():
    '''
    The main body of code that deals with preprocessing data, such as applying cluster memberships, separating vaccinated
    donors and non-vaccinated donors, calculating total values of antigen responses, and binning antigen responses and total
    antigen responses into classes of strengths.
    The code then calls the function "analyse" to produce Correspondence Analysis results.
    The HLA allele excel spreadsheet needs to be named types.xlsx and the response excel spreadsheet needs to be named response.xlsx.
    Both files need to be placed in the spreadsheets directory.
    '''
    os.chdir("{}/spreadsheets".format(os.path.dirname(__file__)))
    main = pd.merge(pd.read_excel("response.xlsx"), pd.read_excel("types.xlsx"), how='inner', on='donor') \
        .set_index('donor') \
        .replace([-777, -888, -999], [np.nan, 0, 0])  # loads the HLA types and responses of donors and merges them

    # sets column pointers for the dataframe
    for column in list(main):
        if "cd4" in column:
            cd4.append(column)
        elif "cd8" in column:
            cd8.append(column)
        elif column in ("A", "B", "C"):
            mhcI.append(column)
        elif column in ("DRB1", "DQB1"):
            mhcII.append(column)

    # calculates summated responses of the CMV antigens, and also saves the column names for reference
    main['cd4_total'] = main[cd4].sum(axis=1)
    main['cd8_total'] = main[cd8].sum(axis=1)
    totals.extend(['cd4_total', 'cd8_total'])
    cd4.append('cd4_total')
    cd8.append('cd8_total')

    os.chdir("{}/output/CA/stats".format(os.path.dirname(__file__)))

    writer = pd.ExcelWriter('stats.xlsx')

    # shows std, mean, variances and interquartile ranges of the responses
    main[cd4+cd8].describe().transpose().to_excel(writer, 'Variation Statistics')

    # bins the summated responses into 5 groups of equal proportions
    main[totals] = main[totals].apply(
        lambda x: pd.qcut(x, 5, labels=["Very Low", "Low", "Moderate", "High", "Very High"]))
    main[totals] = main[totals].astype(str)
    main = main.replace(0, -1)  # replace 0 with -1 for binning of antigen responses(No response class)

    # bins and classifies each CMV antigen responses into No response, Weak, Moderate, Strong
    antigens = (list(filter(lambda x: x not in totals, cd4 + cd8)))
    main[antigens] = main[antigens].apply(lambda x: pd.cut(x, bin, labels=labels, right=False, include_lowest=True))

    # replaces n.t with np.nan and changes values to Yes/No for vaccinated and older patients
    main = main.replace(["n.t.", -1, 1], [np.nan, "No", "Yes"])

    # shows the frequency of No response, Weak, Moderate, Strong for CMV antigens
    # total binned responses were not shown since they were binned according to equal sizes
    main[antigens].apply(pd.value_counts).transpose().to_excel(writer, 'Binned Responses')
    main[['CMVVASC', 'older']].apply(pd.value_counts).to_excel(writer, 'Vaccinated, Older')

    # splits up the two alleles recorded in each HLA type and then places the alleles into clusters they belong to
    # the cluster dictionaries were made using HLA_clusterer.py
    for t in mhcI + mhcII:
        one = '{}.1'.format(t)
        two = '{}.2'.format(t)
        main[one], main[two] = main[t].str.split(', ', 1).str
        main[[one, two]] = main[[one, two]].apply(
            lambda x: t + x)  # places the name of the HLA type in front of the alleles

        main[[one,two]].astype(str).apply(pd.value_counts).to_excel(writer, '{} counts'.format(t))

        # replaces allele names with cluster membership
        if t in mhcI:
            mhcI_ca.extend([one, two])
            main[[one, two]] = main[[one, two]].replace(clusters).astype(str)
        else:
            mhcII_ca.extend([one, two])
            main[[one, two]] = main[[one, two]].replace(clusters2).astype(str)

        # gets rid of the decimals formed when converting int to string
        main[[one, two]] = main[[one, two]].apply(lambda x: x.str[:-2])

        # separates the clustered alleles from vaccinated donors and non-vaccinated donors
        main.loc[main['CMVVASC'] == 'Yes', one] = "V" + main[one]
        main.loc[main['CMVVASC'] == 'Yes', two] = "V" + main[two]
    # drops the columns with the joined HLA alleles
    main = main.drop(columns=mhcI + mhcII).replace(["n", "Vn"],[np.nan, np.nan])

    writer.save()

    os.chdir("{}/output/CA/results".format(os.path.dirname(__file__)))

    # starts Correspondence Analysis
    for col in cd4 + cd8:
        select = set_select(col, main.copy(), mhcI_ca, mhcII_ca)
        analyse(select, main.copy())

def analyse(select, df):
    '''
    The body of code that iterates through each CMV antigen responses and produces results of Correspondence Analysis. The
    results will be saved as .xlsx files for Indexed Residuals and Normalised Residuals(z-score) to show correlations between
    HLA alleles(grouped by similarity) and antigen responses observed by donors
    and how statistically significant the correlations are between them.
    '''
    select = select #columns selected by set_select function
    df = df #main dataframe containing all information
    test = df[select] #selects relevant data for current iteration of antigen responses
    id = select.pop(0) #antigen name
    print("Performing Correspondence Analysis for {}".format(id))

    #creation of contigency table
    test = test.melt(id_vars=id, value_vars=select, value_name="type").dropna()
    test = pd.crosstab(test.type, test[id], margins=True, margins_name="Total")

    #creation of excel spreadsheet and saved contigency table onto the 1st sheet
    writer = pd.ExcelWriter('{}.xlsx'.format(id))
    test.to_excel(writer, 'Contingency Table')

    #uses CA from prince library to perform correspondence analysis and plot graphs
    ca = CA()
    test = test.drop(["Total"]).drop(["Total"], axis=1)
    ca = ca.fit(test)
    ax = ca.plot_coordinates(test, figsize=(12, 12))
    ax.set_title('Clustered alleles vs {} binned responses'.format(id))

    #graph plotted is saved in the ~/CA/graphs folder
    filepath = "{}/output/CA/results/graphs".format(os.path.dirname(__file__))
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    os.chdir(filepath)

    plt.savefig('{}.png'.format(id))
    plt.close()

    #calculation of Indexed Residuals and z-score(Standardised Residuals), which are then saved into the spreadsheet as the 2nd and 3rd sheets
    test = test.apply(lambda x: x / np.sum(np.sum(test)))
    r = ca.row_masses_.values
    c = ca.col_masses_.values
    S = sparse.diags(r) @ (test - np.outer(r, c)) @ sparse.diags(c)
    S2 = sparse.diags(r ** -0.5) @ (test - np.outer(r, c)) @ sparse.diags(c ** -0.5)
    data = S / np.outer(r, c)
    data2 = S2
    data = pd.DataFrame(data=data, index=test.index, columns=test.columns).to_excel(writer, 'Indexed Residuals')
    data2 = pd.DataFrame(data=data2, index=test.index, columns=test.columns).to_excel(writer, 'z-scores')

    filepath = "{}/output/CA/results/{}".format(os.path.dirname(__file__), id[4:])
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    os.chdir(filepath)
    writer.save()
    print("Results saved as {}.xlsx and {}.png".format(id, id))


def set_select(col, df, mhcI_ca, mhcII_ca):
    '''
    The function used to set the dataframe column pointers for the function "analyse".
    '''
    select = []
    df = df.copy() #main dataframe
    if col in list(df): #iterates through columns of dataframe to select relevant columns for data analysis
        select.append(col)
        if "cd4" in col:
            select.extend(mhcII_ca)
        elif "cd8" in col:
            select.extend(mhcI_ca)
        return select

if __name__ == '__main__':
    run()
    print("Correspodence Analysis successful. Please check the results in the \output\CA directory.")
    print("The stats folder contains a spreadsheet with sheets showing information of the data, "
          "such as frequency of alleles and interquartile ranges of responses. Graphs of all"
          " antigen analysis can be found in the graph folder.")
    input('Press ENTER to exit')












