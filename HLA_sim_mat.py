from Bio import pairwise2 as p
from Bio.SubsMat.MatrixInfo import blosum100 as blosum100
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import sys
import multiprocessing as mp

#HLA allele amino acid sequences gotten from HLA_retriever.py
try:
    with open("{}/databases/HLA_alleles.txt".format(os.path.dirname(__file__))) as json_file:
        HLA_dict = json.load(json_file)
except FileNotFoundError:
    print("Please ensure 'HLA_alleles.txt' is in the databases folder.")
    raise

mhcI, mhcI_ca, mhcII, mhcII_ca = ([] for i in range(4)) #dataframe column names used later for easier selection of columns
matrix = blosum100 #substitution matrix used for calculation of similarity between two MHC alleles

def sim_calc(types, i, j):
    """
    calculation of similarity between two alleles, as seen in Alignment.py, will be used in similarity matrix creation of
    MHC I and MHC II alleles
    """
    a1 = HLA_dict[types[i]]
    a2 = HLA_dict[types[j]]
    alignments = p.align.globalds(a1, a2, matrix, -10, -0.5, one_alignment_only=True)[0]
    length = alignments[4]
    max = 17 * length
    sim = 1 - ((max - alignments[2]) / max)
    return (i, j, sim)


def fill(arr, type_list):
    """
    Used to fill in values for similarity matrix creation using similarity value calculated by the function sim_calc
    """
    result = []
    pool = mp.Pool(mp.cpu_count())
    for i in range(len(type_list)):
        for j in range(len(type_list)):
            res = pool.apply_async(sim_calc, (type_list, i, j))
            result.append(res)
    for res in result:
        i, j, sim = res.get()
        arr[i][j] = sim
    pool.terminate()

def run():
    """
    Main method that does dataframe manipulation to build similarity matrices to be later used for spectral embedding and GMM clustering.
    """
    # loads the excel spreadsheet file containing HLA types of donors and replaces missing with nan
    try:
        types = pd.read_excel("{}/spreadsheets/types.xlsx".format(os.path.dirname(__file__))).replace("n.t.", np.nan)
    except FileNotFoundError:
        print("Please ensure 'types.xlsx' is in the spreadsheets folder.")
        raise
    
    #gets the columns to group the MHC I and MHC II alleles separately
    for column in list(types):
        if column in ("A", "B", "C"):
            mhcI.append(column)
        if column in ("DRB1", "DQB1"):
            mhcII.append(column)

    #splits the alleles up for each HLA type in each donor
    for t in mhcI+mhcII:
        one = '{}.1'.format(t)
        two = '{}.2'.format(t)
        if t in mhcI:
            mhcI_ca.extend([one, two])
        else:
            mhcII_ca.extend([one, two])
        types[one], types[two] = types[t].str.split(', ', 1).str
        types[[one,two]] = types[[one,two]].apply(lambda x: t + x)

    #creates two dataframes for similarity matrix calculation, one for MHC I alleles and the other for MHC II alleles
    types = types.drop(columns = mhcI+mhcII)
    types2 = types.copy()

    types = types[mhcI_ca] #MHC I alleles
    types2 = types2[mhcII_ca] #MHC II alleles

    #melts the columns of the dataframes to one column with only unique alleles
    types = types.melt(value_vars= mhcI_ca, value_name="type").dropna()
    types2 = types2.melt(value_vars= mhcII_ca, value_name="type").dropna()

    #changes the dataframes into lists for creation of similarity matrices
    types = types["type"].astype('category').cat.categories.tolist()
    types2 = types2["type"].astype('category').cat.categories.tolist()

    #forms the empty similarity matrices
    sim_matrix = np.zeros((len(types), len(types)))
    sim_matrix2 = np.zeros((len(types2), len(types2)))

    #changes directory to databases, where similarity matrix will be saved for use in HLA_clusterer. The heatmaps of these
    #will also be saved in this directory
    os.chdir("{}/databases/sim_matrix".format(os.path.dirname(__file__)))

    #fills up the similarity matrices using allele similarity calculated by global alignment of the two sequences(scores calculated by blosum100 matrix)
    print("Forming similarity matrix for MHC alleles. The computation will take up at least a few minutes, please wait.")
    fill(sim_matrix, types)
    print("MHC I alleles similarity matrix done...now computing MHC II alleles similarity matrix.")
    fill(sim_matrix2, types2)
    print("MHC II alleles similarity matrix done.")

    #exports to excel spreadsheets
    sim_matrix = pd.DataFrame(data = sim_matrix, index = types, columns = types)
    sim_matrix2 = pd.DataFrame(data = sim_matrix2, index = types2, columns = types2)

    writer = pd.ExcelWriter('sim_matrix.xlsx')

    sim_matrix.to_excel(writer,'MHC I')
    sim_matrix2.to_excel(writer,'MHC II')

    writer.save()

    #creation of MHC allele similarity heatmaps
    datamhcI = pd.read_excel("sim_matrix.xlsx", index_col=0)
    datamhcII = pd.read_excel("sim_matrix.xlsx", 1, index_col=0)

    os.chdir("{}/output/cluster_data".format(os.path.dirname(__file__)))

    #creates and saves heatmaps of similarity matrices into ~/output/cluster_data
    sns.set()
    ax = sns.heatmap(datamhcI.values)
    plt.savefig("mhcI_heatmap.png")
    plt.close()
    ax2 = sns.heatmap(datamhcII.values)
    plt.savefig("mhcII_heatmap.png")
    plt.close()

if __name__ == '__main__':
    run()
    print("Similarity matrices creation successful. They have been saved in the database/sim_matrix folder. Heatmaps of the similarity matrices can also be found in this folder.")
    print("Please proceed to HLA_clusterer.py for clustering of HLA alleles.")
    input('Press ENTER to exit')