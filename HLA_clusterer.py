from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from sklearn.mixture import GaussianMixture as GMM
from sklearn import manifold as mn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json

#opens up HLA database
try:
    with open("{}/databases/HLA_alleles.txt".format(os.path.dirname(__file__))) as json_file:
        HLA_dict = json.load(json_file)
except FileNotFoundError:
    print("Please ensure 'HLA_alleles.txt' is in the databases folder.")
    raise

os.chdir("{}/databases/sim_matrix".format(os.path.dirname(__file__)))

#opens up similarity matrices formed from HLA_sim_mat.py
try:
    datamhcI = pd.read_excel("sim_matrix.xlsx", index_col=0)
    datamhcII = pd.read_excel("sim_matrix.xlsx", 1, index_col=0)
except FileNotFoundError:
    print("Please ensure 'sim_matrix.xlsx' is in the databases folder.")
    raise

def run():
    '''
    The main method that embeds the similarity into a lower dimensional subspace, using laplacian eigenmaps(spectral_embedding). The parameters
    were initialised with 15 dimensions for both MHC I, 16 for MHC II laplacian eigenmaps, as these settings were found by trial and error
    to have a low BIC value with a decent number of clusters that seem to separate HLA-types by locus well(HLA-A clusters will only have
    HLA-A clusters) while keeping a high number of dimensions(retains information).

    The commented out code was used to create plots of the first 2 or 3 eigenvectors with the largest eigenvalues.
    '''
    #removes the labels/names of alleles from the similarity matrices for spectral embedding/laplacian eigenmaps
    data = datamhcI.values
    data2 = datamhcII.values
    

    #laplacian eigenmapping of MHC I similarity matrix, followed by clustering using Gaussian Mixture Models.
    data = mn.spectral_embedding(data, 5, drop_first=True, random_state=11)
    labels = GMM(8, random_state=22).fit(data).predict(data)

    #produces linegraphs showing change in BIC, which suggests the number of clusters for a better model
    #saved in output/cluster_data
    os.chdir("{}/output/cluster_data".format(os.path.dirname(__file__)))

    n_components = np.arange(1, 25)
    models = [GMM(n, covariance_type= 'full', random_state=22).fit(data) for n in n_components]
    plt.plot(n_components, [m.bic(data) for m in models], label='BIC')
    plt.plot(n_components, [m.aic(data) for m in models], label='AIC')
    plt.savefig("BIC_graph_MHC1.png")
    plt.close()

    #earlier steps repeated for MHC II similarity matrix
    data2 = mn.spectral_embedding(data2, 5, drop_first=True, random_state=11)
    labels2 = GMM(7, random_state=22).fit(data2).predict(data2)

    models = [GMM(n, covariance_type= 'full', random_state=22).fit(data2) for n in n_components]
    plt.plot(n_components, [m.bic(data2) for m in models], label='BIC')
    plt.plot(n_components, [m.aic(data2) for m in models], label='AIC')
    plt.savefig("BIC_graph_MHC2.png")
    plt.close()

    #scatter plots showing the 1st dimension against the 2nd, 3rd and 4th dimensions
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(data[:,0], data[:,1], c = labels)
    plt.savefig("1n2 eigenvector MHC 1.png")
    plt.scatter(data[:,0], data[:,2], c = labels)
    plt.savefig("1n3 eigenvector MHC 1.png")
    plt.scatter(data[:,0], data[:,3], c = labels)
    plt.savefig("1n4 eigenvector MHC 1.png")
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(data2[:,0], data2[:,1], c = labels2)
    plt.savefig("1n2 eigenvector MHC 2.png")
    plt.scatter(data2[:,0], data2[:,2], c = labels2)
    plt.savefig("1n3 eigenvector MHC 2.png")
    plt.scatter(data2[:,0], data2[:,3], c = labels2)
    plt.savefig("1n4 eigenvector MHC 2.png")
    plt.close()
    print(data2[:,1])

    #saving of cluster memberships into .json files which can be located in /databases
    dict = pd.DataFrame(data=labels, index=datamhcI.index, columns=['sim_class']).to_dict()["sim_class"]
    dict2 = pd.DataFrame(data=labels2, index=datamhcII.index, columns=['sim_class']).to_dict()["sim_class"]

    os.chdir("{}/databases".format(os.path.dirname(__file__)))

    filename = 'MHCI_clusters.txt'
    with open(filename, 'w') as outfile:
        json.dump(dict, outfile)

    filename = 'MHCII_clusters.txt'
    with open(filename, 'w') as outfile:
        json.dump(dict2, outfile)

    #method used to show cluster members in a more readable format, which is then saved in the ~/output/cluster_data folder.
    def show_cluster_members(dict):
        n = max(dict.values()) + 1
        clusters = {}
        for i in range(n):
            members = ""
            for j in iter(dict.items()):
                if j[1] == i:
                    members += "{} ".format(j[0])
            members = members[:-1]
            clusters[i] = members
        for i in range(n):
            file.write("{}\n".format(list(clusters.items())[i]))
            print(list(clusters.items())[i])


    os.chdir("{}/output/cluster_data".format(os.path.dirname(__file__)))
    print("Clusters for MHC I alleles.\n")
    file = open("MHCI_clust_members.txt", "w")
    show_cluster_members(dict)
    file.close()

    print("Clusters for MHC II alleles.\n")
    file = open("MHCII_clust_members.txt", "w")
    show_cluster_members(dict2)
    file.close()

if __name__ == '__main__':
    run()
    print("Clustering done. Please proceed to CA.py to see how the clustered alleles correlate with response "
          "patterns of specific CMV antigens.")
    print("Keeping this window open will help in reading the correspondence analysis charts later, "
          "though the cluster memberships can be seen in the output/cluster_data folder.")
    input('Press ENTER to exit')
