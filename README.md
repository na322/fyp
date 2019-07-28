These are the python scripts used during my final year project.

CA.py is the main script used for visualising correlations between clustered HLA types and immune responses towards CMV, producing Correspondence Analysis graphs and tables for each antigen of CMV.

HLA_clusterer.py is the script that clusters HLA types into different groups based on the laplacian eigenmaps produced by HLA_sim_mat.py, using GMM clustering.

HLA_sim_mat.py globally aligns each pair of HLA alleles' amino acid sequence found in the spreadsheets to calculates a degree of similarity between them, constructing similarity matrices for MHC class 1 and MHC class 2. High CPU usage as it spawns multiple processes to compute the matrices in parallel.

HLA_typecheck.py checks allele names in the spreadsheets for any errors to avoid key errors when looking up the HLA dictionary. If an allele that is not in the HLA dictionary is found, suggestions will be provided to rename them.

HLA_retriver.py parses the HLA '.txt' files to extract the amino acid sequences and allele names for the formation of the main HLA dictionary used in the other scripts.

Alignment.py is a program with a command line interface that enables users to compare HLA alleles and see differences in their amino acid sequences.

Samples of results for one antigen and cluster memberships can be found in the samples folder.
