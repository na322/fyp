import re
import json
import os
HLA, HLA_g = ({} for i in range(2))

def extract_GP(**string):
    '''
    extract_GP is the method used to extract allele names and the groups they belong to from the text files
    available for download on hla.alleles.org/alleles. Each line is read. and if the current line being read
    contains a G or P group allele, each allele in that group will be stored into dictionary HLA_g, with the
    group(s) it is a member of as its value.

    Alternatively, the alleles and the groups they belong to could have also been extracted by the table shown on the
    webpages themselves, but that approach was not taken.
    '''
    for code, filename in string.items():
        with open(filename) as file:
            alleles = file.read().splitlines()
            for line in alleles:
                if line.endswith(code):
                    hla_class = re.search(r'\w+\*', line).group()
                    group = re.search(r'(\d+:)+\d+G|(\d+:)+\d+P', line).group()
                    groupname = '{}{}'.format(hla_class, group)
                    members = line.replace(group, '').replace(hla_class, '').strip(';').split('/')
                    for member in members:
                        member = '{}{}'.format(hla_class, member)
                        if member in HLA_g:
                            HLA_g[member].append(groupname)
                        else:
                            HLA_g[member] = [groupname]
    with open('HLA_groups.txt', 'w') as outfile:
        json.dump(HLA_g, outfile)


def hla_extract(dict, **string):
    '''
    From the data files Dr. Michael Hallensleben provided, hla_extract extracts the peptide sequences for each group present
    in the HLA_g dictionary, and also for each allele in each HLA class that does not belong in either P or G groups. Only
    HLA-A, HLA-B, HLA-C, HLA_DPB1 and HLA_DQB1 are processed since those are the only ones relevant to the research.

    Then, the HLA dictionary is converted into .json files for use later on.
    '''
    for hla, filename in string.items():
        with open(filename) as file:
            X = file.read().splitlines()
        for line in X:
            x = line.split()
            x[0] = '{}{}'.format(hla+'*', x[0])
            if x[0] in HLA_g:
                for group in HLA_g[x[0]]:
                    dict[group] = x[1]
            else:
                dict[x[0]] = x[1]
        outputfilename = 'HLA_alleles.txt'
        with open(outputfilename, 'w') as outfile:
            json.dump(dict, outfile)


if __name__ == '__main__':
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir('{}/databases'.format(dname))
    extract_GP(G = 'hla_nom_g.txt', P = 'hla_nom_p.txt')
    hla_extract(HLA, A = 'HLA-A.txt', B = 'HLA-B.txt', C = 'HLA-C.txt', DRB1 = 'HLA-DRB1.txt',DQB1 = 'HLA-DQB1.txt')
    print("HLA database prepared. Please proceed to HLA_typecheck to check if your types spreadsheet have alleles were not converted to G/P grouped names, if any.")
    print("You can also proceed to Alignment.py if looking to just compare two different alleles and see their amino acid mismatches.")
    input('Press ENTER to exit')


