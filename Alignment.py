from cmd import Cmd
from Bio import pairwise2 as p
import os
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SubsMat.MatrixInfo import blosum100 as matrix

class Alignment(Cmd):
    '''
    A sequence comparator that is run using a command line interface. Was made to compare HLA mismatches
    between renal transplant data which turned out to be unavailable, though the mismatches between a pair of HLA
    alleles can still be shown.

    Later on, during the development of HLA_sim_mat.py, the similarity calculation function was modified and
    used to create similarity matrices.
    '''
    def __init__(self):
        super().__init__()
        self.ready = False
        self.alleles = ["", ""]
        self.sequences = ["", ""]
        self.alignments = ()
        self.outputname = "alignment"
        self.mismatches = []
        self.match = 0
        self.similarity = 0
        self.info = ""
        self.fig = plt.figure()
        self.data = []
        self.labels = []
        plt.close()
        print(os.path.dirname(__file__))
        with open("{}/databases/HLA_alleles.txt".format(os.path.dirname(__file__))) as json_file:
            self.HLA_dict = json.load(json_file)

    def do_setboth(self, alleles):
        """"Sets first and second alleles to be compared."""
        alleles = alleles.split()
        if not(len(alleles)) == 2:
            print("Please input only two alleles to be compared.")
        else:
            for i in range(len(alleles)):
                if self.check_allele(alleles[i]):
                    self.set_allele(alleles[i], i)
                else:
                    break
        print("""Alleles have been aligned and computed. Use the 'compare' command next to access the similarity measure
        between the two alleles, dissimilarity areas and a heatmap visualisation of the alignment.
        Alternatively, use the 'show' command, input 'help' show for more information.""")

    def do_setfirst(self, allele):
        """Sets first allele to be compared."""
        if self.check_allele(allele):
            self.set_allele(allele, 0)

    def do_setsecond(self, allele):
        """Sets second allele to be compared."""
        if self.check_allele(allele):
            self.set_allele(allele, 1)

    def do_compare(self, args):
        """ Performs global alignment on HLA alleles set up. Afterwards, the program compares the two HLA alleles set,
        producing a detailed '.txt' file that shows where mismatches happened, while also showing a computed degree of
        similarity between them. The program will also produce a heatmap that will visualise both regions of similarity
        and dissimilarity between the two HLA alleles. Warmer regions represent dissimilarity, while cooler regions
        represent similarity. This is in reference to the BLOSUM100 substitution matrix. The heatmap will be saved as a
        '.png' file."""
        if not self.ready:
            print("Both alleles have not been set yet. Please set both alleles to start comparison.")
        else:
            self.compare()

    def do_show(self, args):
        """Shows either only the percent similarity between the HLA alleles set, the mismatch information, or the heatmap
        highlighting the mismatches. If all three are needed to be shown, use 'compare' command instead."""
        if not self.ready:
            print("Both alleles have not been set yet. Please set both alleles to start comparison.")
        elif 'similarity' in args:
            print('HLA sequence similarity between {} and {} is {}%'.format(self.alleles[0], self.alleles[1],
                                                                            self.similarity))
        elif 'mismatch' in args:
            self.show_mismatch_info()
        elif 'heatmap' in args:
            self.create_heatmap()
            plt.show()
        else:
            print("Only accepted arguments are 'similarity', 'mismatch' or 'heatmap'.")

    def do_quit(self,args):
        """Exits the program."""
        print("Quitting.")
        raise SystemExit

    def do_save(self,args):
        """Saves mismatch information and heatmap into 'name.txt' and 'name.png', where name is the argument given.
        If name contains illegal characters for a filename, the will be removed.
        Example usage: 'save compare0' which would result in the files to be saved as 'compare0.txt' and 'compare0.png'
        """
        if not self.ready:
            print("There is nothing to save. Please start a comparison.")
        else:
            self.set_output(args)

    def do_overwrite(self, args):
        """Overwrites files with the same filename selected previously."""
        if not self.ready:
            print("There is nothing to overwrite with. Please start a comparison.")
        elif args == "":
            print("Nothing is inputted as the filename. Please input 'save filename', with filename as desired name for files.")
        elif os.path.exists(self.output_filepath(self.remove_chars(args))):
            self.save_files(self.output_filepath(self.remove_chars(args)))
        else:
            print("There are no such files to overwrite. Please use the 'save' command instead.")

    def check_allele(self, name):
        check = name in self.HLA_dict
        if not check:
            print("Please check input '{}', as this allele cannot be found in the HLA Database present in this program".format(name))
        return check

    def set_allele(self, name, index):
        self.alleles[index] = name
        self.sequences[index] = self.HLA_dict[name]
        print("Allele {} has been set to {}".format(index+1, self.alleles[index]))
        self.ready = (len(self.alleles[0]) > 1 and len(self.alleles[1]) > 1)
        if (self.ready):
            self.align()

    def calc_similarity(self):
        self.match = 0
        self.mismatches = []
        length = self.alignments[4]
        max = 13 * length
        for i in range(length):
            if (self.matrix_refer((self.alignments[0][i], self.alignments[1][i]))) > 0:
                self.match += 1
            else:
                self.mismatches.append(i)
        self.similarity = (max - self.alignments[2])/max

    def write_mismatch_info(self):
        self.info = ''
        self.info = '{} vs {}'.format(self.alleles[0], self.alleles[1])
        self.info = '{}\nMismatches occurred at positions...'.format(self.info)
        for mismatch in self.mismatches:
            self.info = '{}\n{}: {} against {}'.format(self.info, mismatch, self.sequences[0][mismatch], self.sequences[1][mismatch])

    def show_mismatch_info(self):
        print(self.info)

    def save_mismatch_info(self):
        file = open('{}.txt'.format(self.outputname), "w")
        file.write(self.info)
        file.write('\nHLA sequence similarity between {} and {} is {}'.format(self.alleles[0], self.alleles[1], self.similarity))
        file.close()
        print("Mismatch information for {} vs {} has been saved as '{}.txt' in the output\\alignments\{} directory."
              .format(self.alleles[0], self.alleles[1], self.outputname, self.outputname))

    def matrix_refer(self, key):
        if '-' in key:
            value = -10
        elif key not in matrix:
            value = matrix[(key[1], key[0])]
        else:
            value = matrix[key]
        return (value)

    def create_heatmap_data(self):
        self.data = np.zeros([1, len(self.alignments[0])])
        for i in range(1):
            for j in range(self.alignments[4]):
                self.data[i, j] = self.matrix_refer((self.alignments[0][j], self.alignments[1][j]))
        self.labels = []
        for i in range(self.alignments[4]):
            char = "{}\n{}".format(self.alignments[0][i], self.alignments[1][i])
            self.labels.append(char)

    def create_heatmap(self):
        sns.set(font_scale=0.7)
        plt.figure(num="{} vs {}".format(self.alleles[0], self.alleles[1]), figsize=(18, 2))
        plt.title('{} vs {}'.format(self.alleles[0], self.alleles[1]))
        sns.heatmap(self.data, square=True, vmin=-10, vmax=10, yticklabels=False, xticklabels=self.labels,
                    cmap="RdYlBu")
        plt.tight_layout()
        self.fig = plt.gcf()

    def save_heatmap(self):
        self.fig.savefig('{}.png'.format(self.outputname))
        # bbox_inches = 'tight'
        print("Heatmap for {} vs {} saved as '{}.png' in the output\\alignments\{} directory."
              .format(self.alleles[0], self.alleles[1], self.outputname, self.outputname))

    def align(self):
        self.alignments = p.align.globalds(self.sequences[0], self.sequences[1], matrix, -10, -0.5, one_alignment_only = True)[0]
        self.calc_similarity()
        self.write_mismatch_info()
        self.create_heatmap_data()
        self.create_heatmap()
        plt.close()

    def compare(self):
        self.show_mismatch_info()
        print('HLA sequence similarity between {} and {} is {}%'.format(self.alleles[0], self.alleles[1],
                                                                               self.similarity))
        self.create_heatmap()
        plt.show()

    def remove_chars(self, filename):
        filename = filename.translate(dict((ord(char), None) for char in '\/*?:"<>|'))
        return filename

    def output_filepath(self, filename):
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        return '{}\output\{}\{}'.format(dname, "alignments", filename)

    def set_output(self, filename):
        filename = self.remove_chars(filename)
        if(len(filename)==0):
            print('This filename will result in an empty filename, due to only having illegal filename characters. Please'
                  ' choose another filename.')
        else:
            self.outputname = filename
            filepath = self.output_filepath(self.outputname)
            if not os.path.exists(filepath):
                os.makedirs(filepath)
                self.save_files(filepath)
            else:
                print("There are already files with the name '{}', input 'overwrite {}' to confirm overwriting these files."
                      .format(self.outputname, self.outputname))

    def save_files(self, filepath):
        os.chdir(filepath)
        self.save_mismatch_info()
        self.save_heatmap()

if __name__ == '__main__':
    prompt = Alignment()
    prompt.prompt = '> '
    prompt.cmdloop("""
    Set HLA alleles to be compared using the commands 'setboth', 'setfirst' or 'setsecond'. As an example, typing
    "setboth A*24:109 B*39:71" without the quotation marks will set the program to compare A*24:109 and B*39:71.

    To only calculate the similarity between the HLA alleles set, input the command 'show similarity'.

    To only see where mismatches has occurred between the HLA alleles set, input the command 'show mismatch'.

    To see a heatmap representation of the mismatches, input the command 'show heatmap'.

    To do all three at once, use the command 'compare'.

    If you wish to save the mismatch information and the heatmap as files, input the command 'save filename', where
    filename would be the name of your desired files. They can be found in the output\\alignments\\filename directory.
    """)