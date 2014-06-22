#! /usr/bin/env python
try:
    import sys
    import math
    import random
    import argparse
    import os
    from ete2 import Tree
    from nexus import NexusReader
    from summary import partitionparser
    from PTPLLH import lh_ratio_test, exp_distribution, species_setting, exponential_mixture
    from numpy import shape, add, sum, sqrt, argsort, transpose, newaxis
    from numpy.linalg import eigh
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError:
    print("Please install the matplotlib and other dependent package first.")
    print("If your OS is ubuntu or has apt installed, you can try the following:") 
    print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
    sys.exit()


def principal_coordinates_analysis(distance_matrix):
    """Takes a distance matrix and returns principal coordinate results

    point_matrix: each row is an axis and the columns are points within the axis
    eigvals: correspond to the rows and indicate the amount of the variation
        that that the axis in that row accounts for
    NOT NECESSARILY SORTED
    """
    E_matrix = make_E_matrix(distance_matrix)
    F_matrix = make_F_matrix(E_matrix)
    eigvals, eigvecs = run_eig(F_matrix)
    #drop imaginary component, if we got one
    eigvals = eigvals.real
    eigvecs = eigvecs.real
    point_matrix = get_principal_coordinates(eigvals, eigvecs)
    
    vector_order = list(argsort(eigvals))
    vector_order.reverse()
    
    return point_matrix[vector_order], eigvals[vector_order]

def make_E_matrix(dist_matrix):
    """takes a distance matrix (dissimilarity matrix) and returns an E matrix

    input and output matrices are numpy array objects of type Float

    squares and divides by -2 each element in the matrix
    """
    return (dist_matrix * dist_matrix) / -2.0

def make_F_matrix(E_matrix):
    """takes an E matrix and returns an F matrix

    input is output of make_E_matrix

    for each element in matrix subtract mean of corresponding row and 
    column and add the mean of all elements in the matrix
    """
    num_rows, num_cols = shape(E_matrix)
    #make a vector of the means for each row and column
    #column_means = (add.reduce(E_matrix) / num_rows)
    column_means = (add.reduce(E_matrix) / num_rows)[:,newaxis]
    trans_matrix = transpose(E_matrix)
    row_sums = add.reduce(trans_matrix)
    row_means = row_sums / num_cols
    #calculate the mean of the whole matrix
    matrix_mean = sum(row_sums) / (num_rows * num_cols)
    #adjust each element in the E matrix to make the F matrix

    E_matrix -= row_means
    E_matrix -= column_means
    E_matrix += matrix_mean

    #for i, row in enumerate(E_matrix):
    #    for j, val in enumerate(row):
    #        E_matrix[i,j] = E_matrix[i,j] - row_means[i] - \
    #                column_means[j] + matrix_mean
    return E_matrix

def run_eig(F_matrix):
    """takes an F-matrix and returns eigenvalues and eigenvectors"""

    #use eig to get vector of eigenvalues and matrix of eigenvectors
    #these are already normalized such that
    # vi'vi = 1 where vi' is the transpose of eigenvector i
    eigvals, eigvecs = eigh(F_matrix)
    #NOTE: numpy produces transpose of Numeric!

    return eigvals, eigvecs.transpose()

def get_principal_coordinates(eigvals, eigvecs):
    """converts eigvals and eigvecs to point matrix
    
    normalized eigenvectors with eigenvalues"""

    #get the coordinates of the n points on the jth axis of the Euclidean
    #representation as the elements of (sqrt(eigvalj))eigvecj
    #must take the absolute value of the eigvals since they can be negative
    return eigvecs * sqrt(abs(eigvals))[:,newaxis]


class taxa:
    def __init__(self, name, coords = [0.0, 0.0]):
        self.name = name 
        self.coords = coords
        
    def toString(self):
        return repr(self.coords[0])+" "+repr(self.coords[1])+" "+self.name


class Species:
    def __init__(self, species_id = 0):
        self.species_id = species_id
        self.taxon = []

    def addTaxa(self, taxa):
        self.taxon.append(taxa)
    
    def toString(self):
        output = ""
        for taxa in self.taxon:
            output = output + ","+taxa.toString()
        return output[1:]


class phylomap:
    def __init__(self, largetree):
        self.tree = Tree(largetree)
        self.taxaorder = self.tree.get_leaves()
        self.numtaxa = len(self.taxaorder)
        self.dism = np.zeros((self.numtaxa, self.numtaxa))
        self.name_coords = {}
        self.species_list = []
    
    
    def parse_delimitation(self, fin, fout):
        with open(fin) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i]
                if line.startswith("Species"):
                    nextline = lines[i+1].strip()
                    taxas = nextline.split(",")
                    spe = Species()
                    for taxa_name in taxas:
                        Taxa = taxa(taxa_name, self.name_coords[taxa_name])
                        spe.addTaxa(Taxa)
                    self.species_list.append(spe)
                    
        with open(fout, "w") as f2:
            for species in self.species_list:
                f2.write(species.toString() + "\n")
    
    
    def calculate_distancematrix(self):
        for i in range(self.numtaxa):
            taxai = self.taxaorder[i]
            for j in range(i, self.numtaxa):
                taxaj = self.taxaorder[j]
                self.dism[i][j] = taxai.get_distance(taxaj)
                self.dism[j][i] = self.dism[i][j]
        return self.dism
    
    def pcoa(self):
        ppm, ev = principal_coordinates_analysis(self.calculate_distancematrix())
        for i in range(len(self.taxaorder)):
            tname = self.taxaorder[i].name
            self.name_coords[tname] = [ppm[0][i], ppm[1][i]]
        
        




if __name__ == "__main__":
    pm = phylomap(largetree = "/home/zhangje/GIT/SpeciesCounting/example/example.tre")
    pm.pcoa()
    pm.parse_delimitation(fin = "/home/zhangje/GIT/SpeciesCounting/example/exampleout.PTPMLPartition.txt", fout = "/home/zhangje/GIT/SpeciesCounting/example/coords.txt")
    #print(pm.taxaorder)
    #ppm, ev = principal_coordinates_analysis(pm.calculate_distancematrix())
    #print(ppm)
    #print(ev)

