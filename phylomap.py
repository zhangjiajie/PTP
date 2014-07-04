#! /usr/bin/env python
try:
    import sys
    import math
    import random
    import argparse
    import os
    import shutil
    import subprocess
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

def distance(c1, c2):
    dis = 0.0
    for i in range(len(c1)):
        dis = dis + ((c1[i]-c2[i])*(c1[i]-c2[i]))
    return math.sqrt(dis)

def firstDstep(d_tree, d_mds, d_coord_c, d_coord_u):
    d=0.0
    d=((d_tree-d_mds)/(d_tree*d_mds))*(d_coord_u-d_coord_c)
    return d

def secondDstep(d_tree, d_mds, d_coord_c, d_coord_u):
    d=0.0
    d=((d_tree-d_mds)-((d_coord_u-d_coord_c)*(d_coord_u-d_coord_c)/d_mds)*(1+((d_tree-d_mds)/d_mds)))/(d_tree*d_mds)
    return d

def raxmlTreeParser(fin):
    f = open(fin)
    lines = f.readlines()
    f.close()
    trees = []
    for line in lines:
        line = line.strip()
        if not line == "":
            trees.append(line[line.index("("):])
    return trees


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
    def __init__(self, largetree, ptp_result, fout, ftype = "raxml", maxiters = 10000, seed = 2222, debug = False):
        self.ptp = ptp_result
        self.fout = fout
        if ftype == "nexus":
            self.nexus = NexusReader(largetree)
            self.nexus.blocks['trees'].detranslate()
            self.trees = self.nexus.trees.trees
        else:
            self.trees = raxmlTreeParser(largetree)
        self.tree = Tree(self.trees[0])
        self.taxaorder = self.tree.get_leaves()
        self.numtaxa = len(self.taxaorder)
        self.dism = np.zeros((self.numtaxa, self.numtaxa))
        self.name_coords = {}#all taxa name and coords in 2d, the innernodes need to be updated while the exteral nodes and others not in the tree are fixed
        self.species_list = []
        self.represent_taxon = []
        self.all_node_names = []#all node names from the representative tree
        self.inner_node_names = [] 
        self.innernode_dis = {}#inner nodename to 3 connected nodes distances from the tree
        self.innernode_connecting_nodes = {}
        self.innernode_dis_tree_matrix = {} #pair-wise distance matrix for innernodes from tree
        self.innernode_dis_mds_matrix = {}#pair-wise distance matrix for innernodes from mds
        self.rand_nr = random.Random()
        self.rand_nr.seed(seed)
        self.sum_species_tree_length = -0.1
        self.mf = 0.3
        self.min_improvement = 0.0001
        self.maxiters = maxiters
        self.debug = debug
        random.seed(seed)
    
    
    def parse_delimitation(self):
        """step2: must run pcoa first"""
        fin = self.ptp
        fout = self.fout + ".taxa.txt"
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
                self.represent_taxon.append(species.taxon[0].name)
                f2.write(species.toString() + "\n")
    
    
    def _calculate_distancematrix(self):
        for i in range(self.numtaxa):
            taxai = self.taxaorder[i]
            for j in range(i, self.numtaxa):
                taxaj = self.taxaorder[j]
                self.dism[i][j] = taxai.get_distance(taxaj)
                self.dism[j][i] = self.dism[i][j]
        return self.dism
    
    
    def pcoa(self):
        """step1: do pcoa and store first two dimentions in hashmap """
        ppm, ev = principal_coordinates_analysis(self._calculate_distancematrix())
        for i in range(len(self.taxaorder)):
            tname = self.taxaorder[i].name
            self.name_coords[tname] = [ppm[0][i], ppm[1][i]]
        sumev = 0
        for evv in ev:
            sumev += abs(evv)
        return ev[0]/sumev, ev[1]/sumev
    
    
    def extract_species_tree(self):
        """step3:"""
        #prune the big tree to mapping tree
        self.tree.prune(self.represent_taxon, preserve_branch_length=True)
        self.tree.resolve_polytomy(default_dist=0.0000001)
        self.tree.dist = 0
        #name the internal nodes; random init internal nodes coords
        cnt = 0
        for node in self.tree.traverse(strategy="postorder"):
            if not node.is_leaf():
                node.name = "x2j2i2a"+repr(cnt)
                X = [0, 0, 0]
                X[0] = node.dist
                if node.is_root():
                    X[0] = 0.0
                childs = node.get_children()
                X[1] = childs[0].dist
                X[2] = childs[1].dist
                self.innernode_dis[node.name] = X 
                self.name_coords[node.name] = [self.rand_nr.uniform(0.0,1.0), self.rand_nr.uniform(0.0,1.0)]
                self.all_node_names.append(node.name)
                self.inner_node_names.append(node.name)
                cnt = cnt + 1
            else:
                self.all_node_names.append(node.name)
        
        for node in self.tree.traverse(strategy="postorder"):
            if not node.is_leaf():
                N = []
                if node.is_root():
                    childs = node.get_children()
                    N.append(childs[0].name)
                    N.append(childs[1].name)
                else:
                    N.append(node.up.name)
                    childs = node.get_children()
                    N.append(childs[0].name)
                    N.append(childs[1].name)
                self.innernode_connecting_nodes[node.name] = N
        self._calculate_mds_distance()
    
    
    def _calculate_mds_distance(self):
        """generate and update mds pair-wise distance matrix"""
        name_coords_map = self.name_coords
        for node in self.tree.traverse(strategy="postorder"):
            if not node.is_leaf():
                cn = name_coords_map[node.name]
                childs = node.get_children()
                c1 = name_coords_map[childs[0].name]
                c2 = name_coords_map[childs[1].name]
                x1 = distance(cn, c1)
                x2 = distance(cn, c2)
                self.innernode_dis_mds_matrix[node.name+":"+childs[0].name] = x1
                self.innernode_dis_mds_matrix[node.name+":"+childs[1].name] = x2
                self.innernode_dis_mds_matrix[childs[0].name+":"+node.name] = x1
                self.innernode_dis_mds_matrix[childs[1].name+":"+node.name] = x2
    
    
    def calculate_errors(self):
        """compute the mapping errors for every iteration"""
        #generate ground truth, tree space distance matrix
        if self.sum_species_tree_length <=0:
            for node in self.tree.traverse(strategy="postorder"):
                if not node.is_root():
                    self.sum_species_tree_length += node.dist
                if not node.is_leaf():
                    childs = node.get_children()
                    self.innernode_dis_tree_matrix[node.name+":"+childs[0].name] = childs[0].dist
                    self.innernode_dis_tree_matrix[childs[0].name+":"+node.name] = childs[0].dist
                    self.innernode_dis_tree_matrix[node.name+":"+childs[1].name] = childs[1].dist
                    self.innernode_dis_tree_matrix[childs[1].name+":"+node.name] = childs[1].dist
        err = 0.0
        for i in range(len(self.all_node_names)):
            namei = self.all_node_names[i]
            for j in range(i, len(self.all_node_names)):
                namej = self.all_node_names[j]
                d_tree = self.innernode_dis_tree_matrix.get(namei+":"+namej)
                d_mds  = self.innernode_dis_mds_matrix.get(namei+":"+namej)
                if d_tree!=None and d_mds!=None:
                    err=err+(d_tree-d_mds)*(d_tree-d_mds)/d_tree
        return err/self.sum_species_tree_length
    
    
    def update(self, node_name, only_shorter = False):
        coord_upnode = self.name_coords[node_name]
        d1x = d1y = d2x = d2y = 0.00
        connected_nodes = self.innernode_connecting_nodes[node_name]
        changed = False
        for cnode in connected_nodes:
            d_tree = self.innernode_dis_tree_matrix[node_name + ":" +cnode]
            d_mds  = self.innernode_dis_mds_matrix[node_name + ":" +cnode]
            coord_cnode = self.name_coords[cnode]
            
            if only_shorter:
                if d_tree < d_mds:
                    changed = True
                    d1x = d1x+firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0])
                    d2x = d2x+secondDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0])
                    d1y = d1y+firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1])
                    d2y = d2y+secondDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1])
            else:
                d1x = d1x+firstDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0])
                d2x = d2x+secondDstep(d_tree, d_mds, coord_cnode[0], coord_upnode[0])
                d1y = d1y+firstDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1])
                d2y = d2y+secondDstep(d_tree, d_mds, coord_cnode[1], coord_upnode[1])
        
        if only_shorter:
            if changed:
                d2x = abs(d2x)
                d2y = abs(d2y)
                deltaX=-d1x/d2x
                coord_upnode[0]=coord_upnode[0]-self.mf*deltaX
                deltaY=-d1y/d2y
                coord_upnode[1]=coord_upnode[1]-self.mf*deltaY
                self.name_coords[node_name] = coord_upnode
        else:
            d2x = abs(d2x)
            d2y = abs(d2y)
            deltaX=-d1x/d2x
            coord_upnode[0]=coord_upnode[0]-self.mf*deltaX
            deltaY=-d1y/d2y
            coord_upnode[1]=coord_upnode[1]-self.mf*deltaY
            self.name_coords[node_name] = coord_upnode


    def mapping(self):
        ev1, ev2 = self.pcoa()
        print("Variance explained by first axis: {0:.2f} %".format(ev1 * 100))
        print("Variance explained by second axis: {0:.2f} %".format(ev2 * 100))
        self.parse_delimitation()
        self.extract_species_tree()
        mapping_err = self.calculate_errors()
        print("init mapping error = " + repr(mapping_err))
        
        improve_cnt = 0
        last_mapping_err = mapping_err
        index = range(len(self.inner_node_names))
        for i in range(self.maxiters):
            random.shuffle(index)
            for idx in index:
                inodename = self.inner_node_names[idx]
                if idx % 10 == 0:
                    self.update(inodename)
                else:
                    self.update(inodename, only_shorter = True)
                self._calculate_mds_distance()
            mapping_err = self.calculate_errors()
            if abs(last_mapping_err - mapping_err) <= self.min_improvement:
                improve_cnt = improve_cnt + 1
            else:
                improve_cnt = 0
            last_mapping_err = mapping_err
            if improve_cnt >= 200:
                break
            if self.debug:
                print("Mapping error after iteration " + repr(i) + ": " + repr(mapping_err))
        print("Mapping error after iteration " + repr(i) + ": " + repr(mapping_err))
        self.output_branches()


    def output_branches(self):
        outfile = open(self.fout + ".line.txt", "w")
        for node in self.tree.traverse(strategy="postorder"):
            if not node.is_leaf():
                cn = self.name_coords[node.name]
                childs = node.get_children()
                cl = self.name_coords[childs[0].name]
                cr = self.name_coords[childs[1].name]
                l_mds_dis = self.innernode_dis_mds_matrix[node.name + ":" + childs[0].name]
                l_tree_dis = self.innernode_dis_tree_matrix[node.name + ":" + childs[0].name] 
                r_mds_dis = self.innernode_dis_mds_matrix[node.name + ":" + childs[1].name]
                r_tree_dis = self.innernode_dis_tree_matrix[node.name + ":" + childs[1].name] 
                l_ratio = (l_tree_dis - l_mds_dis)/l_mds_dis
                r_ratio = (r_tree_dis - r_mds_dis)/r_mds_dis
                lstroke = math.erf(l_ratio)
                rstroke = math.erf(r_ratio)
                if lstroke < 0:
                    lstroke = 2.0
                else:
                    lstroke = 4.0*lstroke + 2.0
                
                if rstroke < 0:
                    rstroke = 2.0
                else:
                    rstroke = 4.0*rstroke + 2.0 
                
                outfile.write(repr(cn[0])+","+repr(cn[1])+","+repr(cl[0])+","+repr(cl[1])+","+repr(lstroke)+"\n")
                outfile.write(repr(cn[0])+","+repr(cn[1])+","+repr(cr[0])+","+repr(cr[1])+","+repr(rstroke)+"\n")
        outfile.close()


def parse_arguments():
    parser = argparse.ArgumentParser(description="""PhyloMap: an algorithm for visualizing relationships of large sequence data sets.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 
and 
"J. Zhang, A. M. Mamlouk, T. Martinetz, S. Chang, J. Wang, and
R. Hilgenfeld. PhyloMap: an algorithm for visualizing relationships of
large sequence data sets and its application to the influenza A virus
genome. BMC bioinformatics, 12:248, Jan. 2011."


Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com
Visit http://www.exelixis-lab.org/ for more information.

Version 0.1 released by Jiajie Zhang on 28-06-2014.""",
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        prog= "python phylomap.py")
    
    parser.add_argument("-t", dest = "trees", 
                        help = """input comprehensive phylogenetic tree file.""",
                        required = True)

    parser.add_argument("-p", dest = "ptp_result", 
                        help = """PTP delimitation results""",
                        required = True)

    parser.add_argument("-o", dest = "output",
                        help = "output file folder name",
                        required = True)

    parser.add_argument("-m", dest = "max_iters", 
                        help = """max optimization iterations  (default 10000)""",
                        type = int,
                        default = 10000)

    parser.add_argument("-s", dest = "seed", 
                        help = """random seed  (default 222)""",
                        type = int,
                        default = 222)

    parser.add_argument("--debug", 
                        help = """output debug information""",
                        default = False,
                        action="store_true")

    parser.add_argument('--version', action='version', version='%(prog)s 0.1 (28-06-2014)')
    
    return parser.parse_args()


if __name__ == "__main__":
    if len(sys.argv) == 1: 
        sys.argv.append("-h")
    args = parse_arguments()
    
    if not os.path.exists(args.trees):
        print("Input tree file does not exists: %s" % args.trees)
        sys.exit()

    if not os.path.exists(args.ptp_result):
        print("Input PTP results file does not exists: %s" % args.ptp_result)
        sys.exit()

    treetest = open(args.trees)
    l1 = treetest.readline()
    treetest.close()
    inputformat = "nexus"
    if l1.strip() == "#NEXUS":
        inputformat = "nexus"
    else:
        inputformat = "raxml"
    
    basepath = os.path.dirname(os.path.abspath(__file__))
    
    #assert not os.path.isabs(args.output)
    dstdir =  os.path.join(os.path.dirname(args.output), "phylomap")
    dstdir2 = os.path.join(dstdir, "data")
    try:
        os.makedirs(dstdir)
        os.makedirs(dstdir2)
        shutil.copy2(os.path.join(basepath, "phylomap.pde"), dstdir)
    except OSError, e:
        print(e)
        sys.exit()
    
    pm = phylomap(largetree = args.trees, 
                  ptp_result = args.ptp_result, 
                  fout = os.path.join(dstdir2, "data"),
                  ftype = inputformat, 
                  maxiters = args.max_iters,
                  seed = args.seed,
                  debug = args.debug)
    pm.mapping()
    
    subprocess.call([basepath + "/bin/processing-2.2.1/processing-java",
    "--sketch="+dstdir,
    "--output="+dstdir+"/visualapp",
    "--force",
    "--export"])
    #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
    
    print("Java executable processing visualisation app written to: " + dstdir+"/visualapp")
