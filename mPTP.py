#! /usr/bin/env python
try:
	import sys
	import math
	import collections
	import ete2
	import os
	import subprocess
	from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
	from collections import deque
	from scipy import stats
	from numpy import array
	from subprocess import call
	from PTP import *
	from nexus import NexusReader
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()

class mptp:
	def __init__(self, filename, ftype = "nexus"):
		if ftype == "nexus":
			self.nexus = NexusReader(filename)
			self.nexus.blocks['trees'].detranslate()
			self.trees = self.nexus.trees.trees
		else:
			self.trees = self.raxmlTreeParser(filename)
		self.taxa_order = Tree(self.trees[0]).get_leaf_names()
	
	
	def delimit(self, fout, sreroot = False, pvalue = 0.001, weight = 1):
		output = open(fout, "a")
		output.write("#"+self.print_list(self.taxa_order))
		self.partitions = []
		for tree in self.trees:
			me = exponential_mixture(tree= tree, max_iters = 20000, min_br = 0.0001 )
			me.search(reroot = sreroot, strategy = "H0")
			me.count_species(pv = pvalue)
			order, partition = me.output_species(self.taxa_order)
			self.partitions.append(partition)
			output.write(str(weight) + ":" + self.print_list(partition))
			
		output.close()
		return self.partitions
	
	
	def print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + "\t"
		return ss.strip() + "\n"
	
	
	def raxmlTreeParser(self, fin):
		f = open(fin)
		lines = f.readlines()
		f.close()
		trees = []
		for line in lines:
			line = line.strip()
			if not line == "":
				trees.append(line[line.index("("):])
		return trees
	
	
	def summary(self):
		pmap = {}
		



if __name__ == "__main__":
	print("This is PTP - a Poisson tree processes model for species delimitation.")
	print("Version 1.2 released by Jiajie Zhang on 10-11-2013\n")
	print("This program will delimit species on a rooted phylogenetic tree.")
	print("The input tree should be in Newick format (such as the output from RAxML).")
	print("PTP can also infer phylogenetic tree using RAxML, currently only support DNA on GTRGAMMA.\n")
	#print("The program needs ETE(http://ete.cgenomics.org/) package to be installed.\n")
	print("--Please cite: \"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General") 
	print("--Species Delimitation Method with Applications to Phylogenetic Placements. ")
	print("--Bioinformatics (2013), 29 (22): 2869-2876.\" ")
	print("--If you find PTP is useful to your research. \n")
	print("Questions and bug reports, please send to:")
	print("bestzhangjiajie@gmail.com\n")
	
	mp = mptp("/home/zhangje/GIT/SpeciesCounting/example/nex.test")
	mp.delimit(fout = "/home/zhangje/GIT/SpeciesCounting/example/nex.test.out")
