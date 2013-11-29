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
	from GMYC import *
	from nexus import NexusReader
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()


class mgmyc:
	def __init__(self, filename, ftype = "nexus"):
		if ftype == "nexus":
			self.nexus = NexusReader(filename)
			self.nexus.blocks['trees'].detranslate()
			self.trees = self.nexus.trees.trees
		else:
			self.trees = self.raxmlTreeParser(filename)
		self.taxa_order = Tree(self.trees[0]).get_leaf_names()
		self.numtaxa = len(self.taxa_order)
		self.numtrees = len(self.trees)
	
	
	def delimit(self, fout, pvalue = 0.01, weight = 1):
		self.weight = weight
		#output = open(fout, "a")
		#output.write("#"+self.print_list(self.taxa_order))
		self.partitions = []
		for tree in self.trees:
			partition = gmyc_func(tree = tree, taxa_order = self.taxa_order, pv = pvalue)
			self.partitions.append(partition)
			#output.write(str(weight) + ":" + self.print_list(partition))
			
		#output.close()
		besttreeidx = self.summary(fout)
		return self.partitions
	
	
	def print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + "\t"
		return ss.strip() + "\n"
	
	
	def print_2lists(self, l1, l2):
		ss = ""
		for i in range(len(l1)):
			e1 = l1[i]
			e2 = l2[i]
			ss = ss + str(e1)+"|"+str(e2) + "\t"
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
	
	
	def _convert2idx(self, partition):
		a = min(partition)
		b = max(partition) + 1
		par = []
		for i in range(a, b):
			indices = [j for j, x in enumerate(partition) if x == i]
			par.append(tuple(indices))
		return par
	
	
	def summary(self, fout):
		pmap = {}
		idxpars = []
		for partition in self.partitions:
			pars = self._convert2idx(partition)
			idxpars.append(pars)
			for par in pars:
				pmap[par]= pmap.get(par, 0) + self.weight
		self.supports = []
		maxw = 0
		bestpar = None
		bestsupport = None
		output = open(fout, "a")
		output.write("#weight:" + repr(self.weight) + "\n")
		output.write("#taxaorder:"+self.print_list(self.taxa_order))
		for i in range(len(self.partitions)): 
			partition = self.partitions[i]
			pars = idxpars[i]
			support = [1] * self.numtaxa
			sumw = 0.0
			for par in pars:
				w = pmap[par]
				for idx in par:
					support[idx] = float(w)/float(self.numtrees)
					sumw = sumw + float(w)/float(self.numtrees)
			if sumw > maxw:
				maxw = sumw
				bestpar = i
				bestsupport = support
				
			self.supports.append(support)
			output.write(self.print_2lists(partition, support))
		output.write("#best:" + self.print_2lists(self.partitions[bestpar], bestsupport))
		
		bp, bs = self._partition2names(self.partitions[bestpar], bestsupport)
		for i in range(len(bp)):
			pi = bp[i]
			si = bs[i]
			s = si[0]
			output.write("#best: Species " + repr(i) + "-----" + repr(s)+"\n")
			output.write("#best:     " + self.print_list(pi))
		
		output.close()
		return bestpar
	
	
	def _partition2names(self, part, supp):
		nameparts = []
		namesupps = []
		a = min(part)
		b = max(part) + 1
		par = []
		for i in range(a, b):
			onepar = []
			onesup = []
			for j in range(len(part)):
				idfier = part[j]
				sup = supp[j]
				if idfier == i:
					onepar.append(self.taxa_order[j])
					onesup.append(sup)
			nameparts.append(onepar)
			namesupps.append(onesup)
		
		return nameparts, namesupps


def print_options():
		print("usage: python mGMYC.py -t example/umnex.test -o example/umnex.out")
		print("Options:")
		#print("    -a alignment                   Specify the alignment, PTP will build a phylogenetic tree using RAxML, currently only support DNA sequences with GTRGAMMA.\n")
		print("    -t input_tree_file             Specify the input NEXUS file, trees must be ultrametric.")
		#print("                                   if unrooted, please use -r option.\n")
		print("    -o output_name                 Specify output file name.\n")
		#print("    -r                             Rooting the input tree on the longest branch.(default not)\n")
		#print("    -s                             Plot the delimited species on the tree.(default not show)\n")
		#print("    -p                             Print delimited species on the screen.(default not show)\n")
		print("    -pvalue (0-1)                  Set the p-value for likelihood ratio test.(default 0.01)")
		#print("                                   If the test failed, the program will output only one species.")
		#print("                                   Note this could mean there is only one species or all input taxon are different species.\n")


if __name__ == "__main__":
	print("This is mGMYC - GMYC species delimitation using multiple input trees")
	print("Version 0.1 released by Jiajie Zhang on 25-11-2013\n")
	print("This program will delimit species on multiple um-trees.")
	print("The input file should be in NEXUS format (such as the output from MrBayes).")
	#print("PTP can also infer phylogenetic tree using RAxML, currently only support DNA on GTRGAMMA.\n")
	#print("The program needs ETE(http://ete.cgenomics.org/) package to be installed.\n")
	print("--Please cite: \"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General") 
	print("--Species Delimitation Method with Applications to Phylogenetic Placements. ")
	print("--Bioinformatics (2013), 29 (22): 2869-2876.\" ")
	print("--If you find mGMYC is useful to your research. \n")
	print("Questions and bug reports, please send to:")
	print("bestzhangjiajie@gmail.com\n")

	if len(sys.argv) < 3: 
		print_options()
		sys.exit()
	
	salignment = ""
	stree = ""
	sreroot = False
	sstrategy = "H0"
	sprint = False 
	sshow = False
	sscale = 500
	max_iter = 20000
	min_brl = 0.0001
	spe_rate = -1.0
	whiten = False
	pvalue = 0.01
	ptpout = ""
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-t":
			i = i + 1
			stree = sys.argv[i]
		elif sys.argv[i] == "-m":
			i = i + 1
			sstrategy = sys.argv[i]
		elif sys.argv[i] == "-r":
			sreroot = True
		elif sys.argv[i] == "-s":
			sshow = True
		elif sys.argv[i] == "-c":
			i = i + 1
			sscale = int(sys.argv[i])
		elif sys.argv[i] == "-p":
			sprint = True 
		elif sys.argv[i] == "-maxiters":
			i = i + 1
			max_iter = int(sys.argv[i])
		elif sys.argv[i] == "-minbr":
			i = i + 1
			min_brl = float(sys.argv[i])
		elif sys.argv[i] == "-pvalue":
			i = i + 1
			pvalue = float(sys.argv[i])
		elif sys.argv[i] == "-sprate":
			i = i + 1
			spe_rate = float(sys.argv[i])
		elif sys.argv[i] == "-w":
			whiten = True
		elif sys.argv[i] == "-a":
			i = i + 1
			salignment = sys.argv[i]
		elif sys.argv[i] == "-o":
			i = i + 1
			ptpout = sys.argv[i]
		elif i == 0:
			pass
		elif sys.argv[i].startswith("-"):
			print("Unknown options: " + sys.argv[i])
			print_options()
			sys.exit()
	
	
	if stree == "":
		print("The input tree is empty.")
		print_options()
		sys.exit()
		
	if ptpout == "":
		ptpout = stree + ".out"

	try:
		treetest = open(stree)
		l1 = treetest.readline()
		treetest.close()
		inputformat = "nexus"
		if l1.strip() == "#NEXUS":
			inputformat = "nexus"
		else:
			inputformat = "raxml"
		
		mp = mgmyc(filename = stree, ftype = inputformat)
		mp.delimit(fout = ptpout, pvalue = pvalue, weight = 1)

	except ete2.parser.newick.NewickError, e:
		print e
		print("Unexisting tree file or Malformed newick tree structure.")


