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

class partitionparser:
	def __init__(self, fin):
		self.partitions = []
		self.supports = []
		self.weight = 1 #this feature is disabled
		
		with open(fin) as f:
			lines = f.readlines()
			for line in lines:
				if line.startswith("#"):#ignore all other lines starts with #
					if line.startswith("#taxaorder"):
						self.taxa_order = line.strip().split(":")[1].split()
						self.numtaxa = len(self.taxa_order)
				else:
					part = []
					supp = []
					ss = line.strip().split()
					for el in ss:
						els = el.split("|")
						part.append(int(els[0]))
						supp.append(float(els[1]) * self.weight)
					self.partitions.append(part)
					self.supports.append(supp)
		self.numtrees = len(self.partitions)
	
	def print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + "\t"
		return ss.strip() + "\n"
	
	def printout(self):
		outs = ""
		outs = outs + "#taxaorder:" + self.print_list(self.taxa_order)
		for i in range(len(self.partitions)):
			par = self.partitions[i]
			sup = self.supports[i]
			outs = outs + self.print_2lists(par, sup)
		return outs
	
	def get_taxa_order(self):
		return self.taxa_order
	
	def _translate(self, new_taxa_order, old_partition, old_support):
		new_partition = [-1] * len(new_taxa_order)
		new_support   = [-1] * len(new_taxa_order)
		for i in range(len(self.taxa_order)):
			parnum = old_partition[i]
			sup    = old_support[i]
			taxaname = self.taxa_order[i]
			newidx = new_taxa_order.index(taxaname)
			new_partition[newidx] = parnum 
			new_support[newidx] = sup
		return new_partition, new_support
	
	def translate_to(self, p2, fout):
		new_taxa_order = p2.get_taxa_order()
		fo = open(fout, "w")
		fo.write(p2.printout())
		for i in range(len(self.partitions)):
			part = self.partitions[i]
			supp = self.supports[i]
			np, ns = self._translate(new_taxa_order, part, supp)
			outs = self.print_2lists(np, ns)
			fo.write(outs)
		fo.close()
	
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
		for i in range(len(self.partitions)):
			partition = self.partitions[i]
			pars = self._convert2idx(partition)
			idxpars.append(pars)
			for par in pars:
				pmap[par]= pmap.get(par, 0) + 1
		self.supports = []
		maxw = 0
		bestpar = None
		bestsupport = None
		output = open(fout + ".merged_partitions", "a")
		output.write("#taxaorder:"+self.print_list(self.taxa_order))
		for i in range(len(self.partitions)): 
			partition = self.partitions[i]
			pars = idxpars[i]
			support = [0] * self.numtaxa
			sumw = 0.0
			for par in pars:
				w = pmap[par]
				for idx in par:
					support[idx] = float(w)/float(self.numtrees)
					sumw = sumw + w
			if sumw > maxw:
				maxw = sumw
				bestpar = i
				bestsupport = support
				
			self.supports.append(support)
			output.write(self.print_2lists(partition, support))
		
		output.close()
		
		bp, bs = self._partition2names(self.partitions[bestpar], bestsupport)
		print_species(spes = bp, support = bs, fout = fout + ".merged_simpleHeuristics", verbose = False, method = "simple heuristics")
		
		return pmap, maxw
	
	def print_2lists(self, l1, l2):
		ss = ""
		for i in range(len(l1)):
			e1 = l1[i]
			e2 = l2[i]
			ss = ss + str(e1)+"|"+str(e2) + "\t"
		return ss.strip() + "\n"
	
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
			namesupps.append(onesup[0])
		
		return nameparts, namesupps


class pnode:
	def __init__(self, newpartition, lastpnode, all_taxa, pmap, all_partitions, bound, taxa_order, numtrees):
		#all_taxa: set; newpartition: tuple of idices; lastnode: pnode; pmap: dictionary; all_partitions: list; bound: float
		self.numtrees = numtrees
		self.taxa_order = taxa_order
		self.bound = bound
		self.all_partitions = all_partitions
		self.pmap = pmap
		self.all_taxa = all_taxa
		self.all_idx = set(range(len(all_taxa)))
		self.up = lastpnode
		self.partition = newpartition
		self.allprocessed_taxa_idx = set([])
		if not self.up == None: 
			for idx in self.up.allprocessed_taxa_idx:
				self.allprocessed_taxa_idx.add(idx)
		for idx in newpartition:
			self.allprocessed_taxa_idx.add(idx)
		self.next_partition_list = []
		self.support = self.pmap.get(self.partition, 0) * len(self.partition)
		if not self.up == None: 
			self.sumsupport = self.up.sumsupport + self.support
		else:
			self.sumsupport = self.support
		self.isvalid = False
		
	def search(self):
		remains = self.all_idx - self.allprocessed_taxa_idx
		
		if (len(remains) * self.numtrees + self.sumsupport) < self.bound:
			return [] 
		
		if len(remains) > 0:
			next_taxa_to_include_idx = list(remains)[0]
			
			for par in self.all_partitions:
				spar = set(par)
				if (next_taxa_to_include_idx in spar) and (len(spar & self.allprocessed_taxa_idx) == 0):
					self.next_partition_list.append(par)
			
			return self.next_partition_list
		else:
			self.isvalid = True
			return []
	
	def get_support(self):
		return float(self.support) / float(len(self.partition) * self.numtrees)
	
	def __str__(self):
		return repr(self.partition) + " : "+repr(self.support)


def bbsearch(pmap, taxa_order, bound, numtrees):
	#init
	global maxsupport 
	global bestlastnode 
	maxsupport = 0 
	bestlastnode = None
	all_taxa_set = set(taxa_order)
	all_partitions = pmap.keys()
	pnode0 = pnode(newpartition = tuple([]), lastpnode = None, all_taxa = all_taxa_set, pmap = pmap, all_partitions = all_partitions, bound = bound, taxa_order = taxa_order, numtrees = numtrees)
	inipartitions = pnode0.search()
	
	#search for every cases containing the first taxa name
	for par in inipartitions:
		rc_function(newpartition = par , lastpnode = pnode0, all_taxa = all_taxa_set, pmap = pmap, all_partitions = all_partitions, bound = bound, taxa_order = taxa_order, numtrees = numtrees)
	
	#return the results
	spes = []
	support = []
	if not bestlastnode == None:
		print("Max support value: " + repr(maxsupport))
		currnode = bestlastnode
		while not currnode.up == None:
			spe = []
			for idx in currnode.partition:
				spe.append(taxa_order[idx])
			spes.append(spe)
			support.append(currnode.get_support())
			currnode = currnode.up
	else:
		print("Bestlastnode == None")
	
	return spes, support


def print_species(spes, support, fout = "", verbose = True, method = "branch and bound"):
	if not fout == "":
		with open(fout, "a") as fo:
			fo.write("#------------------------------------------------------------------------------------------------------------------------------------------------\n")
			fo.write("#Most supported partitions found by " + method + " search\n")
			for i in range(len(spes)):
				spe = spes[i]
				sup = support[i]
				fo.write("Species " + str(i+1) + " (support = " + repr(sup) + ")\n")
				fo.write("     " + print_list(spe) + "\n")
	
	if verbose:
		print("---------------------------------------------------------------------")
		print("Most supported partitions found by " + method + " search")
		for i in range(len(spes)):
			spe = spes[i]
			sup = support[i]
			print("Species " + str(i+1) + " (support = " + repr(sup) + ")")
			print("    " + print_list(spe))


def print_list(l):
	ss = ""
	for e in l:
		ss = ss + str(e) + ","
	return ss[:-1]


def rc_function(newpartition, lastpnode, all_taxa, pmap, all_partitions, bound, taxa_order, numtrees):
	global maxsupport
	global bestlastnode
	nextnode = pnode(newpartition = newpartition, lastpnode = lastpnode, all_taxa = all_taxa, pmap = pmap, all_partitions = all_partitions, bound = bound, taxa_order = taxa_order, numtrees = numtrees)
	nextpartitions = nextnode.search()
	if len(nextpartitions) == 0:
		if nextnode.isvalid:
			if nextnode.sumsupport > maxsupport:
				maxsupport = nextnode.sumsupport
				bestlastnode = nextnode
	else:
		for par in nextpartitions:
			rc_function(newpartition = par, lastpnode = nextnode, all_taxa = all_taxa, pmap = pmap, all_partitions = all_partitions, bound = bound, taxa_order = taxa_order, numtrees = numtrees)


class mptp:
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
	
	def remove_outgroups(self, ognames, remove = False):
		try:
			if remove:
				for og in ognames:
					self.taxa_order.remove(og)
				self.numtaxa = len(self.taxa_order)
			for i in range(len(self.trees)):
				t = Tree(self.trees[i])
				if len(ognames) < 2:
					t.set_outgroup(ognames[0])
					if remove:
						t.prune(self.taxa_order, preserve_branch_length=True)
				else:
					ancestor = t.get_common_ancestor(ognames)
					if not t == ancestor:
						t.set_outgroup(ancestor)
					if remove:
						t.prune(self.taxa_order, preserve_branch_length=True)
				self.trees[i] = t.write()
		except ValueError, e:
			print(e)
			print("")
			print("")
			print("Somthing is wrong with the input outgroup names")
			print("")
			print("Quiting .....")
			sys.exit()
	
	def delimit(self, fout, sreroot = False, pvalue = 0.001, weight = 1):
		self.weight = weight
		self.partitions = []
		cnt = 1 
		for tree in self.trees:
			print("Delimiting on tree " + repr(cnt) + "........")
			cnt = cnt + 1
			me = exponential_mixture(tree= tree, max_iters = 20000, min_br = 0.0001 )
			me.search(reroot = sreroot, strategy = "H0")
			me.count_species(pv = pvalue)
			order, partition = me.output_species(self.taxa_order)
			self.partitions.append(partition)
			print("")
		pmap, bound = self.summary(fout)
		return pmap, bound
	
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
		output = open(fout + ".mPTP_partitions", "a")
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
					sumw = sumw + w #float(w)/float(self.numtrees)
			if sumw > maxw:
				maxw = sumw
				bestpar = i
				bestsupport = support
				
			self.supports.append(support)
			output.write(self.print_2lists(partition, support))
		output.close()
		
		bp, bs = self._partition2names(self.partitions[bestpar], bestsupport)
		print_species(spes = bp, support = bs, fout = fout + ".mPTP_simpleHeuristics", verbose = False, method = "simple heuristics")
		
		return pmap, maxw

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
			namesupps.append(onesup[0])
		
		return nameparts, namesupps


def print_options():
		print("usage: python mPTP.py -t example/nex.test -o example/nex.out -r")
		print("Options:")
		print("    -t input                       Specify the input NEXUS file, trees can be both rooted or unrooted,")
		print("                                   if unrooted, please use -r option.\n")
		print("    -o output                      Specify output file name.\n")
		print("    -g outgroupnames               t1,t2,t3  commma delimt and no space in between")
		print("    -r                             Rooting the input tree on the longest branch.(default not)\n")
		print("    -d                             Remove outgroups specified by -g. (default not)\n")
		print("    -pvalue (0-1)                  Set the p-value for likelihood ratio test.(default 0.001)")


if __name__ == "__main__":
	print("This is mPTP - a Poisson tree processes model for species delimitation using multiple input trees")
	print("Version 0.1 released by Jiajie Zhang on 25-11-2013\n")
	print("This program will delimit species on multiple phylogenetic trees.")
	print("The input file should be in NEXUS format (such as the output from MrBayes).")
	print("--Please cite: \"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General") 
	print("--Species Delimitation Method with Applications to Phylogenetic Placements. ")
	print("--Bioinformatics (2013), 29 (22): 2869-2876.\" ")
	print("--If you find PTP is useful to your research. \n")
	print("Questions and bug reports, please send to:")
	print("bestzhangjiajie@gmail.com\n")

	if len(sys.argv) < 3: 
		print_options()
		sys.exit()
	
	stree = ""
	sreroot = False
	pvalue = 0.001
	ptpout = ""
	outgroups = []
	rmog = False
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-t":
			i = i + 1
			stree = sys.argv[i]
		elif sys.argv[i] == "-r":
			sreroot = True
		elif sys.argv[i] == "-d":
			rmog = True
		elif sys.argv[i] == "-g":
			i = i + 1
			ogs = sys.argv[i].strip()
			outgroups = ogs.split(",")
		elif sys.argv[i] == "-pvalue":
			i = i + 1
			pvalue = float(sys.argv[i])
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
		mp = mptp(filename = stree, ftype = inputformat)
		if len(outgroups) > 0:
			mp.remove_outgroups(outgroups, remove = rmog)
		
		pmap, b = mp.delimit(fout = ptpout, sreroot = sreroot, pvalue = pvalue, weight = 1)
		
		print("Lower bound support value:" + repr(b))
		
		spes, supports = bbsearch(pmap = pmap, taxa_order = mp.taxa_order, bound = b, numtrees = mp.numtrees)
		
		print_species(spes, supports, fout = ptpout + ".mPTP_bestPartitions", verbose = True)
		
	except ete2.parser.newick.NewickError:
		print("Unexisting tree file or Malformed newick tree structure.")


