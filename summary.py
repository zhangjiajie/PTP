#! /usr/bin/env python
try:
	import sys
	import math
	import numpy
	import collections
	import ete2
	import os
	import random
	from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
	from collections import deque
	from scipy import stats
	from numpy import array
	from nexus import NexusReader
	import matplotlib.pyplot as plt
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()

class partitionparser:
	def __init__(self, fin = None, taxa_order = None, partitions = [], llhs = []):
		self.taxa_order = taxa_order
		self.partitions = []
		self.llhs = []
		if fin != None: 
			with open(fin) as f:
				lines = f.readlines()
				for line in lines:
					if line.startswith("#"):#ignore all other lines starts with #
						if line.startswith("#taxaorder"):
							self.taxa_order = line.strip().split(":")[1].split(",")
							self.numtaxa = len(self.taxa_order)
					else:
						part = line.strip().split(":")
						self.llhs.append(float(part[0]))
						par = part[1].split(",")
						ipar = []
						for p in par:
							ipar.append(int(p))
						self.partitions.append(ipar)
		
		self.numtrees = len(self.partitions)
	
	
	def hpd(self, region = 0.95):
		"""Credible interval / delimitations in our case """
		self.sorted_llhs, self.sorted_partitions = zip(*sorted(zip(self.llhs, self.partitions), reverse = True))
		sumlogl = sum(self.llhs)
		psumlogl = float(sumlogl) * region
		
		idxend = 0
		accumlogl = 0.0
		for i in range(len(self.sorted_partitions)):
			accumlogl = accumlogl + self.sorted_llhs[i]
			if accumlogl >= psumlogl:
				idxend = i 
				break
		self.hpdidx = idxend
		return self.sorted_partitions[:idxend], self.sorted_llhs[:idxend]
	
	
	def hpd_numpartitions(self):
		pmlist = []
		for partition in self.self.sorted_partitions[:self.hpdidx]:
			pmlist.append(len(partition))
		
		return min(pmlist), max(pmlist)
	
	
	def summary(self, fout = "", region = 1.0):
		if region >= 1.0:
			tpartitions = self.partitions
			tllhs = self.llhs
		else:
			tpartitions, tllhs = self.hpd(region = region)
		
		if fout!="":
			fo_partsum = open(fout + ".bPTPPartitonSummary.txt", "w")
			fo_parts   = open(fout + ".bPTPPartitions.txt", "w")
			pmap = {}
			idxpars = []
			for partition in tpartitions:
				pars = self._convert2idx(partition)
				idxpars.append(pars)
				for par in pars:
					pmap[par]= pmap.get(par, 0) + 1
			
			for key, value in sorted(pmap.iteritems(), reverse = True, key=lambda (k,v): (v,k)):
				onespe = ""
				for idx in key:
					onespe = onespe + ", " + self.taxaorder[idx]
				onespe = onespe[1:]
				fo_partsum.write(onespe + ": " + "{0:.3f}".format(float(value)/float(len(tpartitions))) + "\n")
			
			fo_partsum.close()
			
			output = "#taxaorder:"+self._print_list(self.taxaorder)
			for i in range(len(tpartitions)): 
				partition = tpartitions[i]
				output= output + self._print_list(partition)
			fo_parts.write(output)
			fo_parts.close()
			
			plt.plot(tllhs)
			plt.ylabel('Log likelihood')
			plt.xlabel('Iterations')
			plt.savefig(fout + ".llh.png")
			
			self.combine_simple_heuristic(tpartitions = tpartitions, pmap = pmap, fo = fout + ".bPTPhSupportPartition.txt")
	
	
	def combine_simple_heuristic(self, tpartitions, pmap, fo):
		maxw = 0
		bestpar = None
		bestsupport = None
		supports = []
		
		for i in range(len(tpartitions)): 
			partition = tpartitions[i]
			pars = idxpars[i]
			support = [0.0] * self.numtaxa
			sumw = 0.0
			for par in pars:
				w = pmap[par]
				for idx in par:
					support[idx] = float(w)/float(len(tpartitions))
					sumw = sumw + w 
			if sumw > maxw:
				maxw = sumw
				bestpar = i
				bestsupport = support
			supports.append(support)
		
		spes, support = self._partition2names(tpartitions[bestpar], bestsupport)
		
		fo_bestpar = open(fo, "w")
		fo_bestpar.write("# Most supported partition found by simple heuristic search\n")
		for i in range(len(spes)):
			spe = spes[i]
			sup = support[i]
			fo_bestpar.write("Species " + str(i+1) + " (support = " + "{0:.3f}".format(sup) + ")\n")
			fo_bestpar.write("     " + self._print_list(spe) + "\n")
		fo_bestpar.close()
	
	
	def combine_max_MI(self, tpartitions, fout):
		#TODO 
		pass
	
	
	def _print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + ","
		return ss[:-2] + "\n"
	
	
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
