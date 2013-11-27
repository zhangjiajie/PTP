#! /usr/bin/env python
try:
	import sys
	import math
	import os
	import glob
	from ete2 import SeqGroup
	from subprocess import call
	from mPTP import *
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()
	

def crop_species_counting(falin):
	fafa = falin
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/crop", "-i", fafa]) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	foutcrop = open(fafa + ".cluster.list")
	lines = foutcrop.readlines()
	foutcrop.close()
	spes = []
	for line in lines:
		lls = line.strip().split()
		taxas = lls[1]
		taxon = taxas.split(",")
		spes.append(taxon)
	
	jks = glob.glob(fafa + ".*")
	for jk in jks:
		os.remove(jk)
	if os.path.exists("LikelihoodRatio.txt"):
		os.remove("LikelihoodRatio.txt")
	
	return spes

class mcrop:
	def __init__(self, foldname):
		self.partitions = []
		self.supports = []
		self.taxa_order = []
		jks = glob.glob( foldname + "*")
		seq1 = SeqGroup(jks[0], format = "phylip_relaxed")
		for entri in seq1.iter_entries():
			self.taxa_order.append(entri[0])
		self.numtaxa = len(self.taxa_order)
		self.numtrees = len(jks)
		for jk in jks:
			seqs = SeqGroup(jk, format='phylip_relaxed')
			seqs.write(format='fasta', outfile=jk+".fa")
			spes = crop_species_counting(jk+".fa")
			self.partitions.append(self.output_species(self.taxa_order, spes))
			os.remove(jk + ".fa")
			
			
	def output_species(self, taxa_order, species_list):
		"""taxa_order is a list of taxa names, the paritions will be output as the same order"""
		num_taxa = len(taxa_order)
		partion = [-1] * num_taxa
		cnt = 1
		for sp in species_list:
			for leaf in sp:
				idx = taxa_order.index(leaf)
				partion[idx] = cnt
			cnt = cnt + 1
		return partion

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
				pmap[par]= pmap.get(par, 0) + 1
		self.supports = []
		maxw = 0
		bestpar = None
		bestsupport = None
		output = open(fout, "a")
		output.write("#weight:" + repr(1) + "\n")
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


if __name__ == "__main__":
	print("main")
	mc = mcrop("example/bs/")
	mc.summary("example/boot.sum")
