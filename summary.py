#! /usr/bin/env python
try:
	import matplotlib.pyplot as plt
	import sys
	import argparse
	import os
	import numpy
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()

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



def mutual_info(x,y):
	"""Mutual information"""
	N = numpy.double(x.size)
	I = 0.0
	eps = numpy.finfo(float).eps
	for l1 in numpy.unique(x):
		for l2 in numpy.unique(y):
			#Find the intersections
			l1_ids=numpy.nonzero(x==l1)[0]
			l2_ids=numpy.nonzero(y==l2)[0]
			pxy=(numpy.double(numpy.intersect1d(l1_ids,l2_ids).size)/N)+eps
			I+=pxy*numpy.log2(pxy/((l1_ids.size/N)*(l2_ids.size/N)))
	return I



def nmi(x,y):
	"""Normalized mutual information"""
	x = numpy.array(x)
	y = numpy.array(y)
	N = numpy.double(x.size)
	I = mutual_info(x,y)
	Hx = 0
	for l1 in numpy.unique(x):
		l1_count=numpy.nonzero(x==l1)[0].size
		Hx+=-(numpy.double(l1_count)/N)*numpy.log2(numpy.double(l1_count)/N)
	Hy=0
	for l2 in numpy.unique(y):
		l2_count=numpy.nonzero(y==l2)[0].size
		Hy+=-(numpy.double(l2_count)/N)*numpy.log2(numpy.double(l2_count)/N)
	if (Hx+Hy) == 0:
		return 1.0
	else: 
		return I/((Hx+Hy)/2)



class partitionparser:
	def __init__(self, pfin = None, lfin = None, taxa_order = None, partitions = [], llhs = []):
		self.taxaorder = taxa_order
		self.partitions = partitions
		self.llhs = llhs
		if pfin != None:
			with open(pfin) as f:
				lines = f.readlines()
				for line in lines:
					if line.startswith("#"):#ignore all other lines starts with #
						if line.startswith("#taxaorder"):
							self.taxaorder = line.strip().split(":")[1].split(",")
							self.numtaxa = len(self.taxaorder)
					else:
						par = line.strip().split(",")
						ipar = []
						for p in par:
							ipar.append(int(p))
						self.partitions.append(ipar)
		if lfin != None: 
			with open(lfin) as f:
				lines = f.readlines()
				for line in lines:
					llh = float(line.strip())
					self.llhs.append(llh)
		
		self.numtaxa = len(self.taxaorder)
		self.numtrees = len(self.partitions)
		self.hpdidx = len(self.partitions)
		self.sorted_llhs = self.llhs 
		self.sorted_partitions = self.partitions
	
	
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
		idxend = self.hpdidx
		for partition in self.sorted_partitions[:idxend]:
			pmlist.append(max(partition))
		
		return min(pmlist), max(pmlist)
	
	
	def summary(self, fout = "", region = 1.0, bnmi = False):
		if region >= 1.0 or region <=0:
			tpartitions = self.partitions
			tllhs = self.llhs
		else:
			tpartitions, tllhs = self.hpd(region = region)
		
		if fout!="":
			fo_partsum = open(fout + ".PTPPartitonSummary.txt", "w")
			fo_parts   = open(fout + ".PTPPartitions.txt", "w")
			pmap = {}
			idxpars = []
			for partition in tpartitions:
				pars = self._convert2idx(partition)
				idxpars.append(pars)
				for par in pars:
					pmap[par]= pmap.get(par, 0) + 1
			
			"""Output partition summary"""
			for key, value in sorted(pmap.iteritems(), reverse = True, key=lambda (k,v): (v,k)):
				onespe = ""
				for idx in key:
					onespe = onespe + ", " + self.taxaorder[idx]
				onespe = onespe[1:]
				fo_partsum.write(onespe + ": " + "{0:.3f}".format(float(value)/float(len(tpartitions))) + "\n")
			fo_partsum.close()
			
			"""Output all partitions"""
			output = "#taxaorder:"+self._print_list(self.taxaorder)
			for i in range(len(tpartitions)): 
				partition = tpartitions[i]
				output= output + self._print_list(partition)
			fo_parts.write(output)
			fo_parts.close()
			
			"""Output the best partition found"""
			self.combine_simple_heuristic(tpartitions = tpartitions, pmap = pmap, idxpars = idxpars, fo = fout + ".PTPhSupportPartition.txt")
			
			if bnmi:
				self.combine_max_NMI(tpartitions = tpartitions, pmap = pmap, fo = fout + ".PTPhNMIPartition.txt")
			
			"""MCMC LLH"""
			if (region >= 1.0 or region <=0) and len(tllhs)>0:
				plt.plot(tllhs)
				plt.ylabel('Log likelihood')
				plt.xlabel('Iterations')
				plt.savefig(fout + ".llh.png")
				with open(fout + ".PTPllh.txt", "w") as f:
					for llh in tllhs:
						f.write(repr(llh) + "\n")
	
	
	def combine_simple_heuristic(self, tpartitions, pmap, idxpars, fo):
		maxw = 0
		bestpar = None
		bestsupport = None
		
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
		
		spes, support = self._partition2names(tpartitions[bestpar], bestsupport)
		
		fo_bestpar = open(fo, "w")
		fo_bestpar.write("# Most supported partition found by simple heuristic search\n")
		for i in range(len(spes)):
			spe = spes[i]
			sup = support[i]
			fo_bestpar.write("Species " + str(i+1) + " (support = " + "{0:.3f}".format(sup) + ")\n")
			fo_bestpar.write("     " + self._print_list(spe) + "\n")
		fo_bestpar.close()
	
	
	def combine_max_NMI(self, tpartitions, pmap, fo = ""):
		bestpar = None 
		hnmi = 0.0 
		for parx in tpartitions:
			sumnmi = 0 
			for pary in tpartitions: 
				sumnmi += nmi(parx, pary)
			if sumnmi >= hnmi:
				bestpar = parx
				hnmi = sumnmi
		
		idxpar = self._convert2idx(bestpar)
		bestsupport = [0.0] * self.numtaxa
		for par in idxpar:
			w = pmap[par]
			for idx in par:
				bestsupport[idx] = float(w)/float(len(tpartitions))
		
		spes, support = self._partition2names(bestpar, bestsupport)
		
		fo_bestpar = open(fo, "w")
		fo_bestpar.write("# MAX NMI partition found by simple heuristic search\n")
		for i in range(len(spes)):
			spe = spes[i]
			sup = support[i]
			fo_bestpar.write("Species " + str(i+1) + " (support = " + "{0:.3f}".format(sup) + ")\n")
			fo_bestpar.write("     " + self._print_list(spe) + "\n")
		fo_bestpar.close()
	
	
	def _print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + ","
		return ss[:-1] + "\n"
	
	
	def get_taxa_order(self):
		return self.taxaorder
	
	
	def _translate(self, new_taxa_order, old_partition, old_support):
		new_partition = [-1] * len(new_taxa_order)
		new_support   = [-1] * len(new_taxa_order)
		for i in range(len(self.taxaorder)):
			parnum = old_partition[i]
			sup    = old_support[i]
			taxaname = self.taxaorder[i]
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
					onepar.append(self.taxaorder[j])
					onesup.append(sup)
			nameparts.append(onepar)
			namesupps.append(onesup[0])
		
		return nameparts, namesupps



def parse_arguments():
	parser = argparse.ArgumentParser(description="""summary: a helper program to summarize multiple partitions.
Note: it only make sense to use this script if you do Bayesian analysis on a single tree.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 

Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com
Visit http://www.exelixis-lab.org/ for more information.

Version 0.2 released by Jiajie Zhang on 11-02-2014.""",
						formatter_class=argparse.RawDescriptionHelpFormatter,
						prog= "python summary.py")
	
	parser.add_argument("-p", dest = "partitions", 
						help = """Input partitions file""",
						required = True)

	parser.add_argument("-l", dest = "llhs", 
						help = """Input LLH file""",
						required = True)

	parser.add_argument("-o", dest = "output",
						help = "Output file name",
						required = True)

	parser.add_argument("-c", dest = "hpd", 
						help = """Credible interval (or Bayesian confidence interval), must be a value between 0 and 1""",
						type = float,
						required = True,
						default = 1.0)
	
	parser.add_argument("--nmi", 
						help = """Summary mutiple partitions using max NMI, this is very slow for large number of trees""",
						default = False,
						action="store_true")
	
	parser.add_argument('--version', action='version', version='%(prog)s 0.2 (11-02-2014)')
	
	return parser.parse_args()



def check_args(args):
	if not os.path.exists(args.partitions):
		print("Input partitions file does not exists: %s" % args.partitions)
		sys.exit()
	
	if not os.path.exists(args.llhs):
		print("Input LLHS file does not exists: %s" % args.llhs)
		sys.exit()



if __name__ == "__main__":
	if len(sys.argv) == 1: 
		sys.argv.append("-h")
	args = parse_arguments()
	check_args(args)
	pp = partitionparser(pfin = args.partitions, lfin = args.llhs)
	pp.summary(fout = args.output, region = args.hpd, bnmi = args.nmi)
	min_no_p, max_no_p = pp.hpd_numpartitions()
	print("Estimated number of species is between " + repr(min_no_p) + " and " + repr(max_no_p))

