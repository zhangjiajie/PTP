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
	import matplotlib.pyplot as plt
except ImportError:
	print("Please install the matplotlib and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()

class ptpmcmc:
	"""MCMC on a single tree using PTP model"""
	def __init__(self, tree, start_config = None, reroot = False, startmethod = "H0", min_br = 0.0001, seed = 1234, thinning = 100, sampling = 10000, burning = 0.1, taxa_order = []):
		if start_config == None:
			me = exponential_mixture(tree= tree)
			me.search(strategy = startmethod, reroot = reroot)
			me.count_species(print_log = False, pv = 0.0)
			self.tree = me.tree
			self.current_setting = me.max_setting
		else:
			self.current_setting = start_config
			self.tree = Tree(tree, format = 1)
		self.burning = burning
		self.last_setting = self.current_setting
		self.current_logl = self.current_setting.get_log_l()
		self.last_logl = self.last_setting.get_log_l()
		self.min_br = min_br
		self.rand_nr = random.Random()
		self.rand_nr.seed(seed)
		self.thinning = thinning
		self.sampling = sampling
		if taxa_order == []:
			self.taxaorder = self.tree.get_leaf_names()
		else:
			self.taxaorder = taxa_order
		self.numtaxa = len(self.taxaorder)
		self.partitions = []
		self.llhs = []
		self.nsplit = 0
		self.nmerge = 0
		"""remember the ML partition"""
		self.maxllh = self.current_logl
		to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
		self.maxpar = spe
		self.max_setting = self.current_setting
		"""record all delimitation settings for plotting, this could consume a lot of MEM"""
		self.settings = []
	
	
	def split(self, chosen_anode):
		self.nsplit = self.nsplit + 1
		newspenodes = []
		for node in self.current_setting.spe_nodes:
			newspenodes.append(node)
		newspenodes.extend(chosen_anode.get_children())
		self.current_setting = species_setting(spe_nodes = newspenodes, root = self.tree, sp_rate = 0, fix_sp_rate = False, minbr = self.min_br)
		self.current_logl = self.current_setting.get_log_l()
	
	
	def merge(self, chosen_anode):
		self.nmerge = self.nmerge + 1
		mnodes = chosen_anode.get_children()
		newspenodes = []
		for node in self.current_setting.spe_nodes:
			if not node in mnodes:
				newspenodes.append(node)
		self.current_setting = species_setting(spe_nodes = newspenodes, root = self.tree, sp_rate = 0, fix_sp_rate = False, minbr = self.min_br)
		self.current_logl = self.current_setting.get_log_l()
	
	
	def mcmc(self):
		cnt = 0
		accepted = 0
		sample_start = int(self.sampling * self.burning) 
		printinterval = self.thinning * 100
		while cnt < self.sampling:
			cnt = cnt + 1
			if cnt % printinterval == 0:
				print("MCMC generation: " + repr(cnt))
			self.last_setting = self.current_setting
			self.last_logl = self.current_logl
			acceptance = 0.0
			"""proposal"""
			"""First chose to split or merge"""
			rdchoice = self.rand_nr.uniform(0.0,1.0)
			if rdchoice <= 0.5:
				"""split"""
				xinverse = self.current_setting.get_nodes_can_split()
				if xinverse > 0:
					rdidx = self.rand_nr.randint(0, xinverse-1)
					chosen_anode = self.current_setting.node_can_split[rdidx]
					self.split(chosen_anode)
					xpinverse = self.current_setting.get_nodes_can_merge()
					if xpinverse > 0:
						newlogl = self.current_logl
						oldlogl = self.last_logl 
						acceptance = math.exp(newlogl - oldlogl) * float(xinverse)/float(xpinverse)
						if newlogl > self.maxllh:
							self.maxllh = newlogl
							to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
							self.maxpar = spe
							self.max_setting = self.current_setting
			else:
				"""merge"""
				xinverse = self.current_setting.get_nodes_can_merge()
				if xinverse > 0:
					rdidx = self.rand_nr.randint(0, xinverse-1)
					chosen_anode = self.current_setting.node_can_merge[rdidx]
					self.merge(chosen_anode)
					xpinverse = self.current_setting.get_nodes_can_split()
					if xpinverse > 0:
						newlogl = self.current_logl
						oldlogl = self.last_logl  
						acceptance = math.exp(newlogl - oldlogl) * float(xinverse)/float(xpinverse)
						if newlogl > self.maxllh:
							self.maxllh = newlogl
							to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
							self.maxpar = spe
							self.max_setting = self.current_setting
			
			if acceptance > 1.0:
				if cnt % self.thinning == 0 and cnt >= sample_start:
					to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
					self.partitions.append(spe)
					self.llhs.append(newlogl)
					self.settings.append(self.current_setting) 
				accepted = accepted + 1
			else:
				u = self.rand_nr.uniform(0.0,1.0)
				if (u < acceptance):
					if cnt % self.thinning == 0 and cnt >= sample_start:
						to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
						self.partitions.append(spe)
						self.llhs.append(newlogl)
						self.settings.append(self.current_setting) 
					accepted = accepted + 1
				else:
					self.current_setting = self.last_setting
					self.current_logl = self.last_logl
					if cnt % self.thinning == 0 and cnt >= sample_start:
						to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
						self.partitions.append(spe)
						self.llhs.append(self.current_logl)
						self.settings.append(self.current_setting)  
		
		print("Accptance rate: " + repr(float(accepted)/float(cnt)))
		print("Merge: " + repr(self.nmerge))
		print("Split: " + repr(self.nsplit))
		return self.partitions, self.llhs, self.settings



class bayesianptp:
	"""Run MCMC on multiple trees"""
	def __init__(self, filename, ftype = "nexus", reroot = False, method = "H1", seed = 1234, thinning = 100, sampling = 10000, burnin = 0.1, firstktrees = 0, taxa_order = []):
		self.method = method
		self.seed = seed
		self.thinning = thinning 
		self.sampling = sampling
		self.burnin = burnin
		self.firstktrees = firstktrees
		if ftype == "nexus":
			self.nexus = NexusReader(filename)
			self.nexus.blocks['trees'].detranslate()
			self.trees = self.nexus.trees.trees
		else:
			self.trees = self.raxmlTreeParser(filename)
		
		if self.firstktrees > 0 and self.firstktrees <= len(self.trees):
			self.trees = self.trees[:self.firstktrees]
		self.taxa_order = taxa_order
		if len(self.taxa_order) == 0:
			self.taxa_order = Tree(self.trees[0]).get_leaf_names()
		self.numtaxa = len(self.taxa_order)
		self.numtrees = len(self.trees)
		self.reroot = reroot
	
	
	def remove_outgroups(self, ognames, remove = False):
		"""reroot using outgroups and remove them"""
		self.reroot = False
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
	
	
	def delimit(self):
		self.partitions = []
		self.llhs = []
		self.settings = []
		cnt = 1 
		for tree in self.trees:
			print("Running MCMC sampling on tree " + repr(cnt) + ":")
			cnt = cnt + 1
			mcptp = ptpmcmc(tree = tree, reroot = self.reroot, startmethod = self.method, min_br = 0.0001, 
			seed = self.seed, thinning = self.thinning, sampling = self.sampling, burning = self.burnin, taxa_order = self.taxa_order)
			pars, lhs, settings = mcptp.mcmc()
			self.maxhhlpar = mcptp.maxpar
			self.maxhhlsetting = mcptp.max_setting 
			self.partitions.extend(pars)
			self.llhs.extend(lhs)
			self.settings.extend(settings)
			print("")
		return self.partitions, self.llhs, self.settings
	
	
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
	
	
	def get_maxhhl_partition(self):
		return self.maxhhlpar



def parse_arguments():
	parser = argparse.ArgumentParser(description="""bPTP: a Bayesian implementation of the PTP model for species delimitation.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 

Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com
Visit http://www.exelixis-lab.org/ for more information.

Version 0.5 released by Jiajie Zhang on 14-02-2014.""",
						formatter_class=argparse.RawDescriptionHelpFormatter,
						prog= "python bPTP.py")
	
	parser.add_argument("-t", dest = "trees", 
						help = """Input phylogenetic tree file. Trees can be both rooted or unrooted, 
						if unrooted, please use -r option. Supported format: NEXUS (trees without annotation),
						RAxML (simple Newick foramt)""",
						required = True)
	
	parser.add_argument("-o", dest = "output",
						help = "Output file name",
						required = True)
	
	parser.add_argument("-s", dest = "seed", 
						help = """MCMC seed, an integer value""",
						type = int,
						required = True)
	
	parser.add_argument("-r", dest = "reroot",
						help = """Re-rooting the input tree on the longest branch (default not)""",
						default = False,
						action="store_true")
	
	parser.add_argument("-g", dest = "outgroups", 
						nargs='+',
						help = """Outgroup names, seperate by space. If this option is specified, 
						all trees will be rerooted accordingly""")
	
	parser.add_argument("-d", dest = "delete", 
						help = """Remove outgroups specified by -g (default not)""",
						default = False,
						action="store_true")
	
	parser.add_argument("-m", dest = "method", 
						help = """Method for generate the starting partition (H0, H1, H2, H3) (default H1)""",
						choices=["H0", "H1", "H2", "H3"],
						default= "H1")
	
	parser.add_argument("-i", dest = "nmcmc", 
						help = """Number of MCMC iterations (default 10000)""",
						type = int,
						default = 10000)
	
	parser.add_argument("-n", dest = "imcmc", 
						help = """MCMC sampling interval - thinning (default 100)""",
						type = int,
						default = 100)
	
	parser.add_argument("-b", dest = "burnin", 
						help = """MCMC burn-in proportion (default 0.1)""",
						type = float,
						default = 0.1)
	
	parser.add_argument("-k", dest = "num_trees",
						metavar = "NUM-TREES",
						help = """Run bPTP on first k trees (default all trees)""",
						type = int,
						default = 0)
	
	parser.add_argument("--nmi", 
						help = """Summary mutiple partitions using max NMI, this is very slow for large number of trees""",
						default = False,
						action="store_true")
	
	parser.add_argument("--scale", 
						help = """No. pixel per unit of branch length""",
						default = 500,
						type = int)
	
	parser.add_argument('--version', action='version', version='%(prog)s 0.5 (14-02-2014)')
	
	return parser.parse_args()



def print_run_info(args, num_tree):
    print("bPTP finished running with the following parameters:")
    print(" Input tree:.....................%s" % args.trees)
    print(" MCMC iterations:................%d" % args.nmcmc)
    print(" MCMC sampling interval:.........%d" % args.imcmc)
    print(" MCMC burn-in:...................%f" % args.burnin)
    print(" MCMC seed:......................%d" % args.seed)
    print("")
    print(" MCMC samples written to:")
    print("  "+args.output + ".PTPPartitions.txt")
    print("")
    print(" Posterial LLH written to:")
    print("  "+args.output + ".PTPllh.txt")
    print("")
    print(" Posterial LLH plot:")
    print("  "+args.output + ".llh.png")
    print("")
    print(" Posterial Prob. of partitions written to:")
    print("  "+args.output + ".PTPPartitonSummary.txt")
    print("")
    print(" Highest posterial Prob. supported partition written to:")
    print("  "+args.output + ".PTPhSupportPartition.txt")
    print("  Tree plot written to:")
    print("  "+args.output + ".PTPhSupportPartition.txt.png")
    print("  "+args.output + ".PTPhSupportPartition.txt.svg")
    if args.nmi:
        print("")
        print(" MAX NMI partition written to:")
        print("  "+args.output + ".PTPhNMIPartition.txt")
    if num_tree == 1:
        print("")
        print(" Max LLH partition written to:")
        print("  "+args.output + ".PTPMLPartition.txt")
        print("  Tree plot written to:")
        print("  "+args.output + ".PTPMLPartition.txt.png")
        print("  "+args.output + ".PTPMLPartition.txt.svg")



if __name__ == "__main__":
	if len(sys.argv) == 1: 
		sys.argv.append("-h")
	args = parse_arguments()
	
	if not os.path.exists(args.trees):
		print("Input tree file does not exists: %s" % args.trees)
		sys.exit()
	
	treetest = open(args.trees)
	l1 = treetest.readline()
	treetest.close()
	inputformat = "nexus"
	if l1.strip() == "#NEXUS":
		inputformat = "nexus"
	else:
		inputformat = "raxml"
	
	bbptp = bayesianptp(filename = args.trees, ftype = inputformat, 
	reroot = args.reroot, method = args.method, seed = args.seed, 
	thinning = args.imcmc, sampling = args.nmcmc, burnin = args.burnin, 
	firstktrees = args.num_trees)
	
	if args.outgroups!= None and len(args.outgroups) > 0:
		bbptp.remove_outgroups(args.outgroups, remove = args.delete)
	
	pars, llhs, settings = bbptp.delimit()
	
	pp = partitionparser(taxa_order = bbptp.taxa_order, partitions = pars, llhs = llhs, scale = args.scale)
	
	if bbptp.numtrees == 1:
		pp.summary(fout = args.output, bnmi = args.nmi, ML_par = bbptp.get_maxhhl_partition(), 
		ml_spe_setting = bbptp.maxhhlsetting, sp_setting = settings)
	else:
		pp.summary(fout = args.output, bnmi = args.nmi, sp_setting = settings)
	
	min_no_p, max_no_p, mean_no_p = pp.hpd_numpartitions()
	print("Estimated number of species is between " + repr(min_no_p) + " and " + repr(max_no_p))
	print("Mean: " + repr(mean_no_p))
	print("")
	print_run_info(args, bbptp.numtrees)
