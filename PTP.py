#! /usr/bin/env python
try:
	import sys
	import math
	import collections
	import ete2
	import os
	import argparse
	import subprocess
	from ete2 import Tree, SeqGroup
	from subprocess import call
	from nexus import NexusReader
	from PTPLLH import lh_ratio_test, exp_distribution, species_setting, exponential_mixture, showTree
	from summary import partitionparser
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()

class bootstrap_ptp:
	"""Run MCMC on multiple trees"""
	def __init__(self, filename, ftype = "nexus", reroot = False, method = "H1", firstktrees = 0):
		self.method = method
		self.firstktrees = firstktrees
		if ftype == "nexus":
			self.nexus = NexusReader(filename)
			self.nexus.blocks['trees'].detranslate()
			self.trees = self.nexus.trees.trees
		else:
			self.trees = self.raxmlTreeParser(filename)
		
		if self.firstktrees > 0 and self.firstktrees <= len(self.trees):
			self.trees = self.trees[:self.firstktrees]
		
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
			print("\n Somthing is wrong with the input outgroup names \n Quiting ...")
			sys.exit()
	
	
	def delimit(self, args):
		self.partitions = []
		self.settings = []
		cnt = 1
		
		if len(self.trees) == 1:
			tree = self.trees[0]
			me = None 
			if args.spe_rate <= 0:
				me = exponential_mixture(tree= tree, max_iters = args.max_iter, min_br = args.min_brl)
			else:
				me = exponential_mixture(tree= tree, max_iters = args.max_iter, min_br = args.min_brl, sp_rate = args.spe_rate, fix_sp_rate = True)
			
			if args.whiten:
				me.whitening_search(reroot = self.reroot, strategy = args.sstrategy)
			else:
				me.search(reroot = self.reroot, strategy = args.sstrategy)
			
			if args.sprint:
				me.count_species(pv = args.pvalue)
				me.print_species()
			else:
				print("Number of species: " + repr(me.count_species(pv = args.pvalue)))
			
			to, par = me.output_species(taxa_order = self.taxa_order)
			self.partitions.append(par)
			self.settings.append(me.max_setting)
			
		else:
			for tree in self.trees:
				print("Running PTP on tree " + repr(cnt) + " ........")
				cnt = cnt + 1
				me = exponential_mixture(tree= tree, max_iters = args.max_iter, min_br = args.min_brl)
				me.search(reroot = self.reroot, strategy = args.sstrategy)
				me.count_species(pv = args.pvalue, print_log = False)
				to, par = me.output_species(taxa_order = self.taxa_order)
				self.partitions.append(par)
				self.settings.append(me.max_setting)
				print("")
		return self.partitions, self.settings
	
	
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



def build_ref_tree(nfin, num_thread = "2"):
	nfout = "ptptemp"
	nfolder = os.path.dirname(os.path.abspath(nfin)) + "/"
	if os.path.exists(nfolder + nfout + ".tre"):
		print("Using existing reference tree !!")
		return nfolder + nfout + ".tre"
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234", "-T", num_thread, "-w", nfolder] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename(nfolder + "RAxML_bestTree."+nfout, nfolder + nfout + ".tre")
	os.remove(nfolder + "RAxML_info." + nfout)
	os.remove(nfolder + "RAxML_log." + nfout)
	os.remove(nfolder + "RAxML_parsimonyTree." + nfout)
	os.remove(nfolder + "RAxML_result." + nfout)
	return nfolder + nfout + ".tre"



def pick_otu(spe_out, alignment):
	fin = open(spe_out)
	lines = fin.readlines()
	fin.close()
	fout = open(alignment + ".otu", "w")
	aln = SeqGroup(sequences=alignment)
	for i in range(len(lines)):
		line = lines[i]
		if line.startswith("Species"):
			nline = lines[i+1].strip()
			seq = aln.get_seq(nline)
			fout.write(">" + nline + "\n")
			fout.write(seq + "\n")
	fout.close()



def parse_arguments():
	parser = argparse.ArgumentParser(description="""PTP: Poisson Tree Processes model for species delimitation with bootstrap support.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 

Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com.

Version 2.1 released by Jiajie Zhang on 11-02-2014.""",
						formatter_class=argparse.RawDescriptionHelpFormatter,
						prog= "python PTP.py")
	
	parser.add_argument("-t", dest = "stree",
						metavar = "TREE",
						help = """Input phylogenetic tree file. Tree can be both rooted or unrooted, 
						if unrooted, please use -r option. Supported format: NEXUS (trees without annotation),
						RAxML (simple Newick foramt). If the input file contains multiple trees, the program 
						will do bootstrap analysis""",
						required = True)
	
	parser.add_argument("-o", dest = "output",
						help = "Output file name",
						required = True)
	
	parser.add_argument("-a", dest = "salignment",
						metavar = "ALIGNMENT",
						default = "",
						help = """Specify the alignment, PTP will build a phylogenetic tree using RAxML, 
						currently only support DNA sequences with GTRGAMMA""")
	
	parser.add_argument("-u", dest = "ptpout",
						metavar = "REP-SEQUENCES",
						default = "",
						help = """Pick representative sequences and write to this file if combined with -a""")
	
	parser.add_argument("-r", dest = "sreroot",
						help = """Re-rooting the input tree on the longest branch (default not)""",
						default = False,
						action="store_true")
	
	parser.add_argument("-g", dest = "outgroups", 
						nargs='+',
						help = """Outgroup names, seperate by space. If this option is specified, 
						input tree will be rerooted accordingly""")
	
	parser.add_argument("-d", dest = "delete", 
						help = """Remove outgroups specified by -g (default not)""",
						default = False,
						action="store_true")
	
	parser.add_argument("-m", dest = "sstrategy",
						help = """Method for generate the starting partition (H0, H1, H2, H3, Brutal) (default H0)""",
						choices=["H0", "H1", "H2", "H3", "Brutal"],
						default= "H0")
	
	parser.add_argument("-pvalue", dest = "pvalue", 
						help = """Set the p-value for likelihood ratio test.(default 0.001) 
						If the test failed, the program will output only one species.
						Note this could mean there is only one species or all input taxon are different species""",
						type = float,
						default = 0.001)
						
	parser.add_argument("-p", dest = "sprint", 
						help = """Print delimited species on the screen (default not show)""",
						default = False,
						action="store_true")
	
	parser.add_argument("-w", dest = "whiten",
						help = """Specify this option to normalize the No.sequenes of each species 
						from the first run and re-run the program""",
						default = False,
						action="store_true")
	
	parser.add_argument("-minbr", dest = "min_brl", 
						help = """The minimal branch length allowed in tree (default 0.0001)""",
						type = float,
						default = 0.0001)
						
	parser.add_argument("-sprate", dest = "spe_rate", 
						help = """Fix the speciation rate to the input value during model optimization (default not fixed)""",
						type = float,
						default = -1.0)
						
	parser.add_argument("-maxiters", dest = "max_iter", 
						help = """Set the max number of search if using Brutal search (default 20000).
						The program will calculate how many searches are needed for Brutal search,
						if the number of actual search is great than this value, the program will use H0 instead""",
						type = int,
						default = 20000)
	
	parser.add_argument("-c", dest = "sscale", 
						help = """To use with -s option to set how long a branch is displayed in the plot (default 500)""",
						type = int,
						default = 500)
	
	parser.add_argument("-k", dest = "num_trees",
						metavar = "NUM-TREES",
						help = """Run bPTP on first k trees (default all trees)""",
						type = int,
						default = 0)
	
	parser.add_argument("--nmi", 
						help = """Summary mutiple partitions using max NMI, note this is very slow for large number of trees""",
						default = False,
						action="store_true")
	
	parser.add_argument('--version', action='version', version='%(prog)s 2.1 (11-02-2014)')
	
	return parser.parse_args()



def print_run_info(args):
    print("")
    print("PTP finished running with the following parameters:")
    print(" Input tree:.....................%s" % args.stree)
    print(" Search heuristic:...............%s" % args.sstrategy)
    print("")
    print(" Maximal likelihood search results written to:")
    print("  "+args.output + ".PTPPartitions.txt")
    print("")
    print(" Bootstrap values (if input contains multiple trees) of partitions written to:")
    print("  "+args.output + ".PTPPartitonSummary.txt")
    print("")
    print(" Highest bootstrap supported partition (if input contains multiple trees) written to:")
    print("  "+args.output + ".PTPhSupportPartition.txt")
    print("  Tree plot written to:")
    print("  "+args.output + ".PTPhSupportPartition.txt.png")
    print("  "+args.output + ".PTPhSupportPartition.txt.svg")
    if args.nmi:
        print("")
        print(" MAX NMI partition (if input contains multiple trees) written to:")
        print("  "+args.output + ".PTPhNMIPartition.txt")



if __name__ == "__main__":
	if len(sys.argv) == 1: 
		sys.argv.append("-h")
	args = parse_arguments()
	
	if args.ptpout!="" and args.salignment!="":
		pick_otu(spe_out = args.ptpout, alignment = args.salignment)
		sys.exit()
	
	if args.salignment!="":
		basepath = os.path.dirname(os.path.abspath(__file__))
		if not os.path.exists(basepath + "/bin/raxmlHPC-PTHREADS-SSE3"):
			print("The pipeline uses RAxML to infer phylogenetic trees,")
			print("please download the latest source code from: ")
			print("https://github.com/stamatak/standard-RAxML")
			print("Please complie the SSE + PTHREAD version, ")
			print("rename the executable to raxmlHPC-PTHREADS-SSE3 and put it to bin/  \n")
			sys.exit() 
		print("Building phylogenetic tree using RAxML.")
		args.stree = build_ref_tree(nfin = args.salignment, num_thread = "2")
	
	if not os.path.exists(args.stree):
		print("Input tree file does not exists: %s" % args.strees)
		sys.exit()
	
	try:
		tree = args.stree
		treetest = open(args.stree)
		l1 = treetest.readline()
		treetest.close()
		
		inputformat = "nexus"
		if l1.strip() == "#NEXUS":
			inputformat = "nexus"
		else:
			inputformat = "raxml"
		
		bsptp = bootstrap_ptp(filename = args.stree, ftype = inputformat, reroot = args.sreroot, method = args.sstrategy, firstktrees = args.num_trees)
		
		if args.outgroups!= None and len(args.outgroups) > 0:
			bsptp.remove_outgroups(args.outgroups, remove = args.delete)
		
		pars, settings = bsptp.delimit(args=args)
		
		pp = partitionparser(taxa_order = bsptp.taxa_order, partitions = pars, scale = args.sscale)
		pp.summary(fout = args.output, bnmi = args.nmi, sp_setting = settings)
		
		if bsptp.numtrees > 1:
			min_no_p, max_no_p, mean_no_p = pp.hpd_numpartitions()
			print("Estimated number of species is between " + repr(min_no_p) + " and " + repr(max_no_p))
			print("Mean: " + repr(mean_no_p)) 
		
		print_run_info(args = args)
		
	except ete2.parser.newick.NewickError:
		print("Unexisting tree file or Malformed newick tree structure.")


