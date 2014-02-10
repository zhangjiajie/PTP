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
	#from collections import deque
	#from scipy import stats
	#from numpy import array
	from subprocess import call
	from nexus import NexusReader
	from PTPLLH import lh_ratio_test, exp_distribution, species_setting, exponential_mixture, showTree
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	#print(" sudo easy_install -U ete2")
	#print("Otherwise, please go to http://ete.cgenomics.org/ for instructions")
	sys.exit()



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



def remove_outgroups(t, ognames, remove = False):
	"""reroot using outgroups and remove them"""
	taxa_order = t.get_leaf_names()
	try:
		if remove:
			for og in ognames:
				taxa_order.remove(og)
		if len(ognames) < 2:
			t.set_outgroup(ognames[0])
			if remove:
				t.prune(taxa_order, preserve_branch_length=True)
		else:
			ancestor = t.get_common_ancestor(ognames)
			if not t == ancestor:
				t.set_outgroup(ancestor)
			if remove:
				t.prune(taxa_order, preserve_branch_length=True)
		return t.write()
	except ValueError, e:
		print(e)
		print("")
		print("")
		print("Somthing is wrong with the input outgroup names")
		print("")
		print("Quiting .....")
		sys.exit()



def parse_arguments():
	parser = argparse.ArgumentParser(description="""PTP: maximal likelihood search of the Poisson Tree Processes model for species delimitation.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 

Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com.

Version 1.4 released by Jiajie Zhang on 10-02-2014.""",
						formatter_class=argparse.RawDescriptionHelpFormatter,
						prog= "python PTP.py")
	
	parser.add_argument("-t", dest = "stree",
						metavar = "TREE",
						help = """Input phylogenetic tree file. Tree can be both rooted or unrooted, 
						if unrooted, please use -r option. Supported format: NEXUS (trees without annotation),
						RAxML (simple Newick foramt). If the input file contains multiple trees, only the first 
						one will be used """,
						required = True)
	
	parser.add_argument("-a", dest = "salignment",
						metavar = "ALIGNMENT",
						default = "",
						help = """Specify the alignment, PTP will build a phylogenetic tree using RAxML, 
						currently only support DNA sequences with GTRGAMMA.""")
	
	parser.add_argument("-u", dest = "ptpout",
						metavar = "REP-SEQUENCES",
						default = "",
						help = """Pick representative sequences and write to this file if combined with -a""")
	
	parser.add_argument("-r", dest = "sreroot",
						#metavar = "REROOT",
						help = """Re-rooting the input tree on the longest branch (default not).""",
						default = False,
						action="store_true")

	parser.add_argument("-g", dest = "outgroups", 
						nargs='+',
						help = """Outgroup names, seperate by space. If this option is specified, 
						input tree will be rerooted accordingly.""")
	
	parser.add_argument("-d", dest = "delete", 
						help = """Remove outgroups specified by -g (default not).""",
						default = False,
						action="store_true")

	parser.add_argument("-m", dest = "sstrategy",
						#metaval = "METHOD",
						help = """Method for generate the starting partition (H0, H1, H2, H3, Brutal) (default H1).""",
						choices=["H0", "H1", "H2", "H3", "Brutal"],
						default= "H0")

	parser.add_argument("-pvalue", dest = "pvalue", 
						help = """Set the p-value for likelihood ratio test.(default 0.001) 
						If the test failed, the program will output only one species.
						Note this could mean there is only one species or all input taxon are different species.""",
						type = float,
						default = 0.001)
						
	parser.add_argument("-p", dest = "sprint", 
						#metaval = "PRINT",
						help = """Print delimited species on the screen.(default not show)""",
						default = False,
						action="store_true")

	parser.add_argument("-s", dest = "sshow", 
						#metaval = "SHOW",
						help = """Plot delimited species on the tree.(default not show)""",
						default = False,
						action="store_true")

	parser.add_argument("-w", dest = "whiten",
						help = """Specify this option to normalize the No.sequenes of each species 
						from the first run and re-run the program""",
						default = False,
						action="store_true")

	parser.add_argument("-minbr", dest = "min_brl", 
						#metaval = "MIN-BRANCH-LEN",
						help = """The minimal branch length allowed in tree.(default 0.0001)""",
						type = float,
						default = 0.0001)
						
	parser.add_argument("-sprate", dest = "spe_rate", 
						#metaval = "SPECIATION-RATE",
						help = """Fix the speciation rate to the input value during model optimization.(default not fixed)""",
						type = float,
						default = -1.0)
						
	parser.add_argument("-maxiters", dest = "max_iter", 
						#metaval = "MAX-ITERATIONS",
						help = """Set the max number of search if using Brutal search.(default 20000)
						The program will calculate how many searches are needed for Brutal search,
						if the number of actual search is great than this value, the program will use H0 instead.""",
						type = int,
						default = 20000)

	parser.add_argument("-c", dest = "sscale", 
						#metaval = "SCALE",
						help = """To use with -s option to set how long a branch is displayed in the plot. (default 500)""",
						type = int,
						default = 500)
	
	parser.add_argument('--version', action='version', version='%(prog)s 1.4 (08-02-2014)')
	
	return parser.parse_args()



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
		stree = build_ref_tree(nfin = args.salignment, num_thread = "2")
	
	if not os.path.exists(args.stree):
		print("Input tree file does not exists: %s" % args.strees)
		sys.exit()
	
	me = None 
	try:
		tree = args.stree
		treetest = open(args.stree)
		l1 = treetest.readline()
		if l1.strip() == "#NEXUS":
			nexus = NexusReader(args.stree)
			nexus.blocks['trees'].detranslate()
			tree = nexus.trees.trees[0] 
		treetest.close()
		
		if args.outgroups!= None and len(args.outgroups) > 0:
			tree = remove_outgroups(t = Tree(tree), ognames = args.outgroups, remove = args.delete)
		
		if args.spe_rate <= 0:
			me = exponential_mixture(tree= tree, max_iters = args.max_iter, min_br = args.min_brl )
		else:
			me = exponential_mixture(tree= tree, max_iters = args.max_iter, min_br = args.min_brl, sp_rate = args.spe_rate, fix_sp_rate = True)
		
		if args.whiten:
			me.whitening_search(reroot = args.sreroot, strategy = args.sstrategy)
		else:
			me.search(reroot = args.sreroot, strategy = args.sstrategy)
		
		if args.sprint:
			me.count_species(pv = args.pvalue)
			me.print_species()
		else:
			print("Number of species: " + repr(me.count_species(pv = args.pvalue)))
		
		if args.sshow:
			showTree(delimitation = me.max_setting, scale = args.sscale)
		else:
			showTree(delimitation = me.max_setting, scale = args.sscale, render = True, fout = args.stree , form = "pdf")
			showTree(delimitation = me.max_setting, scale = args.sscale, render = True, fout = args.stree , form = "png")
	except ete2.parser.newick.NewickError:
		print("Unexisting tree file or Malformed newick tree structure.")


