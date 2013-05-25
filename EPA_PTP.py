#! /usr/bin/env python
print("This is EPA-PTP pipeline for grouping sequences into species based on a reference database.")
print("Or commonly called open-reference OTU-picking")
print("Version 1.1 released by Jiajie Zhang on 14-05-2013\n")
print("This pipeline will first use EPA to place the query sequences onto the reference tree.")
print("It will then use the PTP model or CROP to count how many species on each branch of the reference tree.")
print("The pipeline needs ETE(http://ete.cgenomics.org/) package to be installed,")
print("and will use USEARCH, RAxML, HAMMER and CROP.\n") 
print("Questions and bug reports, please send to:")
print("bestzhangjiajie@gmail.com\n")

try:
	import sys
	import math
	import collections
	import re
	import json
	import os 
	import glob
	import PTP
	import subprocess
	import random
	import ete2
	from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
	from scipy.optimize import fmin
	from collections import deque
	from scipy import stats
	from numpy import array
	from subprocess import call
except ImportError:
	print("Please install the scipy and ETE package first, and make sure the PTP.py is in the same folder.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	print(" sudo easy_install -U ete2")
	print("Otherwise, please go to http://ete.cgenomics.org/ for instructions")
	sys.exit()


def gen_alignment(seq_names = [], alignment = SeqGroup(), outputfile = "alignment.out"):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup()
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	newalign.write(outfile = outputfile)
	return outputfile, newalign


def gen_alignment2(seq_names = [], alignment = SeqGroup()):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup()
	for taxa in seq_names:
		if taxa.startswith("*R*"):
			seq = alignment.get_seq(taxa[3:])
		elif taxa == "sister":
			continue
		else:
			seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	#newalign.write(outfile = outputfile)
	return newalign


def gen_alignment3(seq_names = [], alignment = SeqGroup()):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup()
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	#newalign.write(outfile = outputfile)
	return newalign


def catalns(refs, alns, sfout):
	cnt = 0
	fout = open(sfout, "w")
	for aln in alns:
		if aln.startswith(">"):
			cnt = cnt + 1
		fout.write(aln)
	#fout.write("\n")
	for ref in refs:
		if ref.startswith(">"):
			cnt = cnt + 1
		fout.write(ref)
	fout.close()
	return sfout, cnt


"""Input ete alignment object, and output representative sequence file; return the logs"""
def count_and_pick_reads(align, outputfile):
	logss = ""
	numreads = 0
	l_reads = None
	max_non_gap = 0
	entries = align.get_entries()
	if len(entries) == 1:
		name = entries[0][0]
		if name.startswith("*R*") or name == "sister":
			return ""
	for entr in align.get_entries():
		name = entr[0]
		seq = entr[1]
		
		if name.startswith("*R*"):
			logss = "R	find reference species: " + name + "\n"
		elif name == "sister":
			pass
		else:
			names = name.split("*")[-1]
			if len(names) > 1:
				try:
					numseqs = int(names[-1])
				except ValueError:
					numseqs = 1
					print("Unable to get the reads count for: " + name)
			else:
				numseqs = 1
			numreads = numreads + numseqs
		if name != "sister":
			seql = count_non_gap(seq)
			if seql > max_non_gap:
				l_reads = entr
				max_non_gap = seql
	if logss == "":
		logss = "D	find new species \n"
	logss = logss + "K	reads number: " + repr(numreads) + "\n"
	if l_reads!=None:
		fout = open(outputfile, "a")
		fout.write(">" + l_reads[0] + "*" + repr(numreads) + "\n")
		fout.write(l_reads[1] + "\n")
		fout.close()
	return logss


def remove_seq_len_smaller_than(f_fasta, min_l):
	fin = open(f_fasta)
	fout = open(f_fasta+".min"+repr(min_l)+".afa", "w")
	line = fin.readline().strip()
	lastname = ""
	while line!="":
		if line.startswith(">"):
			line = line.split()[0]
			lastname = line
			#fout.write(line + "\n")
		else:
			if len(line) >= min_l:
				fout.write(lastname + "\n")
				fout.write(line + "\n")
		line = fin.readline().strip()
	fout.close()


def collapse_identical_seqs(f_fasta):
	if os.path.exists(f_fasta+".collapse.afa"):
		print("Using existing collapsed sequences !!")
		return f_fasta+".collapse.afa"
	
	fin = open(f_fasta)
	fout = open(f_fasta+".collapse.afa", "w")
	seq_list=[]
	line = fin.readline().strip()
	while line!="":
		name = line
		seq = fin.readline().strip()
		flag = True
		for oneseq in seq_list:
			if oneseq[2] == seq:
				oneseq[1] = oneseq[1] + 1
				flag = False
				break
		if flag:
			seq_list.append([name, 1, seq])
		line = fin.readline().strip()
	for oneseq in seq_list:
		fout.write(oneseq[0]+"*"+repr(oneseq[1])+"\n")
		fout.write(oneseq[2]+"\n")
	
	fin.close()
	fout.close()
	return f_fasta+".collapse.afa"


def processHMMseq(seqin):
	newseq = ""
	for s in seqin:
		if s == ".":
			pass
		elif s == "-":
			newseq = newseq + s
		elif s.isupper():
			newseq = newseq + s
	return newseq


def count_non_gap(seqin):
	cnt = 0
	for s in seqin:
		if s!="-":
			cnt = cnt + 1
	return cnt


def parse_HMM(f_stock, l_ref, minl = 50):
	cnt = 0
	fin = open(f_stock)
	line = fin.readline()
	seqs = {}
	while line!="":
		if line.startswith("//"):
			break
		elif line.startswith("#"):
			pass
		elif line.startswith("\n"):
			cnt = cnt + 1
		else:
			line = line.strip()
			if cnt == 1:
				l2 = line.split()
				ss = processHMMseq(l2[1])
				seqs[l2[0]] = ss
			else:
				l2 = line.split()
				seq = seqs[l2[0]]
				ss = processHMMseq(l2[1])
				seqs[l2[0]] = seq + ss 
		line = fin.readline()
	fin.close()
	fout = open(f_stock+".afa", "w")
	for key in seqs.keys():
		if count_non_gap(seqs[key]) >= minl:
			fout.write(">" + key + "\n")
			seq = seqs[key]
			lappd = l_ref - len(seq)
			if lappd > 0:
				appd = "-" * lappd
				seq = seq + appd
			elif lappd < 0:
				print("Warning: query sequence > ref sequence")
			
			fout.write(seq + "\n")
	fout.close()
	return f_stock+".afa"


def chimera_removal(nuseach, nalign, nout, chimeraout):
	align = SeqGroup(nalign)
	newalign = open(nout, "w")
	chalign = open(chimeraout, "w")
	fus = open(nuseach)
	lines = fus.readlines()
	fus.close()
	for line in lines:
		its = line.split()
		c = its[-1]
		sname = its[1]
		if c == "Y" or c =="?":
			seq = align.get_seq(sname)
			chalign.write(">" + sname + "\n")
			chalign.write(seq + "\n")
		else:
			seq = align.get_seq(sname)
			newalign.write(">" + sname + "\n")
			newalign.write(seq + "\n")
	newalign.close()
	chalign.close()


#correct rooting method, this shound run after EPA
def extract_placement(nfin_place, nfin_aln, nfout, min_lw = 0.5, logfile = "spcount.log"):
	if os.path.exists(logfile):
		os.remove(logfile)
	
	if os.path.exists(nfout + "_inode_picked_otus.fasta"):
		os.remove(nfout + "_inode_picked_otus.fasta")
	
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln)
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
	allnodes = root.get_descendants()
	allnodes.append(root)
	
	"""get refseq"""
	refseqset = []
	for leaf in leaves:
		refseqset.append(leaf.name)
	refali = gen_alignment2(seq_names = refseqset, alignment = align_orgin)
	
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			placemap[curredge] = placemap.get(curredge, [])
	
	"""placement quality control"""
	discard_file = open(nfout+".discard.placement.txt", "w")
	"""group taxa to edges"""
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			a = placemap[curredge] 
			a.extend(taxa_names)
			placemap[curredge]  = a
		else:
			discard_file.write(repr(taxa_names) + "\n")
	discard_file.close()
	
	groups = placemap.items()
	cnt_leaf = 0
	cnt_inode = 0
	
	"""check each edge""" 
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		
		"""check if placed on leaf node and find the node being placed on"""
		flag = False
		place_node = None
		for node in allnodes:
			if str(node.B) == str(seqset_name):
				place_node = node
				if node.is_leaf():
					flag = True 
				break
				
		"""find the furthest leaf of the placement node"""
		fnode = place_node.get_farthest_node()[0]
		outgroup_name = fnode.name
		
		"""find sister node"""
		snode = place_node.get_sisters()[0]
		if not snode.is_leaf():
			snode = snode.get_closest_leaf()[0]
		sister_name = snode.name
		
		"""generate aligment"""
		if flag:
			"""process leaf node placement"""
			cnt_leaf = cnt_leaf + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			if len(newalign.get_entries()) < 2:
				#count_and_pick_reads(align = newalign, outputfile = nfout + "_leaf_picked_otus.fasta")
				og_seq = align_orgin.get_seq(outgroup_name)
				sis_seq = align_orgin.get_seq(sister_name)
				newalign.set_seq("sister", sis_seq) #set the sister seqeunce to make 4 taxa
				newalign.set_seq("root_ref", og_seq) #set the outgroup name
				place_seq = align_orgin.get_seq(place_node.name)
				newalign.set_seq("*R*" + place_node.name, place_seq) #set the reference sequence name
				newalign.write(outfile = nfout + "_leaf_"+repr(cnt_leaf) + ".lfa")
			else:
				og_seq = align_orgin.get_seq(outgroup_name)
				newalign.set_seq("root_ref", og_seq) #set the outgroup name
				place_seq = align_orgin.get_seq(place_node.name)
				newalign.set_seq("*R*" + place_node.name, place_seq) #set the reference sequence name
				newalign.write(outfile = nfout + "_leaf_"+repr(cnt_leaf) + ".lfa")
		else:
			"""genrate the newwick string to be inserted into the ref tree"""
			rep = re.compile(r"\{[0-9]*\}")
			multi_fcating = "("
			for seqname in seqset:
				multi_fcating = multi_fcating + seqname + ","
			multi_fcating = multi_fcating[:-1] 
			multi_fcating = "{" + repr(seqset_name) + "}," + multi_fcating + ")"
			mtfc_tree = tree.replace("{" + repr(seqset_name) + "}", multi_fcating)
			mtfc_tree = rep.sub("", mtfc_tree)
			
			cnt_inode = cnt_inode + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			if len(newalign.get_entries()) < 2:
				count_and_pick_reads(align = newalign, outputfile = nfout + "_inode_picked_otus.fasta")
				sp_log(sfout = logfile, logs="I	the palcement is on an internal node \nD	find new species\nK	reads number: 1 \n")
			else:
				#og_seq = align_orgin.get_seq(outgroup_name)
				#newalign.set_seq("root_ref", og_seq)
				for entr in refali.get_entries():
					sname = entr[0]
					seqe = entr[1]
					newalign.set_seq(sname, seq)
				newalign.write(outfile = nfout + "_inode_"+repr(cnt_inode) + ".ifa")
				mtfc_out = open(nfout + "_inode_"+repr(cnt_inode) +  ".mttree", "w")
				mtfc_out.write(mtfc_tree)
				mtfc_out.close()


def extract_placement_crop(nfin_place, nfin_aln, nfout, min_lw = 0.5, logfile = "spcount.log"):
	if os.path.exists(logfile):
		os.remove(logfile)
	
	if os.path.exists(nfout + "_inode_picked_otus.fasta"):
		os.remove(nfout + "_inode_picked_otus.fasta")
	
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln)
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
	allnodes = root.get_descendants()
	allnodes.append(root)
	
	"""get refseq"""
	refseqset = []
	for leaf in leaves:
		refseqset.append(leaf.name)
	refali = gen_alignment2(seq_names = refseqset, alignment = align_orgin)
	
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			placemap[curredge] = placemap.get(curredge, [])
	
	"""placement quality control"""
	discard_file = open(nfout+".discard.placement.txt", "w")
	"""group taxa to edges"""
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			a = placemap[curredge] 
			a.extend(taxa_names)
			placemap[curredge]  = a
		else:
			discard_file.write(repr(taxa_names) + "\n")
	discard_file.close()
	
	groups = placemap.items()
	cnt_leaf = 0
	cnt_inode = 0
	
	"""check each edge""" 
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		
		"""check if placed on leaf node and find the node being placed on"""
		flag = False
		place_node = None
		for node in allnodes:
			if str(node.B) == str(seqset_name):
				place_node = node
				if node.is_leaf():
					flag = True 
				break
		
		"""generate aligment"""
		if flag:
			"""process leaf node placement"""
			cnt_leaf = cnt_leaf + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			place_seq = align_orgin.get_seq(place_node.name)
			newalign.set_seq("*R*" + place_node.name, place_seq) #set the reference sequence name
			newalign.write(outfile = nfout + "_leaf_"+repr(cnt_leaf) + ".lfa")
		else:
			cnt_inode = cnt_inode + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			
			if len(newalign.get_entries()) < 2:
				count_and_pick_reads(align = newalign, outputfile = nfout + "_inode_picked_otus.fasta")
				sp_log(sfout = logfile, logs="I	the palcement is on an internal node \nD	find new species\nK	reads number: 1 \n")
			else:
				#for entr in refali.get_entries():
				#	sname = entr[0]
				#	seqe = entr[1]
				#	newalign.set_seq(sname, seq)
				newalign.write(outfile = nfout + "_inode_"+repr(cnt_inode) + ".ifa")


#build tree with -g
def build_constrain_tree(nsfin, ntfin, nfout, nfolder, num_thread = "2"):
	if os.path.exists(nfolder + nfout + ".tre"):
		print("Using existing constrain tree !!")
		return nfolder + nfout + ".tre"
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n",nfout,"-p", "1234", "-T", num_thread, "-w", nfolder], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename(nfolder + "RAxML_bestTree."+nfout, nfolder + nfout + ".tre")
	os.remove(nfolder + "RAxML_info." + nfout)
	os.remove(nfolder + "RAxML_log." + nfout)
	os.remove(nfolder + "RAxML_result." + nfout)
	return nfolder + nfout + ".tre"


#build a tree
def build_ref_tree(nfin, nfout, nfolder, num_thread = "2"):
	if os.path.exists(nfolder + nfout + ".tre"):
		print("Using existing reference tree !!")
		return nfolder + nfout + ".tre"
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234", "-T", num_thread, "-w", nfolder], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename(nfolder + "RAxML_bestTree."+nfout, nfolder + nfout + ".tre")
	os.remove(nfolder + "RAxML_info." + nfout)
	os.remove(nfolder + "RAxML_log." + nfout)
	os.remove(nfolder + "RAxML_parsimonyTree." + nfout)
	os.remove(nfolder + "RAxML_result." + nfout)
	return nfolder + nfout + ".tre"


def compare_node(node):
	return node.dist


def find_lonest_br(tree):
	node_list = tree.get_descendants()
	node_list.sort(key=compare_node)
	node_list.reverse()
	rootnode = node_list[0]
	return rootnode


"""for leaf node placement: build the phylogenetic tree for ME, extract the subtree - rooting and remove root_ref"""
def raxml_after_epa(nfolder, suf = "lfa", T = "2"):
	naligns = glob.glob(nfolder + "*." + suf)
	cnt = 0
	for aln in naligns:
		print(repr(cnt))
		cnt = cnt + 1
		if os.path.exists(aln.split(".")[0] + ".subtree"):
			pass
		else:
			trename = build_ref_tree(nfin = aln, nfout = "l"+repr(cnt), nfolder = nfolder, num_thread = T)
			full_tree = Tree(trename, format=1)
			rootref = full_tree.get_leaves_by_name("root_ref")[0]
			if rootref.up.is_root():
				newrootnode = rootref.get_farthest_node()[0]
				full_tree.set_outgroup(newrootnode)
			
			rootref = full_tree.get_leaves_by_name("root_ref")[0]
			refroot_brl = rootref.dist
			full_tree.set_outgroup(rootref)
			real_tree = None
			
			for child in full_tree.get_children():
				if not child.is_leaf():
					real_tree = child
					real_tree.up = None
					real_tree.dist = 0.0
					break
			
			lnode = find_lonest_br(real_tree)
			if lnode.dist > refroot_brl:
				real_tree.set_outgroup(lnode)
				real_tree.dist = 0.0
			
			real_tree.write(outfile= aln.split(".")[0] + ".subtree", format=5)
			os.remove(trename)


#build the phylogenetic tree -g option
def raxml_g_after_epa(nfolder, nref_align, suf = "ifa", T = "2"):
	align_orgin = SeqGroup(sequences = nref_align)
	ref_taxa = []
	for entr in align_orgin.get_entries():
		ref_taxa.append(entr[0])
	
	naligns = glob.glob(nfolder + "*." + suf)
	cnt = 0
	for aln in naligns:
		print(repr(cnt))
		cnt = cnt + 1
		if os.path.exists(aln.split(".")[0] + ".subtree"):
			pass
		else:
			mttree = aln.split(".")[0] + ".mttree"
			#raxml constrait search
			trename = build_constrain_tree(nsfin = aln, ntfin = mttree, nfout = "i"+repr(cnt), nfolder = nfolder, num_thread = T)
			#read in the fully resolved tree
			full_tree = Tree(trename, format=1)
			all_taxa = full_tree.get_leaf_names()
			target_taxa = []
			for taxa in all_taxa:
				if taxa in ref_taxa:
					pass
				else:
					target_taxa.append(taxa)
			#the place where the tree can be safely rooted
			ref_node = full_tree.get_leaves_by_name(ref_taxa[0])[0]
			#reroot 
			full_tree.set_outgroup(ref_node)
			#find the common ancestor of the target taxa
			leafA = full_tree.get_leaves_by_name(target_taxa[0])[0]
			leaflist = []
			for n in target_taxa[1:]:
				leaflist.append(full_tree.get_leaves_by_name(n)[0])
			common = leafA.get_common_ancestor(leaflist)
			common.up = None
			common.write(outfile= aln.split(".")[0] + ".subtree", format=5)
			os.remove(trename)
			os.remove(mttree)


def subtrees(nfolder, pref = "RAxML_bestTree"):
	ntrees = glob.glob(nfolder + pref + "*")
	for tree in ntrees:
		#print tree
		if tree.split("/")[-1].startswith("RAxML_bestTree.me_leaf"):
			full_tree = Tree(tree, format=1)
			rootref = full_tree.get_leaves_by_name("root_ref")[0]
			if rootref.up.is_root():
				newrootnode = rootref.get_farthest_node()[0]
				full_tree.set_outgroup(newrootnode)
			
			rootref = full_tree.get_leaves_by_name("root_ref")[0]
			refroot_brl = rootref.dist
			full_tree.set_outgroup(rootref)
			real_tree = None
			
			for child in full_tree.get_children():
				if not child.is_leaf():
					real_tree = child
					real_tree.up = None
					real_tree.dist = 0.0
					break
		
			lnode = find_lonest_br(real_tree)
			if lnode.dist > refroot_brl:
				real_tree.set_outgroup(lnode)
				real_tree.dist = 0.0
			real_tree.write(outfile= nfolder + tree.split(".")[-2] + ".subtree", format=5)


def estimate_ref_exp_rate(nfin):
	ref_model = PTP.exponential_mixture(tree = nfin)
	spe_rate = ref_model.null_model()
	return spe_rate


def sp_log(sfout, logs=""):
	f = open(sfout, "a")
	f.write(logs)
	f.write("\n")
	f.close()


def otu_picking(nfolder, nfout1, nfout2, nref_tree, n_align, suf = "subtree", pvalue = 0.001):
	"""T, tree file; M, search method; N, num cpecies; L, place on leaf; I, place on internal node; R, find reference species; D, find denovo specise; K, read number"""
	trees = glob.glob(nfolder + "*." + suf)
	spe_rate = estimate_ref_exp_rate(nref_tree)
	align = SeqGroup(sequences = n_align)
	for tree in trees:
		print(tree)
		logss = ""
		logss = logss + "T	Searching species in tree: " + tree + ":\n"
		epa_exp = PTP.exponential_mixture(tree, sp_rate = spe_rate, fix_sp_rate = True)
		t = Tree(tree, format = 1)
		tsize = len(t.get_leaves())
		if tsize > 500:
			epa_exp.search(reroot = False, strategy = "H1")
			logss = logss + "M	H1\n"
		elif tsize < 20:
			epa_exp.search(reroot = False, strategy = "Brutal")
			logss = logss + "M	Brutal\n"
		else:
			epa_exp.search(reroot = False, strategy = "H0")
			logss = logss + "M	H0\n"
		num_spe = epa_exp.count_species(print_log = False, pv = pvalue)
		
		logss = logss + "N	find number specices: " + repr(num_spe) + "\n"
		
		idx = tree.find("leaf")
		if idx >= 0:
			logss = logss + "L	the palcement is on a leaf node" + "\n"
		else:
			logss = logss + "I	the palcement is on an internal node" + "\n"
		for spe in epa_exp.species_list:
			newalign = gen_alignment2(seq_names = spe, alignment = align)
			
			if idx >= 0:
				morelog = count_and_pick_reads(newalign, nfout1)
				logss = logss + morelog
			else:
				morelog = count_and_pick_reads(newalign, nfout2)
				logss = logss + morelog
			
		sp_log(sfout = nfolder + "spcount.log", logs = logss)


def crop_otu_picking(nfolder, nfout1, nfout2, nref_tree, n_align):
	lplacement = glob.glob(nfolder + "*.lfa")
	iplacement = glob.glob(nfolder + "*.ifa")
	align = SeqGroup(sequences = n_align)
	for laln in lplacement:
		logss = ""
		logss = logss + "T	Searching species in alignment: " + laln + ":\n"
		spes = crop_species_counting(laln)
		logss = logss + "N	find number specices: " + repr(len(spes)) + "\n"
		logss = logss + "L	the palcement is on a leaf node" + "\n"
		for spe in spes:
			newalign = gen_alignment2(seq_names = spe, alignment = align)
			morelog = count_and_pick_reads(newalign, nfout1)
			logss = logss + morelog
		sp_log(sfout = nfolder + "spcount.log", logs = logss)
	
	for laln in iplacement:
		logss = ""
		logss = logss + "T	Searching species in alignment: " + laln + ":\n"
		spes = crop_species_counting(laln)
		logss = logss + "N	find number specices: " + repr(len(spes)) + "\n"
		logss = logss + "I	the palcement is on an internal node" + "\n"
		for spe in spes:
			newalign = gen_alignment2(seq_names = spe, alignment = align)
			morelog = count_and_pick_reads(newalign, nfout2)
			logss = logss + morelog
		sp_log(sfout = nfolder + "spcount.log", logs = logss)


def auto_stas(sfin):
	f = open(sfin)
	lines = f.readlines()
	f.close()
	min_reads_num_match = sys.maxsize
	
	for k in range(len(lines)):
		line = lines[k]
		if line.startswith("R"):
			nline = lines[k+1]
			numreads = int(nline.split(":")[-1])
			if numreads < min_reads_num_match:
				min_reads_num_match = numreads
	
	"""
	last_min_reads_num_match = -1
	
	min_reads = 0
	
	for i in range(1000):
		if min_reads_num_match < last_min_reads_num_match:
			min_reads = i - 1 
			break
		else:
			last_min_reads_num_match = min_reads_num_match
			min_reads_num_match = 0
		
		for k in range(len(lines)):
			line = lines[k]
			if line.startswith("R"):
				nline = lines[k+1]
				numreads = int(nline.split(":")[-1])
				if numreads > i:
					min_reads_num_match = min_reads_num_match + 1
	"""
	
	match1 = 0 
	match2 = 0
	match3 = 0
	nomatch1 = 0
	otu1 = 0
	
	for k in range(len(lines)):
		line = lines[k]
		if line.startswith("R"):
			nline = lines[k+1]
			numreads = int(nline.split(":")[-1])
			if numreads >= min_reads_num_match:
				match1 = match1 + 1
				otu1 = otu1 + 1
			
			if numreads >= (min_reads_num_match-1):
				match2 = match2 + 1
			
			if numreads >= (min_reads_num_match+1):
				match3 = match3 + 1
			
		if line.startswith("D"):
			nline = lines[k+1]
			numreads = int(nline.split(":")[-1])
			if numreads >= min_reads_num_match:
				nomatch1 = nomatch1 + 1
				otu1 = otu1 + 1
	print("")
	print("The minimal sequence number of reference-matched cluster is: " + repr(min_reads_num_match) + ". Using this value as a cutoff for confident OTUs:")
	print("INFO: >=" + repr(min_reads_num_match - 1) + " sequences match species -   " + repr(match2))
	print("INFO: >=" + repr(min_reads_num_match) + " sequences match species -   " + repr(match1))
	print("INFO: >=" + repr(min_reads_num_match + 1) + " sequences match species -   " + repr(match3))
	print("")
	
	print(">=" + repr(min_reads_num_match) + " sequences species -         " + repr(otu1))
	print(">=" + repr(min_reads_num_match) + " sequences match species -   " + repr(match1))
	print(">=" + repr(min_reads_num_match) + " sequences new species -     " + repr(nomatch1))



def stas(sfin):
	"""T, tree file; M, search method; N, num cpecies; L, place on leaf; I, place on internal node; R, find reference species; D, find denovo specise; K, read number"""
	otu1 = 0
	otu2 = 0 
	otu3 = 0
	otu4 = 0
	otu5 = 0
	match1 = 0
	match2 = 0 
	match3 = 0 
	match4 = 0 
	match5 = 0
	nomatch1 = 0
	nomatch2 = 0
	nomatch3 = 0
	nomatch4 = 0
	nomatch5 =0
	epa = 0
	
	f = open(sfin)
	l = f.readline()
	while l!="":
		if l.startswith("T"):
			epa = epa + 1
		if l.startswith("R"):
			l = f.readline()
			numreads = int(l.split(":")[-1])
			if numreads > 0:
				otu1 = otu1 + 1
				match1 = match1 + 1
			if numreads > 1:
				otu2 = otu2 + 1
				match2 = match2 + 1
			if numreads > 2:
				otu3 = otu3 + 1
				match3 = match3 + 1
			if numreads > 3:
				otu4 = otu4 + 1
				match4 = match4 + 1
			if numreads > 4:
				otu5 = otu5 + 1
				match5 = match5 + 1
		if l.startswith("D"):
			l = f.readline()
			numreads = int(l.split(":")[-1])
			if numreads > 0:
				otu1 = otu1 + 1
				nomatch1 = nomatch1 + 1
			if numreads > 1:
				otu2 = otu2 + 1
				nomatch2 = nomatch2 + 1
			if numreads > 2:
				otu3 = otu3 + 1
				nomatch3 = nomatch3 + 1
			if numreads > 3:
				otu4 = otu4 + 1
				nomatch4 = nomatch4 + 1
			if numreads > 4:
				otu5 = otu5 + 1
				nomatch5 = nomatch5 + 1
		l = f.readline()
	f.close()
	print("Species delimitation summary for: " + sfin + "\n")
	print("No. branches with placements   " + repr(epa))
	print("")
	print(">=5 sequences species -        " + repr(otu5))
	print(">=5 sequences match species:   " + repr(match5))
	print(">=5 sequences new species:     " + repr(nomatch5))
	print("")
	print(">=4 sequences species -        " + repr(otu4))
	print(">=4 sequences match species:   " + repr(match4))
	print(">=4 sequences new species:     " + repr(nomatch4))
	print("")
	print(">=3 sequences species -        " + repr(otu3))
	print(">=3 sequences match species:   " + repr(match3))
	print(">=3 sequences new species:     " + repr(nomatch3))
	print("")
	print(">=2 sequences species -        " + repr(otu2))
	print(">=2 sequences match species:   " + repr(match2))
	print(">=2 sequences new species:     " + repr(nomatch2))
	print("")
	print(">=1 sequences species -        " + repr(otu1))
	print(">=1 sequences match species:   " + repr(match1))
	print(">=1 sequences new species:     " + repr(nomatch1))
	
	auto_stas(sfin)


def random_remove_taxa(falign, num_remove, num_repeat = 1):
	align = SeqGroup(sequences = falign)
	entrs = align.get_entries()
	numseq = len(entrs)
	index = range(numseq)
	namel = []
	
	for i in range(num_repeat):
		newalign = SeqGroup()
		random.shuffle(index)
		idxs = index[num_remove:]
		for idx in idxs:
			newalign.set_seq(entrs[idx][0], entrs[idx][1])
		newalign.write(outfile = falign + "_" + repr(num_remove)+ "_" + repr(i + 1) + ".afa")
		namel.append(falign + "_" + repr(num_remove)+ "_" + repr(i + 1) + ".afa")
	return namel


def count_reads(nfolder, pref = "me_leaf_"):
	cnt = 0
	naligns = glob.glob(nfolder + pref + "*")
	for aln in naligns:
		a = SeqGroup(sequences = aln)
		for ent in a.get_entries():
			name = ent[0]
			if name == "root_ref":
				pass
			elif name.startswith("*R*"):
				pass
			else:
				numread = int(name.split("*")[-1])
				cnt = cnt + numread
	print cnt


def build_hmm_profile(faln, fbase=""):
	#hmmbuild --informat afa refotu.hmm ref_outs_547.fas
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/hmmbuild","--informat", "afa", faln+".hmm", faln]) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	return faln+".hmm"


def hmm_align(fprofile, ffasta, fbase=""):
	#hmmalign -o 454.stock refotu.hmm 454input.fna.min100.fasta
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/hmmalign","-o", ffasta + ".stock", fprofile, ffasta]) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	return ffasta + ".stock"


def trim_refalign_hmm(refaln, hmmprofile):
	sites = []
	hmp = open(hmmprofile)
	l = hmp.readline()
	start = False
	while l!="":
		if l.startswith("//"):
			break
		if start:
			l = l.strip()
			ll = l.split()
			usedsite = int(ll[5])
			sites.append(usedsite)
			l = hmp.readline()
			l = hmp.readline()
		else:
			if l.startswith("HMM "):
				start = True
				l = hmp.readline()
				l = hmp.readline()
				l = hmp.readline()
				l = hmp.readline()
		l = hmp.readline()
	hmp.close()
	align = SeqGroup(refaln)
	fout = open(refaln+".trimed.afa", "w")
	for entr in align.get_entries():
		fout.write(">" + entr[0] + "\n")
		for pos in sites:
			fout.write(entr[1][pos-1])
		fout.write("\n")
	fout.close()
	return refaln+".trimed.afa", len(sites)


def hmm_alignment(ref_align, query, outfolder, lmin = 100, outname = "epa_ready"):
	hmmprofile = build_hmm_profile(faln = ref_align)
	stock = hmm_align(fprofile = hmmprofile, ffasta = query)
	trimaln, ref_len = trim_refalign_hmm(ref_align, hmmprofile)
	afa = parse_HMM(stock, l_ref = ref_len, minl = lmin)
	
	if outname.find("/") >=0:
		os.rename(trimaln, outname + ".ref.afa")
		os.rename(afa, outname + ".query.afa")
	else:
		os.rename(trimaln, outfolder + outname + ".ref.afa")
		os.rename(afa, outfolder + outname + ".query.afa")


def epa_me_species_counting(refaln, queryaln, folder, lw = 0.2, T = "2", pvalue = 0.001):
	"""input reference alignment and alinged query sequences"""
	print("Building refrence tree")
	ref_tree = build_ref_tree(nfin = refaln, nfout = queryaln.split("/")[-1], nfolder = folder, num_thread = T)
	af = open(refaln)
	aln = af.readlines()
	af.close()
	print("Collapsing identical sequences")
	cqali = collapse_identical_seqs(queryaln)
	chimera_free = run_uchime(sref = refaln, squery = cqali)
	bf = open(chimera_free)
	bln = bf.readlines()
	bf.close()
	catalns(bln, aln, queryaln+".epainput")
	
	if os.path.exists(chimera_free):
		os.remove(chimera_free)
	if os.path.exists(cqali):
		os.remove(cqali)
	print("Running epa")
	fplacement = run_epa(query = queryaln+".epainput", reftree = ref_tree, folder = folder, num_thread = T)
	
	extract_placement(nfin_place = fplacement, nfin_aln = queryaln+".epainput", nfout = folder + "me", min_lw = lw, logfile = "spcount.log")
	print("Building trees for epa results:")
	raxml_after_epa(nfolder = folder,suf = "lfa", T = T)
	raxml_g_after_epa(nfolder = folder, nref_align = refaln, suf = "ifa", T = T)
	print("OTU picking:")
	otu_picking(nfolder = folder, nfout1 = folder + "me_leaf_picked_otus.fasta"  , nfout2 = folder + "me_inode_picked_otus.fasta" , nref_tree = ref_tree, n_align = queryaln+".epainput", suf = "subtree", pvalue = pvalue)
	
	clean(sfolder = folder)
	return queryaln+".epainput", ref_tree, fplacement


def crop_species_counting(falin, numseq = -1):
	fafa = falin
	basepath = os.path.dirname(os.path.abspath(__file__))
	if numseq <= 0:
		call([basepath + "/bin/crop", "-i", fafa], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	else:
		b = numseq / 50 
		b = int(b)
		if b < 1:
			b = 1
		call([basepath + "/bin/crop", "-i", fafa, "-b", repr(b)], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
		
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


def epa_crop_species_counting(refaln, queryaln, folder, lw = 0.2, T = "2"):
	"""input reference alignment and alinged query sequences"""
	print("Building refrence tree")
	ref_tree = build_ref_tree(nfin = refaln, nfout = queryaln.split("/")[-1], nfolder = folder, num_thread = T)
	af = open(refaln)
	aln = af.readlines()
	af.close()
	print("Collapsing identical sequences")
	cqali = collapse_identical_seqs(queryaln)
	chimera_free = run_uchime(sref = refaln, squery = cqali)
	bf = open(chimera_free)
	bln = bf.readlines()
	bf.close()
	catalns(bln, aln, queryaln+".epainput")
	os.remove(chimera_free)
	os.remove(cqali)
	print("Running epa")
	fplacement = run_epa(query = queryaln+".epainput", reftree = ref_tree, folder = folder, num_thread = T)
	extract_placement_crop(nfin_place = fplacement, nfin_aln = queryaln+".epainput", nfout = folder + "me", min_lw = lw, logfile = "spcount.log")
	print("OTU picking:")
	crop_otu_picking(nfolder = folder, nfout1 = folder + "me_leaf_picked_otus.fasta"  , nfout2 = folder + "me_inode_picked_otus.fasta" , nref_tree = ref_tree, n_align = queryaln+".epainput")
	clean(sfolder = folder)


def crop_stand_alone(refaln, queryaln, folder):
	af = open(refaln)
	aln = af.readlines()
	af.close()
	aln2 = []
	for line in aln:
		if line.startswith(">"):
			line = ">*R*" + line[1:]
			aln2.append(line)
		else:
			aln2.append(line)
	print("Collapsing identical sequences")
	cqali = collapse_identical_seqs(queryaln)
	chimera_free = run_uchime(sref = refaln, squery = cqali)
	bf = open(chimera_free)
	bln = bf.readlines()
	bf.close()
	cropinput, numseq = catalns(bln, aln2, queryaln+".cropinput")
	os.remove(chimera_free)
	os.remove(cqali)
	print("Running CROP")
	spes = crop_species_counting(cropinput, numseq)
	align = SeqGroup(sequences = cropinput)
	for spe in spes:
		newalign = gen_alignment3(seq_names = spe, alignment = align)
		morelog = count_and_pick_reads(newalign, queryaln + "crop_outs.txt")
		sp_log(sfout = folder + "spcount.log", logs = morelog)


def uchime_ready(sfin):
	fin = open(sfin)
	lines = fin.readlines()
	fin.close()
	fout = open(sfin+".uchime", "w")
	for line in lines:
		line = line.replace("-", "")
		fout.write(line)
	fout.close()
	return sfin+".uchime"


def run_uchime(sref, squery, fbase = ""):
	newref = uchime_ready(sref)
	newquery = uchime_ready(squery)
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/usearch","-uchime_ref", newquery, "-db", newref, "-uchimeout", squery + ".uchimeout", "-strand", "plus"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	chimera_removal(nuseach = squery + ".uchimeout", nalign = squery, nout = squery + ".chimerafree", chimeraout = squery + ".chimera")
	os.remove(newref)
	os.remove(newquery)
	return squery + ".chimerafree"


def run_epa(query, reftree, folder, num_thread = "2", binbase = ""):
	if os.path.exists(folder + query.split("/")[-1] +".jplace"):
		print("Using existing EPA results !!")
		return folder + query.split("/")[-1] +".jplace"
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",query, "-r", reftree, "-n", query.split("/")[-1],"-p", "1234", "-T", num_thread, "-f", "v", "-w", folder], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename(folder + "RAxML_portableTree."+query.split("/")[-1]+".jplace", folder + query.split("/")[-1] + ".jplace")
	os.remove(folder + "RAxML_classification." + query.split("/")[-1])
	os.remove(folder + "RAxML_labelledTree." +query.split("/")[-1])
	os.remove(folder + "RAxML_originalLabelledTree." +query.split("/")[-1])
	os.remove(folder + "RAxML_classificationLikelihoodWeights." +query.split("/")[-1])
	os.remove(folder + "RAxML_entropy." +query.split("/")[-1])
	os.remove(folder + "RAxML_info." +query.split("/")[-1])
	return folder + query.split("/")[-1] +".jplace"


def clean(sfolder):
	jks = glob.glob(sfolder + "*.reduced")
	for jk in jks:
		os.remove(jk)


def print_options():
	print (os.path.dirname(os.path.abspath(__file__)))
	print("usage: ./EPA_PTP.py -step species_counting -refaln example/ref.afa -query example/query.afa\n")
	print("Options example:\n")
	print("Note: the input sequences must be in fasta format.\n")
	print("  ./EPA_PTP.py -step alignment -refaln example/ref.afa  -query example/query.fa -minl 100 -outname epa_ready")
	print("                  This will align the query sequece to the reference sequence,")
	print("                  it requires the reference sequence to be aligned. ")
	print("                  The program will use HAMMER to make the alignment and convert the results to fasta, ")
	print("                  so the executable file hmmbuild and hmmalign must be in the folder bin/ \n")
	print("                  -minl n will remove all the aligned sequences that are short than n excluding gaps. \n")
	print("  ./EPA_PTP.py -step species_counting -refaln example/ref.afa  -query example/query.afa -minlw 0.5 -pv 0.001 -T 2")
	print("                  This example will count how many species in the query data. ")
	print("                  Both reference and query sequences must be aligned, and of the same length.\n")
	print("                  -minlw t(0-1) will set the minimal likelihood weight of the EPA placement,")
	print("                  placements has a lw smaller than t will be removed. Lower t will let more noisy ")
	print("                  sequences to pass and therefore more species will be found. \n")
	print("                  -pv (0-1) set the p-vlaue of likelihood ratio test, higher p-value will lead to more species. \n")
	print("                  -T (2-) number of CPUs. If the query data is from NGS which typically has more than 100,000 reads, ")
	print("                  the compute will be intensive, so we suggest to run this pipeline on a dedicated multi-core server.")
	print("                  100,000 reads on a eight-core Intel i7 server will typically need 24-48 hours to finish.\n")
	print("                  The results will be written to example/spcount.log which will be the input of 'summary' step\n")
	print("  ./EPA_PTP.py -step crop_species_counting -refaln example/ref.afa  -query example/query.afa -minlw 0.5 -T 2")
	print("                  This will use CROP to do open reference OTU-picking with EPA.\n")
	print("  ./EPA_PTP.py -step crop_stand_alone -refaln example/ref.afa  -query example/query.afa ")
	print("                  This will use CROP to do open reference OTU-picking alone.\n")
	print("  ./EPA_PTP.py -step summary -query example/spcount.log")
	print("                  Summarize the species counting results.")


if __name__ == "__main__":
	basepath = os.path.dirname(os.path.abspath(__file__))
	if not os.path.exists(basepath + "/bin/usearch"):
		print("The pipeline uses USEARCH to remove chimera seqeunces,")
		print("please downlaod the programm from:")
		print("http://www.drive5.com/usearch/")
		print("Rename the executable to usearch and put it to bin/  \n")
		sys.exit() 
	
	if not os.path.exists(basepath + "/bin/raxmlHPC-PTHREADS-SSE3"):
		print("The pipeline uses RAxML to infer phylogenetic trees,")
		print("please download the latest source code from: ")
		print("https://github.com/stamatak/standard-RAxML")
		print("Please complie the SSE + PTHREAD version, ")
		print("rename the executable to raxmlHPC-PTHREADS-SSE3 and put it to bin/  \n")
		sys.exit() 
	
	if len(sys.argv) < 3:
		print_options()
		sys.exit() 
		
	sstep = ""
	sfolder = ""
	saln = ""
	sjplace = ""
	sreftree = ""
	sappend = ""
	binbase = ""
	numt = "2"
	squery = ""
	soutname = "epa_ready"
	iminl = 100
	fminlw = 0.5
	pvalue = 0.001
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-step":
			i = i + 1
			sstep = sys.argv[i]
		elif sys.argv[i] == "-outname":
			i = i + 1
			soutname = sys.argv[i]
		#elif sys.argv[i] == "-folder":
		#	i = i + 1
		#	sfolder = sys.argv[i]
		elif sys.argv[i] == "-jplace":
			i = i + 1
			sjplace = sys.argv[i]
		elif sys.argv[i] == "-refaln":
			i = i + 1
			saln = sys.argv[i]
		elif sys.argv[i] == "-reftree":
			i = i + 1
			sreftree = sys.argv[i]
		elif sys.argv[i] == "-minl":
			i = i + 1
			iminl = int(sys.argv[i])
		elif sys.argv[i] == "-minlw":
			i = i + 1
			fminlw = float(sys.argv[i])
		elif sys.argv[i] == "-T":
			i = i + 1
			numt = sys.argv[i]
		elif sys.argv[i] == "-query":
			i = i + 1
			squery = sys.argv[i]
		elif sys.argv[i] == "-pv":
			i = i + 1
			pvalue = float(sys.argv[i])
		elif i == 0:
			pass
		elif sys.argv[i].startswith("-"):
			print("Unknown options: " + sys.argv[i])
			print_options()
			sys.exit()
	
	if squery == "":
		print("Must specify the query sequences file.")
		print_options()
		sys.exit()
	
	sfolder = os.path.dirname(os.path.abspath(squery)) + "/"
	print("Will output results in " + sfolder)
	
	if sstep == "alignment":
		if not os.path.exists(basepath + "/bin/hmmbuild") or not os.path.exists(basepath + "/bin/hmmalign"):
			print("The pipeline uses HAMMER to align the query seqeunces,")
			print("please downlaod the programm from:")
			print("http://hmmer.janelia.org/")
			print("Copy the executables hmmbuild and hmmalign to bin/  \n")
			sys.exit()
		if sfolder == "" or squery == "" or saln == "":
			print("Must specify reference alignment and query sequence files.")
			print_options()
			sys.exit()
		hmm_alignment(ref_align = saln , query = squery,  outfolder = sfolder, lmin = iminl, outname = soutname)
	elif sstep == "species_counting":
		if sfolder == "" or squery == "" or saln == "":
			print("Must specify reference alignment and query alignment files.")
			print_options()
			sys.exit()
		epa_me_species_counting(refaln = saln, queryaln = squery, folder = sfolder, lw = fminlw, T =  numt, pvalue = pvalue)
	elif sstep == "summary":
		if squery == "":
			print("Must specify the species delimitation log file - spcount.log, using -query option")
			print_options()
			sys.exit()
		stas(sfin = squery)
	elif sstep == "reduce_ref":
		random_remove_taxa(falign = saln, num_remove = int(numt), num_repeat = 1)
	elif sstep == "crop_species_counting":
		if not os.path.exists(basepath + "/bin/crop"):
			print("The program CROP does not exist,")
			print("please downlaod the programm from:")
			print("https://code.google.com/p/crop-tingchenlab/")
			print("Rename the executable to crop and put it to bin/  \n")
			sys.exit()
		
		if sfolder == "" or squery == "" or saln == "":
			print("Must specify reference alignment and query alignment files.")
			print_options()
			sys.exit()
		epa_crop_species_counting(refaln = saln, queryaln = squery, folder = sfolder, lw = fminlw, T =  numt)
	elif sstep == "crop_stand_alone":
		if not os.path.exists(basepath + "/bin/crop"):
			print("The program CROP does not exist,")
			print("please downlaod the programm from:")
			print("https://code.google.com/p/crop-tingchenlab/")
			print("Rename the executable to crop and put it to bin/  \n")
			sys.exit()
		
		if sfolder == "" or squery == "" or saln == "":
			print("Must specify the reference alignment and query alignment files.")
			print_options()
			sys.exit()
		crop_stand_alone(refaln = saln, queryaln = squery, folder = sfolder)
	else:
		print("Unknown options: " + sstep)
		print_options()
		sys.exit()
