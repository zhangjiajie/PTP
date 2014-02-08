#! /usr/bin/env python
try:
	import sys
	import math
	import collections
	import ete2
	import os
	import argparse
	import subprocess
	from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
	from collections import deque
	from scipy import stats
	from numpy import array
	from subprocess import call
	from nexus import NexusReader
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	#print(" sudo easy_install -U ete2")
	#print("Otherwise, please go to http://ete.cgenomics.org/ for instructions")
	sys.exit()


class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p


class exp_distribution:
	def __init__(self, data, rate = -1):
		self.data = data
		self.rate = 0.0
		if rate < 0:
			self.estimate_rate()
		else:
			self.rate = rate
		
	def __str__(self):
		return "Exponential distribution with rate = " + repr(self.rate)
	
	def estimate_rate(self):
		sumbr = 0.0
		numbr = len(self.data)
		for br in self.data:
			sumbr = sumbr + br
		if sumbr == 0:
			self.rate = 0
		else:
			self.rate = float(numbr) / sumbr
		
	def log_l(self, x):
		prob = self.rate * math.exp (-1.0 * self.rate * x)
		if prob == 0:
			return float("-inf") 
		else:
			return math.log(prob)
		
	def exp_cdf(self, rate, x):
		prob = 1 - math.exp (-1.0 * rate * x)
		return prob
	
	def sum_log_l(self):
		s = 0.0
		for br in self.data:
			s = s + self.log_l(br)
		return s 
	
	def ks_statistic(self):
		"""
		0.1 :  0.990
		0.05:  1.094
		0.001: 1.308
		"""
		X = self.data
		X.sort()
		s_list = []
		n = float(len(X))
		for i in range(len(X)):
			x = X[i]
			l = abs( (float(i+1)/n) - self.exp_cdf(self.rate, x) )
			r = abs( self.exp_cdf(self.rate, x) - (float(i)/n) )
			s_list.append(l)
			s_list.append(r)
	
		Dn = max(s_list)
		Dtest = ( math.sqrt(n) + 0.26 + (0.5/math.sqrt(n)) ) * ( Dn - (0.2/n) )
		outtest = ""
		if Dtest <= 0.99:
			outtest = "p-value >= 0.1 excellent model fitting"
		elif Dtest <= 1.094:
			outtest = "p-value >= 0.05 good model fitting"
		elif Dtest <=1.308:
			outtest = "p-value >= 0.01 moderate model fitting"
		else:
			outtest = "p-value < 0.01 poor model fitting"
			
		return Dtest, outtest
	
	
	def write_file(self):
		fout  = open(repr(self.rate), "w")
		for br in self.data:
			fout.write(repr(br) + "\n")
		fout.close()


class species_setting:
	def __init__(self, spe_nodes, root, sp_rate = 0, fix_sp_rate = False, minbr = 0.0001):
		self.min_brl = minbr
		self.spe_rate = sp_rate
		self.fix_spe_rate = fix_sp_rate
		self.spe_nodes = spe_nodes
		self.root = root
		self.all_nodes = root.get_descendants()
		self.all_nodes.append(self.root)
		self.coa_nodes = []
		self.logl = 0
		self.spe_list = []
		for node in self.all_nodes:
			if not (node in self.spe_nodes):
				self.coa_nodes.append(node) 
		
		self.active_nodes = []
		for node in self.spe_nodes:
			if node.is_leaf():
				self.active_nodes.append(node)
			else:
				childs = node.get_children()
				flag = False
				for child in childs:
					if not (child in self.spe_nodes):
						flag = True
						break
				if flag:
					self.active_nodes.append(node)
		
		
	def get_log_l(self):
		spe_br = []
		coa_br = []
		for node in self.spe_nodes:
			if node.dist > self.min_brl:
				spe_br.append(node.dist)
		
		for node in self.coa_nodes:
			if node.dist > self.min_brl:
				coa_br.append(node.dist)
				
		self.e1 = exp_distribution(coa_br)
		self.e2 = None
		if self.fix_spe_rate:
			self.e2 = exp_distribution(spe_br, rate = self.spe_rate)
		else:
			self.e2 = exp_distribution(spe_br)
		self.rate1 = self.e1.rate
		self.rate2 = self.e2.rate
		logl = self.e1.sum_log_l() + self.e2.sum_log_l()
		self.logl = logl
		return logl
		
	def count_species(self):
		if len(self.spe_list) != 0:
			return len(self.spe_list), self.spe_list
		else:
			for node in self.active_nodes:
				one_spe = []
				if node.is_leaf():
					one_spe.append(node.name)
				else:
					one_spe.extend(node.get_leaf_names())
				self.spe_list.append(one_spe)
			return len(self.spe_list), self.spe_list
			
	def whiten_species(self):
		if len(self.spe_list) != 0:
			avg_num = 0.0
			for sp in self.spe_list:
				avg_num = avg_num + len(sp)
			avg_num = avg_num / float(len(self.spe_list))
			avg_num = int(math.ceil(avg_num))
			print(avg_num)
			new_spe_list = []
			for sp in self.spe_list:
				if len(sp) > avg_num:
					newsp = sp[:avg_num] 
					new_spe_list.extend(newsp)
				else:
					new_spe_list.extend(sp)
			return new_spe_list


class exponential_mixture:
	"""init(), search() and count_species()"""
	def __init__(self, tree, sp_rate = 0, fix_sp_rate = False, max_iters = 20000, min_br = 0.0001):
		self.min_brl = min_br
		self.tree = Tree(tree, format = 1)
		self.tree.dist = 0.0
		self.fix_spe_rate = fix_sp_rate
		self.fix_spe = sp_rate
		self.max_logl = float("-inf") 
		self.max_setting = None
		self.null_logl = 0.0
		self.null_model()
		self.species_list = None
		self.counter = 0
		self.setting_set = set([])
		self.max_num_search = max_iters


	def null_model(self):
		coa_br = []
		all_nodes = self.tree.get_descendants()
		for node in all_nodes:
			if node.dist > self.min_brl:
				coa_br.append(node.dist)
		e1 = exp_distribution(coa_br)
		self.null_logl = e1.sum_log_l()
		return e1.rate


	def __compare_node(self, node):
		return node.dist


	def re_rooting(self):
		node_list = self.tree.get_descendants()
		node_list.sort(key=self.__compare_node)
		node_list.reverse()
		rootnode = node_list[0]
		self.tree.set_outgroup(rootnode)
		self.tree.dist = 0.0


	def comp_num_comb(self):
		for node in self.tree.traverse(strategy='postorder'):
			if node.is_leaf():
				node.add_feature("cnt", 1.0)
			else:
				acum = 1.0
				for child in node.get_children():
					acum = acum * child.cnt
				acum = acum + 1.0
				node.add_feature("cnt", acum)
		return self.tree.cnt


	def next(self, sp_setting):
		self.setting_set.add(frozenset(sp_setting.spe_nodes))
		logl = sp_setting.get_log_l()
		if logl > self.max_logl:
			self.max_logl = logl
			self.max_setting = sp_setting
		for node in sp_setting.active_nodes:
			if node.is_leaf():
				pass
			else:
				childs = node.get_children()
				sp_nodes = []
				for child in childs:
					sp_nodes.append(child)
				for nod in sp_setting.spe_nodes:
					sp_nodes.append(nod)
				new_sp_setting = species_setting(spe_nodes = sp_nodes, root = sp_setting.root, sp_rate = sp_setting.spe_rate, fix_sp_rate = sp_setting.fix_spe_rate, minbr = self.min_brl)
				if frozenset(sp_nodes) in self.setting_set:
					pass
				else:
					self.next(new_sp_setting)


	def H0(self, reroot = True):
		self.H1(reroot)
		self.H2(reroot = False)
		self.H3(reroot = False)


	def H1(self, reroot = True):
		if reroot:
			self.re_rooting()
			
		#self.init_tree()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
		last_setting = first_setting
		max_logl = last_setting.get_log_l()
		max_setting = last_setting
		
		for node in sorted_node_list:
			if node not in last_setting.spe_nodes:
				curr_sp_nodes = []
				for nod in last_setting.spe_nodes:
					curr_sp_nodes.append(nod)
				
				chosen_branching_node = node.up #find the father of this new node
				if chosen_branching_node in last_setting.spe_nodes:
					for nod in chosen_branching_node.get_children():
						if nod not in curr_sp_nodes:
							curr_sp_nodes.append(nod)
				else:
					for nod in chosen_branching_node.get_children():
						if nod not in curr_sp_nodes:
							curr_sp_nodes.append(nod)
					while not chosen_branching_node.is_root():
						chosen_branching_node = chosen_branching_node.up
						for nod in chosen_branching_node.get_children():
							if nod not in curr_sp_nodes:
								curr_sp_nodes.append(nod)
						if chosen_branching_node in last_setting.spe_nodes:
							break
				new_setting = species_setting(spe_nodes = curr_sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
				new_logl = new_setting.get_log_l()
				if new_logl> max_logl:
					max_logl = new_logl
					max_setting = new_setting 
				last_setting = new_setting
				
			else:
				"""node already is a speciation node, do nothing"""
				pass
		
		if max_logl > self.max_logl:
			self.max_logl = max_logl
			self.max_setting = max_setting


	def H2(self, reroot = True):
		"""Greedy"""
		if reroot:
			self.re_rooting()
			
		#self.init_tree()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
		last_setting = first_setting
		max_logl = last_setting.get_log_l()
		max_setting = last_setting
		contin_flag = True 
		
		
		while contin_flag:
			curr_max_logl = float("-inf") 
			curr_max_setting = None
			contin_flag = False
			for node in last_setting.active_nodes:
				if node.is_leaf():
					pass
				else:
					contin_flag = True 
					childs = node.get_children()
					sp_nodes = []
					for child in childs:
						sp_nodes.append(child)
					for nod in last_setting.spe_nodes:
						sp_nodes.append(nod)
					new_sp_setting = species_setting(spe_nodes = sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
					logl = new_sp_setting.get_log_l()
					if logl > curr_max_logl:
						curr_max_logl = logl
						curr_max_setting = new_sp_setting
			
			if curr_max_logl > max_logl:
				max_setting = curr_max_setting
				max_logl = curr_max_logl
			
			last_setting = curr_max_setting
			
		if max_logl > self.max_logl:
			self.max_logl = max_logl
			self.max_setting = max_setting


	def H3(self, reroot = True):
		if reroot:
			self.re_rooting()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		sorted_br = []
		for node in sorted_node_list:
			sorted_br.append(node.dist)
		maxlogl = float("-inf") 
		maxidx = -1
		for i in range(len(sorted_node_list))[1:]:
			l1 = sorted_br[0:i]
			l2 = sorted_br[i:]
			e1 = exp_distribution(l1)
			e2 = exp_distribution(l2)
			logl = e1.sum_log_l() + e2.sum_log_l()
			if logl > maxlogl:
				maxidx = i
				maxlogl = logl
		
		target_nodes = sorted_node_list[0:maxidx]
		
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
		last_setting = first_setting
		max_logl = last_setting.get_log_l()
		max_setting = last_setting
		contin_flag = True 
		target_node_cnt = 0
		while contin_flag:
			curr_max_logl = float("-inf") 
			curr_max_setting = None
			contin_flag = False
			unchanged_flag = True
			for node in last_setting.active_nodes:
				if node.is_leaf():
					pass
				else:
					contin_flag = True 
					childs = node.get_children()
					sp_nodes = []
					flag = False
					for child in childs:
						if child in target_nodes:
							flag = True
							#target_nodes.remove(child)
					if flag:
						unchanged_flag = False
						for child in childs:
							sp_nodes.append(child)
						for nod in last_setting.spe_nodes:
							sp_nodes.append(nod)
						new_sp_setting = species_setting(spe_nodes = sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
						logl = new_sp_setting.get_log_l()
						if logl > curr_max_logl:
							curr_max_logl = logl
							curr_max_setting = new_sp_setting
			if not unchanged_flag:
				target_node_cnt = target_node_cnt + 1
				if curr_max_logl > max_logl:
					max_setting = curr_max_setting
					max_logl = curr_max_logl
				last_setting = curr_max_setting
			
			if len(target_nodes) == target_node_cnt:
				contin_flag = False
			if contin_flag and unchanged_flag and last_setting!= None:
				for node in last_setting.active_nodes:
					if node.is_leaf():
						pass
					else:
						childs = node.get_children()
						sp_nodes = []
						for child in childs:
							sp_nodes.append(child)
						for nod in last_setting.spe_nodes:
							sp_nodes.append(nod)
						new_sp_setting = species_setting(spe_nodes = sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
						logl = new_sp_setting.get_log_l()
						if logl > curr_max_logl:
							curr_max_logl = logl
							curr_max_setting = new_sp_setting
				if curr_max_logl > max_logl:
					max_setting = curr_max_setting
					max_logl = curr_max_logl
				last_setting = curr_max_setting
				
		if max_logl > self.max_logl:
			self.max_logl = max_logl
			self.max_setting = max_setting


	def Brutal(self, reroot = False):
		if reroot:
			self.re_rooting()
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		num_s = self.comp_num_comb()
		if num_s > self.max_num_search:
			print("Too many search iterations: " + repr(num_s) + ", using H0 instead!!!")
			self.H0(reroot = False)
		else:
			first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate, minbr = self.min_brl)
			self.next(first_setting)


	def search(self, strategy = "H1", reroot = False):
		if strategy == "H1":
			self.H1(reroot)
		elif strategy == "H2":
			self.H2(reroot)
		elif strategy == "H3":
			self.H3(reroot)
		elif strategy == "Brutal":
			self.Brutal(reroot)
		else:# strategy == "H0":
			self.H0(reroot)


	def count_species(self, print_log = True, pv = 0.001):
		lhr = lh_ratio_test(self.null_logl, self.max_logl, 1)
		pvalue = lhr.get_p_value()
		if print_log:
			print("Speciation rate: " + "{0:.3f}".format(self.max_setting.rate2))
			print("Coalesecnt rate: " + "{0:.3f}".format(self.max_setting.rate1))
			print("Null logl: " + "{0:.3f}".format(self.null_logl))
			print("MAX logl: " + "{0:.3f}".format(self.max_logl))
			print("P-value: " + "{0:.3f}".format(pvalue))
			spefit, speaw = self.max_setting.e2.ks_statistic()
			coafit, coaaw = self.max_setting.e1.ks_statistic()
			print("Kolmogorov-Smirnov test for model fitting:")
			print("Speciation: " + "Dtest = {0:.3f}".format(spefit) + " " + speaw)
			print("Coalescent: " + "Dtest = {0:.3f}".format(coafit) + " " + coaaw)
			#self.max_setting.e1.write_file()
			#self.max_setting.e2.write_file()
		if pvalue < pv:
			num_sp, self.species_list = self.max_setting.count_species()
			return num_sp
		else:
			self.species_list = []
			self.species_list.append(self.tree.get_leaf_names()) 
			return 1


	def whitening_search(self, strategy = "H1", reroot = False, pv = 0.001):
		self.search(strategy, reroot, pv)
		num_sp, self.species_list = self.max_setting.count_species()
		spekeep = self.max_setting.whiten_species()
		self.tree.prune(spekeep)
		self.max_logl = float("-inf") 
		self.max_setting = None
		self.null_logl = 0.0
		self.null_model()
		self.species_list = None
		self.counter = 0
		self.setting_set = set([])
		self.search(strategy, reroot, pv)


	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			for leaf in sp:
				print("          " + leaf)
			cnt = cnt + 1


	def output_species(self, taxa_order = []):
		"""taxa_order is a list of taxa names, the paritions will be output as the same order"""
		if len(taxa_order) == 0:
			taxa_order = self.tree.get_leaf_names()
		
		num_taxa = 0
		for sp in self.species_list:
			for leaf in sp:
				num_taxa = num_taxa + 1
		if not len(taxa_order) == num_taxa:
			print("error error, taxa_order != num_taxa!")
			return None, None
		else: 
			partion = [-1] * num_taxa
			cnt = 1
			for sp in self.species_list:
				for leaf in sp:
					idx = taxa_order.index(leaf)
					partion[idx] = cnt
				cnt = cnt + 1
			return taxa_order, partion


	def showTree(self, scale = 500, render = False, fout = "", form = "pdf"):
		style0 = NodeStyle()
		style0["fgcolor"] = "#000000"
		#style2["shape"] = "circle"
		style0["vt_line_color"] = "#0000aa"
		style0["hz_line_color"] = "#0000aa"
		style0["vt_line_width"] = 2
		style0["hz_line_width"] = 2
		style0["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style0["hz_line_type"] = 0
		style0["size"] = 0
		
		for node in self.tree.get_descendants():
			node.set_style(style0)
			node.img_style["size"] = 0
		self.tree.set_style(style0)
		self.tree.img_style["size"] = 0
		
		
		style1 = NodeStyle()
		style1["fgcolor"] = "#000000"
		#style2["shape"] = "circle"
		style1["vt_line_color"] = "#ff0000"
		style1["hz_line_color"] = "#0000aa"
		style1["vt_line_width"] = 2
		style1["hz_line_width"] = 2
		style1["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style1["hz_line_type"] = 0
		style1["size"] = 0
		
		style2 = NodeStyle()
		style2["fgcolor"] = "#0f0f0f"
		#style2["shape"] = "circle"
		style2["vt_line_color"] = "#ff0000"
		style2["hz_line_color"] = "#ff0000"
		style2["vt_line_width"] = 2
		style2["hz_line_width"] = 2
		style2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style2["hz_line_type"] = 0
		style2["size"] = 0
		
		for node in self.max_setting.active_nodes:
			node.set_style(style1)
			node.img_style["size"] = 0
			for des in node.get_descendants():
				des.set_style(style2)
				des.img_style["size"] = 0
		ts = TreeStyle()
		#ts.show_leaf_name = True
		ts.scale =  scale # scale pixels per branch length unit
		if render:
			self.tree.render(fout+"."+form, tree_style=ts)
		else:
			self.tree.show(tree_style=ts)



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
	parser = argparse.ArgumentParser(description="""PTP: maximal likelihood search of the Poisson Tree Processes model for species delimitation.

By using this program, you agree to cite: 
"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General Species 
Delimitation Method with Applications to Phylogenetic Placements.
Bioinformatics (2013), 29 (22): 2869-2876 " 

Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com.

Version 1.3 released by Jiajie Zhang on 08-02-2014.""",
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
			me.showTree(scale = args.sscale)
		else:
			me.showTree(scale = args.sscale, render = True, fout = args.stree , form = "pdf")
			me.showTree(scale = args.sscale, render = True, fout = args.stree , form = "png")
	except ete2.parser.newick.NewickError:
		print("Unexisting tree file or Malformed newick tree structure.")


