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
		self.node_can_split = []
		self.node_can_merge = []
	
	
	def get_nodes_can_split(self):
		if self.node_can_split == []:
			for node in self.active_nodes:
				if not node.is_leaf():
					self.node_can_split.append(node)
			return len(self.node_can_split)
		else:
			return len(self.node_can_split)
	
	
	def get_nodes_can_merge(self):
		if self.node_can_merge == []:
			for node in self.spe_nodes:
				children = node.get_children()
				if len(children) >= 2:
					if (children[0] in self.active_nodes) and (children[1] in self.active_nodes):
						self.node_can_merge.append(node)
			return len(self.node_can_merge)
		else:
			return len(self.node_can_merge)
	
	
	def get_log_l(self):
		if self.logl != 0:
			return self.logl
		else:
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
	
	
	def output_species(self, taxa_order = []):
		"""taxa_order is a list of taxa names, the paritions will be output as the same order"""
		if len(taxa_order) == 0:
			taxa_order = self.root.get_leaf_names()
		
		self.count_species()
		
		num_taxa = 0
		for sp in self.spe_list:
			for leaf in sp:
				num_taxa = num_taxa + 1
		if not len(taxa_order) == num_taxa:
			print("error error, taxa_order != num_taxa!")
			print(repr(num_taxa) + "	" + repr(len(taxa_order)))
			for sp in self.spe_list:
				print sp
			return None, None
		else: 
			partion = [-1] * num_taxa
			cnt = 1
			for sp in self.spe_list:
				for leaf in sp:
					idx = taxa_order.index(leaf)
					partion[idx] = cnt
				cnt = cnt + 1
			return taxa_order, partion



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


	def search(self, strategy = "H1", reroot = False):
		if strategy == "H1":
			self.H1(reroot)
		elif strategy == "H2":
			self.H2(reroot)
		elif strategy == "H3":
			self.H3(reroot)
		elif strategy == "Brutal":
			self.Brutal(reroot)
		else:
			self.H0(reroot)


	def count_species(self, print_log = True, pv = 0.001):
		#lhr = lh_ratio_test(self.null_logl, self.max_logl, 1)
		#pvalue = lhr.get_p_value()
		if print_log:
			print("Speciation rate: " + "{0:.3f}".format(self.max_setting.rate2))
			print("Coalesecnt rate: " + "{0:.3f}".format(self.max_setting.rate1))
			print("Null logl: " + "{0:.3f}".format(self.null_logl))
			print("MAX logl: " + "{0:.3f}".format(self.max_logl))
			#print("P-value: " + "{0:.3f}".format(pvalue))
			#spefit, speaw = self.max_setting.e2.ks_statistic()
			#coafit, coaaw = self.max_setting.e1.ks_statistic()
			#print("Kolmogorov-Smirnov test for model fitting:")
			#print("Speciation: " + "Dtest = {0:.3f}".format(spefit) + " " + speaw)
			#print("Coalescent: " + "Dtest = {0:.3f}".format(coafit) + " " + coaaw)
			#self.max_setting.e1.write_file()
			#self.max_setting.e2.write_file()
		#if pvalue < pv:
		num_sp, self.species_list = self.max_setting.count_species()
		return num_sp
		#else:
		#	self.species_list = []
		#	self.species_list.append(self.tree.get_leaf_names()) 
		#	return 1


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



class ptpmcmc:
	def __init__(self, tree, start_config, min_br = 0.0001, seed = 1234, thinning = 100, sampling = 10000, taxa_order = []):
		self.tree = tree 
		self.current_setting = start_config
		self.last_setting = start_config
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
	
	
	def mcmc(self, sampling = 10000):
		self.sampling = sampling
		cnt = 0
		accepted = 0 
		while cnt < self.sampling:
			cnt = cnt + 1
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
			
			if acceptance > 1.0:
				to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
				self.partitions.append(spe)
				self.llhs.append(newlogl) 
				accepted = accepted + 1
			else:
				u = self.rand_nr.uniform(0.0,1.0)
				if (u < acceptance):
					to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
					self.partitions.append(spe)
					self.llhs.append(newlogl) 
					accepted = accepted + 1
				else:
					self.current_setting = self.last_setting
					self.current_logl = self.last_logl
					to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
					self.partitions.append(spe)
					self.llhs.append(self.current_logl) 
		
		print("Accptance rate: " + repr(float(accepted)/float(cnt)))
		print("Merge: " + repr(self.nmerge))
		print("Split: " + repr(self.nsplit))
	
	
	def summary(self, burning = 0.1, thinning = 100, fout = ""):
		tpartitions = []
		tllhs = []
		sample_start = int(self.sampling * burning)
		for i in range(sample_start, len(self.partitions)):
			if (i % thinning == 0):
				tpartitions.append(self.partitions[i])
				tllhs.append(self.llhs[i])
		
		if fout!="":
			fo = open(fout + ".bPTPresults.txt", "a")
			fo.write("# Mean LLH:  "+str(numpy.mean(tllhs)) + "\n")
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
				fo.write("#" + onespe + ": " + "{0:.3f}".format(float(value)/float(len(tpartitions))) + "\n")
			
			maxw = 0
			bestpar = None
			bestsupport = None
			supports = []
			
			output = "#taxaorder:"+self._print_list(self.taxaorder)
			for i in range(len(tpartitions)): 
				partition = tpartitions[i]
				pars = idxpars[i]
				support = [1] * self.numtaxa
				sumw = 0.0
				for par in pars:
					w = pmap[par]
					for idx in par:
						support[idx] = float(w)/float(len(tpartitions))
						sumw = sumw + w #float(w)/float(self.numtrees)
				if sumw > maxw:
					maxw = sumw
					bestpar = i
					bestsupport = support
				
				supports.append(support)
				output= output + self._print_2lists(partition, support)
			
			spes, support = self._partition2names(tpartitions[bestpar], bestsupport)
			
			fo.write("#------------------------------------------------------------------------------------------------------------------------------------------------\n")
			fo.write("# Most supported partition found by simple heuristic search\n")
			for i in range(len(spes)):
				spe = spes[i]
				sup = support[i]
				fo.write("# Species " + str(i+1) + " (support = " + "{0:.3f}".format(sup) + ")\n")
				fo.write("#     " + self._print_list(spe) + "#\n")
			
			fo.write(output)
			fo.close()
			
			plt.plot(tllhs)
			plt.ylabel('Log likelihood')
			plt.xlabel('Iterations')
			plt.savefig(fout + ".png")
		else:
			return tpartitions, tllhs
	
	
	def _convert2idx(self, partition):
		a = min(partition)
		b = max(partition) + 1
		par = []
		for i in range(a, b):
			indices = [j for j, x in enumerate(partition) if x == i]
			par.append(tuple(indices))
		return par
	
	
	def _print_list(self, l):
		ss = ""
		for e in l:
			ss = ss + str(e) + ", "
		return ss[:-2] + "\n"
	
	
	def _print_2lists(self, l1, l2):
		ss = ""
		for i in range(len(l1)):
			e1 = l1[i]
			e2 = l2[i]
			ss = ss + str(e1)+"|"+"{0:.3f}".format(e2) + "\t"
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
					onepar.append(self.taxaorder[j])
					onesup.append(sup)
			nameparts.append(onepar)
			namesupps.append(onesup[0])
		
		return nameparts, namesupps



class bayesianptp:
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
		print("usage: python bPTP.py -t example/nex.test -o example/nex.out -r")
		print("Options:")
		print("    -t input                       Specify the input NEXUS file, trees can be both rooted or unrooted,")
		print("                                   if unrooted, please use -r option.\n")
		print("    -o output                      Specify output file name.\n")
		print("    -g outgroupnames               t1,t2,t3  commma delimt and no space in between")
		print("    -r                             Rooting the input tree on the longest branch.(default not)\n")
		print("    -d                             Remove outgroups specified by -g. (default not)\n")
		print("    -i                             Number of MCMC iterations. (default 10000)\n")
		print("    -n                             Number of MCMC sampling interval - thinning. (default 10)\n")
		print("    -s                             MCMC seed. (default 22)\n")
		#print("    -pvalue (0-1)                  Set the p-value for likelihood ratio test.(default 0.001)")


if __name__ == "__main__":
	me = exponential_mixture(tree= "/home/zhangje/GIT/SpeciesCounting/example/RAxML_bestTree.s8")
	me.H1(reroot = True)
	me.count_species()
	init_setting = me.max_setting
	bpm = ptpmcmc(tree = me.tree, start_config = init_setting, min_br = 0.0001, seed = 22)
	bpm.mcmc(sampling = 100000)
	bpm.summary(burning = 0.1, thinning = 100, fout = "/home/zhangje/GIT/SpeciesCounting/example/t1")
