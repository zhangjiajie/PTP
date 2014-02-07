#! /usr/bin/env python
try:
	import sys
	import math
	import random
	import argparse
	import os
	from ete2 import Tree
	from nexus import NexusReader
	from summary import *
	import matplotlib.pyplot as plt
except ImportError:
	print("Please install the scipy and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()



class exp_distribution:
	"""Implement exponential distribution"""
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
	"""Store one delimitation"""
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
	"""ML PTP, to use: init(), search() and count_species()"""
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
		if print_log:
			print("Init logl: " + "{0:.3f}".format(self.max_logl))
		num_sp, self.species_list = self.max_setting.count_species()
		return num_sp


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
	"""MCMC on a single tree using PTP model"""
	def __init__(self, tree, start_config = None, reroot = False, startmethod = "H0", min_br = 0.0001, seed = 1234, thinning = 100, sampling = 10000, burning = 0.1, taxa_order = []):
		if start_config == None:
			me = exponential_mixture(tree= tree)
			me.search(strategy = startmethod, reroot = reroot)
			me.count_species()
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
				if cnt % self.thinning == 0 and cnt >= sample_start:
					to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
					self.partitions.append(spe)
					self.llhs.append(newlogl) 
				accepted = accepted + 1
			else:
				u = self.rand_nr.uniform(0.0,1.0)
				if (u < acceptance):
					if cnt % self.thinning == 0 and cnt >= sample_start:
						to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
						self.partitions.append(spe)
						self.llhs.append(newlogl) 
					accepted = accepted + 1
				else:
					self.current_setting = self.last_setting
					self.current_logl = self.last_logl
					if cnt % self.thinning == 0 and cnt >= sample_start:
						to, spe = self.current_setting.output_species(taxa_order = self.taxaorder)
						self.partitions.append(spe)
						self.llhs.append(self.current_logl) 
		
		print("Accptance rate: " + repr(float(accepted)/float(cnt)))
		print("Merge: " + repr(self.nmerge))
		print("Split: " + repr(self.nsplit))
		return self.partitions, self.llhs



class bayesianptp:
	"""Run MCMC on multiple trees"""
	def __init__(self, filename, ftype = "nexus", reroot = False, method = "H1", seed = 1234, thinning = 100, sampling = 10000, burnin = 0.1, firstktrees = 0):
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
		cnt = 1 
		for tree in self.trees:
			print("Running MCMC sampling on tree " + repr(cnt) + "........")
			cnt = cnt + 1
			mcptp = ptpmcmc(tree = tree, reroot = self.reroot, startmethod = self.method, min_br = 0.0001, 
			seed = self.seed, thinning = self.thinning, sampling = self.sampling, burning = self.burnin, taxa_order = self.taxa_order)
			pars, lhs = mcptp.mcmc()
			self.partitions.extend(pars)
			self.llhs.extend(lhs)
			print("")
		return self.partitions, self.llhs
	
	
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



def parse_arguments():
	parser = argparse.ArgumentParser(description="""bPTP: a Bayesian implementation of the PTP model 
									written by Jiajie Zhang.
									Bugs, questions and suggestions please send to bestzhangjiajie@gmail.com""")
	
	parser.add_argument("-t", dest = "trees", 
						help = """Input phylogenetic tree file. Trees can be both rooted or unrooted, 
						if unrooted, please use -r option. Supported format: NEXUS (trees without annotation),
						RAxML (simple Newick foramt).""",
						required = True)
	
	parser.add_argument("-o", dest = "output",
						help = "Output file name",
						required = True)
	
	parser.add_argument("-s", dest = "seed", 
						help = """MCMC seed.""",
						type = int,
						required = True)
	
	parser.add_argument("-r", dest = "reroot",
						help = """Re-rooting the input tree on the longest branch (default not).""",
						default = False,
						action="store_true")
	
	parser.add_argument("-g", dest = "outgroups", 
						nargs='+',
						help = """Outgroup names, seperate by space. If this option is specified, 
						all trees will be rerooted accordingly.""")
	
	parser.add_argument("-d", dest = "delete", 
						help = """Remove outgroups specified by -g (default not).""",
						default = False,
						action="store_true")
	
	parser.add_argument("-m", dest = "method", 
						help = """Method for generate the starting partition (H0, H1, H2, H3) (default H1).""",
						choices=["H0", "H1", "H2", "H3"],
						default= "H1")
	
	parser.add_argument("-i", dest = "nmcmc", 
						help = """Number of MCMC iterations (default 10000).""",
						type = int,
						default = 10000)
	
	parser.add_argument("-n", dest = "imcmc", 
						help = """MCMC sampling interval - thinning (default 100).""",
						type = int,
						default = 100)
	
	parser.add_argument("-b", dest = "burnin", 
						help = """MCMC burn-in proportion (default 0.1).""",
						type = float,
						default = 0.1)
	
	parser.add_argument("-k", dest = "num_trees",
						help = """Run bPTP on first k trees (default all trees)""",
						type = int,
						default = 0)
	
	parser.add_argument('--version', action='version', version='%(prog)s 0.1 (07-02-2014)')
	
	return parser.parse_args()



def print_run_info(args):
    print("bPTP finished running with the following parameters:")
    print(" Input tree:.....................%s" % args.trees)
    print(" MCMC iterations:................%d" % args.nmcmc)
    print(" MCMC sampling interval:.........%d" % args.imcmc)
    print(" MCMC burn-in:...................%f" % args.burnin)
    print(" MCMC seed:......................%d" % args.seed)
    print("")
    print(" MCMC samples written to:")
    print("  "+args.output + ".bPTPPartitions.txt")
    print("")
    print(" Posterial LLH written to:")
    print("  "+args.output + ".bPTPllh.txt")
    print("")
    print(" Posterial LLH plot:")
    print("  "+args.output + ".llh.png")
    print("")
    print(" Posterial Prob. of partitions written to:")
    print("  "+args.output + ".bPTPPartitonSummary.txt")
    print("")
    print(" Highest posterial Prob. supported partition written to:")
    print("  "+args.output + ".bPTPPartitions.txt")


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
	
	pars, llhs = bbptp.delimit()
	pp = partitionparser(taxa_order = bbptp.taxa_order, partitions = pars, llhs = llhs)
	pp.summary(fout = args.output)
	
	print_run_info(args)
