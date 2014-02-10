#! /usr/bin/env python
try:
	import math
	import random
	import sys
	from ete2 import Tree, NodeStyle, TreeStyle
	from scipy import stats
except ImportError:
	print("Please install the matplotlib and other dependent package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	sys.exit()


class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p



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
			print("#taxa_order != num_taxa!")
			print(repr(num_taxa) + "	" + repr(len(taxa_order)))
			for sp in self.spe_list:
				print sp
			sys.exit()
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
	"""ML search PTP, to use: __init__(), search() and count_species()"""
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
		else:
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



def showTree(delimitation, scale = 500, render = False, fout = "", form = "pdf"):
	"""delimitation: species_setting class"""
	tree = delimitation.root
	style0 = NodeStyle()
	style0["fgcolor"] = "#000000"
	style0["vt_line_color"] = "#0000aa"
	style0["hz_line_color"] = "#0000aa"
	style0["vt_line_width"] = 2
	style0["hz_line_width"] = 2
	style0["vt_line_type"] = 0 
	style0["hz_line_type"] = 0
	style0["size"] = 0
	
	for node in tree.get_descendants():
		node.set_style(style0)
		node.img_style["size"] = 0
	
	tree.set_style(style0)
	tree.img_style["size"] = 0
	
	style1 = NodeStyle()
	style1["fgcolor"] = "#000000"
	style1["vt_line_color"] = "#ff0000"
	style1["hz_line_color"] = "#0000aa"
	style1["vt_line_width"] = 2
	style1["hz_line_width"] = 2
	style1["vt_line_type"] = 0 
	style1["hz_line_type"] = 0
	style1["size"] = 0
	
	style2 = NodeStyle()
	style2["fgcolor"] = "#0f0f0f"
	style2["vt_line_color"] = "#ff0000"
	style2["hz_line_color"] = "#ff0000"
	style2["vt_line_width"] = 2
	style2["hz_line_width"] = 2
	style2["vt_line_type"] = 0 
	style2["hz_line_type"] = 0
	style2["size"] = 0
	
	for node in delimitation.active_nodes:
		node.set_style(style1)
		node.img_style["size"] = 0
		for des in node.get_descendants():
			des.set_style(style2)
			des.img_style["size"] = 0
	ts = TreeStyle()
	#ts.show_leaf_name = True
	"""scale pixels per branch length unit"""
	ts.scale =  scale 
	if render:
		tree.render(fout+"."+form, tree_style=ts)
	else:
		tree.show(tree_style=ts)
