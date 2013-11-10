#! /usr/bin/env python
try:
	import sys
	import math
	import ete2
	import os
	import subprocess
	from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
	from subprocess import call
	from scipy.optimize import fmin
	from scipy.optimize import fmin_powell
	from scipy.optimize import fmin_l_bfgs_b
	from scipy.optimize import fmin_tnc
	from scipy import stats
	import matplotlib.pyplot as plt
	from subprocess import call
except ImportError:
	print("Please install the scipy, matplotlib package first.")
	print("If your OS is ubuntu or has apt installed, you can try the following:") 
	print(" sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib")
	#print(" sudo easy_install -U ete2")
	#print("Otherwise, please go to http://ete.cgenomics.org/ for instructions")
	sys.exit()

class um_tree:
	def __init__(self, tree):
		self.tree = Tree(tree, format = 1)
		self.tree.resolve_polytomy(default_dist=0.000001, recursive=True)
		self.tree.dist = 0
		self.tree.add_feature("age", 0)
		self.nodes = self.tree.get_descendants()
		internal_node = []
		cnt = 0
		for n in self.nodes:
			node_age = n.get_distance(self.tree)
			n.add_feature("age", node_age)
			if not n.is_leaf():
				n.add_feature("id", cnt)
				cnt = cnt + 1
				internal_node.append(n)
		self.nodes = internal_node
		one_leaf = self.tree.get_farthest_node()[0]
		one_leaf.add_feature("id", cnt+1)
		if one_leaf.is_leaf():
			self.nodes.append(one_leaf)
		self.nodes.sort(key=self.__compare_node)
		self.species_list = []
		self.coa_roots = None


	def __compare_node(self, node):
		return node.age


	def get_waiting_times(self, threshold_node = None, threshold_node_idx = 0):
		wt_list = []
		reach_t = False
		curr_age = 0.0
		curr_spe = 2
		curr_num_coa = 0
		coa_roots = []
		min_brl = 1000
		num_spe = -1
		
		if threshold_node == None:
			threshold_node = self.nodes[threshold_node_idx]
		
		last_coa_num = 0
		tcnt = 0 
		for node in self.nodes:
			num_children = len(node.get_children())
			wt = None
			times = node.age - curr_age
			if times >= 0:
				if times < min_brl and times > 0:
					min_brl = times
				curr_age = node.age
				assert curr_spe >=0
				 
				if reach_t:
					if tcnt == 0:
						last_coa_num = 2
					fnode = node.up
					coa_root = None
					
					idx = 0
					while not fnode.is_root():
						idx = 0 
						for coa_r in coa_roots:
							if coa_r.id == fnode.id:
								coa_root = coa_r
								break
							idx = idx + 1
						
						if coa_root!=None:
							break
						else:
							fnode = fnode.up
							
					wt = waiting_time(length = times, num_coas =curr_num_coa, num_lines = curr_spe)
					
					for coa_r in coa_roots:
						coa = coalescent(num_individual = coa_r.curr_n)
						wt.coas.add_coalescent(coa)
					
					wt.coas.coas_idx = last_coa_num
					wt.num_curr_coa = last_coa_num
					if coa_root == None: #here can be modified to use multiple T
						curr_spe = curr_spe - 1
						curr_num_coa = curr_num_coa + 1
						node.add_feature("curr_n", 2)
						coa_roots.append(node)
						last_coa_num = 2
					else:
						curr_n = coa_root.curr_n
						coa_root.add_feature("curr_n", curr_n + 1)
						last_coa_num = curr_n + 1
					tcnt = tcnt + 1
				else:
					if node.id == threshold_node.id:
						reach_t = True
						tcnt = 0 
						wt = waiting_time(length = times, num_coas = 0, num_lines = curr_spe)
						num_spe = curr_spe
						curr_spe = curr_spe - 1
						curr_num_coa = 2
						node.add_feature("curr_n", 2)
						coa_roots.append(node)
					else:
						wt = waiting_time(length = times, num_coas = 0, num_lines = curr_spe)
						curr_spe = curr_spe + 1
				if times > 0.00000001:
					wt_list.append(wt)
		
		
		for wt in wt_list:
			wt.count_num_lines()
		
		self.species_list = []
		all_coa_leaves = []
		self.coa_roots = coa_roots
		for coa_r in coa_roots:
			leaves = coa_r.get_leaves()
			all_coa_leaves.extend(leaves)
			self.species_list.append(leaves)
		
		all_leaves = self.tree.get_leaves()
		for leaf in all_leaves:
			if leaf not in all_coa_leaves:
				self.species_list.append([leaf])
		
		return wt_list, num_spe


	def show(self, wt_list):
		cnt = 1
		for wt in wt_list:
			print("Waitting interval "+ repr(cnt))
			print(wt)
			cnt = cnt + 1


	def get_species(self):
		sp_list = []
		for sp in self.species_list:
			spe = []
			for taxa in sp:
				spe.append(taxa.name)
			sp_list.append(spe)
		
		all_taxa_name = []
		
		#self.tree.convert_to_ultrametric(tree_length = 1.0, strategy='balanced')
		
		for leaf in self.tree.get_leaves():
			all_taxa_name.append(leaf.name)
		
		
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
		
		for node in self.coa_roots:
			node.set_style(style1)
			node.img_style["size"] = 0
			for des in node.get_descendants():
				des.set_style(style2)
				des.img_style["size"] = 0
		
		return [all_taxa_name], sp_list


	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			cnt = cnt + 1
			taxas = ""
			for taxa in sp:
				taxas = taxas + taxa.name + ", "
			print("	" + taxas[:-1])


	def num_lineages(self, wt_list):
		nl_list = []
		times = []
		last_time = 0.0
		for wt in wt_list:
			nl_list.append(wt.get_num_branches())
			times.append(last_time)
			last_time = wt.length + last_time
		
		plt.plot(times, nl_list)
		plt.ylabel('Number of lineages')
		plt.xlabel('Time')
		plt.savefig("Time_Lines")
		plt.show()


class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p


class speciation:
	def __init__(self, num_lineage, rate = 0 , p = 1.0):
		self.num_lineages = num_lineage #n
		self.rate = rate #speciation rate 
		self.p = p
		
	def __str__(self):
		s = "spe_event with n=" + repr(self.num_lineages) +  ", rate= " + repr(self.rate) + "\n" 
		return s
	
	def update(self,spe_rate, spe_p):
		self.rate = spe_rate
		self.p = spe_p
	
	def update_p(self, spe_p):
		self.p = spe_p
	
	def update_rate(self, spe_rate):
		self.rate = spe_rate
	
	def getBirthRate(self): # this is b
		if self.num_lineages <= 0:
			return 0
		else:
			rv = 0
			try:
				rv = self.rate * math.pow(self.num_lineages, self.p)
			except ValueError:
				rv = 0
			return rv
		
	def getNumDivEvent(self):
		return self.num_lineages-1
		
	def bprime(self):
		if self.num_lineages <= 0:
			return 0.0
		else:
			bp = 0
			try:
				bp = self.p * self.rate * math.pow(self.num_lineages, (self.p - 1.0))
			except ValueError:
				bp = 0
			
			return bp


class coalescent:
	def __init__(self, num_individual, rate = 0, p = 1.0):
		 self.num_individual = num_individual #n
		 self.rate = rate # 1/2N
		 self.p = p
	
	def getCoalesecntRate(self):
		rv = 0
		try:
			rv = self.rate * math.pow(self.num_individual * (self.num_individual - 1.0), self.p)
		except ValueError:
			rv = 0 
		return rv
		
	def getNumDivEvent(self):
		return self.num_individual -1
	
	def bprime(self):
		bp = 0
		try:
			bp = self.p * self.rate * math.pow(self.num_individual * (self.num_individual - 1.0), (self.p - 1.0))
		except ValueError:
			bp = 0
		return bp


class coalescents:
	def __init__(self, num_coalescent, rate = 0, p = 1, const_rate = True, const_p = True):
		self.coa_list = []
		self.num_coa =  num_coalescent
		self.rate = rate 
		self.p = p
		self.const_rate = const_rate
		self.const_p = const_p
		self.coas_idx = -1
	
	def __str__(self):
		s = ""
		cnt = 1
		for coa in self.coa_list:
			s = s + "coa_event" + repr(cnt) + ", with n= "+repr(coa.num_individual) + ", rate="+ repr(coa.rate)+ "\n" 
			cnt = cnt + 1
		return s
	
	def add_coalescent(self,coa):
		coa.rate = self.rate 
		coa.p = self.p
		self.coa_list.append(coa)
		
	def getSumCoaRate(self): #this is b
		sr = 0
		for coa in self.coa_list:
			sr = sr + coa.getCoalesecntRate()
		return sr
	
	def update(self, rate = 0, p = 1):
		self.p = p
		self.rate = rate
		for coa in self.coa_list:
			coa.rate = rate
			coa.p = p
	
	def update_p(self, p):
		self.p = p
		for coa in self.coa_list:
				coa.p = p
				
	
	def update_rate(self, rate):
		self.rate = rate
		for coa in self.coa_list:
				coa.rate = rate
				
				
	def getNumDivEvent(self):
		ndiv = 0
		for coa in self.coa_list:
			ndiv = ndiv + coa.getNumDivEvent()
		return ndiv
		
	def bprime(self):
		bp = 0.0
		for coa in self.coa_list:
			bp = bp + coa.bprime()
		return bp


class waiting_time:
	def __init__(self, length, num_coas = 0, num_lines = 0):
		self.length = length #x
		self.spe = speciation(num_lineage = num_lines)
		self.coas  =  coalescents(num_coalescent = num_coas)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		self.spe_rate = 0
		self.spe_p = 1
		self.spe_n = num_lines
		self.coa_rate = 0
		self.coa_p = 1
		self.coa_n = num_coas
		self.num_lines = None
		self.num_curr_coa = 0
	
	def count_num_lines(self):
		"""This is used for the null model only!!!"""
		self.num_lines = self.spe_n
		for coa in self.coas.coa_list:
			self.num_lines = self.num_lines + coa.num_individual
		self.num_lines = self.num_lines * (self.num_lines - 1)
		
	def get_num_branches(self):
		nbr = 0
		for coa in self.coas.coa_list:
			nbr = nbr + coa.num_individual
		nbr = nbr + self.spe.num_lineages
		return nbr
	
	def __str__(self):
		s = "Waitting time = " + repr(self.length) + "\n"
		s = s + str(self.spe) + "\n"
		s = s + str(self.coas)# + "\n"
		s = s + "Numer curr coa lines:"+ str(self.num_curr_coa) + "\n"
		s = s + "---------------------------------------------------------\n"
		return s 
	
	def update(self, sp_rate, sp_p, coa_rate, coa_p):
		self.spe_rate = sp_rate 
		self.spe_p = sp_p
		self.coa_rate = coa_rate 
		self.coa_p = coa_p
		self.spe.update(sp_rate, sp_p)
		self.coas.update(coa_rate, coa_p)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate() 
	
	def update_p(self, sp_p, coa_p):
		self.spe_p = sp_p
		self.coa_p = coa_p
		self.spe.update_p(sp_p)
		self.coas.update_p(coa_p)
	
	def update_rate(self, sp_rate, coa_rate):
		self.spe_rate = sp_rate 
		self.coa_rate = coa_rate 
		self.spe.update_rate(sp_rate)
		self.coas.update_rate(coa_rate)
	
	def logl(self):
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		prob = self.b * math.exp (-1.0 * self.b * self.length)
		if prob <=0:
			print("branch lenghth:" + repr(self.length))
			print("b:" + repr(self.b))
			return None
		else:
			return math.log(prob)
	
	def scaleSpeBranchL(self):
		if self.spe_n < 0:
			self.spe_n = 0.000001
		
		rvalue = 0
		try:
			rvalue = math.pow(self.spe_n, self.spe_p) * self.length
		except ValueError:
			rvalue = 0
		return rvalue
	
	def scaleCoaBranchL(self):
		br = 0
		for coa in self.coas.coa_list:
			try:
				br = br + math.pow(coa.num_individual * (coa.num_individual-1), self.coa_p) * self.length 
			except ValueError:
				pass
				
		return br
		
	def bprime_spe(self):
		bp = self.spe.bprime()
		bps = bp * math.exp(-1.0 * self.b * self.length) + self.b * math.exp (-1.0 * self.b * self.length) * -1.0 * self.length * bp
		return bps
		
	def bprime_coa(self):
		bp = self.coas.bprime()
		bps = bp * math.exp(-1.0 * self.b * self.length) + self.b * math.exp (-1.0 * self.b * self.length) * -1.0 * self.length * bp
		return bps


class tree_time:
	def __init__(self, wtimes, num_spe):
		self.w_time_list = wtimes
		self.llh = 0
		self.spe_rate = 0.001
		self.spe_p = 1
		self.coa_rate = 0.001
		self.coa_p = 1 
		self.num_species = num_spe
		self.numSpeEvent = self.num_species - 1
		last_wc = self.w_time_list[-1]
		self.numCoaEvent = 0
		for coa in last_wc.coas.coa_list:
			self.numCoaEvent = self.numCoaEvent + coa.getNumDivEvent() 
		spe_rate_dn = 0
		coa_rate_dn = 0 
		for w_time in self.w_time_list:
			spe_rate_dn = spe_rate_dn + w_time.scaleSpeBranchL()
			coa_rate_dn = coa_rate_dn + w_time.scaleCoaBranchL()
		
		if spe_rate_dn == 0:
			self.spe_rate = 0
		else:
			self.spe_rate = self.numSpeEvent/spe_rate_dn
			
		if coa_rate_dn == 0:
			self.coa_rate = 0
		else:
			self.coa_rate = self.numCoaEvent/coa_rate_dn
		
	def show(self):
		print("This is tree_time with spe event: " + repr(self.numSpeEvent) + ", coa event: " + repr(self.numCoaEvent)) 
	
	def sum_llh(self):
		self.llh = 0.0
		for w_time in self.w_time_list:
			logl = w_time.logl()
			if logl == None:
				self.llh = -sys.float_info.max
				print("wtime logl infinity!!")
				break
			else:
				self.llh = self.llh + logl
		return self.llh
	
	def update(self, spe_p, coa_p):
		self.spe_p = spe_p
		self.coa_p = coa_p
		
		for w_time in self.w_time_list:
			w_time.update_p(spe_p, coa_p)
		spe_rate_dn = 0
		coa_rate_dn = 0 
		for w_time in self.w_time_list:
			spe_rate_dn = spe_rate_dn + w_time.scaleSpeBranchL()
			coa_rate_dn = coa_rate_dn + w_time.scaleCoaBranchL()
		
		if spe_rate_dn == 0:
			self.spe_rate = 0
		else:
			self.spe_rate = self.numSpeEvent/spe_rate_dn
			
		if coa_rate_dn == 0:
			self.coa_rate =0
		else:
			self.coa_rate = self.numCoaEvent/coa_rate_dn
		
		
		for w_time in self.w_time_list:
			w_time.update_rate(self.spe_rate, self.coa_rate)
	
	def bprime_spe(self):
		bp = 0.0
		for wt in self.w_time_list:
			bp = bp + wt.bprime_spe()
		return bp
		
	def bprime_coa(self):
		bp = 0.0
		for wt in self.w_time_list:
			bp = bp + wt.bprime_coa()
		return bp


class null_model:
	"""The null model using Coalescent"""
	def __init__(self, wt_list, tree):
		nodes = tree.get_leaves()
		self.num_speEvent = len(nodes) - 1
		self.wt_list = wt_list
		self.p = 1.0
		self.rate = 0.0
	
	def logl(self, p = 1.0):
		self.p = p
		br_de = 0.0
		for wt in self.wt_list:
			if wt.num_lines < 0:
				wt.num_lines = 0
			try:
				br_de = br_de + wt.length * math.pow(wt.num_lines, self.p)
			except ValueError:
				pass
				
		self.rate = self.num_speEvent / br_de
		logl = 0
		for wt in self.wt_list:
			try:
				logl = logl + math.log(self.rate * math.pow(wt.num_lines, self.p) * math.exp(self.rate * math.pow(wt.num_lines, self.p) * wt.length * -1.0))
			except ValueError:
				pass
				
		return logl


def tar_fun(x, *args):
	"""args[0] is tree_time"""
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return args[0].sum_llh() * (-1.0)


def tar_fun_null(x, *args):
	return args[0].logl(p = x[0]) * (-1.0)


def prime_fun(x, *args):
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return [args[0].bprime_spe() * (-1.0) , args[0].bprime_coa() * (-1.0)]


def optimize_null_model(umtree):
	min_change = 0.1
	max_iters = 100
	wt_list, num_spe = umtree.get_waiting_times(threshold_node_idx = 0)
	nm = null_model(wt_list, umtree.tree)
	last_llh = float("-inf")
	change = float("inf")
	cnt = 0
	while change > min_change and cnt < max_iters:
		cnt = cnt + 1
		para, nn, cc = fmin_l_bfgs_b(tar_fun_null, [1], args = [nm], bounds = [[0, 10]], approx_grad = True)
		curr_logl = nm.logl(p = para[0])
		change = abs(curr_logl - last_llh)
		last_llh = curr_logl
	return last_llh


def gmyc(tree, print_detail = False, show_tree = False, show_llh = False, show_lineages = False, print_species = False, pv = 0.01):
	llh_list = []
	min_change = 0.1
	max_iters = 100
	best_llh = float("-inf")
	best_num_spe = -1
	best_node = None
	utree = um_tree(tree)
	for tnode in utree.nodes:
		wt_list, num_spe = utree.get_waiting_times(threshold_node = tnode)
		tt = tree_time(wt_list, num_spe)
		last_llh = float("-inf")
		change = float("inf")
		cnt = 0
		
		while change > min_change and cnt < max_iters:
			cnt = cnt + 1
			para, nn, cc = fmin_l_bfgs_b(tar_fun, [1, 1], args = [tt], bounds = [[0, 10], [0, 10]], approx_grad = True)
			#para, nn, cc = fmin_tnc(tar_fun, [0, 0], args = [tt], disp = False, bounds = [[0, 10], [0, 10]], approx_grad = True)
			tt.update(para[0], para[1])
			logl = tt.sum_llh()
			change = abs(logl - last_llh)
			last_llh = logl
		
		if print_detail:
			print("Num spe:" + repr(num_spe) + ": " + repr(tt.sum_llh()))
			print("spe_lambda:" + repr(tt.spe_rate))
			print("coa_lambda:" + repr(tt.coa_rate))
			print("spe_p:" + repr(tt.spe_p))
			print("coa_p:" + repr(tt.coa_p))
			print("-----------------------------------------------------")
		final_llh = tt.sum_llh()
		if final_llh > best_llh:
			best_llh = final_llh
			best_num_spe = num_spe
			best_node = tnode
		llh_list.append(final_llh)
	
	null_logl = optimize_null_model(utree)
	
	wt_list, num_spe = utree.get_waiting_times(threshold_node = best_node)
	one_spe, spes = utree.get_species()
	lrt = lh_ratio_test(null_llh = null_logl, llh = best_llh, df = 2)
	
	print("Highest llh:" + repr(best_llh))
	print("Num spe:" + repr(best_num_spe))
	print("Null llh:" + repr(null_logl))
	print("P-value:" + repr(lrt.get_p_value()))
	
	if show_lineages:
		utree.num_lineages(wt_list)
	
	if show_llh:
		plt.plot(llh_list)
		plt.ylabel('Log likelihood')
		plt.xlabel('Time')
		plt.savefig("Likelihood")
		plt.show()
	
	if print_species:
		utree.print_species()
	
	if show_tree:
		utree.tree.show()
	
	if lrt.get_p_value() >= pv:
		return one_spe
	else:
		return spes


def call_upgma(fin):
	basepath = os.path.dirname(os.path.abspath(__file__))
	call([basepath + "/bin/FastTree","-nt",fin], stdout=open(fin+".upgmaout", "w"), stderr=subprocess.STDOUT)
	fall = open(fin+".upgmaout")
	ss = fall.readlines()
	fall.close()
	fout = open(fin+".upgma", "w")
	fout.write(ss[2])
	fout.close()
	os.remove(fin+".upgmaout")
	return fin+".upgma"


def print_options():
	print("usage: ./GMYC.py -t example/gmyc_example.tre -st")
	print("usage: ./GMYC.py -a example/query.afa -st")
	print("Options:")
	print("    -a alignment                     Specify the alignment, GMYC.py will build a UPGMA tree using FastTree.\n")
	print("    -t input_umtree_file             Specify the input ultrametric tree.\n")
	print("    -pvalue (0-1)                    Set the p-value for likelihood ratio test with TWO degrees of freedom.(default 0.01)")
	print("                                     Note the old split package used three degrees of freedom which has been proven wrong.\n")
	print("    -st                              Plot the delimited species on the tree.(default not show)\n")
	print("    -ps                              Print delimited species on the screen.(default not show)\n")
	print("    -pd                              Print optimization details.(default not)\n")
	print("    -sl                              Show the log likelihood value plot.(default not)\n")
	print("    -sn                              Show lineages through time plot. (default not)\n")


if __name__ == "__main__":
	print("This is pGMYC - a python implementation of GMYC model for species delimitation.")
	print("Version 1.1 released by Jiajie Zhang on 10-11-2013\n")
	print("This program will delimit species on a rooted ultrametric tree, ")
	print("using single threshold GMYC model.")
	print("The input tree should be in Newick format and must be ultrametric.")
	print("Some common programs to infer ultrametric tree are: BEAST, DPPDIV and r8s." )
	print("pGMYC needs scipy and matplotlib packages to be installed.\n")
	print("**This new version experimentally support multifurcating tree, which ")
	print("**is quite common for many ultrametric tree inference programs.")
	print("**Note: pGMYC does not check the ultrametricity of the input tree! \n")
	print("--Please cite: \"J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: A General") 
	print("--Species Delimitation Method with Applications to Phylogenetic Placements. ")
	print("--Bioinformatics (2013), 29 (22): 2869-2876.\" ")
	print("--If you found pGMYC is useful to your research. \n")
	print("Questions and bug reports, please send to:")
	print("bestzhangjiajie@gmail.com\n")
	if len(sys.argv) < 3: 
		print_options()
		sys.exit()
	stree = ""
	sprint_detail = False 
	sshow_tree = False 
	sshow_llh = False 
	sshow_lineages = False 
	sprint_species = False
	p_value = 0.01
	salignment = ""
	
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-t":
			i = i + 1
			stree = sys.argv[i]
		if sys.argv[i] == "-pvalue":
			i = i + 1
			p_value = float(sys.argv[i])
		elif sys.argv[i] == "-pd":
			sprint_detail = True
		elif sys.argv[i] == "-st":
			sshow_tree = True
		elif sys.argv[i] == "-sl":
			sshow_llh = True
		elif sys.argv[i] == "-sn":
			sshow_lineages = True
		elif sys.argv[i] == "-ps":
			sprint_species = True
		elif sys.argv[i] == "-a":
			i = i + 1
			salignment = sys.argv[i]
		elif i == 0:
			pass
		elif sys.argv[i].startswith("-"):
			print("Unknown options: " + sys.argv[i])
			print_options()
			sys.exit()
	
	if salignment!="":
		basepath = os.path.dirname(os.path.abspath(__file__))
		if not os.path.exists(basepath + "/bin/FastTree"):
			print("GMYC.py uses FastTreeUPGMA to infer ultramatric trees,")
			print("please download the latest source code from: ")
			print("http://meta.microbesonline.org/fasttree/FastTreeUPGMA.c")
			print("Please complie with gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTreeUPGMA.c -lm, ")
			print("and put FastTree it to bin/  \n")
			sys.exit() 
		print("Building UPGMA tree using FastTree.")
		stree = call_upgma(salignment)
	
	if stree == "":
		print("Input tree is empty.")
		print_options()
		sys.exit()
	
	try:
		sp = gmyc(tree = stree, print_detail = sprint_detail, show_tree = sshow_tree, show_llh = sshow_llh, show_lineages = sshow_lineages, print_species = sprint_species, pv = p_value)
		print("Final number of estimated species by GMYC: " +  repr(len(sp)) )
	except ete2.parser.newick.NewickError:
		print("Unexisting tree file or Malformed newick tree structure.")


