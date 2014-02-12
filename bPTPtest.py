#! /usr/bin/env python
import json
import os 
import glob
import sys
import subprocess
import re
from ete2 import Tree
from ete2 import SeqGroup
from subprocess import call
from collections import deque
import numpy
from numpy import *
from scipy import *
from bPTP import bayesianptp
from summary import *


def rm_redudent_seqs_m(nfin, nfolder):
	call(["bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-f","c","-n","ck", "-T", "2", "-w", nfolder], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	if os.path.exists(nfin + ".reduced"):
		os.remove(nfolder + "RAxML_info." + "ck")
		return nfin +".reduced"
	else:
		return nfin


def build_tree(nfin, nfout, nfolder, num_thread = "2"):
	call(["bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234", "-T", num_thread, "-w", nfolder], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename(nfolder + "RAxML_bestTree."+nfout, nfolder + nfout + ".tre")
	os.remove(nfolder + "RAxML_info." + nfout)
	os.remove(nfolder + "RAxML_log." + nfout)
	os.remove(nfolder + "RAxML_parsimonyTree." + nfout)
	os.remove(nfolder + "RAxML_result." + nfout)
	return nfolder + nfout + ".tre"


def batch_bPTP(folder="./", suf = "phy", t = "2"):
	phyl = glob.glob(folder + "*." + suf)
	print(folder + "*." + suf) 
	rt_correct = 0
	rt_num = 0
	rt_nmi_b = 0.0
	rt_nmi_m = 0.0
	nmi_blist = []
	nmi_mlist = []
	
	for phy in phyl:
		fin1 = rm_redudent_seqs_m(nfin=phy, nfolder = folder)
		truth = ground_truth(fin1)
		taxaorder = truth.get_taxa_order()
		fin2 = build_tree(nfin = fin1, nfout = "temp", nfolder = folder , num_thread = t)
		
		bbptp = bayesianptp(filename = fin2, ftype = "raxml", 
		reroot = True, method = "H0", seed = 222, 
		thinning = 100, sampling = 10000, burnin = 0.1, taxa_order = taxaorder)
		
		pars, llhs = bbptp.delimit()
		pp = partitionparser(taxa_order = bbptp.taxa_order, partitions = pars, llhs = llhs)
		newpar_m = bbptp.get_maxhhl_partition()
		newpar_b = pp.summary(fout = folder + "tmp", region = 0.95)
		
		nmi_b = truth.get_nmi(array(newpar_b))
		nmi_m = truth.get_nmi(array(newpar_m))
		
		nmi_blist.append(nmi_b)
		nmi_mlist.append(nmi_m)
		os.remove(fin2)
		print("NMI_m: " + repr(nmi_m))
		print("NMI_b: " + repr(nmi_b))
		rt_nmi_b = rt_nmi_b + nmi_b
		rt_nmi_m = rt_nmi_m + nmi_m
	
	print("Average NMI bayesian: "  +  repr(rt_nmi_b/float(len(nmi_blist))))
	print("Average NMI ml: "  +  repr(rt_nmi_m/float(len(nmi_mlist))))


class ground_truth:
	def __init__(self, refaln, type = ""):
		if type == "fasta":
			self.aln = SeqGroup(sequences=refaln)
		else:
			self.aln = SeqGroup(sequences=refaln, format='phylip_relaxed')
		self.true_spe = {}
		self._get_truth()
		self._get_cluster_label()
		
	
	def _get_truth(self):
		for entr in self.aln.get_entries():
			name = entr[0]
			gid = name.split(".")[0]
			self.true_spe[gid] = []
		
		for entr in self.aln.get_entries():
			name = entr[0]
			gid = name.split(".")[0]
			group = self.true_spe[gid]
			group.append(name)
			self.true_spe[gid] = group
	
	def _get_cluster_label(self):
		self.seq_list = []
		self.seq_cid_list = [] 
		for entr in self.aln.get_entries():
			seq_name = entr[0]
			cid = int(seq_name.split(".")[0])
			self.seq_list.append(seq_name)
			self.seq_cid_list.append(cid)
		self.C0 = array(self.seq_cid_list) 
		
	
	def get_taxa_order(self):
		return self.seq_list
	
		
	def set_new_cluster_label(self, new_cid_list, seq_list, newid):
		if len(new_cid_list) == 0:
			for i in range(len(self.seq_list)):
				new_cid_list.append(-1)
		
		for i in range(len(self.seq_list)):
			name = self.seq_list[i]
			if name in seq_list:
				new_cid_list[i] = newid
		return new_cid_list
	
	#Mutual information
	def mutual_info(self,x,y):
		N=float(len(x))
		I=0.0
		eps = numpy.finfo(float).eps
		for l1 in numpy.unique(x):
			for l2 in numpy.unique(y):
				#Find the intersections
				l1_ids=nonzero(x==l1)[0]
				l2_ids=nonzero(y==l2)[0]
				pxy=(double(intersect1d(l1_ids,l2_ids).size)/N)+eps
				I+=pxy*log2(pxy/((l1_ids.size/N)*(l2_ids.size/N)))
		return I

	#Normalized mutual information
	def nmi(self,x,y):
		N=x.size
		I=self.mutual_info(x,y)
		Hx=0
		for l1 in unique(x):
			l1_count=nonzero(x==l1)[0].size
			Hx+=-(double(l1_count)/N)*log2(double(l1_count)/N)
		Hy=0
		for l2 in unique(y):
			l2_count=nonzero(y==l2)[0].size
			Hy+=-(double(l2_count)/N)*log2(double(l2_count)/N)
		if (Hx+Hy) == 0:
			return 1.0
		else: 
			return I/((Hx+Hy)/2)
	
	def get_seq_list(self):
		return self.seq_list
	
	def get_nmi(self, new_cluster_labels):
		return self.nmi(self.C0, new_cluster_labels)
	
	def is_correct(self,names):
		#*R*
		newnames = []
		for name in names:
			if name.startswith("*R*"):
				pass
			else:
				newnames.append(name)
			
		names_set = set(newnames)
		for key in self.true_spe.keys():
			sps = self.true_spe[key]
			sps_set = set(sps)
			if names_set == sps_set:
				return True
		return False
	
	def get_num_species(self):
		return len(self.true_spe.keys())


def phy2fasta(fin):
	f = open(fin)
	fout = open(fin+".afa", "w")
	line = f.readline()
	line = f.readline()
	while line!="":
		lls = line.split()
		fout.write(">" + lls[0] + "\n")
		fout.write(lls[1] + "\n")
		line = f.readline().strip()
	f.close()
	fout.close()
	return fin + ".afa"


if __name__ == "__main__":
	
	#batch_bPTP(folder="/home/zhangje/Research/PTP/bPTP/simulation/even/set5/l1000/", suf = "phy", t = "2")
	
	if len(sys.argv) < 3: 
		print("usage: ./bPTPtest.py -folder <folder contain test files>")
		sys.exit()
	
	sfolder = "./"
	ssuf = "phy"
	snum_spe = 10
	method = "gmyc"
	sfout = "ptp_raxml_log.txt"
	sT = "1"
	pvflag = False
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-folder":
			i = i + 1
			sfolder = sys.argv[i]
		elif sys.argv[i] == "-suf":
			i = i + 1
			ssuf = sys.argv[i]
		elif sys.argv[i] == "-num_spe":
			i = i + 1
			snum_spe = int(sys.argv[i])
		elif sys.argv[i] == "-m":
			i = i + 1
			method = sys.argv[i]
		elif sys.argv[i] == "-o":
			i = i + 1
			sfout = sys.argv[i]
		elif sys.argv[i] == "-T":
			i = i + 1
			sT = sys.argv[i]
		elif sys.argv[i] == "-pvflag":
			pvflag = True
	
	batch_bPTP(folder=sfolder, suf = "phy", t = "2")


