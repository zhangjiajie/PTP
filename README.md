Species delimitation using a Poisson tree processes model
           by Jiajie Zhang 03-04-2013
           bestzhangjiajie@gmail.com
           
before reading the following text, please check:
https://github.com/zhangjiajie/SpeciesCounting
for latest updates
========================================================================================================================================
(1) What's in the package?
    This package contains serveral programs written in Python that can give species delimitation hypothesis 
    based on a gene tree inferred from molecular sequences.
    
    GMYC.py    Inplements the single threshold general mixed Yule coalescent (GMYC) model which was first proposed by Pons et al 
               (Sequence-Based Species Delimitation for the DNA Taxonomy of Undescribed Insects. Systematic Biology, 55(4), 595–609)
               This model requires the input tree to be time calibrated ultrametric tree, in other words, the branch length of the 
               input ultrametric tree should reprent time. 
               The most commonly used programs for getting an ultrametric tree are BEAST, DPPDIV and r8s.
               There is also an R inplementation of this model called "splits" by Tomochika Fujisawa
               (http://barralab.bio.ic.ac.uk/downloads.html)
               To find out how to use it, type ./GMYC.py
             
    PTP.py     This is a newly introduced model we call it Poisson tree processes(PTP) model. In PTP, we model speciations or branching 
               events in terms of number of mutations. So it only requires a phylogenetic input tree, for example the output of RAxML. 
               To be more clear, the branch lengths should represent number of mutations. Our numerous tests show PTP results are comparable 
               to GMYC in most cases, and still give reasonable results when GMYC fails. However, PTP is much easier to use, since it can 
               use the phylogenetic tree directly without needing the difficult and error prone prcedures of time calibration required by GMYC.
               To find out how to use it, type ./PTP.py  
             
    EPA_PTP.py This a pipeline that uses evolutionary placement algorithm (EPA) and PTP to count species number when reference data is 
               avaliable. For details of EPA, please read this paper: Performance, accuracy, and Web server for evolutionary placement of short 
               sequence reads under maximum likelihood. Systematic biology, 60(3), 291–302.
               The pipeline will first run USEARCH to remove the chimera seqeunces, the it will use EPA to place the query reads to optimal 
               position on the reference tree inferred from the reference alignment. PTP will then be applied to the reads been placed on each 
               branch, with a fixed speciation rate inferred from the reference data.  
               Similar analysis used in bacterial metagenomics studies are called OTU-picking. For discussions about OTU-picking and EPA_PTP 
               species counting, please have a look at our paper. 
 

(2) Which operating system is required?
    I wrote and tested all the python code under Ubuntu Linux. So everything should run well under Linux if you follow the instructions below or 
    the output of the python programs. I do not have the time and chance to test them under windows or mac yet, however, I think GMYC.py and PTP.py 
    should be able to run on windows and mac if you have properly installed the dependent python packages (see below). EPA_PTP.py was desinged to 
    run with NGS data, which means the calculations might be intense if you have say 10,000 reads. So ideally it shoud be run on a multi-core 
    Linux server such that it can speedup using the PTHREADS version of RAxML. The biodiversity soup data in our paper, for example, will need 
    24-48 hour to finish on our 8-core i7 server. If you encounter any problems to run the program under Linux, simply drop me an e-mail.  
               


(3) Install dependent python packages
    1.Dependent
    sudo apt-get install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml python-matplotlib

    2.ETE2
    sudo easy_install -U ete2


(4) Download and complie required programs for EPA_PTP pipeline



(5) Important notes on the input data
    a. All trees must be in Newwick format.
    b. All sequences must be in Fasta format.
    c. The input tree to GMYC must be strictely bifercating and ultrametric.
    d. The input tree to PTP should ideally be rooted with some outgroups, if an unrooted tree is used, pleae specify the -r option.
    e. The input to EPA_PTP.py should be two alignments, one should be the query sequences that you want to know how many species are there;
       the other should be the reference alignment, which should contain ONLY ONE sequence for each known species. The reference data represents our 
       knowledge about the speciation history, so EPA_PTP will estimate the speciation rate from it. If multiple sequences exist for a single
       species, then the speciation rate will be over-estimated.
    f. If your query sequences are not aligned, EPA_PTP.py can align them using HMMER, however, the reference sequences must be perperly aligned.
       Both reference and query alignment should be of the same length, if not, please use the same option in EPA_PTP.py to align the query sequences 
       to the reference alignment.  
    

(6) Examples



