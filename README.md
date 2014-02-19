A Bayesian implementation of the Poisson tree processes (PTP) model for 
species delimitation.

By Jiajie Zhang 19-02-2014.
Questions and bugs report please sent to bestzhangjiajie[at]gmail[dot]com.
Before reading the following text, please check: 
https://github.com/zhangjiajie/SpeciesCounting for latest updates

=========================================================================
(0) What's new
    
    I implemented a Markov chain Monte Carlo sampler for the PTP model - 
    bPTP.py, it can generate posterior probability of species delimitations 
    on phylogenetic trees. Since I give a flat prior to all possible 
    delimitations, bPTP.py now replaces the old maximal likelihood search 
    program PTP.py. bPTP.py searches a much large space and gives Bayesian
    support values to delimitations. PTP.py is now used for bootstrap
    analysis if you have trees from bootstrap analysis. The new programs 
    also feature various improvements, such as NEXUS format support, 
    rooting on outgroups, SVG tree output and code optimizations. 


(1) What's in the package?

    This package contains several programs written in Python that can give
    species delimitation hypothesis based on a gene trees inferred from 
    molecular sequences. We introduce a new model called Poisson tree 
    processes(PTP) model. In PTP, we model speciations or branching events 
    in terms of number of mutations. So it only requires a phylogenetic 
    input tree, for example the output of RAxML. To be more clear, the 
    branch lengths should represent number of mutations. Our numerous tests
    show PTP outperforms GMYC. Furthermore, PTP is much easier to use, 
    since it can use the phylogenetic tree directly without needing the 
    difficult and error prone procedures of time calibration required by 
    GMYC.
    
    bPTP.py    A Bayesian implementation of the PTP model for species 
               delimitation. It uses a Markov chain Monte Carlo sampler 
               to produce posterior probability of species delimitations 
               on phylogenetic trees. bPTP.py can work both on a single 
               maximal likelihood phylogenetic tree, and multiple trees 
               from Bayesian phylogenetic analysis. Using trees from 
               Bayesian phylogenetic analysis can account for uncertainties 
               in phylogenetics inference, but with the cost of much longer
               run time. 
               To find out how to use it, type: python bPTP.py  
    
    PTP.py     Maximal likelihood search of the PTP model by heuristics, 
               it can also provide bootstrap support values of species 
               delimitation if given multiple phylogenetic trees derived 
               from bootstrap analysis.
               To find out how to use it, type: python PTP.py  
               
    GMYC.py    Implements the single threshold general mixed Yule coalescent 
               (GMYC) model which was first proposed by Pons et al 
               (Sequence-Based Species Delimitation for the DNA Taxonomy 
               of Undescribed Insects. Systematic Biology, 55(4), 595–609)
               This model requires the input tree to be time calibrated 
               ultrametric tree, in other words, the branch length of the 
               input ultrametric tree should represent time. 
               The most commonly used programs for getting an ultrametric
               tree are BEAST, DPPDIV and r8s. There is also an R 
               implementation of this model called "splits" by Tomochika 
               Fujisawa (http://barralab.bio.ic.ac.uk/downloads.html)
               To find out how to use it, type ./GMYC.py
             
    EPA_PTP.py (Broken now, I will fix it when I have time)This is a pipeline
               that uses evolutionary placement algorithm (EPA) and PTP to 
               count species number when reference data is available. For 
               details of EPA, please read this paper: Performance, accuracy, 
               and Web server for evolutionary placement of short sequence 
               reads under maximum likelihood. Systematic biology, 60(3), 
               291–302. The pipeline will first run USEARCH to remove the 
               chimera sequences, then it will use EPA to place the query 
               reads to optimal position on the reference tree inferred 
               from the reference alignment. PTP will then be applied to 
               the reads been placed on each branch, with a fixed speciation
               rate inferred from the reference data.  
               Similar analysis used in bacterial metagenomics studies are 
               called OTU-picking. For discussions about OTU-picking and EPA_PTP 
               species counting, please have a look at our paper. 
               Type ./EPA_PTP.py for help and instructions.
 

(2) Which operating system is required?

    I wrote and tested all the python code under Ubuntu Linux. So everything
    should run well under Linux if you follow the instructions below or 
    the output of the python programs. I did not have the time and chance 
    to test them under windows or mac yet, however, I think bPTP.py, PTP.py
    and GMYC.py should be able to run on windows and mac if you have properly
    installed the dependent python packages (see below). EPA_PTP.py was 
    designed to run with NGS data, which means the calculations might be 
    intense if you have say 10,000 reads. So ideally it should be run on a 
    multi-core Linux server such that it can speedup using the PTHREADS 
    version of RAxML. The biodiversity soup data in our paper, for example, 
    will need 24-48 hour to finish on our 8-core i7 server. If you encounter
    any problems to run the program under Linux, simply drop me an e-mail.  
               


(3) Install dependent python packages

    The programs used ETE package (http://ete.cgenomics.org/) for tree 
    manipulations, and some functions from scipy and matplotlib. I included 
    a copy of ETE package, so there is no need for seperate installation, 
    however, ETE is dependent on some python packages, The following 
    python packages are needed:  python-setuptools python-numpy python-qt4
    python-scipy python-mysqldb python-lxml python-matplotlib 
    
    If you are running Ubuntu, or Debian GNU/Linux distribution, 
    you can try the following:
    
    sudo apt-get install python-setuptools python-numpy python-qt4 
    python-scipy python-mysqldb python-lxml python-matplotlib


(4) Prepare input data

    bPTP.py: accept phylogenetic trees in Newick format and NEXUS format. 
             The trees should not be annotated, if bPTP.py can not parse 
             your tree, you can try to use FigTree to convert the tree 
             format. If you input NEXUS format, your trees will be automatically 
             detranslated if possible. However, if taxon are named using 
             only numbers, this step is will cause errors and you should 
             rename your taxa. 
             If a single tree is used as input, this should ideally be a 
             maximal likelihood tree inferred by, for example RAxML or 
             PhlML. Alternatively, this single tree can also be the consensus 
             tree from Bayesian or bootstrap phylogenetic analysis. 
             If multiple trees are used, those trees should come from 
             Bayesian phylogenetic analysis programs such as MrBayes.
    
    PTP.py: accept the same tree format as bPTP.py, both single and multiple 
            trees can be used as input. Multiple trees should come from 
            bootstrap phylogenetic analysis.
    
    GMYC.py:input tree to GMYC must be strictly ultrametric (I do not check 
            for this!) and in Newick format only.


(5) Output



(6) Examples

    python bPTP.py -t example/ptp_example.tre -o example/myoutput -s 1234 -r


(7) How to cite

    If you find PTP and bPTP useful to your research, please cite: 
    J. Zhang, P. Kapli, P. Pavlidis, A. Stamatakis: "A General Species 
    Delimitation Method with Applications to Phylogenetic Placements". 
    In Bioinformatics (2013), 29 (22): 2869-2876.


(8) Web server

    There is also a web server avaliable for bPTP:
    http://species.h-its.org/ptp/
    The server version only accept a single input tree and has limitation
    on the number of MCMC iterations.  


(4) Download and compile required programs for EPA_PTP pipeline

    The EPA_PTP pipeline requires the following three programs to run, 
    you can download, compile and put them in the included bin folder. 
    I included the binary executable of RAxML and HMMMER for 64-bit Linux 
    in the bin folder, However, USEARCH does NOT allow for redistribution. 
    USEARCH provides a free binary 32-bit version that will also run on 
    64-bit platform, but must be requested per e-mail. 
    
    a. USEARCH: http://www.drive5.com/usearch/ please rename the executable
       file to "usearch", and copy to bin/ folder
    b. HMMER: http://hmmer.janelia.org/ please copy "hmmbuild" and 
       "hmmalign" to bin/ folder
    c. RAxML: https://github.com/stamatak/standard-RAxML please compile 
       the PTHREADS version, rename to "raxmlHPC-PTHREADS-SSE3" and copy to bin/
