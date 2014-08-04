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
               phylogenetic tree, and multiple trees from Bayesian 
               phylogenetic analysis. Using trees from Bayesian phylogenetic
               analysis can account for uncertainties in phylogenetics 
               inference, but with the cost of much longer run time. 
               To find out how to use it, type: python bPTP.py  
    
    PTP.py     Maximal likelihood search of the PTP model by heuristics, 
               it can also provide bootstrap support values of species 
               delimitation if given multiple phylogenetic trees derived 
               from bootstrap analysis.
               To find out how to use it, type: python PTP.py  
               

 

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
    
    MAC: for mac users, first you have to make sure you have python installed, 
    go to the terminal and type python. Once you have python, you can try to use 
    easy_install to install required packages:
    
    

(4) Prepare input data

    bPTP.py: accept phylogenetic trees in Newick format and NEXUS format. 
             The trees should NOT be annotated, if bPTP.py can not parse 
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


(5) Output

    a) outputname.PTPPartitions.txt: bPTP - MCMC samples of delimitations 
                                     after thinning, all posterial probabilities 
                                     are computed based this file. 
                                     PTP - bootstrap delimitations.
                                     
    b) outputname.PTPllh.txt: Posterial Log likelihood trace file.
    
    c) outputname.llh.pdf: Posterial Log likelihood trace plot, visual 
                           check for convergence. 
                           
    d) outputname.PTPPartitonSummary.txt: bPTP - summary of Posterial Prob. of 
                                          delimited species.
                                          PTP - bootstrap values of delimited 
                                          species.
                                          
    e) outputname.PTPhSupportPartition.txt: bPTP - highest posterial Prob. 
                                            supported delimitation.
                                            PTP - highest bootstrap
                                            supported delimitation.
                                            
    f) outputname.PTPhSupportPartition.txt.png/svg: Tree plot from e
    
    g) outputname.PTPMLPartition.txt: Maximal likelihood species delimitation
    
    h) outputname.PTPMLPartition.txt.png/svg: Tree plot from g
    
    
    bPTP.py:
    Support values shown on the tree plot are computed as the "number of 
    occurrence of all the descendants under this node"/ "number of samples
    from MCMC sampling". They are the posterial probabilities of those taxa 
    form one species under the PTP model and a flat prior. From tests on 
    simulated data, support values are strongly correlated with the accuracy 
    of the delimitation, r = 0.91.  
    If input a single tree, output will be a - h; if using multiple 
    trees, output will be a - e 
    
    PTP.py: 
    Support values shown on the tree plot are the bootstrap support of 
    all taxa under this node form one species under the PTP model.
    Output will always be: a, d, e and f
     

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

