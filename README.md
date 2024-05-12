OrientAGraph
============

OrientAGraph implements **Maximum Likelihood Network Orientation (MNLO)** within [TreeMix](https://doi.org/10.1371/journal.pgen.1002967), a popular package for estimating admixture graphs from f-statistics (and related quantities).
OrientAGraph can be used to find the MLNO of a user-provided graph (option: `-gf <vertex file> <edge file> -score mlno`) or incorporated into TreeMix's heursitic search for the best fitting admixture graph (option:  `-mlno` ).
In an experimental study, we found that MLNO improved (or else did not impact) the accuracy of the original TreeMix search heuristic.
The current implementation exhaustively searches for the MLNO, and thus, we expect it to be computationally intensive on very large admixture graphs; in this case,  MLNO could be run only after the addition of the first two admixture edges (option:  `-mlno 1,2` ).
To learn more, check out [this example](example/README.md) and [this paper](https://doi.org/10.1093/bioinformatics/btab267) with Arun Durvasula and Sriram Sankararaman.


Acknowledgements
----------------
OrientAGraph is built from the [TreeMix code](https://bitbucket.org/nygcresearch/treemix/src/master/) by J.K. Pickrell and J.K. Pritchard, and like the TreeMix code, is provided under the [GNU General Public License v3.0](LICENSE). TreeMix is presented in [Pickrell and Pritchard (2012)](https://doi.org/10.1371/journal.pgen.1002967) and in [Pickrell et al. (2012)](https://doi.org/10.1038/ncomms2140).

OrientAGraph implements algorithms / utilizes theoretical results from [Huber et al. (2019)](https://arxiv.org/abs/1906.07430) and [Francis and Steel (2015)](https://doi.org/10.1093/sysbio/syv037).


Installation
------------
OrientAGraph has only been tested on Mac (using Apple clang version 12.0.0) and Linux (using gcc versions 4.8.5 and 4.9.3), both with GSL version 2.6. 

1. If your system does not have GSL installed, then you will need to install it.
This can also be done manually using the following commands:
```
cd $HOME
mkdir gsl-2.6-local-install
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
cd gsl-2.6
./configure --prefix=$HOME/gsl-2.6-local-install
make
make check
make install
```
2. After GSL is installed, export the related environmental variables.
If GSL was installed manually using the commands above, then export:
```
export INCLUDE_PATH="$HOME/gsl-2.6-local-install/include"
export LIBRARY_PATH="$HOME/gsl-2.6-local-install/lib"
```
If you are using Mac systems, you could alternatively perform step 1 using [homebrew](https://brew.sh): 
```
brew install gsl
```
When this README was created, homebrew installed GSL with the following paths:
```
export INCLUDE_PATH="/usr/local/Cellar/gsl/2.6/include"
export LIBRARY_PATH="/usr/local/Cellar/gsl/2.6/lib"
```
Now it installs GSL with the following paths:
```
export INCLUDE_PATH=/opt/homebrew/Cellar/gsl/2.6/include
export LIBRARY_PATH=/opt/homebrew/Cellar/gsl/2.6/lib
```
although note that you may also need to update 2.6 depending on the GSL version installed by Homebrew.
As of the updating of this README, you also need to install boost when using homebrew:
```
brew install boost
```
3. Now download and build OrientAGraph.
```
cd ..
git clone https://github.com/ekmolloy/OrientAGraph.git
cd OrientAGraph
./configure CPPFLAGS=-I${INCLUDE_PATH} LDFLAGS=-L${LIBRARY_PATH}
make
```
If you used homebrew, additionally include the flag:
```
--with-boost="/opt/homebrew/Cellar/boost/1.82.0_1
```
Note that you need to include the flag `LDFLAGS="-static"` option to the configure commands for the binaries to be compiled statically.
4. Update `~/.bash_profile` to include 
```
export LD_LIBRARY_PATH="$HOME/gsl-2.6-local-install/lib:$LD_LIBRARY_PATH"
```
if you installed GSL from source.
Additionally, include 
```
export PATH="$HOME/OrientAGraph/src:$PATH"
```
if you want to be able to use orientagraph in any directory on your system.

5. If everything has gone well, typing
```
source ~/.bash_profile
orientagraph
```
will produce the help message:
```
OrientAGraph 1.1

OrientAGraph is built from TreeMix v1.13 Revision 231 by
J.K. Pickrell and J.K. Pritchard and implements several new
features, including Maximum Likelihood Network Orientation
(MLNO), which can be used as a graph search heuristic.

Contact: Erin Molloy (ekmolloy@umd.edu)

COMMAND:  ./src/orientagraph -h

TreeMix Options:
-h Display this help
-i [file] Input file (e.g. containing allele frequencies)
-o [stem] Output prefix (i.e. output will be [stem].treeout.gz,
    [stem].cov.gz, [stem].modelcov.gz, etc.)
-k [int] Number of SNPs per block for estimation of covariance matrix (1)
-global Do a round of global rearrangements after adding all populations
-tf [newick file] Read tree from a file, rather than estimating it
-m [int] Number of migration edges to add (default: 0)
-root [string] Comma-delimited list of populations to put on one side of root
-gf [vertices file] [edges file] Read graph from files (e.g. [stem].vertices.gz
    and [stem].edges.gz from a previous TreeMix run)
-se Calculate standard errors of migration weights (computationally expensive)
-micro Input is microsatellite data
-bootstrap Perform a single bootstrap replicate
-cor_mig [file] List of known migration events to include (also use -climb)
-noss Turn off sample size correction
-seed [int] Set the seed for random number generation
-n_warn [int] Display first N warnings

Options added for OrientAGraph:
-freq2stat Estimate covariances or f2-statistics from allele frequencies
    and then exit;
    the resulting files can be given as input using the -givenmat option
-givenmat [matrix file] Allows user to input matrix (e.g. [stem].cov.gz)
    with the -i flag, the file after this flag should contain the standard
    error (e.g. [stem].covse.gz)
-refit Refit model parameters on starting tree (-tf) or graph (-gf)
-score [string] Score input tree (-tf) or graph (-gf) and then exit:
    'asis' = score 'as is' i.e. without refitting,
    'rfit' = score after refitting (default),
    'mlbt' = score each base tree (with refitting) and return best,
    'mlno' = score each network orientation (with refitting) and return best
-mlno [string] Comma-delimited list of integers, indicating when to run
    maximum likelihood network orientation (MLNO) as part of heuristic search
    (e.g. '1,2' means run MLNO only after adding the first two migration edges
    and no string means run MLNO after adding each migration edge)
-allmigs [string] Comma-delimited list of integers, indicating when to run
    evaluate all legal ways of adding migration edge to base tree instead of
    using heuristic
-popaddorder [population list file] Order to add populations when building
    starting tree
-checkpoint Write checkpoint files
```

Recommended Usage
------------
+ To execute OrientAGraph as described in our paper, you must include the `-allmigs` and `-mlno` flags to the end of your command. These options expand the maximum likelihood search!
+ To determine how to set the other options, we recommend reading the [TreeMix manual](https://bitbucket.org/nygcresearch/treemix/downloads/). Based on our experience, it is important to use the `-root <outgroup>` flag. This option roots the starting tree at the user-specified outgroup before starting the network search.
