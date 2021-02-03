OrientAGraph
============

OrientAGraph is a built from [TreeMix](https://bitbucket.org/nygcresearch/treemix/src/master/) version 1.13 revision 231. Like the TreeMix code from J.K. Pickrell and J.K. Pritchard, OrientAGraph is provided under the [GNU General Public License v3.0](LICENSE).

OrientAGraph implements several new features on top of TreeMix, most notable of which is  **Maximum Likelihood Network Orientation (MNLO)**. When the `-mlno` option is provided, the MLNO is performed after the addition of each migration edge. In our experimental study, we found that this improved (or else did not impact) the accuracy of the original TreeMix method.

TreeMix is described in [Pickrell and Pritchard  (2012)](https://doi.org/10.1371/journal.pgen.1002967).

OrientAGraph is described in [this bioRxiv preprint](https://doi.org/10.1101/2021.02.02.429467).

# Installation
OrientAGraph has only been tested on Linux using gcc versions 4.8.5 and 4.9.3 and GSL version 2.6. 

1. If your system does not have GSL installed, then you will need to install it, for example with the following commands:
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
2. Export the environmental variables:
```
export INCLUDE_PATH="$HOME/gsl-2.6-local-install/include
export LIBRARY_PATH="$HOME/gsl-2.6-local-install/lib"
```
3. Then, download and build OrientAGraph.
```
cd ..
git clone https://github.com/ekmolloy/OrientAGraph.git
cd treemix-mlno
./configure CPPFLAGS=-I${INCLUDE_PATH} LDFLAGS=-L${LIBRARY_PATH}
make
```
4. Update `~/.base_profile` to contain the following lines:
```
export LD_LIBRARY_PATH="$HOME/gsl-2.6-local-install/lib:$LD_LIBRARY_PATH"
export PATH="$HOME/OrientAGraph/src:$PATH"
```
5. If everything has gone well, then typing
```
source ~/.bash_profile
treemix
```
should produce the help message:
```
OrientAGraph v1.0

OrientAGraph is built from TreeMix v1.13 Revision 231
by J.K. Pickrell and J.K. Pritchard and has several new
features, including the option to run Maximum Likelihood
Network Orientation (MLNO) as part of the admixture graph
search heuristic.

Contact: Erin Molloy (ekmolloy@cs.ucla.edu)

TreeMix Options:
-h display this help
-i [file name] input file
-o [stem] output stem (will be [stem].treeout.gz, [stem].cov.gz, [stem].modelcov.gz)
-k [int] number of SNPs per block for estimation of covariance matrix (1)
-global Do a round of global rearrangements after adding all populations
-tf [file name] Read the tree topology from a file, rather than estimating it
-m [int] number of migration edges to add (0)
-root [string] comma-delimited list of populations to set on one side of the root (for migration)
-g [vertices file name] [edges file name] read the graph from a previous TreeMix run
-se Calculate standard errors of migration weights (computationally expensive)
-micro microsatellite data
-bootstrap Perform a single bootstrap replicate
-cor_mig [file] list of known migration events to include (also use -climb)
-noss Turn off sample size correction
-seed [int] Set the seed for random number generation
-n_warn [int] Display first N warnings

OrientAGraph Options:
-mlno Run maximum likelihood network orientation subroutine as part of search heuristic
-allmigs Try all legal ways of adding a migration edge instead of using the minimum residual heuristic
-popaddorder [file with list of populations] Specify the order to add populations when building the starting tree
-givenmat [se matrix file] Allows user to input matrix (e.g. [stem].cov) with the -i flag, the matrix after this option should contain the standard error (e.g. [stem].covse); if no matrix is provided after this option, then 0.0001 is used.\n"
```
