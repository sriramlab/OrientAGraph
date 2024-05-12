OrientAGraph
============

OrientAGraph implements **Maximum Likelihood Network Orientation (MNLO)** within [TreeMix](https://doi.org/10.1371/journal.pgen.1002967), a popular package for estimating admixture graphs from f-statistics (and related quantities). In our experimental study, we found that MLNO either improved or else did not impact the accuracy of the original TreeMix search heuristic. To learn more, check out [this paper](https://doi.org/10.1093/bioinformatics/btab267) with Arun Durvasula and Sriram Sankararaman. 

**Starting in version 1.2:** By default, OrientAGraph searches for the MLNO only after each of the first two admixture edges are added (equivalent to using the flag: `-mlno 1,2`). To get started, see the installation instructions below as well as [this example](example/arctic-data/README.md).

Other Useages
-------------
+ To execute OrientAGraph as described in our original paper, you must include the `-allmigs` and `-mlno` flags to the end of your command. These options expand the maximum likelihood search but may be prohibitively expensive for large populations!
+ Based on our experience, it is important to use the `-root <outgroup>` flag. This option roots the starting tree at the user-specified outgroup before starting the network search.
+ To determine how to set the other options, we recommend reading the [TreeMix manual](https://bitbucket.org/nygcresearch/treemix/downloads/). 
+ There are also some options specific to OrientAGraph. As an example, OrientAGraph can be used to find the MLNO of a user-provided graph (option: `-gf <vertex file> <edge file> -score mlno`).


Acknowledgements
----------------
OrientAGraph is built from the [TreeMix code](https://bitbucket.org/nygcresearch/treemix/src/master/) by J.K. Pickrell and J.K. Pritchard, and like the TreeMix code, is provided under the [GNU General Public License v3.0](LICENSE). 

TreeMix is presented in [Pickrell and Pritchard (2012)](https://doi.org/10.1371/journal.pgen.1002967) and in [Pickrell et al. (2012)](https://doi.org/10.1038/ncomms2140).

OrientAGraph implements algorithms / utilizes theoretical results from [Huber et al. (2019)](https://arxiv.org/abs/1906.07430) and [Francis and Steel (2015)](https://doi.org/10.1093/sysbio/syv037).

OrientAGraph has only been tested on Mac (using Apple clang version 12.0.0) and Linux (using gcc versions 4.8.5 and 4.9.3), both with GSL version 2.6.


Installation on MAC OS with [homebrew](https://brew.sh)
-------------------------------------------------------

1. Install GSL (we used version 2.7.1)
```
brew install gsl
```

2. Export the path information for GSL.
```
GSL_VERSION=$(ls /opt/homebrew/Cellar/gsl/)
export INCLUDE_PATH="/opt/homebrew/Cellar/gsl/${GSL_VERSION}/include"
export LIBRARY_PATH="/opt/homebrew/Cellar/gsl/${GSL_VERSION}/lib"
```
Now check if these paths are working by typing
```
ls $INCLUDE_PATH
```
which should return `gsl` and typing
```
ls $LIBRARY_PATH
```
which should return `libgsl.a` among other files.

3. Install BOOST (we used version 1.83.0).
```
brew install boost
```

4. Export the path information for BOOST.
```
BOOST_VERSION=$(ls /opt/homebrew/Cellar/boost/)
export BOOST_PATH="/opt/homebrew/Cellar/boost/${BOOST_VERSION}"
```
Now check to see if these paths are working by typing
`ls $BOOST_PATH`
which should return `include` and `lib`, among other files.

5. Now download and build OrientAGraph.
```
git clone https://github.com/ekmolloy/OrientAGraph.git
cd OrientAGraph
LDFLAGS="-static"
./configure CPPFLAGS=-I${INCLUDE_PATH} LDFLAGS="-L${LIBRARY_PATH}" --with-boost="${BOOST_PATH}"
make
```

6. If everything has gone well, typing
```
source ~/.bash_profile
orientagraph
```
will produce the help message:
```
OrientAGraph 1.2

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
    e.g. '1,2' means run MLNO only after adding the first two migration edges (default)
         '0' means do NOT run MLNO
         '' (no string) means run MLNO after adding each migration edge
-allmigs [string] Comma-delimited list of integers, indicating when to run
    evaluate all legal ways of adding migration edge to base tree instead of
    using heuristic (similar to -mlno but default is -allmigs 0)
-popaddorder [population list file] Order to add populations when building
    starting tree
-checkpoint Write checkpoint files
```

You can also check the compile by typing
```
file src/orientagraph 
```
which should produce something like
```
src/orientagraph: Mach-O 64-bit executable arm64
```

7. To use orientagraph from any directory on your system, you can add its path to your profile. Check the path with the `pwd` command, then add the following line 
```
export PATH=<result of pwd>/src:$PATH"
```
is for the bash profile. Similar commands exist for other profiles.

Installation on Linux 
----------------------
A similar installation can be done for linux but you will need to use apt-get (or some other tool) instead of homebrew. If you are using a Linux system, you can ask your system admin about how to access GSL and BOOST; however, you may not be able to compile static binaries then.

Installation from source
------------------------
1. Build and install BOOST.
```
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
tar -zxvf boost_1_82_0.tar.gz
cd boost_1_82_0
mkdir install
BOOST_PATH="$(pwd)/install"
./bootstrap.sh --prefix="$BOOST_PATH"
./b2 install
```
Check there are static libs
```
ls ${BOOST_PATH}/lib/*a
```
should return a bunch of files.

2. Build and install GSL.
```
wget https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz
tar -zxvf gsl-2.7.1.tar.gz
cd gsl-2.7.1
mkdir install
GSL_PATH="$(pwd)/install"
export INCLUDE_PATH="${GSL_PATH}/include"
export LIBRARY_PATH="${GSL_PATH}/lib"
./configure --prefix="$GSL_PATH" 
make
make check
make install
```
Check there are static libs
```
ls $LIBRARY_PATH/*a
```
should return `libgsl.a` and `libgslcblas.a`.

3. Download and build OrientAGraph
```
git clone https://github.com/ekmolloy/OrientAGraph.git
cd OrientAGraph
LDFLAGS="-static"
./configure CPPFLAGS=-I${INCLUDE_PATH} LDFLAGS="-L${LIBRARY_PATH} -static" --with-boost="${BOOST_PATH}"
make
```

4. Check if the build is static with command
```
ldd ./src/orientagraph
```
should return `not a dynamic executable`
