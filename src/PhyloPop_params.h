/*
 * PhyloPop_params.h
 *
 *  Created on: Jul 1, 2011
 *      Author: pickrell
 */

#ifndef PHYLOPOP_PARAMS_H_
#define PHYLOPOP_PARAMS_H_

#include "Settings.hpp"

class PhyloPop_params{
public:
	PhyloPop_params();
	bool bias_correct, global, readtree;
	string treefile;
	int window_size;
	int alfreq_scaling; // 0 = no scaling, 1 = asin(sqrt(f))

	//optimization of weights
	double tau; //stopping criterion
	double minweight, maxweight;
	int maxit, maxit2; //maximum number of iterations when maximizing weights

	int nmig; //number of migration edges to add

	int m_neigh; //"neighborhood" size for addition of migration events

	// setting root
	bool set_root;
	string root;

	// read from previous run
	bool read_graph;
	string vfile;
	string efile;

	//quick optimization of weights
	bool quick;
	double min_migw;

	bool nofrac;


	bool smooth_lik;
	double smooth_scale;

	int nrand;

	bool print_hzy;
	//only use a certain number of populations
	bool restrict_pop;
	int pops2use;

	//
	bool sample_size_correct;
	bool calc_se;

	bool f2;

	//penalty for negative branch lengths
	double neg_penalty;

	set<string> hold;

	bool snpinfo;
	double epsilon; //for determining whether the likelihood is being increased

	int nresid; // number of maximum residuals to search when adding migration
	double search_delta; // in golden section search, do "backup search" in window of size search_delta around current weight

	bool end_mig; //optimize migration weights at end

	bool dotarget;
	string target; //population to target

	bool bootstrap; // bootstrap from covariance matrices

	pair<string, string> mig_pops;
	bool forcemig;

	pair<int, int> mig_index;
	bool forcemig_index;

	bool cov_snp;
	int which_cov_snp;

	//read heterozygosity from a file (rather than estimate it)
	bool read_hzy;
	string hzyfile;

	//for correcting f2 stats [IGNORE]
	bool cor_f2;
	map<string, map<string, double> > migfracs;
	string f2_corpop;
	double f2_mixdist;
	void read_migfracs(string);

	//same but for forcing migration during tree building
	bool cor_mig;
	bool fitmig;
	//string corpop;
	bool climb; //do hill climbing

	//microsats
	bool micro;

	bool flip;
	string flipstring;
	// seed random number generator
	unsigned int seed;

	// warning suppression - define number of warnings to show (-1 indicates show all)
	int num_warnings;
};

#endif /* PHYLOPOP_PARAMS_H_ */
