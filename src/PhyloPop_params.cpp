
/*
 * PhyloPop_params.cpp
 *
 *  Created on: Jul 1, 2011
 *      Author: pickrell
 */

#include "PhyloPop_params.h"

PhyloPop_params::PhyloPop_params(){
	bias_correct = true;
	window_size = 1;
	alfreq_scaling =0;
	global = false;
	readtree = false;
	treefile = "NA";
	tau = 0.001;
	minweight = -15;
	maxweight = 10;
	nmig = 0;
	m_neigh = 3;
	maxit  = 100;
	maxit2 = 20;
	set_root = false;
	root = "NA";
	readgraph = false;
	vfile = "NA";
	efile = "NA";
	quick = false;
	min_migw = 0.001;
	nofrac = false;
	smooth_lik = true;
	smooth_scale = 1;
	nrand = 0;
	restrict_pop = false;
	pops2use = 0;
	sample_size_correct = true;
	calc_se = false;
	f2 = false;
	neg_penalty = 100;
	snpinfo = false;
	epsilon = 1e-3;
	nresid = 4;
	search_delta = 0.1;
	end_mig = false;
	target ="NA";
	dotarget = false;
	bootstrap = false;
	forcemig = false;
	forcemig_index = false;
	cov_snp = false;
	hzyfile = "";
	read_hzy = false;
	cor_f2 = false;
	cor_mig = false;
	climb = false;
	micro = false;
	fitmig = true;
	flip = false;
	print_hzy = false;
	seed = 0;
	num_warnings = -1; // warning suppression - number of warnings to show (-1 indicates show all)

	// Start of parameters added by EKM
	freq2stat = false;
	givenmat = false;
	matfile = "NA";
	refit = false;
	doscore = false;
	scoremethod = "rfit";
	givenpopaddorder = false;
	popaddorderfile = "NA";
	checkpoint = false;
	// End of parameters added by EKM
}

void PhyloPop_params::read_migfracs(string infile){
	migfracs.clear();
	ifstream in(infile.c_str());
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;
    intStat = stat(infile.c_str(), &stFileInfo);
     if (intStat !=0){
             std::cerr<< "ERROR: cannot open file " << infile << "\n";
             exit(1);
     }

    while(getline(in, st)){
    		// Start of addition by EKM
            int back = st.size() - 1;
            if (st[back] == '\r') {
                    cout << "ERROR: " << infile << " generated using Windows. Unable to read!\n";
                    exit(1);
            }
            // End of addition by EKM

            buf.clear();
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            string sourcep = line[0];
            string p = line[1];
            float f = atof(line[2].c_str());
            if (migfracs.find(p) == migfracs.end()){
            	map<string, double> tmpmap;
            	tmpmap.insert(make_pair(sourcep, f));
            	migfracs.insert(make_pair(p, tmpmap));
            }
            else{
            	migfracs[p].insert(make_pair(sourcep, f));
            	//migfracs.insert(make_pair(p, f));
            }
            hold.insert(p);
    }
}
