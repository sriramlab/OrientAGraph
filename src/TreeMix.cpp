/*
 * TreeMix.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "CountData.h"
#include "GraphState2.h"
#include "PhyloPop_params.h"

string infile;
string outstem = "TreeMix";

void printv() {
    cout << "OrientAGraph 1.2\n\n"
         << "OrientAGraph is built from TreeMix v1.13 Revision 231 by\n"
         << "J.K. Pickrell and J.K. Pritchard and implements several new\n"
         << "features, including Maximum Likelihood Network Orientation\n"
         << "(MLNO), which can be used as a graph search heuristic.\n\n"
         << "Contact: Erin Molloy (ekmolloy@umd.edu)\n\n";
}

void printopts() {
    cout << "TreeMix Options:\n";
    cout << "-h Display this help\n";
    cout << "-i [file] Input file (e.g. containing allele frequencies)\n";
    cout << "-o [stem] Output prefix (i.e. output will be [stem].treeout.gz,\n"
         << "    [stem].cov.gz, [stem].modelcov.gz, etc.)\n";
    cout << "-k [int] Number of SNPs per block for estimation of covariance matrix (1)\n";
    cout << "-global Do a round of global rearrangements after adding all populations\n";
    cout << "-tf [newick file] Read tree from a file, rather than estimating it\n";
    cout << "-m [int] Number of migration edges to add (default: 0)\n";
    cout << "-root [string] Comma-delimited list of populations to put on one side of root\n";
    cout << "-gf [vertices file] [edges file] Read graph from files (e.g. [stem].vertices.gz\n"
         << "    and [stem].edges.gz from a previous TreeMix run)\n";
    cout << "-se Calculate standard errors of migration weights (computationally expensive)\n";
    cout << "-micro Input is microsatellite data\n";
    cout << "-bootstrap Perform a single bootstrap replicate\n";
    cout << "-cor_mig [file] List of known migration events to include (also use -climb)\n";
    cout << "-noss Turn off sample size correction\n";
    cout << "-seed [int] Set the seed for random number generation\n";
    cout << "-n_warn [int] Display first N warnings\n"; 
    // Start of additions by EKM
    cout << "\nOptions added for OrientAGraph:\n";
    cout << "-freq2stat Estimate covariances or f2-statistics from allele frequencies\n"
         << "    and then exit;\n"
         << "    the resulting files can be given as input using the -givenmat option\n";
    cout << "-givenmat [matrix file] Allows user to input matrix (e.g. [stem].cov.gz)\n"
         << "    with the -i flag, the file after this flag should contain the standard\n"
         << "    error (e.g. [stem].covse.gz)\n";
    cout << "-refit Refit model parameters on starting tree (-tf) or graph (-gf)\n";
    cout << "-score [string] Score input tree (-tf) or graph (-gf) and then exit:\n"
         << "    'asis' = score 'as is' i.e. without refitting,\n"
         << "    'rfit' = score after refitting (default),\n"
         << "    'mlbt' = score each base tree (with refitting) and return best,\n"
         << "    'mlno' = score each network orientation (with refitting) and return best\n";
    cout << "-mlno [string] Comma-delimited list of integers, indicating when to run\n"
         << "    maximum likelihood network orientation (MLNO) as part of heuristic search\n"
         << "    e.g. '1,2' means run MLNO only after adding the first two migration edges (default)\n"
         << "         '0' means do NOT run MLNO\n"
         << "         '' (no string) means run MLNO after adding each migration edge\n";
    cout << "-allmigs [string] Comma-delimited list of integers, indicating when to run\n"
         << "    evaluate all legal ways of adding migration edge to base tree instead of\n"
         << "    using heuristic (similar to -mlno but default is -allmigs 0)\n";
    cout << "-popaddorder [population list file] Order to add populations when building\n"
         << "    starting tree\n";
    cout << "-checkpoint Write checkpoint files\n";
    // End of additions by EKM
    cout << "\n";
}


void process_list_of_ints(string listofints, set<int> &amigs, int nmigs) {
    if (listofints.empty()) {
        for (int i = 0; i < nmigs; i++) {
            amigs.insert(i + 1);
        }
    } else {
        istringstream ss(listofints);
        string word;
        while(getline(ss, word, ',')) {
            amigs.insert(atoi(word.c_str()));
        }
    }
}

void write_check_point(string outprefix, GraphState2 &state, CountData &counts, PhyloPop_params &p) {
    string cov_se_cp = outprefix + ".covse.gz";
    string llik_cp = outprefix + ".llik"; 

    // Print graph (vertices and edges)
    state.tree->print(outprefix);

    // Print Cov SE matrix
    counts.print_cov_var(cov_se_cp);

    // Print likelihood
    ofstream likout(llik_cp.c_str());
    likout << state.llik() << "\n";
}

void write_graph(string outprefix, GraphState2 &state, PhyloPop_params &p) {
    // Print graph (vertices and edges)
    if (p.sample_size_correct == true) {
        map<Graph::edge_descriptor, double> maxw = state.get_edge_maxw();
        state.tree->print(outprefix, maxw);
    } else state.print_trimmed(outprefix);
}

int main(int argc, char *argv[]){
    printv();

    // Start of addition by EKM
    bool skipmlno = false;

    cout << "COMMAND: ";
    for (int i = 0; i < argc; i++) {
        printf(" %s", argv[i]);
    }
    cout << "\n\n";
    // End of addition by EKM

    CCmdLine cmdline;
    PhyloPop_params p;
    if (cmdline.SplitLine(argc, argv) < 1){
    	printopts();
    	exit(0);
    }
    if (cmdline.HasSwitch("-h")) {
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-o"))	outstem = cmdline.GetArgument("-o", 0).c_str();
    if (cmdline.HasSwitch("-tf"))	{
    	p.treefile = cmdline.GetArgument("-tf", 0).c_str();
    	p.readtree = true;
    }
    if (cmdline.HasSwitch("-gf"))	{
      	p.vfile = cmdline.GetArgument("-gf", 0);
      	p.efile = cmdline.GetArgument("-gf", 1);
      	p.readgraph = true;
    }
    if (cmdline.HasSwitch("-noss")) p.sample_size_correct = false;
    if (cmdline.HasSwitch("-printhzy")) p.print_hzy = true;
    if (cmdline.HasSwitch("-arcsin")) p.alfreq_scaling = 1;
    if (cmdline.HasSwitch("-nofrac")) p.nofrac = true;
    if (cmdline.HasSwitch("-scale")) p.alfreq_scaling = 3;
    if (cmdline.HasSwitch("-nothing")) p.alfreq_scaling = 4;
    if (cmdline.HasSwitch("-quick")) p.quick = true;
    if (cmdline.HasSwitch("-global")) p.global = true;
    if (cmdline.HasSwitch("-micro")) p.micro = true;
    if (cmdline.HasSwitch("-penalty")) p.neg_penalty = atof(cmdline.GetArgument("-penalty", 0).c_str());;
    if (cmdline.HasSwitch("-se")) p.calc_se = true;
    if (cmdline.HasSwitch("-emig")) p.end_mig = true;
    if (cmdline.HasSwitch("-k"))	p.window_size = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-m"))	p.nmig = atoi(cmdline.GetArgument("-m", 0).c_str());
    if (cmdline.HasSwitch("-r"))	p.nrand = atoi(cmdline.GetArgument("-r", 0).c_str());
    if (cmdline.HasSwitch("-bootstrap"))	p.bootstrap = true;
    if (cmdline.HasSwitch("-climb"))	p.climb = true;
    if (cmdline.HasSwitch("-flip"))	{
    	p.flip = true;
    	p.flipstring = cmdline.GetArgument("-flip", 0);
    }
    if (cmdline.HasSwitch("-hzy")){
    	p.hzyfile = cmdline.GetArgument("-hzy", 0);
    	p.read_hzy = true;
    }
    if (cmdline.HasSwitch("-force"))	{
    	p.forcemig = true;
    	p.mig_pops.first = cmdline.GetArgument("-force", 0);
    	p.mig_pops.second = cmdline.GetArgument("-force", 1);
    }
    if (cmdline.HasSwitch("-forcei"))	{
     	p.forcemig_index = true;
     	p.mig_index.first = atoi(cmdline.GetArgument("-forcei", 0).c_str());
     	p.mig_index.second = atoi(cmdline.GetArgument("-forcei", 1).c_str());
     }
    if (cmdline.HasSwitch("-target"))	{
    	p.dotarget = true;
    	p.target = cmdline.GetArgument("-target", 0);
    }
    if (cmdline.HasSwitch("-hold")){
    	string tmp = cmdline.GetArgument("-hold", 0);
    	vector<string> strs;
    	boost::split(strs, tmp, boost::is_any_of(","));
    	for(vector<string>::iterator it = strs.begin(); it!= strs.end(); it++) 	p.hold.insert(*it);

    }
    if (cmdline.HasSwitch("-f2"))	{
    	p.f2 = true;
    	p.alfreq_scaling = 4;
    }
    if (cmdline.HasSwitch("-covsnp"))	{
    	p.cov_snp = true;
    	p.which_cov_snp = atoi(cmdline.GetArgument("-covsnp", 0).c_str());
    }
    if (cmdline.HasSwitch("-root")) {
    	p.set_root = true;
    	p.root = cmdline.GetArgument("-root", 0);
        size_t found = p.root.find(',');
        if (found != string::npos) {
            cerr << "ERROR: Currently outgroup must be a single population when running MLNO!\n";
            exit(1);
        }
    }
    if (cmdline.HasSwitch("-f2_cor")){
    	p.cor_f2 = true;
    	p.f2_corpop = cmdline.GetArgument("-f2_cor", 0);
    	p.read_migfracs(cmdline.GetArgument("-f2_cor", 1) );
    	p.f2_mixdist = atof(cmdline.GetArgument("-f2_cor", 2).c_str());
    }
    if (cmdline.HasSwitch("-cor_mig")){
     	p.cor_mig = true;
     	//p.corpop = cmdline.GetArgument("-cor_mig", 0);
     	p.read_migfracs(cmdline.GetArgument("-cor_mig", 0) );
    }
    if (cmdline.HasSwitch("-seed")){
    	p.seed = atoi(cmdline.GetArgument("-seed", 0).c_str());
    }
    else p.seed = unsigned( time(NULL));

    if (cmdline.HasSwitch("-n_warn")) p.num_warnings = atoi(cmdline.GetArgument("-n_warn", 0).c_str());

    // Start of parameters added by EKM
    if (cmdline.HasSwitch("-freq2stat")) {
        p.freq2stat = true;
    }
    if (cmdline.HasSwitch("-givenmat")) {
        p.givenmat = true;
        p.matfile = cmdline.GetArgument("-givenmat", 0);
    }
    if (cmdline.HasSwitch("-refit")) {
        p.refit = true;
    }
    if (cmdline.HasSwitch("-score")) {
        p.doscore = true;
        string smthd = cmdline.GetArgument("-score", 0).c_str();
        if (!smthd.empty()) p.scoremethod = smthd;
    }
    if (cmdline.HasSwitch("-mlno")) {
        string listofints = cmdline.GetArgument("-mlno", 0);
        //cout << "listofints = " << listofints << "\n";
        if (listofints.compare("0") == 0) {
            skipmlno = true;
        } else {
            process_list_of_ints(listofints, p.domlno, p.nmig);
        }
    }
    if (cmdline.HasSwitch("-allmigs")) {
        string listofints = cmdline.GetArgument("-allmigs", 0);
        //cout << "listofints = " << listofints << "\n";
        if (listofints.compare("0") != 0) {
            process_list_of_ints(listofints, p.tryallmigsbt, p.nmig);
        }
    }
    if (cmdline.HasSwitch("-popaddorder")) {
        p.givenpopaddorder = true;
        p.popaddorderfile = cmdline.GetArgument("-popaddorder", 0);
    }
    if (cmdline.HasSwitch("-checkpoint")) {
        p.checkpoint = true;
    }

    if ((p.domlno.size() == 0) and (skipmlno != true)) {
        if (p.nmig >= 2) {
            string listofints = "1,2";
            process_list_of_ints(listofints, p.domlno, 2);
        }
        else if (p.nmig == 1) {
            string listofints = "1";
            process_list_of_ints(listofints, p.domlno, 1);
        }
    }
    cout << "Network search will include " << p.nmig << " admixture edge additions\n";
    std::set<int>::iterator itr;
    cout << "Exhaustive search will be performed for the following admixture edge additions: ";
    for (itr = p.tryallmigsbt.begin(); itr != p.tryallmigsbt.end(); itr++) {
        cout << *itr << " ";
    }
    cout << "\n";
    cout << "MLNO search will be performed after each of the following admixture edge additions: ";
    for (itr = p.domlno.begin(); itr != p.domlno.end(); itr++) {
        cout << *itr << " ";
    }
    cout << "\n";
    // End of parameters added by EKM

    cout.precision(8);

    // Set-up random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs2;
    r = gsl_rng_alloc(T);
    int seed = (int) time(0);
    int start_nmig;
    gsl_rng_set(r, p.seed);

    // Output files
    string covfile = outstem + ".cov.gz";
    string modelcovfile = outstem + ".modelcov.gz";
    string cov_sefile = outstem + ".covse.gz";
    string llikfile = outstem + ".llik";

    // Read input data
    CountData counts(infile, &p);
    if (p.bootstrap) counts.set_cov_bootstrap(r);
    if (p.cov_snp) counts.set_cov_singlesnp(p.which_cov_snp);
    counts.print_cov(covfile);
    counts.print_cov_var(cov_sefile);
    if (p.freq2stat) return 0;

    ofstream likout(llikfile.c_str());

    if (p.smooth_lik) p.smooth_scale = 1; //sqrt( (double) counts.nsnp / (double) p.window_size);
    GraphState2 state(&counts, &p);

    if (p.readtree) state.set_graph_from_file(p.treefile);
    else if (p.readgraph) {
        state.set_graph_from_file(p.vfile, p.efile);
        //while (state.current_llik <= -DBL_MAX) {
        //	p.smooth_scale = p.smooth_scale *2;
        //	state.current_llik = state.llik();
        //}
    }

    start_nmig = state.get_nmig();

    // Start of addition by EKM - Score input graph and exit
    if (p.doscore) {
        if (p.scoremethod.compare("asis") == 0) {
            cout << "Scoring input tree or graph\n";

            state.print_sigma_cor(modelcovfile);

            likout << setprecision(12) << "Input ln(likelihood) "
                   << state.llik() << " with "
                   << start_nmig << " migration events\n";
            likout.close();

            cout << "Final Admixture"; state.mlno_print_graph_w_params();
            cout << "Log-likelihood = " << state.current_llik << "\n";
            cout << "DONE.\n";

            return 0;
        } 

        if (p.scoremethod.compare("rfit") == 0) {
            cout << "Refitting parameters and scoring input graph\n";
            state.mlno_fit_graph();
        } else if (p.scoremethod.compare("mlbt") == 0) {
            cout << "Finding ML base tree for input graph\n";
            state.mlno_fit_graph_mlbt();
        } else if (p.scoremethod.compare("mlno") == 0) {
            cout << "Finding MLNO for input graph\n";
            state.mlno_fit_graph_mlno();
        } else {
            cerr << "ERROR: Invalid string given for -score option!\n";
            exit(1);
        }

        state.print_treeout(outstem);

        write_graph(outstem, state, p);

        state.print_sigma_cor(modelcovfile);

        likout << setprecision(12) << "Input ln(likelihood) "
               << state.llik() << " with "
               << start_nmig << " migration events\n";
        likout.close();

        cout << "Final Admixture"; state.mlno_print_graph_w_params();
        cout << "Log-likelihood = " << state.current_llik << "\n";
        cout << "DONE.\n";
        return 0;
    }

    if (p.refit) state.mlno_fit_graph();
    // End of addition by EKM

    if (start_nmig > 0) {
        if (p.set_root) state.mlno_reroot_at_outgroup(); // Added by EKM
    } else {
        // add remaining populations
        while (counts.npop > state.current_npops){
            state.add_pop();
            //state.iterate_hillclimb();
            if (p.cor_mig) {
                p.fitmig = false;
                state.iterate_local_hillclimb_wmig_all();
                p.fitmig = true;
            }
            else state.iterate_hillclimb();
            cout << "ln(likelihood): " << state.current_llik << " \n";
            cout << state.tree->get_newick_format() << "\n";
        }

        //do global rearrangements
        if (p.global) {
            cout << "Testing global rearrangements\n"; cout.flush();
            state.iterate_global_hillclimb();
            if (p.f2) state.set_branches_ls_f2();
            else state.set_branches_ls();
        }

        //place the root
        if (p.set_root) {
            state.place_root(p.root);
        }
    }

    if (p.checkpoint) write_check_point(outstem + "-checkpoint-" + to_string(start_nmig), state, counts, p);

    //print the starting likelihood (after tree building)
    likout << "Starting ln(likelihood) with " << start_nmig << " migration events: " << state.llik() << " \n";
    if (p.dotarget)	state.target_pop();
    if (p.climb) state.iterate_all_hillclimb();
    for (int i = start_nmig; i < p.nmig; i++){

    	state.current_llik = state.llik();
       	while (state.current_llik <= -DBL_MAX){
       		cout << "RESCALING\n"; cout.flush();
       		p.smooth_scale = p.smooth_scale *2;
       		state.current_llik = state.llik();
       	}
    	double current_nsum = state.negsum;
    	pair<bool, pair<int, int> > add;

        if (p.tryallmigsbt.find(i + 1) != p.tryallmigsbt.end()) {
            // Start of addition by EKM
            cout << "Performing exhaustive search to add migration edge to base tree\n";
            add = state.mlno_add_mig_to_base_tree_exhaustive();
            // End of addition by EKM
        } else {
            cout << "Performing targeted search to add migration edge to base tree\n";
            if (p.f2) add = state.add_mig_targeted_f2();
            else add = state.add_mig_targeted();
        }

        if (p.checkpoint) write_check_point(outstem + "-checkpoint-" + to_string(i+1), state, counts, p);

        if (add.first == true) {
            cout << "Migration edge #" << i + 1 << " added\n";  // Added by EKM
            state.iterate_mig_hillclimb_and_optimweight(add.second, current_nsum);

            // Start of addition by EKM
            if (p.domlno.find(i + 1) != p.domlno.end()) {
                // Check that you should be doing an add
                cout << "Performing exhaustive search for MLNO\n";
                bool is_reoriented = state.mlno_fit_graph_mlno();
                if (is_reoriented) {
                    cout << "Found orientation with higher likelihood!\n";
                    if (p.checkpoint) write_check_point(outstem + "-checkpoint-" + to_string(i+1) + "-reoriented", state, counts, p);
                }
            }
            // End of addition by EKM
    	}

        state.optimize_weights();
        if (p.f2) state.set_branches_ls_f2();
        else state.set_branches_ls();
        state.flip_mig();
        cout << "ln(likelihood):" << state.llik() << " \n";

        // Start of addition by EKM
        if (add.first == false) {
            cout << "Failed to add migration edge so won't continue trying\n";
            break;
        }
        // End of addition by EKM
    }

    if (p.end_mig) state.optimize_weights();
    if (p.forcemig) {
    	state.add_mig(p.mig_pops.first, p.mig_pops.second);
        cout << "ln(likelihood):" << state.llik() << " \n";
    }
    if (p.forcemig_index){
    	state.add_mig(p.mig_index.first, p.mig_index.second);
        cout << "ln(likelihood):" << state.llik() << " \n";
    }
    if (!p.cor_mig && !p.flip && p.nmig > 0) state.flip_mig();
    if (p.flip) state.flip_mig(p.flipstring);

    // if (p.set_root) state.mlno_reroot_at_outgroup();  // Added by EKM

    state.print_treeout(outstem);  // Change by EKM

    write_graph(outstem, state, p);

    state.print_sigma_cor(modelcovfile);

    //print the likelihood and number of migration events begin exiting
    likout << "Exiting ln(likelihood) with " << state.get_nmig()
           << " migration events: " << state.llik() << "\n";
    likout.close();

    cout << "Final Admixture"; state.mlno_print_graph_w_params();
    cout << "Log-likelihood = " << state.current_llik << "\n";
    cout << "DONE.\n";
    return 0;
}
