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
#include "../config.h"

string infile;
string outstem = "OrientAGraph";


void printv(){
	cout << PACKAGE_NAME << " " << PACKAGE_VERSION << "\n\n";
	cout << "OrientAGraph is built from TreeMix v1.13 Revision 231\n";
	cout << "by J.K. Pickrell and J.K. Pritchard and has several new\n";
	cout << "features, including the option to run Maximum Likelihood\n";
	cout << "Network Orientation (MLNO) as part of the admixture graph\n";
	cout << "search heuristic.\n\n";
	cout << "Contact: Erin Molloy (ekmolloy@cs.ucla.edu)\n\n";
}


void printopts(){
    cout << "Options:\n";
    cout << "-h display this help\n";
    cout << "-i [file name] input file\n";
    cout << "-o [stem] output stem (will be [stem].treeout.gz, [stem].cov.gz, [stem].modelcov.gz)\n";
    cout << "-k [int] number of SNPs per block for estimation of covariance matrix (1)\n";
    cout << "-global Do a round of global rearrangements after adding all populations\n";
    cout << "-tf [file name] Read the tree topology from a file, rather than estimating it\n";
    cout << "-m [int] number of migration edges to add (0)\n";
    cout << "-root [string] comma-delimited list of populations to set on one side of the root (for migration)\n";
    cout << "-g [vertices file name] [edges file name] read the graph from a previous TreeMix run\n";
    cout << "-se Calculate standard errors of migration weights (computationally expensive)\n";
    cout << "-micro microsatellite data\n";
    cout << "-bootstrap Perform a single bootstrap replicate\n";
    cout << "-cor_mig [file] list of known migration events to include (also use -climb)\n";
    cout << "-noss Turn off sample size correction\n";
    cout << "-seed [int] Set the seed for random number generation\n";
    cout << "-n_warn [int] Display first N warnings\n"; 
    // Start of additions by EKM
    cout << "\nOrientAGraph Options::\n";
    cout << "-mlno Run maximum likelihood network orientation subroutine as part of search heuristic\n";
    cout << "-allmigs Try all legal ways of adding migration edge to base tree instead of using the minimum residual heuristic\n";
    cout << "-popaddorder [file with list of populations] Specify the order to add populations when building the starting tree\n";
    cout << "-givenmat [se matrix file] Allows user to input matrix (e.g. [stem].cov) with the -i flag,\n"
         << "    the matrix after this option should contain the standard error (e.g. [stem].covse);\n"
         << "    if no matrix is provided after this option, then 0.0001 is used.\n";
    cout << "-refit Allows user to (re)fit model parameters on starting tree (-tf) or graph (-g)\n";
    cout << "-score [1, 2, 3, 4] Score input tree (-tf) or graph (-g) without refitting (1), with refitting (2),\n"
         << "    evaluating each base tree and returning the best (3),\n"
         << "    evaluating each network orientation and returning the best (4)\n";
    // End of additions by EKM
    cout << "\n";
}


void process_command_line(CCmdLine &cmdline, PhyloPop_params &p);
void print_treelike_output(ogzstream &treeout, GraphState2 &state, PhyloPop_params &p);


int main(int argc, char *argv[]){
    printv();

    // Added by EKM for record keeping
    cout << "COMMAND: ";
    for(int i=0; i <argc; i++) {
        printf(" %s", argv[i]);
    }
    cout << "\n\n";

    // Process parameters options from command line
    CCmdLine cmdline;
    PhyloPop_params p;
    if (cmdline.SplitLine(argc, argv) < 1){
    	printopts();
    	exit(0);
    }
    process_command_line(cmdline, p);

    // Initialize random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs2;
    r = gsl_rng_alloc(T);
    int seed = (int) time(0);
    gsl_rng_set(r, p.seed);

    // Set-up output files
    string treefile = outstem + ".treeout.gz";
    string covfile = outstem + ".cov.gz";
    string modelcovfile = outstem + ".modelcov.gz";
    string cov_sefile = outstem + ".covse.gz";
    string llikfile = outstem + ".llik";
    string mltreefile = outstem + ".mltree";  // Added by EKM
    ofstream mltout(mltreefile.c_str());  // Added by EKM
    ofstream likout(llikfile.c_str());
    ogzstream treeout(treefile.c_str());
    cout.precision(12);

    // Process the input data - either counts or matrix
    CountData counts(infile, &p);
    cout << "Read input data\n";
    if (p.bootstrap) counts.set_cov_bootstrap(r);
    if (p.cov_snp) counts.set_cov_singlesnp(p.which_cov_snp);
    counts.print_cov(covfile);
    counts.print_cov_var(cov_sefile);
    //counts.print_cov_samp("test.gz");
    if (p.smooth_lik) p.smooth_scale = 1; //sqrt((double) counts.nsnp / (double) p.window_size);

    // Initialize a graph state, giving it the input data
    // and the user-specified parameters
    GraphState2 state(&counts, &p);

    // Read input tree or graph from file
    if (p.readtree) {
	state.set_graph_from_file(p.treefile);
	cout << "Read input tree " << state.tree->get_newick_format() << "\n";
    } else if (p.read_graph){
    	state.set_graph(p.vfile, p.efile);
    	cout << "Read graph with base tree " << state.tree->get_newick_format() << "\n";
    }

    // Start of addition by EKM
    if (p.score) {
        cout << "Scoring input tree or graph\n";

        if (p.score == 2) state.mlno_fit_graph();
        else if (p.score == 3) state.mlno_score_graph();
        else if (p.score == 4) state.mlno_doit();

        state.current_llik = state.llik();
        cout << "ln(likelihood): " << state.current_llik << "\n";
        likout << setprecision(12) << "Input ln(likelihood) "
               << state.llik() << " with " 
               << state.get_nmig() << " migration events\n";
        likout.close();

        print_treelike_output(treeout, state, p);
        if (p.sample_size_correct == true) {
            map<Graph::edge_descriptor, double> maxw = state.get_edge_maxw();
            state.tree->print(outstem, maxw);
        } else {
            state.print_trimmed(outstem);
        }
        state.print_sigma_cor(modelcovfile);

        cout << "Final Admixture Graph:\n"; state.mlno_print_graph_w_params();
        cout << "Log-likelihood = " << state.current_llik << "\n";
        cout << "DONE.\n";

        return 0;
    }

    // Process starting tree
    if (p.readtree || p.read_graph) {
         cout << "Refitting parameters on input tree or graph\n";
         if (p.refit) state.mlno_fit_graph();
    } else {
        cout << "Building starting tree via population addition\n";
        while (counts.npop > state.current_npops) {
    	    state.add_pop();
    	    //state.iterate_hillclimb();
    	    if (p.cor_mig) {
                p.fitmig = false;
    		state.iterate_local_hillclimb_wmig_all();
    		p.fitmig = true;
    	    } else {
                state.iterate_hillclimb();
            }
    	    cout << "ln(likelihood): " << state.current_llik << " \n";
    	    cout << state.tree->get_newick_format() << "\n";
        }
        mltout << "[starting " << state.current_llik << "] " 
               << state.tree->get_newick_format() << "\n";
    }

    // Perform global rearrangements on starting tree
    if ((p.global) && (!p.read_graph)) {
    	cout << "Testing global rearrangements\n"; cout.flush();
    	state.iterate_global_hillclimb();

    	if (p.f2) state.set_branches_ls_f2();
    	else state.set_branches_ls();

        mltout << "[global " << state.llik() << "] " 
               << state.tree->get_newick_format() << "\n";
    }

    // Place root at outgroup
    if ((p.set_root) && !p.read_graph) {
        state.place_root(p.root);
        mltout << "[root " << state.llik() << "] " 
               << state.tree->get_newick_format() << "\n";
    }
    mltout.close();

    // Print the starting likelihood, after tree building
    likout << setprecision(12) << "Starting ln(likelihood) "
           << state.llik() << " with " 
           << state.get_nmig() << " migration events\n";

    // Try to add migration edge around population, specified by user
    if (p.dotarget) state.target_pop();
    if (p.climb) state.iterate_all_hillclimb();

    // Add migration edges
    for (int i = 0; i < p.nmig; i++){
        // Store log-likelihood
    	state.current_llik = state.llik();
       	while (state.current_llik <= -DBL_MAX){
       		cout << "RESCALING\n"; cout.flush();
       		p.smooth_scale = p.smooth_scale * 2;
       		state.current_llik = state.llik();
       	}
        
        // Try to add migration edge
    	double current_nsum = state.negsum;
    	pair<bool, pair<int, int> > add;
	if (p.tryallmigs) {
		// Start of addition by EKM
		cout << "Performing exhaustive search to add migration edge to base tree\n";
		add = state.mlno_add_mig_exhaustive();
		// End of addition by EKM
	} else {
                cout << "Performing targeted search to add migration edge to base tree\n";
		if (p.f2) add = state.add_mig_targeted_f2();
    		else add = state.add_mig_targeted();
	}

        // If successful, do some hillclimbing etc
	if (add.first == true) {
    		state.iterate_mig_hillclimb_and_optimweight(add.second, current_nsum);
		// Start of addition by EKM
		if (p.domlno) {
			cout << "Evaluating network orientations\n";
			bool is_reoriented = state.mlno_doit();
		}
		// TODO: Add global hill climbing here!
		// End of addition by EKM
	}

        // Re-optimize parameters
    	state.optimize_weights();
    	if (p.f2) state.set_branches_ls_f2();
    	else state.set_branches_ls();

        // Evaluate flipping edge included in base tree
        state.flip_mig();

        cout << "ln(likelihood):" << state.current_llik << " \n";
        likout << setprecision(12) << "Current ln(likelihood) "
               << state.current_llik << " with " 
               << state.get_nmig() << " migration events\n";
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

    if (p.set_root) state.mlno_reroot_at_outgroup();  // Added by EKM

    // Write output files
    print_treelike_output(treeout, state, p);

    if (p.sample_size_correct == true) {
        map<Graph::edge_descriptor, double> maxw = state.get_edge_maxw();
        state.tree->print(outstem, maxw);
    } else {
        state.print_trimmed(outstem);
    }

    state.print_sigma_cor(modelcovfile);

    state.current_llik = state.llik(); 
    likout << setprecision(12) << "Exiting ln(likelihood) "
           << state.current_llik << " with "
           << state.get_nmig() << " migration events\n";
    likout.close();

    cout << "Final Admixture Graph:\n"; state.mlno_print_graph_w_params();
    cout << "Log-likelihood = " << state.current_llik << "\n";
    cout << "DONE.\n";

    return 0;
}


void process_command_line(CCmdLine &cmdline, PhyloPop_params &p) {
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
    if (cmdline.HasSwitch("-g"))	{
      	p.vfile = cmdline.GetArgument("-g", 0);
      	p.efile = cmdline.GetArgument("-g", 1);
      	p.read_graph = true;
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
    if (cmdline.HasSwitch("-givenmat")) {
         p.givenmat = true;
         p.matfile = cmdline.GetArgument("-givenmat", 0);
    }
    if (cmdline.HasSwitch("-refit")) {
	// Allows user to refit parameters on the starting tree (-tf) or graph (-g)
	p.refit = true;
    }
    if (cmdline.HasSwitch("-score")) {
	// Allows user to score input tree (-tf) or graph (-g) and exit
	// 1 = score as is
	// 2 = score after refitting
	// 3 = score each base tree; return best
	// 4 = score each network orientation; return best
	p.score = atoi(cmdline.GetArgument("-score", 0).c_str());
    }
    if (cmdline.HasSwitch("-popaddorder")) {
	// Allows user to specify the order of populations when building the starting tree
        p.givenpopaddorder = true;
        p.popaddorderfile = cmdline.GetArgument("-popaddorder", 0);
    }
    if (cmdline.HasSwitch("-allmigs")) {
        // Allows users to test all possibilities for adding an edge
        p.tryallmigs = true;
    }
    if (cmdline.HasSwitch("-mlno")) {
        // Allows user to include MLNO after each edge addition
        p.domlno = true;
    }
    if (cmdline.HasSwitch("-driftunits")) {
        // Allows user to convert branch lengths to genetic drift units
        // TODO: Change units using technique used by MixMapper
        p.driftunits = true;
    }

    if (p.refit || p.score) {
        if (!(p.readtree || p.read_graph)) {
            cout << "Must use -tf or -g option with -refit or -score!\n";
            exit(1);
        }
    }
    // End of parameters added by EKM
}


void print_treelike_output(ogzstream &treeout, GraphState2 &state, PhyloPop_params &p) {
    treeout << state.tree->get_newick_format() << "\n";

    if (p.sample_size_correct == false) treeout << state.get_trimmed_newick() << "\n";

    pair<Graph::edge_iterator, Graph::edge_iterator> eds = edges(state.tree->g);
    Graph::edge_iterator it = eds.first;
    p.smooth_lik = false;
    while (it != eds.second){
        if (state.tree->g[*it].is_mig) {
     	    double w = state.tree->g[*it].weight;

            treeout << state.tree->g[*it].weight<< " ";
            if (p.calc_se){
                Graph::vertex_descriptor p1 = source( *it, state.tree->g);
                p1 = state.tree->get_child_node_mig(p1);
                Graph::vertex_descriptor p2 = target(*it, state.tree->g);
                cout << state.tree->get_newick_format(p1) << " ";
     	        cout << state.tree->get_newick_format(p2) << "\n"; cout.flush();
                p.neg_penalty = 0;
                pair<double, double> se = state.calculate_se(*it);
                treeout << se.first << " "<< se.second << " ";
                double test = se.first/ se.second;
                double pval = 1-gsl_cdf_gaussian_P(test, 1);
                if (pval < DBL_MIN){
     	            pval = DBL_MIN;
                    treeout << "<"<< pval << " ";
                } else {
                    treeout << pval << " ";
                }
            } else {
                treeout << "NA NA NA ";
            }

            state.tree->g[*it].weight = w;
            if (p.f2) state.set_branches_ls_f2();
            else state.set_branches_ls();

            Graph::vertex_descriptor p1 = source( *it, state.tree->g);
            p1 = state.tree->get_child_node_mig(p1);
            Graph::vertex_descriptor p2 = target(*it, state.tree->g);
            treeout << state.tree->get_newick_format(p1) << " ";
            treeout << state.tree->get_newick_format(p2) << "\n";
        }
        it++;
    }
}

