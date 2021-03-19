/*
 * GraphState2.h
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#ifndef GRAPHSTATE2_H_
#define GRAPHSTATE2_H_

#include "Settings.hpp"
#include "PopGraph.h"
#include "CountData.h"
#include "PhyloPop_params.h"
#include "nnls.h"

#include <iomanip>          // Added by EKM
#include <limits>           // Added by EKM
#include <list>             // Added by EKM
#include <queue>            // Added by EKM
#include <string>           // Added by EKM
//#include <unordered_set>  // Added by EKM
#include <vector>           // Added by EKM

class GraphState2{
public:
	GraphState2();
	GraphState2(CountData*, PhyloPop_params*);

	PhyloPop_params* params; //paramters for run
	PopGraph* tree;
	PopGraph* tree_bk;
	PopGraph* tree_bk2;
	PopGraph* tree_bk3;

	gsl_matrix *sigma;
	gsl_matrix *sigma_cor;

	PopGraph* scratch_tree;
	gsl_matrix *scratch_sigma_cor;

	CountData* countdata; //pointer to the data
	vector<string> allpopnames; //names of populations, will be added one at a time after the first 3
	map<string, int> popname2index; //go from the population name to the index in the above vector
	int current_npops; //current total number of populations
	double current_llik, current_llik_w;
	double scatter_det, scatter_gamma;
	gsl_matrix *scatter; //current scatter matrix
	double phi, resphi;

	//set the graph structure to a Newick string
	void set_graph(string);
	void set_graph_from_file(string);
	void set_graph_from_string(string);
	void set_countdata(CountData*);
	//set graph from files with the vertices and edges
	void set_graph(string, string);
	void set_graph(GraphState2*);

	//set the root to a given clade
	void place_root(string);

	//print
	void print_trimmed(string);

	//covariance matrix
	void compute_sigma();
	void set_sigmacor_from_sigma();
	void print_sigma();
	void print_sigma_cor(string);

	//local hill-climbing
	int local_hillclimb(int);
	int local_hillclimb_root();
	int many_local_hillclimb();
	void iterate_hillclimb();

	//global hill-climbing
	int global_hillclimb(int);
	int many_global_hillclimb();
	void iterate_global_hillclimb();

	// Start of functions added by EKM
	PopGraph* tree_binary;
	PopGraph* tree_maxllik;
	double mlno_epsilon;
	bool mlno_doit();
	int mlno_treemix_to_binary_graph();
	void mlno_binary_to_treemix_graph();
	int mlno_check_if_nonbinary();
	void mlno_get_root_einds(vector<pair<int, int> > &root_einds);
	void mlno_get_admixture_vind_combos(int &nmig, vector<set<int> > &admixture_vind_combos);
	bool mlno_reorient_huber2019(pair<int, int> &root_eind, set<int> &admixture_vinds);
	void mlno_move_root(pair<int, int> &root_eind);
	void mlno_fit_graph();
	double mlno_score_graph();
	bool mlno_reroot_at_outgroup();
        pair<int, int> mlno_get_root_eind_to_outgroup();
	bool mlno_label_topology_treechild();
	bool mlno_is_labeled_treebased();
	bool mlno_label_topology_treebased();
	bool mlno_label_and_propagate_treebased(Graph::edge_descriptor e, bool is_mig);
	bool mlno_apply_c1_treebased(Graph::edge_descriptor e, bool is_mig);
	bool mlno_apply_c2_treebased(Graph::edge_descriptor e, bool is_mig);
	void mlno_print_graph_w_params();
	void mlno_print_edge_w_params(Graph::edge_descriptor e);
	void mlno_print_graph();
	void mlno_print_edge(Graph::edge_descriptor e);
	void mlno_print_vertex(Graph::vertex_descriptor v);
	void mlno_print_edges(vector<Graph::edge_descriptor> &edges);
	pair<bool, pair<int, int> > mlno_add_mig_exhaustive();
	void mlno_get_vertices_below(Graph::vertex_descriptor v, set<Graph::vertex_descriptor> &vset);
	// End of functions added by EKM

	//add a new population
	void add_pop();
	void add_pop(string, map<string, double>);
	void add_pop(string, string);
	void process_scatter();

	//under normal model, get the max lik branch lengths
	// for a given topology by least squares
	void set_branches_ls();
	void set_branches_ls_wmig();
	void set_branches_ls_wmig_estmig();
	void set_branches_ls_f2();
	void set_branches_ls_f2_nnls();
	void set_branches_ls_f2_nnls_precompute();
	void set_branches_ls_f2_precompute_old();
	void set_branches_ls_f2_old();
	void set_branch_coefs(gsl_matrix*, gsl_vector*, map<Graph::edge_descriptor, int>*, map<Graph::edge_descriptor, double>*);
	void set_branch_coefs_f2(gsl_matrix*, gsl_vector*, map<Graph::edge_descriptor, int>*, map<Graph::edge_descriptor, double>*);
	void set_branch_coefs_f2_nnls(double *, double *, int, int, map<Graph::edge_descriptor, int>*, map<Graph::edge_descriptor, double>*);
	void set_branch_coefs_nnls(double *, double *, int, int, map<Graph::edge_descriptor, int>*, map<Graph::edge_descriptor, double>*);
	//functions used by the above least squares fitting
	map<Graph::vertex_descriptor, int> get_v2index();

	//maximize the weights on the branches. This will be iterative on each individual weight
	void optimize_weights();
	void optimize_weights_quick();
	void optimize_weights_quick(set<int>);
	void optimize_weights(set<int>);
	void optimize_weights_wish();
	void optimize_weights(Graph::edge_descriptor);
	void optimize_weight(Graph::edge_descriptor);
	void optimize_weight_quick(Graph::edge_descriptor);
	void optimize_weight_wish(Graph::edge_descriptor);
	void quick_optimize_weight(Graph::edge_descriptor);
	int golden_section_weight(Graph::edge_descriptor, double, double, double, double, int*);
	int golden_section_weight_quick(Graph::edge_descriptor, double, double, double, double, int*);
	int golden_section_weight_noexp(Graph::edge_descriptor, double, double, double, double, int*);
	int golden_section_weight_noexp_quick(Graph::edge_descriptor, double, double, double, double, int*);
	int golden_section_weight_wish(Graph::edge_descriptor, double, double, double, double);
	void optimize_fracs();
	void optimize_fracs(Graph::edge_descriptor);
	void optimize_fracs_wish();
	void optimize_frac(Graph::edge_descriptor);
	void optimize_frac_wish(Graph::edge_descriptor);
	void quick_optimize_frac(Graph::edge_descriptor);
	int golden_section_frac(Graph::edge_descriptor, double, double, double, double);
	int golden_section_frac_wish(Graph::edge_descriptor, double, double, double, double);

	//do golden section on an edge length
	int golden_section_edge(Graph::edge_descriptor, double, double, double, double);

	//likelihoods
	double llik(bool w = false);
	double llik_normal();

	double llik_wishart();
	double llik_mvn();

	//migration
	pair<string, string> get_max_resid();
	bool try_mig(Graph::vertex_descriptor, Graph::vertex_descriptor, gsl_matrix*);
	void add_mig();
	void add_mig(string, string);
	pair<bool, pair<int, int> > add_mig_targeted();
	pair<bool, pair<int, int> > add_mig_targeted_f2();
	pair< pair<bool, bool>, pair<double, pair<int, int> > > add_mig_targeted(string, string);
	double get_mig_targeted(string, string, set<pair<int, int> >*);
	pair<set<int>, set<int> > get_neighborhood(int);
	pair<set<int>, set<int> > get_neighborhood(int, int);
	void get_neighborhood(Graph::vertex_descriptor, int, pair< set<int>, set<int> >*);
	int local_hillclimb_wmig(int);
	int local_hillclimb_wmig_all(int);
	int local_hillclimb_wmig(int, set<int>);
	int iterate_local_hillclimb_wmig(pair< set<int>, set<int> >);
	int iterate_local_hillclimb_wmig_all();
	void iterate_mig_hillclimb_and_optimweight(pair<int, int>, double);
	void iterate_all_hillclimb();
	int many_local_hillclimb_wmig(pair<set<int>, set<int> >);
	int many_local_hillclimb_wmig_all();

	int try_changedir(Graph::edge_descriptor); //try changing the direction of a migration event
	int all_try_changedir();

	void flip_mig();
	void flip_mig(string);
	void trim_mig();
	//get newick string with trimmed terminal branch lengths
	string get_trimmed_newick();
	int iterate_movemig(int);
	pair<bool, int> movemig(int);
	pair<bool, int> movemig_limit(int);
	int all_try_movemig();
	int all_try_movemig_limit();
	pair<double, double> calculate_se(Graph::edge_descriptor); //get the SE of the weight on an edge (jackknife on blocks)
	pair<double, double> calculate_se_fromsamp(Graph::edge_descriptor); //get the SE of the weight on an edge (from each block)
	pair<double, double> calculate_se_bootstrap(gsl_rng *, Graph::edge_descriptor); //bootstrap
	//alterations to the tree
	pair<bool, Graph::edge_descriptor> add_mig(int, int);
	Graph::edge_descriptor add_mig_noopt(int, int);
	void rearrange(int, int);
	bool has_loop();
	void clean_negedge();

	//For rapid estimation of migration rates, need to keep the paths to the root, the least square design matrix, etc., in memory
	gsl_matrix *X_current;
	gsl_vector *y_current;
	map<Graph::edge_descriptor, int> e2index;
	map<Graph::edge_descriptor, double> e2frac;
	map<Graph::edge_descriptor, set<int> > e2tips; // the tips affected by each edge
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > popname2paths;
	map<string, map<string, int> > popnames2index;
	void update_mig(Graph::edge_descriptor, double);
	void initialize_migupdate();
	void set_branches_ls_f2_precompute();
	void print_X();
	double negsum;


	//targeting population
	void target_pop();

	//
	map<Graph::edge_descriptor, double> get_edge_maxw();

	//get the number of migration edges
	int get_nmig();
};

double lndgauss(double, double);
#endif /* GRAPHSTATE2_H_ */
