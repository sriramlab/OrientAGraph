/*
 * GraphState2.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#include "GraphState2.h"

GraphState2::GraphState2(){
}

GraphState2::GraphState2(CountData* counts, PhyloPop_params* pa) {
	params = pa;
	countdata = counts;
	allpopnames = counts->list_pops();
	cout << "SEED: "<< pa->seed << "\n";
	srand ( pa->seed );

	if (pa->givenpopaddorder) {
		// Start of addition by EKM
		allpopnames.clear();
		
		ifstream fptr(pa->popaddorderfile.c_str());
		if (!fptr.is_open()) {
			cerr << "Unable to open " << pa->popaddorderfile << "\n";
			exit(1);
		}

		cout << "Populations will be added in the following order:";
		string line;
		int back;
		int i = 0;
		while (getline(fptr, line)) {
			back = line.size() - 1;
			if (line[back] == '\r') line.erase(back);  // Windows

			allpopnames.push_back(line);
			cout << " " << allpopnames[i];
			i++;
		}
		cout << "\n";
		fptr.close();
		// End of addition by EKM
	} else {
		random_shuffle(allpopnames.begin(), allpopnames.end() );
	}

	vector<string> tmppopnames;
	vector<string> tmphold;
	for (int i = 0; i < allpopnames.size(); i++){
		string tmp = allpopnames[i];
		if (params->hold.find(tmp) == params->hold.end()) tmppopnames.push_back(tmp);
		else tmphold.push_back(tmp);
	}
	for (int i = 0; i < tmphold.size(); i++) tmppopnames.push_back(tmphold[i]);
	for (int i = 0; i <  allpopnames.size(); i++)	{
		allpopnames[i] = tmppopnames[i];
		popname2index.insert(make_pair( allpopnames[i], i));
	}
	vector<string> startpops;
	startpops.push_back(allpopnames[0]); startpops.push_back(allpopnames[1]); startpops.push_back(allpopnames[2]);

	tree = new PopGraph(startpops);
	tree_bk = new PopGraph(startpops);
	tree_bk2 =new PopGraph(startpops);
	tree_bk3 = new PopGraph(startpops);

	// Start of addition by EKM
	tree_binary = new PopGraph(startpops);
	tree_maxllik = new PopGraph(startpops);
	// End of addition by EKM

	current_npops = 3;
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);

	scatter = gsl_matrix_alloc(current_npops, current_npops);
	scatter_det = counts->scatter_det;
	scatter_gamma = counts->scatter_gamma;
	X_current = gsl_matrix_alloc(current_npops, current_npops);
	y_current = gsl_vector_alloc(current_npops);
	gsl_matrix_set_zero(sigma);
	negsum = 0;

	if (pa->f2) set_branches_ls_f2();
	else set_branches_ls();
	current_llik = llik();

	//vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(3);

	local_hillclimb(3);
	if (!params->readtree && !params->readgraph)
		cout << "Starting from\n" << tree->get_newick_format() << "\n";
	phi = (1+sqrt(5))/2;
	resphi = 2-phi;
}

void GraphState2::set_countdata(CountData* counts){
	countdata = counts;
	allpopnames = counts->list_pops();
}
void GraphState2::print_sigma(){
	for (int i = 0; i < current_npops; i++){
		cout << allpopnames[i]<< " ";
	}
	cout << "\n";
	for(int i = 0; i < current_npops; i++){
		cout << allpopnames[i] << " ";
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
	for(int i = 0; i < current_npops; i++){
		cout << allpopnames[i] << " ";
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma_cor, i, j) << " ";
		}
		cout << "\n";
	}
}

void GraphState2::set_graph(string newick){
	cerr << "WARNING: Check function is correct before using!\n";  // Added by EKM

	tree->set_graph(newick);
	current_npops = allpopnames.size();
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

    set_branches_ls_wmig();

    current_llik = llik();
}

void GraphState2::set_graph(GraphState2* g){
	cerr << "WARNING: Check function is correct before using!\n";  // Added by EKM

	tree->copy( g->tree);
	current_npops = allpopnames.size();
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

    set_branches_ls_wmig();
    cout << llik()<< "\n";
    current_llik = llik();

}


void GraphState2::set_graph_from_file(string vfile, string efile) {
	// Define variables
	struct stat stFileInfo;

	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	map<string, Graph::vertex_descriptor>::iterator mi;
	vector<string>::iterator si;

	map<string, Graph::vertex_descriptor> tips;
	map<int, Graph::vertex_descriptor> i2v;
	queue<Graph::vertex_descriptor> vqueue;
	set<Graph::vertex_descriptor> vfound;
	vector<string> line;

	Graph::edge_descriptor e;
	Graph::vertex_descriptor v, sv, tv;

	string st;
	float len, w, mig_frac, totalw, maxf;
	int nv, index, index1, index2, intStat, back;

	// Set-up
	tree->g.clear();
	current_npops = 0;
	tree->popnames.clear();

	// Process vertex file
	intStat = stat(vfile.c_str(), &stFileInfo);
	if (intStat != 0) {
		std::cerr << "ERROR: cannot open file " << vfile << "\n";
		exit(1);
	}
	igzstream vin(vfile.c_str());

	nv = 0;
	while(getline(vin, st)){
		back = st.size() - 1;
		if (st[back] == '\r') st.erase(back);  // Windows

		string buf;
		stringstream ss(st);
		line.clear();
		while (ss >> buf){
			line.push_back(buf);
		}

		index = atoi(line[0].c_str());
		if (index >= tree->indexcounter) tree->indexcounter = index+1;

		string name = line[1];
		string root = line[2];
		string mig = line[3];
		string tip = line[4];

		v = add_vertex(tree->g);

		tree->g[v].index = index;
		tree->g[v].name = name;

		if (root == "ROOT") {
			tree->g[v].is_root = true;
			tree->root = v;
		} else {
			tree->g[v].is_root = false;
		}
		if (mig == "MIG") {
			tree->g[v].is_mig = true;
			tree->g[v].mig_frac = 0.5;
		} else {
			tree->g[v].is_mig = false;
			tree->g[v].mig_frac = 0.0;  // EKM: Changed 1 to 0 based on PopGraph
		}
		if (tip == "TIP") {
			tree->g[v].is_tip = true;
			current_npops++;
			tree->popnames.push_back(name);
		} else {
			tree->g[v].is_tip = false;
		}

		tree->g[v].height = 0.0;           // Added by EKM based on PopGraph
		tree->g[v].desired_indegree = -1;  // Added by EKM
		tree->g[v].current_indegree = -1;  // Added by EKM
		tree->g[v].rev = false;

		++nv;
	}
	vin.close();
	i2v = tree->index2vertex();

	// Process edge file
	intStat = stat(efile.c_str(), &stFileInfo);
	if (intStat != 0) {
		std::cerr << "ERROR: cannot open file " << efile << "\n";
		exit(1);
	}
	igzstream ein(efile.c_str());

	while(getline(ein, st)) {
		back = st.size() - 1;
		if (st[back] == '\r') st.erase(back);  // Windows

		string buf;
		stringstream ss(st);
		line.clear();
		while (ss >> buf){
			line.push_back(buf);
		}

		index1 = atoi(line[0].c_str());
		index2 = atoi(line[1].c_str());
		len = atof(line[2].c_str());
		w = atof(line[3].c_str());
		string mig = line[4];
		mig_frac = atof(line[5].c_str());  // EKM: Doesn't get used?
		maxf = atof(line[6].c_str());      // EKM: Doesn't get used?

		if (len < 0) {
			// Note from EKM: TreeMix can produce negative branch lengths!
			cerr << "WARNING: Input graph has negative branch lengths or weights!\n";
		}
		if (w < 0) {
			cerr << "ERROR: Input graph has negative weights!\n";
			exit(1);
		}

		sv = i2v[index1];
		tv = i2v[index2];
		e = add_edge(sv, tv, tree->g).first;

		tree->g[e].weight = w;
		if (mig == "MIG") {
			tree->g[e].is_mig = true;
			len = 0.0;
		} else {
			tree->g[e].is_mig = false;
			tree->g[e].len = len;
		}
		tree->g[e].is_labeled = false;   // Added by EKM
		tree->g[e].is_oriented = false;  // Added by EKM

		//if (tree->g[e].is_mig && !tree->g[sv].is_mig) {
		//	cerr << "ERROR: Input graph is broken!\n";
		//	exit(1);
		//}
	}
	ein.close();

	// Start of addition by EKM
	// Check leaves of graph match input data
	tips = tree->get_tips(tree->root);
	for (mi = tips.begin(); mi != tips.end(); mi++) {
		si = std::find(allpopnames.begin(), allpopnames.end(), mi->first);
		if (si == allpopnames.end()) {
			cerr << "ERROR: Population " << mi->first << " found in input graph but not input data!\n";
			exit(1);
		}
	}
	for (si = allpopnames.begin(); si != allpopnames.end(); si++) {
		if (tips.find(*si) == tips.end()) {
			cerr << "ERROR: Population " << *si << " found in input data but not graph!\n";
			exit(1);
		}
	}

	// Check that input weights on in-edges sum to one
//	viter = vertices(tree->g);
//	for (vi = viter.first; vi != viter.second; ++vi) {
//		if (!tree->g[*vi].is_root) { 
//			totalw = 0;
//			in_eiter = in_edges(*vi, tree->g);
//			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei)
//				totalw += tree->g[*in_ei].weight;
//
	//		if ((totalw < 0.999) || (totalw > 1.001)) {
	//			cerr << "ERROR: Weights in input graph are broken!\n";
	//			exit(1);
	//		}
	//	}
	//}
	// TODO: Check all vertices are reachable from the root, etc.
	// TODO: Check it is 'tree-based enough' TreeMix and then warn user if not...
	// End of addition by EKM

	// Start of comment by EKM -- didn't look like this was doing anything...
	//igzstream ein2(efile.c_str());
	//while(getline(ein2, st)) {
     	//	string buf;
 	//	stringstream ss(st);
 	//	line.clear();
 	//	while (ss>> buf){
 	//		line.push_back(buf);
 	//	}
 	//	sv = i2v[atoi(line[0].c_str())];
 	//	tv = i2v[atoi(line[1].c_str())];
 	//	mig_frac = atof(line[5].c_str());
 	//	e = edge(sv, tv, tree->g).first;	
 	//	//if (tree->g[e].is_mig) tree->set_mig_frac(e, mig_frac);
     	//}
     	//ein2.close();
	// End of comment by EKM

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

	mlno_print_graph_w_params();

	mlno_compute_sigma_cor();  // Change by EKM

	current_llik = llik();
}

void GraphState2::set_graph_from_file(string infile) {
	// Define variables
	struct stat stFileInfo;

	std::map<string, Graph::vertex_descriptor>::iterator mi;
	vector<string>::iterator si;

	map<string, Graph::vertex_descriptor> tips;
	vector<string> line;

	string st, buf;
	int intStat;

	intStat = stat(infile.c_str(), &stFileInfo);
	if (intStat != 0) {
		std::cerr<< "ERROR: cannot open file " << infile << "\n";
		exit(1);
	}
	igzstream in(infile.c_str());

	getline(in, st);
	int back = st.size() - 1;
	if (st[back] == '\r') st.erase(back);  // Windows

	current_npops = allpopnames.size();

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);

	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

	tree->set_graph(st);
	tips = tree->get_tips(tree->root);

	// Start of change by EKM
	for (mi = tips.begin(); mi != tips.end(); mi++) {
		si = std::find(allpopnames.begin(), allpopnames.end(), mi->first);
		if (si == allpopnames.end()) {
			cerr << "ERROR: Population " << mi->first << " found in input tree but not input data!\n";
			exit(1);
		}
	}
	for (si = allpopnames.begin(); si != allpopnames.end(); si++){
		if (tips.find(*si) == tips.end()) {
			cerr << "ERROR: Population " << *si << " found in input data but not tree!\n";
			exit(1);
		}
	}

	// TODO: 
	// Enable user to provide a tree on a subset of populations, so that the remaining
	// populations are added. This requires updating the ordering of the list allpopnames
	// as well as the matrix files, and so on.

	current_npops = tips.size();
	mlno_compute_sigma_cor();
	// End of change by EKM

	current_llik = llik();
}


void GraphState2::set_graph_from_string(string newick){
	cerr << "WARNING: Check function is correct before using!\n";  // Added by EKM

	set_graph(newick);
	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);
	allpopnames.clear();
	for (map<string, Graph::vertex_descriptor>::iterator it = tips.begin(); it != tips.end(); it++){
		allpopnames.push_back(it->first);
		if (countdata->pop2id.find(it->first) == countdata->pop2id.end()){
			cerr << "ERROR: cannot find population " << it->first << "\n";
			exit(1);
		}
	}

    current_npops = allpopnames.size();
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

	// Start of change by EKM
	//set_branches_ls();
	mlno_compute_sigma_cor();
	// End of change by EKM

	current_llik = llik();
}

map<Graph::vertex_descriptor, int> GraphState2::get_v2index(){

	map<Graph::vertex_descriptor, int> vertex2index;

	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;
	//get all the migration nodes
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}
	return vertex2index;

}
void GraphState2::set_branches_ls_wmig(){

	/* one parameter for each non-migration node (minus 1 for the root and 1 for the unidentifiable branch length next to the root), one for each migration nodde
	 *
	 *
    */

	map<Graph::vertex_descriptor, int> vertex2index;
	map<Graph::vertex_descriptor, float> vertex2frac;
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			vertex2frac.insert(make_pair(i_nodes[i], 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;
		//get all the migration nodes
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}

	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3 + mig_nodes.size(); // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);


	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations

	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);

	//get all paths to the root for all tips
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::vertex_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::vertex_descriptor> > > tmpset = tree->get_paths_to_root(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}


	// Set up the matrices from the tree
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y, index, empirical_cov);

			set<pair<double, set<Graph::vertex_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if (tree->g[*it3].is_root) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = vertex2index[*it3];
							double add = it1->first*it2->first;
							double inv2add = it1->first*it2->first*inv_total2;
							if ( root_adj.find(*it3) != root_adj.end()) {
									add = add/2;
									inv2add = inv2add/2;
							}
							gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
							for (int z = 0; z < n ; z++){
								gsl_matrix_set(X, z, addindex, gsl_matrix_get(X, z, addindex) + inv2add);
							}
						}
					}
				}
			}
			index++;
		}
	}

	// Now add the bias terms
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_2 = name2paths[p2];

			for (int k = 0; k < current_npops; k++){
				string p3 = allpopnames[k];
				set<pair<double, set<Graph::vertex_descriptor> > > paths_3 = name2paths[p3];
				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
					for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_root) continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = vertex2index[*it3];
								double weight = it1->first*it2->first;
								double invadd = weight*inv_total;
								if ( root_adj.find(*it3) != root_adj.end()) invadd = invadd/2;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}

				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_root) continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = vertex2index[*it3];
								double weight = it1->first*it2->first;
								double invadd = weight*inv_total;
								if ( root_adj.find(*it3) != root_adj.end()) invadd = invadd/2;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}
			}
			index++;
		}

	}

/*
	for(int i = 0; i < n ; i ++){
		cout << gsl_vector_get(y, i) << " ";
		for (int j = 0; j < p; j++){
			cout << gsl_matrix_get(X, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
*/

	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	//and put in the solutions in the graph
	for( int i = 0; i < i_nodes2.size(); i++){
		Graph::vertex_descriptor v = i_nodes2[i];
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, tree->g).first;
		double l = gsl_vector_get(c, i);
		//if (l < 0) l = 1E-8;
		tree->g[*in_i].len = l;
	}
	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++){
		Graph::vertex_descriptor v = *it;
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, tree->g).first;
		double l = gsl_vector_get(c, joint_index);
		//if (l < 0) l = 1E-8;
		tree->g[*in_i].len = l/2;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X, index, k) * gsl_vector_get(c, k);
			}
			//cout << i << " "<< j << " "<< pred << "\n";
			//cout << index << " "<< i << " "<< j << " "<< pred << "\n";
			gsl_matrix_set(sigma_cor, i, j, pred);
			//cout << "not here\n";
			index++;
		}

	}


	//free memory
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}

void GraphState2::set_branch_coefs(gsl_matrix* X, gsl_vector* y, map<Graph::edge_descriptor, int>* edge2index, map<Graph::edge_descriptor, double>* edge2frac){

	//
	// y = Xc
	//
	// y is the observed matrix, X contains the contribution of each branch length to each entry in y
	//
	//

	gsl_matrix_set_zero(X);
	// get all the paths to the root from each tip

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}

	double inv_total = 1.0/ (double) countdata->npop;
	double inv_total2 = 1.0/ ( (double) countdata->npop * (double) countdata->npop);

	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y, index, empirical_cov);

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;
							//if ( i != j ) add = add*2;
							double inv2add = add*inv_total2;
							gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
							for (int z = 0; z < current_npops*current_npops ; z++){
								gsl_matrix_set(X, z, addindex, gsl_matrix_get(X, z, addindex) + inv2add);

							}
						}
					}
				}
			}
			index++;
		}
	}


	for(map<Graph::edge_descriptor, int>::iterator it = edge2index->begin(); it != edge2index->end(); it++){
		cout << tree->g[source(it->first, tree->g)].index << " "<<  tree->g[target(it->first, tree->g)].index << " "<< it->second << " "<< edge2frac->find(it->first)->second << "\n";

	}
	cout <<  "\n";

	// Now add the bias terms
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];

			for (int k = 0; k < current_npops; k++){
				string p3 = allpopnames[k];
				set<pair<double, set<Graph::edge_descriptor> > > paths_3 = name2paths[p3];
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_mig)continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first*frac;
								double invadd = weight*inv_total;
								//if (i != j) invadd = invadd*2;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}

				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (it2->second.find(*it3) != it2->second.end()){
								if (tree->g[*it3].is_mig)continue;
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first *frac;
								double invadd = weight*inv_total;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}
			}
			index++;
		}

	}



}


void GraphState2::set_branch_coefs_f2(gsl_matrix* X, gsl_vector* y, map<Graph::edge_descriptor, int>* edge2index, map<Graph::edge_descriptor, double>* edge2frac){

	//
	// y = Xc
	//
	// y is the observed matrix, X contains the contribution of each branch length to each entry in y
	//
	//

	gsl_matrix_set_zero(X);
	gsl_vector_set_zero(y);
	// get all the paths to the root from each tip

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}

	map<string, int> pair2index;
	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			if (i == j) {
				index++;
				continue;
			}
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			string both = p1+p2;
			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y, index, empirical_cov);

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
					if ( tree->g[*it3].is_mig) continue;
					double frac = edge2frac->find(*it3)->second;
					int addindex = edge2index->find(*it3)->second;
					double add = it1->first * it1->first *frac;
					//if (p1 == "pop3" && p2 == "pop10") cout << index << " "<< addindex << " "<< add << "\n";
					//if (index  == 12 && addindex   == 0) cout << "addind for pop "<< p2 << " "<< add << "\n";
					gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
				}
			}

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						double frac = edge2frac->find(*it3)->second;
						int addindex = edge2index->find(*it3)->second;
						double add = it1->first * it1->first *frac;
						//if (p1 == "pop3" && p2 == "pop10") cout << index << " "<< addindex << " "<< add << "\n";
						//if (index  == 12 && addindex   == 0) cout << "addind for pop "<< p2 << " "<< add << "\n";
						gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
					}
				}

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;
							add = -2*add;
							//if (index  == 12 && addindex   == 0) cout << "addind for overlap "<< add << "\n";
							gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
						}
					}
				}
			}
			pair2index.insert(make_pair(p1+p2, index));
			pair2index.insert(make_pair(p2+p1, index));
			index++;
		}
	}
}


void GraphState2::set_branch_coefs_f2_nnls(double * A, double * b, int n, int p, map<Graph::edge_descriptor, int>* edge2index, map<Graph::edge_descriptor, double>* edge2frac){

	//
	// y = Xc
	//
	// y is the observed matrix, X contains the contribution of each branch length to each entry in y
	//
	//

	for (int i = 0; i < n; i++){
		for (int j = 0; j < p ; j++) A[j*n+i] = 0.0;

	}

	for (int i = 0; i < p; i ++) b[i] = 0.0;
	// get all the paths to the root from each tip

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}

	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			if (i == j) {
				index++;
				continue;
			}
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			string both = p1+p2;
			double empirical_cov = countdata->get_cov(p1, p2);
			b[index] = empirical_cov;

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];

			//Variances of population 1
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
					if ( tree->g[*it3].is_mig) continue;
					double frac = edge2frac->find(*it3)->second;
					int addindex = edge2index->find(*it3)->second;
					double add = it1->first * it1->first *frac;
					set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
					it4++;
					while (it4 != paths_1.end()){
						if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
						it4++;
					}
					A[addindex* n + index]+=add;
				}
			}

			//Variances of population 2
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						double frac = edge2frac->find(*it3)->second;
						int addindex = edge2index->find(*it3)->second;
						double add = it1->first * it1->first *frac;
						set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
						it4++;
						while (it4 != paths_2.end()){
							if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
							it4++;
						}
						A[addindex* n + index]+=add;
					}
				}

			//Covariances
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;
							add = -2*add;
							A[addindex* n + index]+=add;
						}
					}
				}
			}
			index++;
		}
	}
}


void GraphState2::set_branch_coefs_nnls(double * A, double * b, int n, int p, map<Graph::edge_descriptor, int>* edge2index, map<Graph::edge_descriptor, double>* edge2frac){

	//
	// y = Xc
	//
	// y is the observed matrix, X contains the contribution of each branch length to each entry in y
	//
	//
	gsl_matrix * tmpA = gsl_matrix_alloc(n, p);
	gsl_matrix_set_zero(tmpA);
	gsl_matrix *tmpN = gsl_matrix_alloc(current_npops, p);
	gsl_matrix_set_zero(tmpN);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < p ; j++) A[j*n+i] = 0.0;

	}

	for (int i = 0; i < p; i ++) b[i] = 0.0;
	// get all the paths to the root from each tip

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}
	// trying to speed up by calculating all covariances once
	//
	//
	//
	double inv_total = 1.0/ (double) countdata->npop;
	double inv_total2 = 1.0/ ( (double) countdata->npop * (double) countdata->npop);
	gsl_vector * pbias = gsl_vector_alloc(p);
	gsl_vector_set_zero(pbias);
	int indextmp =0;
	map<pair<int, int>, int> pop_pair2index; //hold the index for i and j

	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			pop_pair2index.insert(make_pair( make_pair(i, j), indextmp));
			pop_pair2index.insert(make_pair( make_pair(j, i), indextmp));
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;
							gsl_matrix_set(tmpA, indextmp, addindex, gsl_matrix_get(tmpA, indextmp, addindex)+add);

							//negative biases
							double invadd = add*inv_total;
							//cout << invadd << " "<< frac << " "<< add << "\n";
							gsl_matrix_set(tmpN, i, addindex, gsl_matrix_get(tmpN, i, addindex)+invadd);
							if (i != j) gsl_matrix_set(tmpN, j, addindex, gsl_matrix_get(tmpN, j, addindex)+invadd);


							//positive bias
							if (i !=j) add = add*2;
							double inv2add = add*inv_total2;
							gsl_vector_set(pbias, addindex, gsl_vector_get(pbias, addindex)+inv2add);
						}
					}
				}
			}
			indextmp++;
		}
	}

	//
	// /TEST
	//

	/*
	cout << "tmpN\n";
	for(int i = 0; i < current_npops; i++){
		for (int j = 0; j < p; j++){
			cout <<gsl_matrix_get(tmpN, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n\n";
*/

	// TEST
	int index =0;



	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			b[index] = empirical_cov;

			for (int k = 0 ; k < p; k++){
				//1. Covariance term
				A[k* n + index]+=gsl_matrix_get(tmpA, index, k);

				//2. Positive bias
				A[k* n + index]+=gsl_vector_get(pbias, k);

				//3. Negative biases
				A[k* n + index] -= gsl_matrix_get(tmpN, i, k);
				A[k* n + index] -= gsl_matrix_get(tmpN, j, k);
			}

			// 2. Now the negative bias terms
			//for (int k = 0; k < current_npops; k++){
			//	int index1 = pop_pair2index[ make_pair(i, k)];
			//	int index2 = pop_pair2index[ make_pair(j, k)];
			//	for (int l = 0; l < p; l++){
			//		A[l* n + index]-= inv_total*gsl_matrix_get(tmpA, index1, l);
			//		A[l* n + index]-= inv_total*gsl_matrix_get(tmpA, index2, l);
			//	}
			//}
			index++;
		}
	}


	// /TEST

	/*
	// Set up the matrices from the tree
	int index = 0;
	int zmax = current_npops*(current_npops-1)/2 +current_npops;
	for( int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			b[index] = empirical_cov;

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;

							A[addindex* n + index]+=add;
							if (i !=j) add = add*2;
							double inv2add = add*inv_total2;
							int index2 = 0;
							//for (int y = 0; y < current_npops; y++){
							//	for (int z = y; z< current_npops; z++)
							//		if ( z == y)
							//}
							for (int z = 0; z < zmax ; z++){
								A[addindex* n + z]+=inv2add;
							}
						}
					}
				}
			}
			index++;
		}
	}


	// Now add the bias terms
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];

			for (int k = 0; k < current_npops; k++){
				string p3 = allpopnames[k];
				set<pair<double, set<Graph::edge_descriptor> > > paths_3 = name2paths[p3];
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_mig)continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first*frac;
								//if (i != k && i !=j) weight = weight*2;
								double invadd = weight*inv_total;
								A[addindex* n + index]-=invadd;
							}
						}
					}
				}

				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (it2->second.find(*it3) != it2->second.end()){
								if (tree->g[*it3].is_mig)continue;
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first *frac;
								//if (j != k && i !=j) weight = weight*2;
								double invadd = weight*inv_total;
								A[addindex* n + index]-=invadd;
							}
						}
					}
				}
			}
			index++;
		}

	}
	 */

	/*
	int zmax = current_npops*(current_npops-1)/2 +current_npops;
	for (int i = 0; i < zmax; i++){
		cout << b[i] << " ";
		for(int j = 0; j < p ; j++){
			cout << A[j* n + i] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
	tree->print();
*/
	gsl_matrix_free(tmpA);
	gsl_matrix_free(tmpN);
	gsl_vector_free(pbias);
}


map<Graph::edge_descriptor, double> GraphState2::get_edge_maxw(){
	map<Graph::edge_descriptor, double> toreturn;

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}
	// trying to speed up by calculating all covariances once
	//
	//
	//


	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							double add = it1->first * it2->first;
							if (toreturn.find(*it3) == toreturn.end()){
								toreturn.insert(make_pair(*it3, add));
							}
							else if (toreturn[*it3]< add) toreturn[*it3] = add;
						}
					}
				}
			}
		}
	}


	return toreturn;
}


void GraphState2::optimize_weights(){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik;
	//cout << current_llik << " "<< start_llik << "\n";
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++) optimize_weight(*it);
		if (current_llik < start_llik+params->epsilon) done = true;
		else start_llik = current_llik;
		nit++;
	}

}


void GraphState2::optimize_weights_quick(){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */
	//cout << "in quick\n"; cout.flush();
	if (!params->f2) optimize_weights();
	else{
		initialize_migupdate();

		double start_llik = current_llik;
		//cout << current_llik << " "<< start_llik << "\n";
		bool done = false;
		int nit = 0;
		while(!done){
			for (map<Graph::edge_descriptor, set<int> >::iterator it = e2tips.begin(); it != e2tips.end(); it++) 	{
			//cout << tree->g[source(it->first, tree->g)].index << "\n";
			//cout << "quick\n"; cout.flush();
				optimize_weight_quick(it->first);
			}
			if (current_llik < start_llik+params->epsilon) done = true;
			else start_llik = current_llik;
		//cout << nit << " "<< current_llik << " "<< start_llik << "\n"; cout.flush();
			nit++;
		}
	}
	//cout << "done\n";
}



void GraphState2::optimize_weights_quick(set<int> n){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */
	if (!params->f2) optimize_weights(n);
	initialize_migupdate();
	double start_llik = current_llik;
	bool done = false;
	int nit = 0;
	while(!done && nit < params->maxit2){
		for (map<Graph::edge_descriptor, set<int> >::iterator it = e2tips.begin(); it != e2tips.end(); it++) 	{
			int test = tree->g[source(it->first, tree->g)].index;
			if (n.find(test) == n.end()) continue;
			optimize_weight_quick(it->first);
		}
		//cout << nit << " "<< current_llik << " "<< start_llik << "\n";
		if ( (current_llik < start_llik+0.1) || n.size() < 2) done = true;
		else start_llik = current_llik;
		nit++;
	}
	//cout << "done\n";
}



void GraphState2::optimize_weights(set<int> n){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */
	//initialize_migupdate();
	double start_llik = current_llik;
	bool done = false;
	int nit = 0;
	while(!done && nit < params->maxit2){
		vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++) 	{
			int test = tree->g[source(*it, tree->g)].index;
			if (n.find(test) == n.end()) continue;
			optimize_weight_quick(*it);
		}
		if ( (current_llik < start_llik+0.1) || n.size() < 2) done = true;
		else start_llik = current_llik;
		nit++;
	}
}

void GraphState2::optimize_weights(Graph::edge_descriptor e){
	/*
	 *  Go through each migration edge except the one in the argument, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik;
	//cout << current_llik << " "<< start_llik << "\n";
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++){
			if (*it == e ) continue;
			double min, max, guess;
			guess = tree->g[*it].weight;
			guess = log(guess/ (1-guess));
			min = params->minweight;
			max = params->maxweight;
			int nit = 0;
			golden_section_weight( *it, min, guess, max, params->tau, &nit);
			cout << tree->g[*it].weight << "\n";
			}
		//cout << current_llik << " " <<start_llik << "\n";
		if (current_llik < start_llik+0.1) done = true;
		else start_llik = current_llik;
		nit++;
	}

}



void GraphState2::optimize_fracs(){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik;
	//cout << start_llik << " "<< current_llik << " llik\n";
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++){
			double min, max, guess;
			guess = tree->g[source(*it, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac( *it, min, guess, max, params->tau);
			cout << tree->g[source(*it, tree->g)].mig_frac << " frac\n";
			}
		//cout << start_llik << " "<< current_llik << " llik\n";
		if (current_llik < start_llik+0.1) done = true;
		else start_llik = current_llik;
		nit++;
	}

}



void GraphState2::optimize_fracs_wish(){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik_w;
	//cout << start_llik << " "<< current_llik_w << "  w llik\n";
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++){
			double min, max, guess;
			guess = tree->g[source(*it, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac_wish( *it, min, guess, max, params->tau);
			cout << tree->g[source(*it, tree->g)].mig_frac << " w frac\n";
			}
		//cout << start_llik << " "<< current_llik << " llik\n";
		if (current_llik_w < start_llik+0.1) done = true;
		else start_llik = current_llik_w;
		nit++;
	}

}



void GraphState2::optimize_fracs(Graph::edge_descriptor e){
	/*
	 *  Go through each migration edge (except one in argument), optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik;
	cout << start_llik << " "<< current_llik << " llik\n";
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++){
			if (*it == e) continue;
			double min, max, guess;
			guess = tree->g[source(*it, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac( *it, min, guess, max, params->tau);
			cout << tree->g[source(*it, tree->g)].mig_frac << " frac\n";
			}
		cout << start_llik << " "<< current_llik << " llik\n";
		if (current_llik < start_llik+0.1) done = true;
		else start_llik = current_llik;
		nit++;
	}

}

void GraphState2::quick_optimize_weight(Graph::edge_descriptor e){
	tree->g[e].weight = 0.1;
	set_branches_ls_wmig();
	double best = 0.1;
	double start_llik = llik();
	for (double w = 0.3; w < 1; w+=0.2){
		tree->g[e].weight = w;
		set_branches_ls_wmig();
		double test_llik = llik();
		if (test_llik > start_llik){
			start_llik  = test_llik;
			best = w;
		}
	}
	tree->g[e].weight = best;
	set_branches_ls_wmig();
	current_llik = llik();
}



void GraphState2::optimize_weight(Graph::edge_descriptor e){


	double start_llik = current_llik;

	double min, max, guess;
	guess = tree->g[e].weight;
	guess = log(guess/ (1-guess));
	min = params->minweight;
	max = params->maxweight;
	int nit = 0;
	golden_section_weight(e, min, guess, max, params->tau, &nit);


}


void GraphState2::optimize_weight_quick(Graph::edge_descriptor e){
	if (!params->f2) optimize_weight(e);
	else{
		double start_llik = current_llik;
		double min, max, guess, start_weight;
		guess = tree->g[e].weight;
		start_weight = guess;
		guess =log(guess/ (1-guess));

		min = params->minweight;
		max = params->maxweight;
		int nit = 0;
		golden_section_weight_quick(e, min, guess, max, params->tau, &nit);

		if (start_llik > current_llik){
		//cout << "optim_weight did not improve "<< start_llik << " "<< current_llik << " "<< start_weight << " "<< tree->g[e].weight << "\n"; cout.flush();
			update_mig(e, start_weight);
			current_llik = start_llik;
		}
	}

}


void GraphState2::optimize_frac(Graph::edge_descriptor e){


		double start_llik = current_llik;
		bool done = false;
		int nit = 0;
		while(!done && nit < params->maxit){
			//cout << nit << "\n"; cout.flush();
			double min, max, guess;
			guess = tree->g[source(e, tree->g)].mig_frac;
			guess = log (guess/(1-guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac(e, min, guess, max, params->tau);
			if (current_llik < start_llik+0.01) done = true;
			else start_llik = current_llik;
			//cout << guess << " "<< start_llik << " "<< current_llik << "\n";
			nit++;
		}


}



void GraphState2::optimize_frac_wish(Graph::edge_descriptor e){


		double start_llik = current_llik_w;
		bool done = false;
		int nit = 0;
		while(!done && nit < params->maxit){
			//cout << nit << "\n"; cout.flush();
			double min, max, guess;
			guess = tree->g[source(e, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac_wish(e, min, guess, max, params->tau);
			if (current_llik_w < start_llik+0.01) done = true;
			else start_llik = current_llik_w;
			//cout << guess << " "<< start_llik << " "<< current_llik << "\n";
			nit++;
		}


}

int GraphState2::golden_section_weight(Graph::edge_descriptor e, double min, double guess, double max, double tau, int* nit){
	double x;

	//cout << guess << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max)) || *nit > params->maxit) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));
		tree->g[e].weight = neww;
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();

		current_llik = llik();
		return 0;
	}
	*nit = *nit+1;
	double w = 1/(1+exp(-x));
	tree->g[e].weight = w;
	//cout << "here\n"; cout.flush();
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();

	//cout << "not here\n"; cout.flush();
	double f_x = -llik();
	//cout << w << " "<< -llik() << "\n";
	w = 1/(1+exp(-guess));
	tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik();
	//cout << guess << " "<< -llik() << "\n";

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_weight(e, guess, x, max, tau, nit);

		else return golden_section_weight(e, min, x, guess, tau, nit);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight(e, min, guess, x, tau, nit);
		else return golden_section_weight(e, x, guess, max, tau, nit);
	}
}

int GraphState2::golden_section_weight_quick(Graph::edge_descriptor e, double min, double guess, double max, double tau, int* nit){
	double x;


	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max)) || *nit > params->maxit) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));

		update_mig(e, neww);
		current_llik = llik();
		if (current_llik < DBL_MIN){
			neww = 1/ (1+exp(-guess));
			update_mig(e, neww);
			current_llik = llik();
		}
		//cout << "returning "<< neww << " "<< current_llik<< "\n";
		return 0;
	}
	//cout << *nit << "\n";
	*nit = *nit+1;
	double w = 1/(1+exp(-x));
	update_mig(e, w);
	double f_x = -llik();
	//cout << w << " x " << f_x << "\n";
	w = 1/(1+exp(-guess));
	update_mig(e, w);
	double f_guess = -llik();
	//cout << w << " guess " << f_guess << "\n";

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_weight_quick(e, guess, x, max, tau, nit);

		else return golden_section_weight_quick(e, min, x, guess, tau, nit);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight_quick(e, min, guess, x, tau, nit);
		else return golden_section_weight_quick(e, x, guess, max, tau, nit);
	}
}


int GraphState2::golden_section_weight_noexp(Graph::edge_descriptor e, double min, double guess, double max, double tau, int *nit){
	double x;

	//cout << min << " "<< guess << " "<< max << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max)) || *nit > params->maxit) {
		//cout << "here "<< min << " "<< max << "\n";
		double new_logweight = (min+max)/2;
		double neww = new_logweight;
		tree->g[e].weight = neww;
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		current_llik = llik();
		return 0;
	}
	*nit = *nit+1;
	double w = x;
	tree->g[e].weight = w;
	//cout << "here\n"; cout.flush();
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();

	//cout << "not here\n"; cout.flush();
	double f_x = -llik();

	w = guess;
	tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik();
	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_weight_noexp(e, guess, x, max, tau, nit);

		else return golden_section_weight_noexp(e, min, x, guess, tau, nit);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight_noexp(e, min, guess, x, tau, nit);
		else return golden_section_weight_noexp(e, x, guess, max, tau, nit);
	}
}


int GraphState2::golden_section_weight_noexp_quick(Graph::edge_descriptor e, double min, double guess, double max, double tau, int *nit){
	double x;

	//cout << min << " "<< guess << " "<< max << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max)) || *nit > params->maxit) {
		//cout << "here "<< min << " "<< max << "\n";
		double new_logweight = (min+max)/2;
		double neww = new_logweight;
		update_mig(e, neww);

		current_llik = llik();
		return 0;
	}
	*nit = *nit+1;
	double w = x;
	update_mig(e, w);
	//cout << "not here\n"; cout.flush();
	double f_x = -llik();

	w = guess;
	tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	update_mig(e, w);
;
	double f_guess = -llik();
	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_weight_noexp_quick(e, guess, x, max, tau, nit);

		else return golden_section_weight_noexp_quick(e, min, x, guess, tau, nit);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight_noexp_quick(e, min, guess, x, tau, nit);
		else return golden_section_weight_noexp_quick(e, x, guess, max, tau, nit);
	}
}

int GraphState2::golden_section_frac(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;

	//cout << guess << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));
		tree->set_mig_frac(e, neww);
		set_branches_ls_wmig();
		current_llik = llik();
		return 0;
	}
	double w = 1/(1+exp(-x));
	tree->set_mig_frac(e, w);
	//cout << "here\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here\n"; cout.flush();
	double f_x = -llik();

	w = 1/(1+exp(-guess));
	tree->set_mig_frac(e, w);
	//tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik();

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_frac(e, guess, x, max, tau);

		else return golden_section_frac(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_frac(e, min, guess, x, tau);
		else return golden_section_frac(e, x, guess, max, tau);
	}
}


int GraphState2::golden_section_edge(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;

	//cout << min << " "<< guess << " "<< max << " "<< tau << " "<< resphi << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_loglen = (min+max)/2;
		double neww = exp(new_loglen);
		tree->g[e].len = neww;
		compute_sigma();
		set_sigmacor_from_sigma();
		current_llik = llik();
		return 0;
	}
	double w = exp(x);
	tree->g[e].len = w;
	compute_sigma();
	set_sigmacor_from_sigma();

	double f_x = -llik();
	//cout << x << " " << f_x << "\n";
	w = exp(guess);
	tree->g[e].len = w;
	compute_sigma();
	set_sigmacor_from_sigma();

	double f_guess = -llik();
	//cout << f_guess << " "<< f_x << "\n";
	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_edge(e, guess, x, max, tau);

		else return golden_section_edge(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_edge(e, min, guess, x, tau);
		else return golden_section_edge(e, x, guess, max, tau);
	}
}


int GraphState2::golden_section_frac_wish(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;

	//cout << guess << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));
		tree->set_mig_frac(e, neww);
		set_branches_ls_wmig();
		current_llik_w = llik(true);
		return 0;
	}
	double w = 1/(1+exp(-x));
	tree->set_mig_frac(e, w);
	//cout << "here\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here\n"; cout.flush();
	double f_x = -llik(true);

	w = 1/(1+exp(-guess));
	tree->set_mig_frac(e, w);
	//tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik(true);

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_frac(e, guess, x, max, tau);

		else return golden_section_frac(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_frac(e, min, guess, x, tau);
		else return golden_section_frac(e, x, guess, max, tau);
	}
}

void GraphState2::set_branches_ls_wmig_estmig(){

	/* one parameter for each non-migration node (minus 1 for the root and 1 for the unidentifiable branch length next to the root)

	   Complication when doing this: paths to root in terms of edges (migration coming into nodes makes nodes not possible).
	   Many edge lengths are not identifiable, so have a single parameter which is their sum. Need to figure out which edge goes with which parameter, how to weight them.
    */


	map<Graph::vertex_descriptor, int> vertex2index;
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	//cout << "here1.1\n"; cout.flush();
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	//cout << "here1.2\n"; cout.flush();
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;
	//cout << "here1.1\n"; cout.flush();
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}
	//cout << "here2\n"; cout.flush();
	// now get the edge to index and edge to fraction maps
	map<Graph::edge_descriptor, int> edge2index;
	map<Graph::edge_descriptor, double> edge2frac;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		Graph::vertex_descriptor index_vertex, t;
		double f = 1.0;
		double f2 = 1.0;
		t = target(*it, tree->g);
		Graph::vertex_descriptor t2 = source(*it, tree->g);
		//cout << tree->g[t].index << "\n";
		if (tree->g[*it].is_mig) continue;
		else if ( tree->g[t].is_mig) {
			index_vertex = tree->get_child_node_mig(t);
		}
		else index_vertex = t;

		if (tree->g[t].is_mig || tree->g[t2].is_mig){
			if (tree->g[t].is_mig && tree->g[t2].is_mig)	{
				//cout << "here?\n"; cout.flush();
				f = tree->g[t].mig_frac - tree->g[t2].mig_frac;
			}
			else if (tree->g[t].is_mig) {
				//cout << "here2 "<< tree->g[t].mig_frac << "\n";
				f = tree->g[t].mig_frac;
			}
			else if (tree->g[t2].is_mig) {
				//cout << "here3 "<< tree->g[t2].mig_frac << "\n";
				f = 1-tree->g[t2].mig_frac;
			}
		}

		//f = tree->g[*it].len / tree->get_parent_node(index_vertex).second;

		//cout << tree->g[t2].index << " "<< tree->g[t].index << " "<< f<< " "<< tree->g[*it].len << "\n";
		//cout << tree->g[source(*it, tree->g)].index << " "<< tree->g[t].index << " "<< f << " "<< f2 << "\n";
		//cout << tree->g[t].index << "\n";
		if (root_adj.find(index_vertex) != root_adj.end()) f = f/2;
		int i;
		if (vertex2index.find(index_vertex) == vertex2index.end()){
			cerr << "ERROR: Error in least squares estimation: vertex "<< tree->g[index_vertex].index << " not found in the list of vertices\n";
			exit(1);
		}
		else i = vertex2index[index_vertex];
		//cout << tree->g[source(*it, tree->g)].index<< " "<< tree->g[target(*it, tree->g)].index << " "<<  tree->g[index_vertex].index << " "<< f<< "\n";
		edge2index.insert(make_pair(*it, i));
		edge2frac.insert(make_pair(*it, f));
	}
	//cout <<  "\n";
	//cout << "here3\n"; cout.flush();
	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3; // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);

	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations

	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);
	//cout << "here\n";
	set_branch_coefs(X, y, &edge2index, &edge2frac);
	//cout << "here2\n";
	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);
	//and put in the solutions in the graph
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = gsl_vector_get(c, i);
		tree->g[*it].len = l*frac;
	}
	//cout << "here3.1\n";
	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X, index, k) * gsl_vector_get(c, k);
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			index++;
		}

	}

	//cout << "here4\n";
	//free memory
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}


void GraphState2::set_branches_ls(){

	/*
	 * Also estimate branch lengths to/from migration nodes
    */

	//cout << "HERE\n"; cout.flush();
	map<Graph::edge_descriptor, int> edge2index;
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;
	map<Graph::edge_descriptor, double> edge2frac;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			edge2index.insert(make_pair(*it, index));
			edge2frac.insert(make_pair(*it, 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		edge2index[*it] = joint_index;
		edge2frac.insert(make_pair(*it, 0.5));
	}
	index++;

	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2 + current_npops; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated

	//set up the workspace
	// solve Ax = b  with x >=0
	// A = weighted branches
	// x = branch lengths (to be solved)
	// b = observed f_2 matrix
	//cout << p << " P \n";
	//cout << n << " N \n";
	NNLS_SOLVER nnls(n, p);


	double * A = new double[n*p]; //A holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root.

	double * b  = new double[n];  // y contains the entries of the empirical covariance matrix
	double * x = new double[p];   // x will be estimated, contains the fitted branch lengths
	double rNorm;

	set_branch_coefs_nnls(A, b, n, p, &edge2index, &edge2frac);

	//fit NNLS
	bool converged = nnls.solve(A, p, b, x, rNorm);
	//cout << converged << " conv\n";


	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = x[i];
		tree->g[*it].len = l*frac;
		//if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += A[ k * n + index ] * x[k];
			}
			//cout << i << " "<< j << " "<< pred << "\n";
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//free memory

	delete [] A;
	delete [] b;
	delete [] x;

}



void GraphState2::set_branches_ls_f2_old(){

	/*
	 * Also estimate branch lengths to/from migration nodes
    */

	//cout << "HERE\n"; cout.flush();
	map<Graph::edge_descriptor, int> edge2index;
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;
	map<Graph::edge_descriptor, double> edge2frac;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			edge2index.insert(make_pair(*it, index));
			edge2frac.insert(make_pair(*it, 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		edge2index[*it] = joint_index;
		edge2frac.insert(make_pair(*it, 0.5));
	}
	index++;

	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated

	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations

	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	//cout << "here\n";
	set_branch_coefs_f2(X, y, &edge2index, &edge2frac);

	/*for (int i = 0 ; i < n ; i++){
		cout << gsl_vector_get(y, i) << " " << "y"<<" ";
		for (int j = 0; j < p ; j++){
			cout << gsl_matrix_get(X, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n\n";
	*///cout << "here2\n";
	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);
	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = gsl_vector_get(c, i);
		tree->g[*it].len = l*frac;
		if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
		if (l < 0) 	negsum += -l;


	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X, index, k) * gsl_vector_get(c, k);
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//cout << "here4\n";
	//free memory
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}



void GraphState2::set_branches_ls_f2(){

	/*
	 * Also estimate branch lengths to/from migration nodes
    */

	//cout << "HERE\n"; cout.flush();
	map<Graph::edge_descriptor, int> edge2index;
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;
	map<Graph::edge_descriptor, double> edge2frac;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			edge2index.insert(make_pair(*it, index));
			edge2frac.insert(make_pair(*it, 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		edge2index[*it] = joint_index;
		edge2frac.insert(make_pair(*it, 0.5));
	}
	index++;

	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated

	//set up the workspace
	// solve Ax = b  with x >=0
	// A = weighted branches
	// x = branch lengths (to be solved)
	// b = observed f_2 matrix

	NNLS_SOLVER nnls(n, p);


	double * A = new double[n*p]; //A holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root.

	double * b  = new double[n];  // y contains the entries of the empirical covariance matrix
	double * x = new double[p];   // x will be estimated, contains the fitted branch lengths
	double rNorm;

	set_branch_coefs_f2_nnls(A, b, n, p, &edge2index, &edge2frac);

	//fit NNLS
	bool converged = nnls.solve(A, p, b, x, rNorm);
	//cout << converged << " conv\n";

	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = x[i];
		tree->g[*it].len = l*frac;
		//if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += A[ k * n + index ] * x[k];
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//free memory

	delete [] A;
	delete [] b;
	delete [] x;

}



double GraphState2::llik_normal(){
	double toreturn = 0;
	for (int i = 0; i < current_npops; i++){
		int stj = i;
		if (params->f2) stj++;
		for (int j = stj; j < current_npops; j++){

			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			double pred = gsl_matrix_get(sigma_cor, i, j);
			double obs = countdata->get_cov(p1, p2);
			double se = countdata->get_cov_var(p1, p2);

			double dif = obs-pred;
			//double scale = params->smooth_scale;
			double toadd = lndgauss(dif, se);
			toreturn += toadd;
			//double toadd = gsl_ran_gaussian_pdf(dif, se * scale);
			//toreturn+= log(toadd);

			//}
		}
	}
	return toreturn;
}

double lndgauss(double dif, double se){
	double toreturn = 0;
	toreturn += -log (se * sqrt(2.0*M_PI));
	toreturn += -(dif*dif) / (2*se*se);
	return toreturn;
}
int GraphState2::local_hillclimb(int inorder_index){
	// if there was a rearrangement, return 1. otw 0.
	//

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	if (params->cor_mig){
		Graph::vertex_descriptor v = inorder[inorder_index];
		Graph::vertex_descriptor p = tree->get_parent_node(v).first;
		if (tree->g[p].is_root || tree->g[v].is_root) {
			return 0;
		}
	}
	if ( tree->g[ inorder[inorder_index]].is_root) return local_hillclimb_root();

	double llik1, llik2, llik3;

	tree_bk->copy(tree);


	llik1 = current_llik;

	tree->local_rearrange(inorder[inorder_index], 1);

	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	llik2 =  llik();
	tree_bk2->copy(tree);


	tree->copy(tree_bk);
	inorder = tree->get_inorder_traversal(current_npops);
	tree->local_rearrange(inorder[inorder_index], 2);

	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	llik3 =  llik();

	tree_bk3->copy(tree);
	tree->copy(tree_bk);

	if (llik1 < llik2 || llik1 < llik3){
		if (llik2 > llik3){
			tree->copy(tree_bk2);
			current_llik = llik2;
			return 1;
		}
		else{
			tree->copy(tree_bk3);
			current_llik = llik3;
			return 1;
		}
	}
	return 0;
}

int GraphState2::local_hillclimb_root(){
	double best_llik = current_llik;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);

	int toreturn = 0;
	if (params->cor_mig) return toreturn;
	for (int i = 1; i <= 4; i++){
		tree->move_root(i);
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		double tmplik = llik();
		if (tmplik > best_llik){
			toreturn = 1;
			tree_bk->copy(tree);
			best_llik = tmplik;
		}
		tree->copy(tree_bk2);
	}
	if (toreturn == 1){
		tree->copy(tree_bk);
		current_llik = best_llik;
	}
	return toreturn;
}

int GraphState2::global_hillclimb(int inorder_index){
	// take the node in inorder_index, try putting it above all the other nodes
	// return 1 if there is an improvement, 0 otw
	double max = current_llik;
	gsl_matrix* tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);

	int maxindex = inorder_index;
	tree_bk ->copy(tree);
	int toreturn = 0;
	for (int i = 0; i < 2*current_npops-1; i++){
		vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
		Graph::vertex_descriptor v1 = inorder[inorder_index];
		Graph::vertex_descriptor v2 = inorder[i];

		//make sure its a reasonable move (can't attach a clade within itself)
		if (i == inorder_index) continue;
		if ( tree->get_parent_node( v1 ).first == tree->get_parent_node( v2 ).first) continue;
		if ( tree->get_parent_node(v1).first == v2) continue;

		set<Graph::vertex_descriptor> path = tree->get_path_to_root( v2 );
		if ( path.find(v1) != path.end() ) continue;
		if (!try_mig(v1, v2, tmpfitted)) continue;
		cout << tree->get_newick_format(v1)<< " "<< tree->get_newick_format(v2)<< "\n";
		tree->global_rearrange(v1, v2);
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		double lk =  llik();
		cout << lk << " "<<  max << "\n";
		if ( lk > max ){
			max = lk;
			maxindex = i;
			toreturn = 1;
		}
		tree->copy(tree_bk);
	}
	if (toreturn ==1){
		vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
		tree->global_rearrange( inorder[ inorder_index], inorder[maxindex] );
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		double lk =  llik();
		current_llik = lk;
		cout << tree->get_newick_format() <<"\n";
		cout << "ln(lk): "<< lk <<"\n"; cout.flush();
	}
	gsl_matrix_free(tmpfitted);
	return toreturn;

}

int GraphState2::iterate_local_hillclimb_wmig(pair< set<int>, set<int> > indices){
	double start_lik = current_llik;
	int moving = many_local_hillclimb_wmig(indices);
	//cout << moving <<" first iter\n";
	if (moving  == 0 || current_llik < start_lik + params->epsilon) return 0;
	start_lik = current_llik;
	int nit = 0;
	while (moving > 0 && nit < params->maxit2) {
		moving = many_local_hillclimb_wmig(indices);
		nit++;
		if ( current_llik < start_lik +params->epsilon) moving = 0;
		else start_lik = current_llik;
	}

	if (nit >= params->maxit2) return 0;
	return 1;
}


int GraphState2::iterate_local_hillclimb_wmig_all(){
	double start_lik = current_llik;
	int moving = many_local_hillclimb_wmig_all();
	if (moving  == 0 || current_llik < start_lik+ params->epsilon) return 0;
	start_lik = current_llik;
	while (moving > 0) {
		moving = many_local_hillclimb_wmig_all();
		if (current_llik < start_lik +params->epsilon) moving = 0;
		start_lik = current_llik;
	}
	return 1;
}

void GraphState2::iterate_mig_hillclimb_and_optimweight(pair<int, int> indices, double start_nsum){
	pair< set<int>, set<int> > p1 = get_neighborhood(indices.first);
	pair< set<int>, set<int> > p2 = get_neighborhood(indices.second);
	//tree->print("incoming");
	int moving1 = iterate_local_hillclimb_wmig(p1);
	cout << "Local updates around node 1: "<< moving1 << " ln(lk):"<< current_llik << " \n"; cout.flush();
	if (moving1> 0) p1 = get_neighborhood(indices.first);

	//tree->print("after_1");

	int moving2 = iterate_local_hillclimb_wmig(p2);
	cout << "Local updates around node 2: "<< moving2 << " ln(lk):"<< current_llik << " \n"; cout.flush();
	if (moving2 > 0) 	p2 = get_neighborhood(indices.second);

	//tree->print("after_2");

	int moving4 = all_try_movemig();
	cout << "Updates in migration position: "<< moving4 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();
	if (moving4 > 0) {
		p1 = get_neighborhood(indices.first);
		p2 = get_neighborhood(indices.second);
	}
	//tree->print("after_move");
	int moving3 = all_try_changedir();
	cout << "Switches in migration direction: "<< moving3 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();
	if (moving3 > 0){
		p2 = get_neighborhood(indices.first);
		p1 = get_neighborhood(indices.second);
	}

	int moving = moving1+moving2+moving3+moving4;
	int moving_local = moving;
	while (moving > 0){

			moving1 = iterate_local_hillclimb_wmig(p1);
			cout << "Local updates around node 1: "<< moving1 << " ln(lk):"<< current_llik << " \n"; cout.flush();
			if (moving1> 0) p1 = get_neighborhood(indices.first);

			///tree->print("after1");
			moving2 = iterate_local_hillclimb_wmig(p2);
			if (moving2 > 0) p2 = get_neighborhood(indices.second);

			cout << "Local updates around node 2: "<< moving2 << " ln(lk):"<< current_llik << " \n"; cout.flush();
			//tree->print("after2");
			int moving3 = all_try_changedir();
			cout << "Switches in migration direction: "<< moving3 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();
			if (moving3 > 0){
				p2 = get_neighborhood(indices.first);
				p1 = get_neighborhood(indices.second);

			}

			int moving4 = all_try_movemig();
			cout << "Updates in migration position: "<< moving4 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();
			if (moving4 > 0) {
				p1 = get_neighborhood(indices.first);
				p2= get_neighborhood(indices.second);
			}
			//tree->print("after4");
			moving = moving1+moving2+moving3+moving4;

			int moving5 = 0;
			if (moving1 == 0 and moving2 == 0 and moving3 ==0 and moving4==0){
				cout << "Trying all local rearrangements\n"; cout.flush();
				moving5 = iterate_local_hillclimb_wmig_all();
				cout << "Local rearrangements: "<< moving5 << " ln(lk):" << current_llik << " \n";cout.flush();
				moving = moving+moving5;
			}
	}


}

void GraphState2::iterate_all_hillclimb(){

	//tree->print("start");
	optimize_weights();

	cout << "Trying all local rearrangements\n"; cout.flush();
	int moving5 = iterate_local_hillclimb_wmig_all();
	cout << "Local rearrangements: "<< moving5 << " ln(lk):" << current_llik << " \n";cout.flush();

	int moving4 = all_try_movemig_limit();
	cout << "Updates in migration position: "<< moving4 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();

	//tree->print("after_move0");
//	int moving3 = all_try_changedir();
	//cout << "Switches in migration direction: "<< moving3 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();





	int moving = moving4+moving5;
	int iter = 1;
	while (moving > 0){

		cout << "Trying all local rearrangements\n"; cout.flush();
		moving5 = iterate_local_hillclimb_wmig_all();
		cout << "Local rearrangements: "<< moving5 << " ln(lk):" << current_llik << " \n";cout.flush();
			//moving3 = all_try_changedir();
			//cout << "Switches in migration direction: "<< moving3 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();

		moving4 = all_try_movemig_limit();
		cout << "Updates in migration position: "<< moving4 <<  " ln(lk):"<< current_llik << " \n"; cout.flush();
		stringstream ss;
		ss << "after_move";
		ss << iter;
		string tmp;
		tmp = ss.str();
		//tree->print(tmp);
		moving = moving4+moving5;
		iter++;
	}
}

int GraphState2::many_local_hillclimb_wmig(pair<set<int>, set<int> > indices){
	int toreturn = 0;
	for (set<int>::iterator it = indices.first.begin(); it != indices.first.end(); it++){
		toreturn+= local_hillclimb_wmig(*it, indices.second);
	}
	return toreturn;
}


int GraphState2::many_local_hillclimb_wmig_all(){
	int toreturn = 0;
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	vector<int> allindex;
	for (vector<Graph::vertex_descriptor>::iterator it = inorder.begin(); it != inorder.end(); it++) {
		if (tree->g[*it].is_root) continue;
		allindex.push_back(tree->g[*it].index);
	}

	for (vector<int>::iterator it = allindex.begin(); it != allindex.end(); it++){
		toreturn+= local_hillclimb_wmig_all(*it);
	}
	return toreturn;
}

int GraphState2::all_try_changedir(){
	int toreturn =0;
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	int i = 0;
	while (i < mig_edges.size()){
		//cout << "here\n"; cout.flush();
		//cout << tree->g[source(mig_edges[i], tree->g)].index << "\n";
		//tree->print("test0");
		int test = try_changedir(mig_edges[i]);
		//cout <<"not here\n"; cout.flush();
		mig_edges = tree->get_mig_edges();
		if (test ==1) {
			toreturn++;
			i=0;
		}
		else i++;
	}
	return toreturn;
}


int GraphState2::all_try_movemig(){
	int toreturn =0;
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	vector<int> mig_index;

	int i = 0;
	while (i < mig_edges.size() && toreturn < params->maxit2){
		//cout << i << " " << current_llik << " "<< llik() << " "<< tree->g[source(mig_edges[i], tree->g)].index <<  "\n";
		bool test = movemig(tree->g[source(mig_edges[i], tree->g)].index).first;
		mig_edges = tree->get_mig_edges();
		if (test) {
			//tree->print("moved");
			toreturn++;
			i=0;
		}
		else i++;
	}
	if ( toreturn >= params->maxit2) return 0;
	return toreturn;
}


int GraphState2::all_try_movemig_limit(){
	int toreturn =0;
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	vector<int> mig_index;

	int i = 0;
	while (i < mig_edges.size() && toreturn < params->maxit2){
		//cout << i << " " << current_llik << " "<< llik() << " "<< tree->g[source(mig_edges[i], tree->g)].index <<  "\n";
		bool test = movemig_limit(tree->g[source(mig_edges[i], tree->g)].index).first;
		mig_edges = tree->get_mig_edges();
		if (test) {
			//tree->print("moved");
			toreturn++;
			i=0;
		}
		else i++;
	}
	if ( toreturn >= params->maxit2) return 0;
	return toreturn;
}

int GraphState2::try_changedir(Graph::edge_descriptor e){
	double max = current_llik;
	double lik_bk = current_llik;
	double negsum_bk = negsum;
	double max_negsum = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	int toreturn = 0;
	//tree->print("changedir");
	Graph::vertex_descriptor s = source(e, tree->g);
	Graph::vertex_descriptor t = target(e, tree->g);
	//cout << tree->g[s].index << " "<<tree->g[t].index << "\n";
	Graph::vertex_descriptor newt = tree->get_child_node_mig(s);
	tree->remove_mig_edge(e);
	if (!tree->is_legal_migration(t, newt)){
		tree->copy(tree_bk);
		return toreturn;
	}
	Graph::edge_descriptor e2 = tree->add_mig_edge(t, newt);

	if (params->f2){
		initialize_migupdate();
		optimize_weight_quick(e2);
	}
	else optimize_weight(e2);
	double lk = llik();
	//cout << "here3\n";
	if (lk > max){
		max = lk;
		max_negsum = negsum;
		//cout << i << "\n";
		toreturn = 1;
		tree_bk2->copy(tree);
		gsl_matrix_memcpy( tmpfitted, sigma_cor);
	}
	tree->copy(tree_bk);
	if (toreturn == 1){
		tree->copy(tree_bk2);
		current_llik = max;
		negsum = max_negsum;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	else{
		tree->copy(tree_bk);
		current_llik = lik_bk;
		negsum = negsum_bk;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	gsl_matrix_free(tmpfitted);
	return toreturn;

}

int GraphState2::local_hillclimb_wmig(int index){
	double max = current_llik;
	double lik_bk = current_llik;
	double max_negsum = negsum;
	double negsum_bk = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	int toreturn = 0;

	for (int i = 1; i <=4; i++){
		map<int, Graph::vertex_descriptor> index2v = tree->index2vertex();
		Graph::vertex_descriptor v = index2v[index];
		bool rearr = tree->local_rearrange_wmig(v, i);

		if ( has_loop() ) {
			tree->copy(tree_bk);
			continue;
		}


		if (rearr){
			optimize_weights_quick();
			//cout << "here2\n";
			double lk = llik();
			if (lk > max){
				max = lk;
				max_negsum = negsum;
				//cout << i << "\n";
				toreturn = 1;
				tree_bk2->copy(tree);
				gsl_matrix_memcpy( tmpfitted, sigma_cor);
			}
			tree->copy(tree_bk);
		}
	}
	if (toreturn == 1){
		tree->copy(tree_bk2);
		current_llik = max;
		negsum = max_negsum;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	else{
		tree->copy(tree_bk);
		current_llik = lik_bk;
		negsum = negsum_bk;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	gsl_matrix_free(tmpfitted);
	return toreturn;
}


int GraphState2::local_hillclimb_wmig_all(int index){
	double max = current_llik;
	double lik_bk = current_llik;
	double max_negsum = negsum;
	double negsum_bk = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	int toreturn = 0;
	//tree->print("test0");

	for (int i = 1; i <=2; i++){
		map<int, Graph::vertex_descriptor> index2v = tree->index2vertex();
		Graph::vertex_descriptor v = index2v[index];
		pair<set<int>, set<int> > ne = get_neighborhood(index, 2);
		if (tree->g[v].is_root || tree->g[v].is_tip) return toreturn;
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch = tree->get_child_nodes(v);
		Graph::vertex_descriptor p = tree->get_parent_node(v).first;
		if (tree->g[p].is_root) return toreturn;
		if (i == 1){
			bool rearr = tree->local_rearrange_wmig(ch.first, 4);
			//tree->print("test0");
			if ( has_loop() ) {
				//cout << "has loop!\n";
				tree->copy(tree_bk);
				continue;
			}
			if (rearr){

				if (params->f2){
					set_branches_ls_f2();
					optimize_weights_quick(ne.second);
				}
				else{
					set_branches_ls();
					if (params->fitmig) optimize_weights(ne.second);
				}
				//tree->print("test0");
				double lk = llik();
				//cout << i << " "<< lk << " "<< max << "\n";
				if (lk > max){
					max = lk;
					max_negsum = negsum;

					toreturn = 1;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted, sigma_cor);
				}
				tree->copy(tree_bk);
			}
		}

		else if (i == 2){
			bool rearr = tree->local_rearrange_wmig(ch.second, 4);
			if ( has_loop() ) {

				tree->copy(tree_bk);
				continue;
			}
			if (rearr){

				if (params->f2){
					set_branches_ls_f2();
					optimize_weights_quick(ne.second);
				}
				else{
					set_branches_ls();
					if (params->fitmig) optimize_weights(ne.second);
				}

				double lk = llik();

				//cout << i << " "<< lk << " "<< max << "\n";
				if (lk > max){
					max = lk;
					max_negsum = negsum;
					toreturn = 1;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted, sigma_cor);
				}
				tree->copy(tree_bk);
			}
		}
	}
	if (toreturn == 1){
		tree->copy(tree_bk2);
		current_llik = max;
		negsum = max_negsum;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	else{
		tree->copy(tree_bk);
		current_llik = lik_bk;
		negsum = negsum_bk;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	gsl_matrix_free(tmpfitted);
	return toreturn;
}



int GraphState2::local_hillclimb_wmig(int index, set<int> n){
	double max = current_llik;
	double lik_bk = current_llik;
	double max_negsum = negsum;
	double negsum_bk = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	int toreturn = 0;

	for (int i = 1; i <=2; i++){
		map<int, Graph::vertex_descriptor> index2v = tree->index2vertex();
		Graph::vertex_descriptor v = index2v[index];
		if (tree->g[v].is_root || tree->g[v].is_tip) return toreturn;
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch = tree->get_child_nodes(v);
		Graph::vertex_descriptor p = tree->get_parent_node(v).first;
		if (tree->g[p].is_root) return toreturn;
		if (i == 1){
			bool rearr = tree->local_rearrange_wmig(ch.first, 4);
			//tree->print("test0");
			if ( has_loop() ) {
				//cout << "has loop!\n";
				tree->copy(tree_bk);
				continue;
			}
			if (rearr){
				if (params->f2){
					set_branches_ls_f2();
					optimize_weights_quick(n);
				}
				else optimize_weights(n);
				//tree->print("test0");
				double lk = llik();
				//cout << i << " "<< lk << " "<< max << "\n";
				if (lk > max){
					max = lk;
					max_negsum = negsum;

					toreturn = 1;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted, sigma_cor);
				}
				tree->copy(tree_bk);
			}
		}

		else if (i == 2){
			bool rearr = tree->local_rearrange_wmig(ch.second, 4);
			if ( has_loop() ) {

				tree->copy(tree_bk);
				continue;
			}
			if (rearr){

				if (params->f2){
					set_branches_ls_f2();
					optimize_weights_quick(n);
				}
				else optimize_weights(n);

				//tree->print("test1");

				double lk = llik();

				//cout << i << " "<< lk << " "<< max << "\n";
				if (lk > max){
					max = lk;
					max_negsum = negsum;
					toreturn = 1;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted, sigma_cor);
				}
				tree->copy(tree_bk);
			}
		}
	}

	if (toreturn == 1){
		tree->copy(tree_bk2);
		current_llik = max;
		negsum = max_negsum;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	else{
		tree->copy(tree_bk);
		current_llik = lik_bk;
		negsum = negsum_bk;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	gsl_matrix_free(tmpfitted);
	return toreturn;
}
int GraphState2::many_local_hillclimb(){
	int leninorder = 2*current_npops -1;
	int toreturn = 0;
	for (int i = 1; i < leninorder; i+=2){
		toreturn += local_hillclimb(i);
	}
	return toreturn;
}

int GraphState2::many_global_hillclimb(){
	int leninorder = 2*current_npops -1;
	int toreturn = 0;
	for (int i = 0; i < leninorder; i++){
		toreturn += global_hillclimb(i);
	}
	return toreturn;
}

void GraphState2::iterate_hillclimb(){
	//cout << "here\n"; cout.flush();
	int changes = many_local_hillclimb();
	//cout << "not here1\n"; cout.flush();
	//tree->print("test1");
	//cout << "Hill climbing "<< changes << " changes\n";
	while (changes > 0) {
		changes = many_local_hillclimb();
		//cout << "Hill climbing "<< changes << " changes\n";
	}
	//cout << "not here2\n"; cout.flush();
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	//cout << "nothere\n\n";
}


void GraphState2::iterate_global_hillclimb(){
	int changes = many_global_hillclimb();
	while (changes > 0) changes = many_global_hillclimb();

}

void GraphState2::add_pop(){

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma_cor);

	//cout << "here?\n"; cout.flush();
	//current_npops++;
	//process_scatter();
	//current_npops--;
	//cout << "not here?\n"; cout.flush();

	string toadd = allpopnames[current_npops];

	cout << "Adding "<< toadd << " ["<< current_npops+1 << "/" << allpopnames.size() <<"]\n";

	//If adding with a migration event
	if (params->cor_mig && params->migfracs.find(toadd) != params->migfracs.end()){
		//If this is the first time and the root isn't yet set, set it

		if (params->migfracs.find(allpopnames[current_npops-1]) == params->migfracs.end() && params->set_root){
			//tree->print("before");
			place_root(params->root);
		}

		add_pop(toadd, params->migfracs[toadd]);
		//tree->print("after2");
	}
	else{
		vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
		int max_index;
		double max_llik;
		tree_bk->copy(tree);

		for (int i = 0; i < inorder.size(); i++){

			Graph::vertex_descriptor tmp = tree->add_tip(inorder[i], toadd);
			current_npops++;
			if (params->f2) set_branches_ls_f2();
			else set_branches_ls();
			double llk = llik();
			if (i == 0){
				max_index = i;
				max_llik = llk;
			}
			else if (llk > max_llik){
				max_index = i;
				max_llik = llk;
			}
			tree->copy(tree_bk);
			current_npops--;
			inorder = tree->get_inorder_traversal(current_npops);
		}
		if (max_llik <= -DBL_MAX){
			cerr <<"ERROR: numerical problem in likelihood [ln(lk) = "<< max_llik << "]\n";
			cerr <<"ERROR: Please report this bug.\n";
			exit(1);
			/*
			cout << "RESCALING\n"; cout.flush();
			params->smooth_scale = params->smooth_scale *2;
			for (int i = 0; i < inorder.size(); i++){

				Graph::vertex_descriptor tmp = tree->add_tip(inorder[i], toadd);
				current_npops++;
				if (params->f2) set_branches_ls_f2();
				else set_branches_ls();
				double llk = llik();
				if (i == 0){
					max_index = i;
					max_llik = llk;
				}
				else if (llk > max_llik){
					max_index = i;
					max_llik = llk;
				}
				tree->copy(tree_bk);
				current_npops--;
				inorder = tree->get_inorder_traversal(current_npops);
			}
			*/
		}
		Graph::vertex_descriptor tmp = tree->add_tip(inorder[max_index], toadd);
		current_npops++;

		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		current_llik = max_llik;
	}
}

void GraphState2::add_pop(string name, map<string, double> migfracs){
	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);

	for (map<string, double>::iterator it = migfracs.begin(); it != migfracs.end(); it++){
		if (tips.find(it->first) == tips.end()){
			cerr << "ERROR: cannot find population "<< it->first << "\n";
			exit(1);
		}
	}

	string toadd = name;
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	int max_index;
	double max_llik;
	tree_bk->copy(tree);
	int first = 0;
	for (int i = 0; i < inorder.size(); i++){

		if ( params->set_root && params->root ==tree->g[ inorder[i] ].name) { //
			if (i == first) first++;
			continue;
		}
		if ( params->set_root && tree->g[ inorder[i] ].is_root) { //
			if (i == first) first++;
			continue;
		}
		Graph::vertex_descriptor tmp = tree->add_tip(inorder[i], toadd);
		map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();

		bool stillok = true;
		for (map<string, double>::iterator it = migfracs.begin(); it != migfracs.end(); it++){
			int sourceindex = tree->g[tips[it->first]].index;
			Graph::vertex_descriptor sourcev = i2v[sourceindex];

			if (tree->is_legal_migration(sourcev, tmp)){
				//cout << "legal!\n";
				Graph::edge_descriptor e = tree->add_mig_edge(sourcev, tmp);
				tree->g[e].weight = it->second;
			}
			else {
				//cout << "nope\n";
				stillok = false;
			}
		}
		if (!stillok){
			if (i == first) first++;
			tree->copy(tree_bk);
			inorder = tree->get_inorder_traversal(current_npops);
			tips = tree->get_tips(tree->root);
			continue;
		}

		current_npops++;
		if (params->f2) set_branches_ls_f2();
		else set_branches_ls();
		double llk = llik();
		if (i == first){
			max_index = i;
			max_llik = llk;
		}
		else if (llk > max_llik){
			max_index = i;
			max_llik = llk;
		}
		tree->copy(tree_bk);
		current_npops--;
		inorder = tree->get_inorder_traversal(current_npops);
		tips = tree->get_tips(tree->root);
	}
	if (max_llik <= -DBL_MAX){
		cerr <<"ERROR: numerical problem in likelihood [ln(lk) = "<< max_llik << "]\n";
		cerr <<"ERROR: Please report this bug.\n";
		exit(1);
		/*
		cout << "RESCALING\n"; cout.flush();
		params->smooth_scale = params->smooth_scale *2;
		for (int i = 0; i < inorder.size(); i++){

				if ( params->set_root && params->root ==tree->g[ inorder[i] ].name) { //
					if (i == first) first++;
					continue;
				}
				if ( params->set_root && tree->g[ inorder[i] ].is_root) { //
					if (i == first) first++;
					continue;
				}
				Graph::vertex_descriptor tmp = tree->add_tip(inorder[i], toadd);
				map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();

				bool stillok = true;
				for (map<string, double>::iterator it = migfracs.begin(); it != migfracs.end(); it++){
					int sourceindex = tree->g[tips[it->first]].index;
					Graph::vertex_descriptor sourcev = i2v[sourceindex];

					if (tree->is_legal_migration(sourcev, tmp)){
						//cout << "legal!\n";
						Graph::edge_descriptor e = tree->add_mig_edge(sourcev, tmp);
						tree->g[e].weight = it->second;
					}
					else {
						//cout << "nope\n";
						stillok = false;
					}
				}
				if (!stillok){
					if (i == first) first++;
					tree->copy(tree_bk);
					inorder = tree->get_inorder_traversal(current_npops);
					tips = tree->get_tips(tree->root);
					continue;
				}

				current_npops++;
				if (params->f2) set_branches_ls_f2();
				else set_branches_ls();
				double llk = llik();
				if (i == first){
					max_index = i;
					max_llik = llk;
				}
				else if (llk > max_llik){
					max_index = i;
					max_llik = llk;
				}
				tree->copy(tree_bk);
				current_npops--;
				inorder = tree->get_inorder_traversal(current_npops);
				tips = tree->get_tips(tree->root);
		}
		*/
	}
	//cout << tree->g[inorder[max_index]].index << "\n";
	Graph::vertex_descriptor tmp = tree->add_tip(inorder[max_index], toadd);

	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	for (map<string, double>::iterator it = migfracs.begin(); it != migfracs.end(); it++){
		int sourceindex = tree->g[tips[it->first]].index;
		Graph::vertex_descriptor sourcev = i2v[sourceindex];
		Graph::edge_descriptor e = tree->add_mig_edge(sourcev, tmp);
		tree->g[e].weight = it->second;
	}
	current_npops++;

	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	current_llik = max_llik;

}


void GraphState2::add_pop(string name, string popname){

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma_cor);

	string toadd = popname;
	map<string, Graph::vertex_descriptor> p2v = tree->get_tips(tree->root);
	Graph::vertex_descriptor t = tree->add_tip(p2v[name], toadd);
	current_npops++;
	compute_sigma();

	set_sigmacor_from_sigma();
	//print_sigma();
	current_llik = llik();
}


double GraphState2::llik( bool w){
	if (!w) return llik_normal();
	else return llik_wishart();
}

double GraphState2::llik_wishart(){
	// density of the wishart distribution with covariance matrix sigma, n = number of snps-1, p = number of populations
	// 		scatter matrix has been stored in countdata->scatter, the ln(determinant) is in countdata->scatter_det
	// 		and the ln of the relevant multiariate gamma is in scatter_gamma
	//
	// density is ( [n-p-1]/2 * scatter_det - [1/2] trace [sigma^-1 * scatter] - [np/2] ln(2) - [n/2] ln(det(sigma)) - scatter_gamma
	gsl_matrix_free(scatter);
	scatter = gsl_matrix_alloc(current_npops, current_npops);
	double scale = (double) countdata->ne/ (double) countdata->nsnp;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			//cout << allpopnames[i] << " "<< allpopnames[j] << " "<< countdata->get_scatter(allpopnames[i], allpopnames[j]) << "\n";
			gsl_matrix_set(scatter, i, j, countdata->get_scatter( allpopnames[i], allpopnames[j] ) *scale);
			//cout << countdata->get_scatter( allpopnames[i], allpopnames[j]) << " "<< scale << " "<< gsl_matrix_get(scatter, i, j) << "\n";
		}
	}
/*
	ofstream tmpout("scatter_tmp");

	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			tmpout << gsl_matrix_get(scatter, i, j) << " ";
		}
		tmpout << "\n";
	}
*/
	double toreturn = 0;
	int s;
	int p = current_npops;
	int n = countdata->ne;

	//set scatter gamma
	scatter_gamma = ( (double) (p-1) * ( (double)  (p-1)-1.0) /4.0) * log (M_PI);
	for (int i = 1; i <= p-1; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);

	gsl_matrix * U = gsl_matrix_alloc(p-1, p);
	gsl_matrix *scatter_prime = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix * W_prime = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix * W_inv = gsl_matrix_alloc(p-1, p-1);

	// 1. Take SVD of scatter
	gsl_matrix * A = gsl_matrix_alloc(p,p);
	gsl_matrix * VT = gsl_matrix_alloc(p,p);
	gsl_vector * S = gsl_vector_alloc(p);
	gsl_vector * work = gsl_vector_alloc(p);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);


	// Now copy the first npop-1 eigenvectors to U

	for (int i = 0; i < p-1; i++){
		for(int j = 0; j < p; j++){
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}
/*
	cout << "\n";
	for (int i = 0; i < current_npops-1; i++){
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(U, i, j) << " ";
		}
		cout << "\n";
	}

	cout << "\n";
*/
	// 2. transform scatter into m-1 space
	// S' = U S U^T

	gsl_matrix * US = gsl_matrix_alloc(p-1, p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, scatter, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U, 0.0, scatter_prime);

	/*
	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(scatter_prime, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";
*/

	// 3. Same thing on predicted covariance matrix
	gsl_matrix * UW = gsl_matrix_alloc(p-1, p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, sigma_cor, 0.0, UW);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, UW, U, 0.0, W_prime);


	// 4. Do LU decomposition and get determinant of scatter_prime
	gsl_matrix_free(A);
	A = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix_memcpy( A, scatter_prime );
	gsl_permutation * perm = gsl_permutation_alloc(p-1);
	gsl_linalg_LU_decomp( A, perm, &s );
	scatter_det = gsl_linalg_LU_lndet( A );

	// 5. Do LU decomposition, get inverse of W_prime
	gsl_matrix_free(A);
	gsl_permutation_free(perm);
	A = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix_memcpy( A, W_prime );
	perm = gsl_permutation_alloc(p-1);
	gsl_linalg_LU_decomp( A, perm, &s );
	gsl_linalg_LU_invert( A, perm, W_inv );
	double w_det = gsl_linalg_LU_lndet( A );
/*
	cout << "scatter_prime_det "<< scatter_det << " w_det "<< w_det <<"\n";
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma_cor, i, j) << " ";
		}
		cout << "\n";
	}

	cout<<"\n";
*/
/*
	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(W_prime, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";
	*/
/*
	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(W_inv, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";

	*/
	//multiply inverse of cov by scatter, get trace
	gsl_matrix * ViU = gsl_matrix_alloc(p-1, p-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W_inv, scatter_prime, 0.0, ViU);

	double trace = 0;
	for (int i = 0; i < p-1 ; i++) trace+= gsl_matrix_get(ViU, i, i);

	toreturn+= ( (double) n- (double) (p-1)-1.0 )/2.0 * scatter_det - trace/2.0;
	toreturn+= -( (double) n* (double) (p-1)/2.0) * log(2.0);
	toreturn += -((double) n/2.0)*w_det;
	toreturn+= -scatter_gamma;

	gsl_matrix_free( U );
	gsl_matrix_free(scatter_prime);
	gsl_matrix_free(W_prime);
	gsl_matrix_free(W_inv);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(ViU);
	gsl_matrix_free( US );
	gsl_matrix_free( UW );
	gsl_permutation_free(perm);

	//cout << "llik: "<< toreturn <<" "<< n << " \n";
	return toreturn;
}

/*
double GraphState2::llik_mvn(){
	double toreturn;
	int n = countdata->ncomp;
	gsl_matrix * D = gsl_matrix_alloc(countdata->nblock, countdata->ncomp);
	gsl_matrix_memcpy( D, countdata->cov_samp);
	gsl_matrix * V = gsl_matrix_alloc(countdata->ncomp, countdata->ncomp);
	gsl_matrix_memcpy( V, countdata->cov_cov);

	//cout << "here\n";
	gsl_vector * means = gsl_vector_alloc(n);
	int index = 0;
	for (int i = 0; i < countdata->npop; i++){
		for (int j = i; j < countdata->npop; j++){
			string s1 = countdata->id2pop[i];
			string s2 = countdata->id2pop[j];
			int index1 = popname2index[s1];
			int index2 = popname2index[s2];

			double m= gsl_matrix_get(sigma_cor, index1, index2);
			//cout << s1 <<" "<< s2 << " "<< index1 << " "<< index2 << "\n";
			gsl_vector_set(means, index, m);
			index++;
		}
	}


	gsl_matrix * S = gsl_matrix_alloc(countdata->ncomp, countdata->ncomp);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, D, D, 0.0, S);
	gsl_matrix * A = gsl_matrix_alloc(n,n);
	gsl_matrix * VT = gsl_matrix_alloc(n,n);
	gsl_vector * Sv = gsl_vector_alloc(n);
	gsl_vector * work = gsl_vector_alloc(n);
	gsl_matrix_memcpy( A, S );
	gsl_linalg_SV_decomp(A, VT, Sv, work);
	int is = 0;
	while ( is < n && gsl_vector_get(Sv, is)  > 1e-10) is++; //this is the number of eigenvectors to use
	//cout << is << "\n";
	//is = 5;
	gsl_matrix * U = gsl_matrix_alloc(is, n);
	for (int i = 0; i < is ; i++){
		for (int j = 0; j < n; j++){
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}
	gsl_matrix * Dnew = gsl_matrix_alloc(countdata->nblock, is);
	gsl_matrix *Dnew_t = gsl_matrix_alloc(is, countdata->nblock);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, D, 0.0, Dnew_t);
	gsl_matrix_transpose_memcpy(Dnew, Dnew_t);

	gsl_vector * mean_new = gsl_vector_alloc(is);
	for(int i = 0; i < is; i++){
		double tmp = 0;
		for(int j = 0; j < n; j++){
			tmp += gsl_matrix_get(U, i, j)* gsl_vector_get(means, j);
		}
		gsl_vector_set(mean_new, i, tmp);
	}
	gsl_matrix * Vnew = gsl_matrix_alloc(is, is);
	gsl_matrix * UV = gsl_matrix_alloc(is, n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, V, 0.0, UV);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, UV, U, 0.0, Vnew);

	ofstream tmpm("means_new");
	for (int i = 0; i < is ; i++){
		tmpm << gsl_vector_get(mean_new, i) << "\n";
		//cout << gsl_vector_get(mean_new, i) << "\n";
	}
	ofstream tmpd("data_new");
	for (int i = 0; i < countdata->nblock; i++){
		for (int j = 0; j< is; j++){
			tmpd << gsl_matrix_get(Dnew, i, j)<< " ";
		}
		tmpd << "\n";
	}
	ofstream tmpv("var_new");
	for (int i = 0; i < is; i++){
		for (int j = 0; j< is; j++){
			tmpv << gsl_matrix_get(Vnew, i, j)<< " ";
		}
		tmpv << "\n";
	}

	for (int i = 0; i < countdata->nblock; i++){
	//for (int i = 0; i < 1; i++){
		//gsl_vector * tmpvec  = gsl_vector_alloc(n);
		//for (int j = 0; j < n; j++) gsl_vector_set(tmpvec, j, gsl_matrix_get(D, i, j));
		//double d = dmvnorm(n, tmpvec, means, V);
		gsl_vector * tmpvec  = gsl_vector_alloc(is);
		for (int j = 0; j < is; j++) gsl_vector_set(tmpvec, j, gsl_matrix_get(Dnew, i, j));
		double d = dmvnorm(is, tmpvec, mean_new, Vnew);
		toreturn+= log(d);
		//cout << d << " "<< log(d)<< "\n";
		gsl_vector_free(tmpvec);
	}

	gsl_matrix_free(S);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_matrix_free(D);
	gsl_matrix_free(V);
	gsl_matrix_free(Dnew);
	gsl_matrix_free(Dnew_t);
	gsl_matrix_free(Vnew);
	gsl_matrix_free(UV);
	gsl_vector_free(means);
	gsl_vector_free(Sv);
	gsl_vector_free(work);
	//cout << toreturn << "\n";
	return toreturn;
}
*/
void GraphState2::process_scatter(){
	scatter_det = 0;
	scatter_gamma = 0;

	int n = countdata->nsnp-1;
	size_t pop = current_npops;
	gsl_matrix_free(scatter);
	scatter = gsl_matrix_alloc(pop, pop);
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			//cout << i <<  " "<< j << "\n";
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			gsl_matrix_set(scatter, i, j, countdata->get_scatter(p1, p2));
		}
	}

	int s;

	gsl_matrix * work = gsl_matrix_alloc(pop,pop);
	gsl_permutation * p = gsl_permutation_alloc(pop);



	//do LU decomposition and get determinant
	gsl_linalg_LU_decomp( work, p, &s );
	scatter_det = gsl_linalg_LU_lndet( work );
	//get the log sum of the gammas
	scatter_gamma = ( (double) current_npops * ( (double)  current_npops-1.0) /4.0) * log (M_PI);
	for (int i = 1; i <= current_npops; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);
	//cout << "scatter_gamma "<< scatter_gamma << "\n";
}


void GraphState2::print_sigma_cor(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < current_npops; i++) out << allpopnames.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < current_npops; i++){
		out << allpopnames.at(i);
		for(int j = 0; j < current_npops; j++)	 out << " "<< gsl_matrix_get(sigma_cor, i, j);
		out << "\n";
	}
}

pair<bool, Graph::edge_descriptor> GraphState2::add_mig(int index1, int index2){
	//pair<bool, Graph::vertex_descriptor> toreturn;
	Graph::edge_descriptor e;
	bool added = false;
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	if ( tree->is_legal_migration( i2v[index1], i2v[index2])){
		added = true;
		e = tree->add_mig_edge(i2v[index1], i2v[index2]);

		//initialize_migupdate();
		//optimize_weight_quick(e);
		optimize_weight(e);
		//cout << tree->g[e].weight <<"\n";

	}
	current_llik = llik();
	return(make_pair(added, e));
}


Graph::edge_descriptor GraphState2::add_mig_noopt(int index1, int index2){
	Graph::edge_descriptor e;
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	if ( tree->is_legal_migration( i2v[index1], i2v[index2])){
		e = tree->add_mig_edge(i2v[index1], i2v[index2]);
		Graph::vertex_descriptor v = source(e, tree->g);
		tree->g[e].weight =0;
	}
	else{
		cerr << "ERROR: not a legal migration between index " << index1 << " and "<< index2 << "\n";
		exit(1);
	}
	current_llik = llik();
	return(e);
}

void GraphState2::rearrange(int index1, int index2){
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	if (i2v.find(index1) == i2v.end()){
		cout << "ERROR: no such index "<< index1 << "\n";
		exit(1);
	}
	if (i2v.find(index2) == i2v.end()){
		cout << "ERROR: no such index "<< index1 << "\n";
		exit(1);
	}
	tree->global_rearrange(i2v[index1], i2v[index2]);
	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	current_llik = llik();
}

void GraphState2::add_mig(){

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	int maxst;
	int maxsp;
	int size = inorder.size();
	tree_bk->copy(tree);
	double maxllk = current_llik;

	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){

			inorder = tree->get_inorder_traversal(current_npops);
			if (tree->is_legal_migration(inorder[i], inorder[j])){
				Graph::edge_descriptor e = tree->add_mig_edge( inorder[i], inorder[j]);
				optimize_weights();

				if (current_llik > maxllk){
					maxst = i;
					maxsp = j;
					maxllk = current_llik;
				}

			}
			tree->copy(tree_bk);
		}
	}
	cout << "here2\n"; cout.flush();
	cout << maxst <<  " "<< maxsp << " here\n";
	inorder = tree->get_inorder_traversal(current_npops);
	cout << tree->g[ inorder[maxst]].index <<  " "<< tree->g[ inorder[maxsp]].index<< " here2\n";
	tree->add_mig_edge(inorder[maxst], inorder[maxsp]);
	optimize_weights();
}

pair<string, string> GraphState2::get_max_resid(){
	string pop1, pop2;
	double max = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			if( i == j) continue;
			double cov = countdata->get_cov( allpopnames[i], allpopnames[j] );
			double fitted = gsl_matrix_get(sigma_cor, i, j);
			double diff = cov-fitted;
			if (diff > max){
				pop1 = allpopnames[i];
				pop2 = allpopnames[j];
				max = diff;
			}
		}
	}
	return make_pair(pop1, pop2);
}



pair<bool, pair<int, int> > GraphState2::add_mig_targeted_f2(){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	// tmpfitted will have a backup of the current covariance matrix
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);

	//tmpfitted2 will hold the covariance matrix at the tree with the max likelihood
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);

	// resid holds the current residuals
	gsl_matrix *resid = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(resid);


	double llik_bk = current_llik;
	double negsum_bk = negsum;
	double max_negsum = negsum;

	pair<bool, pair<int, int> > toreturn;
	toreturn.first = false;

	//1. set the residual matrix

	for ( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double cov = countdata->get_cov( allpopnames[i], allpopnames[j] );
			double fitted = gsl_matrix_get(sigma_cor, i, j);
			double diff = cov-fitted;
			gsl_matrix_set(resid, i, j, diff);
		}
	}

	double max;
	//2. Get the minimum nresid residuals
	set<pair<string, string> > minresids;
	cout << "Targeting migration to vicinity of:\n";
	for (int i = 0; i < params->nresid; i++){
		size_t min_i, min_j;
		gsl_matrix_min_index(resid, &min_i, &min_j);
		string p1 = allpopnames[min_i];
		string p2 = allpopnames[min_j];
		if (i == 0) max = gsl_matrix_get(resid, min_i, min_j);
		minresids.insert(make_pair(p1, p2));
		gsl_matrix_set(resid, min_i, min_j, 0);
		gsl_matrix_set(resid, min_j, min_i, 0);
		cout << p1 << " "<< p2 << "\n";
	}
	cout << "\n";
	//3. get populations in the neighborhood of those to target

	set<pair<int, int> > tested; //hold the pairs of vertices that have been tested
	set<pair<string, string> > pops2test;


	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();

	for (set<pair<string, string> >::iterator it = minresids.begin(); it != minresids.end(); it++){
		string pop1 = it->first;
		string pop2 = it->second;
		pair<set<int>, set<int> > n1 = get_neighborhood( tree->g[ tips[pop1] ].index, 4);
		pair<set<int>, set<int> > n3 = get_neighborhood( tree->g[ tips[pop2] ].index, 4);

		for (set<int>::iterator it = n1.first.begin();it != n1.first.end(); it++){
			Graph::vertex_descriptor v1 = i2v[*it];
			if ( ! tree->g[v1].is_tip) continue;
			for (set<int>::iterator it2 = n3.first.begin(); it2 != n3.first.end(); it2++){
				Graph::vertex_descriptor v2 = i2v[*it2];
				if ( ! tree->g[v2].is_tip) continue;
				string p1 = tree->g[v1].name;
				string p2 = tree->g[v2].name;
				double cov = countdata->get_cov(p1, p2);
				double fitted = gsl_matrix_get(sigma_cor, popname2index[p1], popname2index[p2]);
				double diff = cov-fitted;
				if (diff < max * 0.3){
					pops2test.insert(make_pair(p1, p2));
					cout << p1 << " "<< p2 << "\n";
				}
			}
		}
	}



	//3. try migration to all the pairwise combinations

	double max_llik = current_llik;
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);
	//tree_bk holds a backup of the current tree; bk_2 the best tree
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	pair<int, int> best_edge;
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj();
	for (set<pair<string, string> >::iterator pit = pops2test.begin(); pit != pops2test.end(); pit++){
		cout << "Trying "<< pit->first<< " "<< pit->second << "\n"; cout.flush();
		set<Graph::vertex_descriptor> p1_s = tree->get_path_to_root(p2node[pit->first]);
		set<Graph::vertex_descriptor> p2_s = tree->get_path_to_root(p2node[pit->second]);
		for (set<Graph::vertex_descriptor>::iterator it = p1_s.begin(); it != p1_s.end(); it++){
			for (set<Graph::vertex_descriptor>::iterator it2 = p2_s.begin(); it2 != p2_s.end(); it2++){
				if (!try_mig(*it, *it2, tmpfitted)) continue;

				cout << tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
				pair<int, int> totest = make_pair( tree->g[*it].index, tree->g[*it2].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it, *it2)){


					Graph::edge_descriptor e = tree->add_mig_edge( *it, *it2);
					//set_branches_ls_f2();
					//optimize_weight(e);
					initialize_migupdate();
					optimize_weight_quick(e);

					cout << "1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";

					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[*it].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[*it].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[*it2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it2, *it)){

					Graph::edge_descriptor e = tree->add_mig_edge( *it2, *it);
					//set_branches_ls_f2();
					//optimize_weight(e);
					initialize_migupdate();
					optimize_weight_quick(e);

					cout << "2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[*it2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);

					tested.insert(make_pair( tree->g[*it2].index, tree->g[*it].index));
				}

				Graph::vertex_descriptor p1 = tree->get_parent_node(*it).first;
				Graph::vertex_descriptor p2 = tree->get_parent_node(*it2).first;

				totest = make_pair( tree->g[p1].index, tree->g[*it2].index);
				if (tested.find(totest) == tested.end() && tree->is_legal_migration(p1, *it2)){
					//cout << "trying "<< tree->get_newick_format(p1)<< " "<< tree->get_newick_format(*it2)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p1, *it2);
					//set_branches_ls_f2();
					//optimize_weight(e);
					initialize_migupdate();
					optimize_weight_quick(e);

					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[p1].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p1].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[p2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(p2, *it)){
					//cout << "trying "<< tree->get_newick_format(p2)<< " "<< tree->get_newick_format(*it)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p2, *it);

					//set_branches_ls_f2();
					//optimize_weight(e);
					initialize_migupdate();

					optimize_weight_quick(e);

					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[p2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p2].index, tree->g[*it].index));
				}
			}
		}
	}
	if (toreturn.first == true)	{
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max_llik;
		negsum = max_negsum;
		toreturn.second.first = best_edge.first;
		toreturn.second.second = best_edge.second;

	}
	else {
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = llik_bk;
		negsum = negsum_bk;
		tree->copy(tree_bk);
	}

	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	gsl_matrix_free(resid);
	return toreturn;

}

pair<bool, pair<int, int> > GraphState2::add_mig_targeted(){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	// tmpfitted will have a backup of the current covariance matrix
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);

	//tmpfitted2 will hold the covariance matrix at the tree with the max likelihood
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);

	// resid holds the current residuals
	gsl_matrix *resid = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(resid);


	double llik_bk = current_llik;
	double negsum_bk = negsum;
	double max_negsum = negsum;

	pair<bool, pair<int, int> > toreturn;
	toreturn.first = false;

	//1. set the residual matrix

	for ( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double cov = countdata->get_cov( allpopnames[i], allpopnames[j] );
			double fitted = gsl_matrix_get(sigma_cor, i, j);
			double diff = cov-fitted;
			gsl_matrix_set(resid, i, j, diff);
		}
	}

	double max;
	//2. Get the minimum nresid residuals
	set<pair<string, string> > minresids;
	cout << "Targeting migration to vicinity of:\n";
	for (int i = 0; i < params->nresid; i++){
		size_t min_i, min_j;
		gsl_matrix_max_index(resid, &min_i, &min_j);
		string p1 = allpopnames[min_i];
		string p2 = allpopnames[min_j];
		if (i == 0) max = gsl_matrix_get(resid, min_i, min_j);
		minresids.insert(make_pair(p1, p2));
		gsl_matrix_set(resid, min_i, min_j, 0);
		gsl_matrix_set(resid, min_j, min_i, 0);
		cout << p1 << " "<< p2 << "\n";
	}
	cout << "\n";
	//3. get populations in the neighborhood of those to target

	set<pair<int, int> > tested; //hold the pairs of vertices that have been tested
	set<pair<string, string> > pops2test;


	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();

	for (set<pair<string, string> >::iterator it = minresids.begin(); it != minresids.end(); it++){
		string pop1 = it->first;
		string pop2 = it->second;
		pair<set<int>, set<int> > n1 = get_neighborhood( tree->g[ tips[pop1] ].index, 3);
		pair<set<int>, set<int> > n3 = get_neighborhood( tree->g[ tips[pop2] ].index, 3);

		for (set<int>::iterator it = n1.first.begin();it != n1.first.end(); it++){
			Graph::vertex_descriptor v1 = i2v[*it];
			if ( ! tree->g[v1].is_tip) continue;
			for (set<int>::iterator it2 = n3.first.begin(); it2 != n3.first.end(); it2++){
				Graph::vertex_descriptor v2 = i2v[*it2];
				if ( ! tree->g[v2].is_tip) continue;
				string p1 = tree->g[v1].name;
				string p2 = tree->g[v2].name;
				double cov = countdata->get_cov(p1, p2);
				double fitted = gsl_matrix_get(sigma_cor, popname2index[p1], popname2index[p2]);
				double diff = cov-fitted;
				if (diff > max * 0.3){
					pops2test.insert(make_pair(p1, p2));
					cout << p1 << " "<< p2 << "\n";
				}
			}
		}
	}



	//3. try migration to all the pairwise combinations

	double max_llik = current_llik;
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);
	//tree_bk holds a backup of the current tree; bk_2 the best tree
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	pair<int, int> best_edge;
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj();
	for (set<pair<string, string> >::iterator pit = pops2test.begin(); pit != pops2test.end(); pit++){
		cout << "Trying "<< pit->first<< " "<< pit->second << "\n"; cout.flush();
		set<Graph::vertex_descriptor> p1_s = tree->get_path_to_root(p2node[pit->first]);
		set<Graph::vertex_descriptor> p2_s = tree->get_path_to_root(p2node[pit->second]);
		for (set<Graph::vertex_descriptor>::iterator it = p1_s.begin(); it != p1_s.end(); it++){
			for (set<Graph::vertex_descriptor>::iterator it2 = p2_s.begin(); it2 != p2_s.end(); it2++){
				if (!try_mig(*it, *it2, tmpfitted)) continue;

				cout << tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
				pair<int, int> totest = make_pair( tree->g[*it].index, tree->g[*it2].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it, *it2)){


					Graph::edge_descriptor e = tree->add_mig_edge( *it, *it2);
					//set_branches_ls_f2();
					optimize_weight(e);


					cout << "1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n"; cout.flush();

					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[*it].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[*it].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[*it2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it2, *it)){

					Graph::edge_descriptor e = tree->add_mig_edge( *it2, *it);
					//set_branches_ls_f2();
					optimize_weight(e);


					cout << "2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n"; cout.flush();
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[*it2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);

					tested.insert(make_pair( tree->g[*it2].index, tree->g[*it].index));
				}

				Graph::vertex_descriptor p1 = tree->get_parent_node(*it).first;
				Graph::vertex_descriptor p2 = tree->get_parent_node(*it2).first;

				totest = make_pair( tree->g[p1].index, tree->g[*it2].index);
				if (tested.find(totest) == tested.end() && tree->is_legal_migration(p1, *it2)){
					//cout << "trying "<< tree->get_newick_format(p1)<< " "<< tree->get_newick_format(*it2)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p1, *it2);
					//set_branches_ls_f2();
					optimize_weight(e);


					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n"; cout.flush();
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[p1].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p1].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[p2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(p2, *it)){
					//cout << "trying "<< tree->get_newick_format(p2)<< " "<< tree->get_newick_format(*it)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p2, *it);

					//set_branches_ls_f2();
					optimize_weight(e);


					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n"; cout.flush();
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						max_negsum = negsum;
						best_edge.first = tree->g[p2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p2].index, tree->g[*it].index));
				}
			}
		}
	}
	if (toreturn.first == true)	{
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max_llik;
		negsum = max_negsum;
		toreturn.second.first = best_edge.first;
		toreturn.second.second = best_edge.second;

	}
	else {
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = llik_bk;
		negsum = negsum_bk;
		tree->copy(tree_bk);
	}

	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	gsl_matrix_free(resid);
	return toreturn;

}


pair< pair<bool, bool>, pair<double, pair<int, int> > > GraphState2::add_mig_targeted(string p1, string p2){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	//2. Get the paths to the root for both
	pair< pair<bool, bool>, pair<double, pair<int, int> > > toreturn;
	toreturn.first.first = false;
	toreturn.first.second = false;
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);

	vector<Graph::vertex_descriptor> p1_s = tree->get_path_to_root_vec(p2node[p1]);
	vector<Graph::vertex_descriptor> p2_s = tree->get_path_to_root_vec(p2node[p2]);

	//3. try migration to all the pairwise combinations
	double max_llik = current_llik;
	tree_bk->copy(tree);
	double migw = 0;
	pair<int, int> best_edge;
	for(int i = 0; i < p1_s.size(); i++){
		for (int j = 0; j < p2_s.size(); j++){

			if ( tree->is_legal_migration(p1_s[i], p2_s[j])){
				Graph::edge_descriptor e = tree->add_mig_edge( p1_s[i] , p2_s[j]);

				optimize_weights();

				if (current_llik > max_llik){
					tree_bk->copy(tree);
					max_llik = current_llik;
					best_edge.first = i;
					best_edge.second = j;
					toreturn.first.first = true;
					migw = tree->g[e].weight;
				}

				tree->remove_mig_edge(e);
			}
			if ( tree->is_legal_migration(p2_s[j], p1_s[i])){
				Graph::edge_descriptor e = tree->add_mig_edge( p2_s[j], p1_s[i]);

				optimize_weights();
				if (current_llik > max_llik){
					tree_bk->copy(tree);
					max_llik = current_llik;
					best_edge.first = i;
					best_edge.second = j;
					toreturn.first.first = true;
					toreturn.first.second = true;
					migw = tree->g[e].weight;
				}
				tree->remove_mig_edge(e);
			}
		}
	}
	if (toreturn.first.first == true)	{

		tree->copy(tree_bk);
		set_branches_ls_wmig();
		toreturn.second.first = migw;
		toreturn.second.second = best_edge;

	}
	return toreturn;

}

string GraphState2::get_trimmed_newick() {
	if (params->givenmat) {
		cerr << "ERROR: Cannot use options -givenmat and -noss!\n";
		exit(1);
	}

	map<string, double> trim;
	string toreturn;

	for (map<string, int>::iterator it = countdata->pop2id.begin(); it!= countdata->pop2id.end(); it++) {
		int id = it->second;
		string pop = it->first;
		double meanhzy = countdata->mean_hzy.find(id)->second;
		double mean_n = countdata->mean_ninds.find(id)->second;
		double t = meanhzy / (4.0 * mean_n);
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	toreturn = tree->get_newick_format(&trim);
	return toreturn;
}


void GraphState2::print_trimmed(string stem) {
	if (params->givenmat) {
		cerr << "ERROR: Cannot use options -givenmat and -noss!\n";
		exit(1);
	}

	map<string, double> trim;
	string toreturn;

	for ( map<string, int>::iterator it = countdata->pop2id.begin(); it!= countdata->pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = countdata->mean_hzy.find(id)->second;
		double mean_n = countdata->mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	tree->print(stem, &trim);
}

pair< set<int>, set<int> > GraphState2::get_neighborhood(int index){
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	pair<set<int>, set<int> > toreturn;
	toreturn.first.insert(index);
	Graph::vertex_descriptor v = i2v[index];
	//cout << params->m_neigh << " neigh\n";
	get_neighborhood(v, params->m_neigh, &toreturn);
	return toreturn;
}



pair< set<int>, set<int> > GraphState2::get_neighborhood(int index, int number){
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	pair<set<int>, set<int> > toreturn;
	toreturn.first.insert(index);
	Graph::vertex_descriptor v = i2v[index];
	//cout << params->m_neigh << " neigh\n";
	get_neighborhood(v, number, &toreturn);
	return toreturn;
}


void GraphState2::get_neighborhood(Graph::vertex_descriptor v, int dist, pair<set<int>, set<int> >* alln){
	Graph::vertex_descriptor p, c1, c2;
	if (!tree->g[v].is_root){
		p = tree->get_parent_node_wmig(v).first;
		while (tree->g[p].is_mig){
			//cout << tree->g[p].index << "\n";
			alln->second.insert(tree->g[p].index);
			p = tree->get_parent_node_wmig(p).first;
		}
		int pindex = tree->g[p].index;
		//cout << pindex << "\n";
		set<Graph::edge_descriptor> tmpset = tree->get_in_mig_edges(p);
		for (set<Graph::edge_descriptor>::iterator it = tmpset.begin(); it != tmpset.end(); it++) alln->second.insert(tree->g[source(*it, tree->g)].index);

		if (alln->first.find(pindex) == alln->first.end()){
			alln->first.insert(pindex);
			if (dist > 1) get_neighborhood(p, dist-1, alln);
		}
	}
	if (!tree->g[v].is_tip){
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch = tree->get_child_nodes_wmig(v);
		c1 = ch.first;
		c2 = ch.second;
		while (tree->g[c1].is_mig){
			//cout << "c1\n"; cout.flush();
			alln->second.insert(tree->g[c1].index);
			c1 = tree->get_child_node_mig(c1);
		}

		while (tree->g[c2].is_mig){
			//cout << "c1\n"; cout.flush();
			alln->second.insert(tree->g[c2].index);
			c2 = tree->get_child_node_mig(c2);
		}

		int c1index = tree->g[c1].index;
		int c2index = tree->g[c2].index;
		set<Graph::edge_descriptor> tmpset = tree->get_in_mig_edges(c1);
		for (set<Graph::edge_descriptor>::iterator it = tmpset.begin(); it != tmpset.end(); it++) alln->second.insert(tree->g[source(*it, tree->g)].index);
		tmpset = tree->get_in_mig_edges(c2);
		for (set<Graph::edge_descriptor>::iterator it = tmpset.begin(); it != tmpset.end(); it++) alln->second.insert(tree->g[source(*it, tree->g)].index);

		if (alln->first.find(c1index) == alln->first.end()){
			alln->first.insert(c1index);
			if (dist> 1) get_neighborhood(c1, dist-1, alln);
		}
		if (alln->first.find(c2index) == alln->first.end()){
			alln->first.insert(c2index);
			if (dist> 1) get_neighborhood(c2, dist-1, alln);
		}
	}
}



bool GraphState2::try_mig(Graph::vertex_descriptor v1, Graph::vertex_descriptor v2, gsl_matrix * fitted){
	map<string, Graph::vertex_descriptor> tips1 = tree->get_tips(v1);
	map<string, Graph::vertex_descriptor> tips2 = tree->get_tips(v2);
	string test1 = tips1.begin()->first;
	if (tips2.find(test1) == tips2.end()){
		int totalcount = 0;
		int totalneg = 0;
		for (map<string, Graph::vertex_descriptor>::iterator it1 = tips1.begin(); it1 != tips1.end(); it1++){
			string p1 = it1->first;
			int index1 = popname2index[p1];
			for (map<string, Graph::vertex_descriptor>::iterator it2 = tips2.begin(); it2 != tips2.end(); it2++){
				string p2 = it2->first;
				int index2 = popname2index[p2];
				if (p1==p2) continue;
				double resid = countdata->get_cov(p1, p2) -  gsl_matrix_get(fitted, index1, index2);
				if (resid < 0 ) totalneg++;
				totalcount++;

			}
		}
		double frac = (double) totalneg/ (double) totalcount;
		if ( !params->f2 && frac > 0 ) return false;
		else if ( params->f2 && frac < 1 ) return false;

	}
	else{
		int totalcount = 0;
		int totalneg = 0;
		for (map<string, Graph::vertex_descriptor>::iterator it1 = tips1.begin(); it1 != tips1.end(); it1++){
			string p1 = it1->first;
			int index1 = popname2index[p1];

			for (map<string, Graph::vertex_descriptor>::iterator it2 = tips2.begin(); it2 != tips2.end(); it2++){
				string p2 = it2->first;
				int index2 = popname2index[p2];
				if (p1==p2) continue;
				double resid = countdata->get_cov(p1, p2) -  gsl_matrix_get(fitted, index1, index2);
				if ( resid < 0) totalneg++;
				totalcount++;
			}
		}
		double frac = (double) totalneg/ (double) totalcount;
		if (!params->f2 && frac < 1) return false;
		else if (params->f2 && frac > 0) return false;
	}

	return true;
}

void GraphState2::place_root(string r){
	tree->place_root(r);
	cout << "Set root above "<< r << "\n"; cout.flush();
	//set_branches_ls();
	//current_llik = llik();
	//cout << tree->get_newick_format()<< "\n"; cout.flush();
	//cout << "ln(lk) = "<< current_llik << "\n";
}

int GraphState2::get_nmig(){
	return tree->get_mig_edges().size();
}

void GraphState2::flip_mig(){
	vector<Graph::edge_descriptor> m = tree->get_mig_edges();

	int i = 0;
	while (i < m.size()){

		Graph::vertex_descriptor v = target(m[i], tree->g);
		set<Graph::edge_descriptor> inm = tree->get_in_mig_edges(v);
		double w = 0;
		double max = 0;
		Graph::edge_descriptor e;
		for (set<Graph::edge_descriptor>::iterator it3 = inm.begin(); it3!= inm.end(); it3++ ){
			w += tree->g[*it3].weight;
			if (tree->g[*it3].weight > max){
				max = tree->g[*it3].weight;
				e = *it3;
			}
		}
		double treew = 1-w;
		if (max > treew){
			// t holds the target
			// and s the source
			Graph::vertex_descriptor t = target(e, tree->g);
			Graph::vertex_descriptor s = source(e, tree->g);
			cout << "Flipping migration edge\n"; cout.flush();
			if ( tree->g[tree->get_parent_node(t).first].is_root) {
				i++;
				continue;
			}
			if ( tree->g[ tree->get_parent_node_wmig(t).first].is_mig ) {
				i++;
				continue;
			}

			//get remaining weight
			Graph::edge_descriptor inc;
			pair<Graph::in_edge_iterator, Graph::in_edge_iterator> ine = in_edges(t, tree->g);
			while (ine.first != ine.second){
				if (tree->g[*ine.first].is_mig == false) inc = *ine.first;
					ine.first++;
			}

			Graph::vertex_descriptor newm = source(inc, tree->g);

			set<Graph::edge_descriptor> inm2 = tree->get_in_mig_edges(newm);
			if (inm2.size() > 0 ) {
				i++;
				continue;
			}
			tree->g[newm].is_mig = true;
			tree->g[newm].mig_frac = 0.5;

			tree->g[s].is_mig = false;
			tree->g[s].mig_frac = 0;

			tree->g[inc].is_mig = true;
			tree->g[inc].weight = treew;
			tree->g[inc].len = 0;


			tree->g[e].is_mig = false;
			tree->g[e].len = 1;

			if (params->f2) set_branches_ls_f2();
			else set_branches_ls();

			optimize_weight(inc);
			current_llik = llik();
			//iterate_movemig( tree->g[newm].index);
			//iterate_local_hillclimb_wmig(tree->g[t].index);
			m = tree->get_mig_edges();
			i = 0;
		}
		else i++;
	}

}


void GraphState2::flip_mig(string pop){
	vector<Graph::edge_descriptor> m = tree->get_mig_edges();

	int i = 0;
	while (i < m.size()){

		Graph::vertex_descriptor v = target(m[i], tree->g);
		set<Graph::edge_descriptor> inm = tree->get_in_mig_edges(v);
		double w = 0;
		double max = 0;
		Graph::edge_descriptor e;
		for (set<Graph::edge_descriptor>::iterator it3 = inm.begin(); it3!= inm.end(); it3++ ){
			w += tree->g[*it3].weight;
			if (tree->g[*it3].weight > max){
				max = tree->g[*it3].weight;
				e = *it3;
			}
		}
		double treew = 1-w;
		//cout << tree->g[target(m[i], tree->g)].name << " "<< pop << "\n";
		if (tree->g[target(m[i], tree->g)].name == pop){
			// t holds the target
			// and s the source
			Graph::vertex_descriptor t = target(e, tree->g);
			Graph::vertex_descriptor s = source(e, tree->g);
			cout << "Flipping migration edge\n"; cout.flush();
			if ( tree->g[tree->get_parent_node(t).first].is_root) {
				i++;
				continue;
			}
			while ( tree->g[ tree->get_parent_node_wmig(t).first].is_mig ) {
				Graph::vertex_descriptor mn = tree->get_parent_node_wmig(t).first;
				pair<Graph::out_edge_iterator, Graph::out_edge_iterator> oute = out_edges(mn, tree->g);
				Graph::out_edge_iterator ei = oute.first;
				while (!tree->g[*ei].is_mig) ei++;
				Graph::vertex_descriptor pt = target(*ei, tree->g);
				Graph::vertex_descriptor p = tree->get_parent_node(t).first;
				tree->remove_mig_edge(*ei);
				add_mig(tree->g[p].index, tree->g[pt].index);
			}
			//cout << "Still good\n";
			//get remaining weight
			Graph::edge_descriptor inc;
			pair<Graph::in_edge_iterator, Graph::in_edge_iterator> ine = in_edges(t, tree->g);
			while (ine.first != ine.second){
				if (tree->g[*ine.first].is_mig == false) inc = *ine.first;
					ine.first++;
			}

			Graph::vertex_descriptor newm = source(inc, tree->g);

			set<Graph::edge_descriptor> inm2 = tree->get_in_mig_edges(newm);
			if (inm2.size() > 0 ) {
				i++;
				continue;
			}
			tree->g[newm].is_mig = true;
			tree->g[newm].mig_frac = 0.5;

			tree->g[s].is_mig = false;
			tree->g[s].mig_frac = 0;

			tree->g[inc].is_mig = true;
			tree->g[inc].weight = treew;
			tree->g[inc].len = 0;


			tree->g[e].is_mig = false;
			tree->g[e].len = 1;

			if (params->f2) set_branches_ls_f2();
			else set_branches_ls();

			optimize_weight(inc);
			current_llik = llik();
			//iterate_movemig( tree->g[newm].index);
			//iterate_local_hillclimb_wmig(tree->g[t].index);
			m = tree->get_mig_edges();
			i = m.size();
		}
		else i++;
	}

}

void GraphState2::trim_mig(){

	vector<Graph::edge_descriptor> m = tree->get_mig_edges();
	int i = 0;
	while (i < m.size()){
			//cout << i <<" in trim\n";
			double w = tree->g[ m[i] ].weight;
			Graph::vertex_descriptor v = target( m[i], tree->g);
			if (w < params->min_migw){
				tree->remove_mig_edge(m[i]);
				if (params->f2) set_branches_ls_f2();
				else set_branches_ls_wmig();
				m = tree->get_mig_edges();
				i = 0;
				continue;
			}
			i++;
	}
}



bool GraphState2::has_loop(){
	//cout << "in has loop\n"; cout.flush();
	bool toreturn = false;
	tree_bk3->copy(tree);
	set<pair<Graph::vertex_descriptor, Graph::vertex_descriptor> > vmig;
	vector<Graph::edge_descriptor> migedges = tree_bk3->get_mig_edges();
	for (vector<Graph::edge_descriptor>::iterator it = migedges.begin(); it != migedges.end(); it++){
		Graph::vertex_descriptor t1 = source(*it, tree_bk3->g);
		t1 = tree_bk3->get_child_node_mig(t1);
		Graph::vertex_descriptor t2 = target(*it, tree_bk3->g);
		vmig.insert(make_pair(t1, t2));
	}
	for (vector<Graph::edge_descriptor>::iterator it = migedges.begin(); it != migedges.end(); it++) tree_bk3->remove_mig_edge(*it);
	for (set<pair<Graph::vertex_descriptor, Graph::vertex_descriptor> >::iterator it = vmig.begin(); it != vmig.end(); it++){
		if (!tree_bk3->is_legal_migration(it->first, it->second) )	return true;
		tree_bk3->add_mig_edge(it->first, it->second);
	}

	return false;
}

int GraphState2::iterate_movemig(int index){
	int toreturn = 0;
	pair<bool, int> moving = movemig(index);
	while( moving.first){
		toreturn++;
		moving = movemig(moving.second);
	}
	return toreturn;

}

pair<bool, int> GraphState2::movemig( int index ){

	pair<bool, int> toreturn = make_pair(false, 0);
	double max = current_llik;
	double lik_bk = current_llik;
	double max_negsum = negsum;
	double negsum_bk = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);
	//tree->print("test");
	for (int i = 0; i < 7 ; i++){
		//cout << "moving "<< i <<"\n";
		map<int, Graph::vertex_descriptor> vindex = tree->index2vertex();
		Graph::vertex_descriptor s = vindex[index];
		Graph::vertex_descriptor p = tree->get_parent_node(s).first;
		Graph::vertex_descriptor ch = tree->get_child_node_mig(s);
		Graph::edge_descriptor e = tree->get_out_mig_edge(s);
		Graph::vertex_descriptor t = target(e, tree->g);
		Graph::vertex_descriptor tp = tree->get_parent_node(t).first;
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> t_ch = tree->get_child_nodes(t);
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch2 = tree->get_child_nodes(ch);
		//tree->print("test");
		if ( i == 0){
			if (tree->g[p].is_root){
				continue;
				//cout << "here\n";
				//pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				//if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				//else p = ch3.first;
				//ch3 = tree->get_child_nodes(p);
				//p = ch3.first;
				//cout.flush();
				//cout << tree->g[p].index << "\n";
			}

			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(p, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(p, t);
				if (params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test0");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max+params->epsilon){
					max = current_llik;
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;


				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 4){
			tree->remove_mig_edge(e);
			if (!tree->g[tp].is_root && tree->is_legal_migration(ch, tp)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch, tp);

					if (params->f2){
						initialize_migupdate();
						optimize_weight_quick(e2);
					}
					else optimize_weight(e2);
					//tree->print("test0");
					//cout << i << " "<< current_llik << " "<< max << "\n";
					if ( current_llik > max +params->epsilon){
						max = current_llik;

						max_negsum = negsum;
						tree->g[source(e2, tree->g)].index = index;
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						toreturn.first = true;
						tree->g[source(e2, tree->g)].index = index;
						toreturn.second = tree->g[source(e2, tree->g)].index;
					}

			}
			tree->copy(tree_bk);
		}
		else if (i == 5){
			tree->remove_mig_edge(e);
			if (!tree->g[t].is_tip && tree->is_legal_migration(ch, t_ch.first)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch, t_ch.first);

				if (params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else	optimize_weight(e2);
					//cout << "here\n"; cout.flush();
					//tree->print("test1");
					//cout << i << " "<< current_llik << " "<< max << "\n";
					if ( current_llik > max +params->epsilon){
						max = current_llik;
						//cout << "move5 "<< current_llik << " "<< llik() << "\n";
						max_negsum = negsum;
						tree->g[source(e2, tree->g)].index = index;
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						toreturn.first = true;
						tree->g[source(e2, tree->g)].index = index;
						toreturn.second = tree->g[source(e2, tree->g)].index;
					}

			}
			tree->copy(tree_bk);
		}

		else if (i == 6){
			tree->remove_mig_edge(e);
			if (!tree->g[t].is_tip && tree->is_legal_migration(ch, t_ch.second)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch, t_ch.second);
				if(params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
					//tree->print("test2");
					//cout << i << " "<< current_llik << " "<< max << "\n";
					if ( current_llik > max +params->epsilon){
						max = current_llik;
						//cout << "move6 "<< current_llik << " "<< llik() << "\n";
						max_negsum = negsum;
						tree->g[source(e2, tree->g)].index = index;
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						toreturn.first = true;
						tree->g[source(e2, tree->g)].index = index;
						toreturn.second = tree->g[source(e2, tree->g)].index;
					}

			}
			tree->copy(tree_bk);
		}
		else if ( i == 3){
			if (tree->g[p].is_root){
				pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				else p = ch3.first;
				ch3 = tree->get_child_nodes(p);
				p = ch3.second;
				//cout << tree->g[p].index << "\n";
			}
			else{
				pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				else p = ch3.first;
			}
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(p, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(p, t);

				if (params->f2){
				initialize_migupdate();
				optimize_weight_quick(e2);
				}
				else optimize_weight(e2);


				//tree->print("test3");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){
					max = current_llik;
					//cout << "move3 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);

					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 1){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.first, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.first, t);
				if (params->f2){
				initialize_migupdate();
				optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test1");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){
					max = current_llik;
					//cout << "move1 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 2){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.second, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.second, t);
				//initialize_migupdate();
				if (params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test2");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){

					max = current_llik;
					//cout << "move2 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}


	}
	if (toreturn.first){
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max;
		negsum = max_negsum;
		//tree->print("m0_2");
	}
	else{
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = lik_bk;
		negsum = negsum_bk;
		tree->copy(tree_bk);
	}
	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	return toreturn;
}

pair<bool, int> GraphState2::movemig_limit( int index ){

	pair<bool, int> toreturn = make_pair(false, 0);
	double max = current_llik;
	double lik_bk = current_llik;
	double max_negsum = negsum;
	double negsum_bk = negsum;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);
	//tree->print("test");
	for (int i = 0; i < 4 ; i++){
		//cout << "moving "<< i <<"\n";
		map<int, Graph::vertex_descriptor> vindex = tree->index2vertex();
		Graph::vertex_descriptor s = vindex[index];
		Graph::vertex_descriptor p = tree->get_parent_node(s).first;
		Graph::vertex_descriptor ch = tree->get_child_node_mig(s);
		Graph::edge_descriptor e = tree->get_out_mig_edge(s);
		Graph::vertex_descriptor t = target(e, tree->g);
		Graph::vertex_descriptor tp = tree->get_parent_node(t).first;
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> t_ch = tree->get_child_nodes(t);
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch2 = tree->get_child_nodes(ch);
		//tree->print("test");
		if ( i == 0){
			if (tree->g[p].is_root){
				continue;
				//cout << "here\n";
				//pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				//if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				//else p = ch3.first;
				//ch3 = tree->get_child_nodes(p);
				//p = ch3.first;
				//cout.flush();
				//cout << tree->g[p].index << "\n";
			}

			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(p, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(p, t);
				if (params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test0");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max+params->epsilon){
					max = current_llik;
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;


				}
			}
			tree->copy(tree_bk);
		}

		else if ( i == 3){
			if (tree->g[p].is_root){
				pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				else p = ch3.first;
				ch3 = tree->get_child_nodes(p);
				p = ch3.second;
				//cout << tree->g[p].index << "\n";
			}
			else{
				pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				else p = ch3.first;
			}
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(p, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(p, t);

				if (params->f2){
				initialize_migupdate();
				optimize_weight_quick(e2);
				}
				else optimize_weight(e2);


				//tree->print("test3");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){
					max = current_llik;
					//cout << "move3 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);

					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 1){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.first, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.first, t);
				if (params->f2){
				initialize_migupdate();
				optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test1");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){
					max = current_llik;
					//cout << "move1 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 2){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.second, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.second, t);
				//initialize_migupdate();
				if (params->f2){
					initialize_migupdate();
					optimize_weight_quick(e2);
				}
				else optimize_weight(e2);
				//tree->print("test2");
				//cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max +params->epsilon){

					max = current_llik;
					//cout << "move2 "<< current_llik << " "<< llik() << "\n";
					max_negsum = negsum;
					tree->g[source(e2, tree->g)].index = index;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}


	}
	if (toreturn.first){
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max;
		negsum = max_negsum;
		//tree->print("m0_2");
	}
	else{
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = lik_bk;
		negsum = negsum_bk;
		tree->copy(tree_bk);
	}
	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	return toreturn;
}


void GraphState2::compute_sigma(){

	gsl_matrix_set_zero(sigma);
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}

	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							double add = it1->first * it2->first * tree->g[*it3].len;
							gsl_matrix_set(sigma, i, j, gsl_matrix_get(sigma, i, j)+add);
						}
					}
				}
			}
		}
	}
}


void GraphState2::set_sigmacor_from_sigma(){

	gsl_matrix_set_zero(sigma_cor);
	double c1 = 1.0 / (double) current_npops;
	double c2 = c1*c1;
	double shared = 0;
	for(int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			shared += gsl_matrix_get(sigma, i, j);
		}
	}

	shared = shared * c2;
	//cout << "here "<< shared << " "<< current_npops << " \n";
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			//cout << i << " "<< j;
			double vij = gsl_matrix_get(sigma, i, j);
			double sum_i = 0;
			double sum_j = 0;
			for (int k = 0; k < current_npops; k++){
				sum_i += gsl_matrix_get(sigma, k, i);
				sum_j += gsl_matrix_get(sigma, j, k);
			}
			sum_i = sum_i * c1;
			sum_j = sum_j * c1;
			double w_ij = vij - sum_i - sum_j + shared;
			//cout << " "<< vij << " "<< sum_i << " "<< sum_j << " "<< w_ij << "\n";
			gsl_matrix_set(sigma_cor, i, j, w_ij);
			gsl_matrix_set(sigma_cor, j, i, w_ij);
		}
	}
}


pair<double, double> GraphState2::calculate_se(Graph::edge_descriptor e){
	vector<double> samps;
	double oldweight = tree->g[e].weight;
	//for (int i = 0; i < 1; i++){
	for (int i = 0; i < countdata->nblock; i++){
		countdata->set_cov_jackknife(i);
		double min, max, guess;
		guess = oldweight;
		min = params->minweight;
		max = params->maxweight;
		int nit = 0;

		if (params->f2){
			initialize_migupdate();
			golden_section_weight_noexp_quick(e, guess -0.02, guess, guess+0.02, 0.005, &nit);
		}
		else golden_section_weight_noexp(e, guess -0.02, guess, guess+0.02, 0.005, &nit);

		double tmpllk = llik();
		double tmpw = tree->g[e].weight;
		int tmpnit = nit;
		nit = 0;
		if (params->f2) golden_section_weight_noexp_quick(e, -1, 0, 1, 0.005, &nit);
		else golden_section_weight_noexp(e, -1, 0, 1, 0.005, &nit);
		double newllk = llik();
		//cout <<i << " "<< tree->g[e].weight << " " <<  llik() << " "<< tmpw << " "<< tmpllk << " "<< nit << " "<< tmpnit << "\n"; cout.flush();
		//for (int j = 0; j < 90; j++){
		//	double tj = (double) j/ 100.0;
		//	update_mig(e, tj);
		//	cout << i << " "<< j << " "<< llik() << "\n";
		//}
		if (tmpllk > newllk) samps.push_back(tmpw);
		else samps.push_back(tree->g[e].weight);

		//stringstream ss;
		//ss << "cov" << i << ".gz";
		//string tmp = ss.str();
		//countdata->print_cov(tmp);

	}
	double mean = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) mean += *it;
	mean = mean/ (double) countdata->nblock;
	double sum = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) {
		double toadd = (*it - mean);
		toadd = toadd*toadd;
		sum += toadd;
	}
	double se = ( (double) (countdata->nblock -1)/ (double) countdata->nblock) * sum;
	se = sqrt(se);
	tree->g[e].weight = oldweight;
	//cout << mean << " " << se << " se\n";
	return make_pair(mean, se);
}


pair<double, double> GraphState2::calculate_se_bootstrap(gsl_rng *r, Graph::edge_descriptor e){
	vector<double> samps;
	double oldweight = tree->g[e].weight;
	for (int i = 0; i < countdata->nblock; i++){
		countdata->set_cov_bootstrap(r);
		double min, max, guess;
		guess = oldweight;
		min = params->minweight;
		max = params->maxweight;
		int nit = 0;
		golden_section_weight_noexp(e, -1, guess, 1, 0.005, &nit);
		samps.push_back(tree->g[e].weight);
		cout<< i << " "<< tree->g[e].weight << "\n";

	}
	double mean = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) mean += *it;
	mean = mean/ (double) countdata->nblock;
	double sum = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) {
		double toadd = (*it - mean);
		toadd = toadd*toadd;
		sum += toadd;
	}

	double se = sqrt(sum) / ((double) countdata->nblock - 1.0);
	tree->g[e].weight = oldweight;
	//cout << se << " se\n";
	return make_pair(mean, se);
}


pair<double, double> GraphState2::calculate_se_fromsamp(Graph::edge_descriptor e){
	vector<double> samps;
	double oldweight = tree->g[e].weight;
	//for (int i = 0; i < 1; i++){
	for (int i = 0; i < countdata->nblock; i++){
		countdata->set_cov_fromsamp(i);
		double min, max, guess;
		guess = oldweight;
		min = params->minweight;
		max = params->maxweight;
		int nit = 0;
		golden_section_weight_noexp(e, -1, 0, 1, 0.0005, &nit);
		samps.push_back(tree->g[e].weight);
		cout << i << " "<< tree->g[e].weight << "\n";
	}
	double mean = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) mean += *it;
	mean = mean/ (double) countdata->nblock;
	double sum = 0;
	for (vector<double>::iterator it = samps.begin(); it!= samps.end(); it++) {
		double toadd = (*it - mean);
		toadd = toadd*toadd;
		sum += toadd;
	}
	double se = sqrt(sum) / (double) countdata->nblock;
	tree->g[e].weight = oldweight;
	return make_pair(mean, se);
}

void GraphState2::initialize_migupdate(){

	e2index.clear();
	e2frac.clear();
	e2tips.clear();
	popnames2index.clear();
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root
	int index = 0;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			e2index.insert(make_pair(*it, index));
			e2frac.insert(make_pair(*it, 1));
			index++;
		}
		else if (tree->g[*it].is_mig == true){
			Graph::vertex_descriptor t = target(*it, tree->g);
			Graph::vertex_descriptor s = source(*it, tree->g);
			Graph::vertex_descriptor s2 = tree->get_child_node_mig(s);
			map<string, Graph::vertex_descriptor> t1 = tree->get_tips(s2);
			map<string, Graph::vertex_descriptor> t2 = tree->get_tips(t);
			set<int> tmpnames;
			for(map<string, Graph::vertex_descriptor>::iterator it2 = t1.begin(); it2 != t1.end(); it2++) {
				//cout << tree->g[it2->second].name << " "<< tree->g[it2->second].index << "\n"; cout.flush();
				tmpnames.insert(tree->g[it2->second].index);
			}
			for(map<string, Graph::vertex_descriptor>::iterator it2 = t2.begin(); it2 != t2.end(); it2++) {
				//cout << tree->g[it2->second].name << " "<< tree->g[it2->second].index <<"\n"; cout.flush();
				tmpnames.insert(tree->g[it2->second].index);
			}
			e2tips.insert(make_pair(*it, tmpnames));
		}
	}
	//cout << "here1\n"; cout.flush();
	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		e2index[*it] = joint_index;
		e2frac.insert(make_pair(*it, 0.5));
	}
	index++;
	//cout << "here2\n"; cout.flush();
	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated
	gsl_matrix_free(X_current);
	X_current = gsl_matrix_alloc(n, p);
	gsl_vector_free(y_current);
	y_current = gsl_vector_alloc(n);
	gsl_matrix_set_zero(X_current);
	//cout << "here2.01\n"; cout.flush();
	popname2paths.clear();
	//cout << "here2.1\n"; cout.flush();
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		popname2paths.insert(make_pair(it->first, tmpset));
	}
	//cout << "here3\n"; cout.flush();

	map<string, int> pair2index;
	// Set up the matrices from the tree
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y_current, index, empirical_cov);

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = popname2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = popname2paths[p2];

			//variances of 1
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
					if ( tree->g[*it3].is_mig) continue;
					double frac = e2frac.find(*it3)->second;
					int addindex = e2index.find(*it3)->second;

					double add = it1->first * it1->first *frac;
					set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
					it4++;
					while (it4 != paths_1.end()){
						if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
						it4++;
					}
					gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
				}
			}

			//variances of 2
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
				for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
					if ( tree->g[*it3].is_mig) continue;
					double frac = e2frac.find(*it3)->second;
					int addindex = e2index.find(*it3)->second;
					double add = it1->first * it1->first *frac;

					set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
					it4++;
					while (it4 != paths_2.end()){
						if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
						it4++;
					}
						//if (p1 == "pop3" && p2 == "pop10") cout << index << " "<< addindex << " "<< add << "\n";
						//if (index  == 12 && addindex   == 0) cout << "addind for pop "<< p2 << " "<< add << "\n";
						//if (addindex ==1){
						//	cout << add << " "<< it1->first << " "<< frac << " "<< p2 << "\n";
						//}
					gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
				}
			}


			//covariances
			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = e2index.find(*it3)->second;
							double frac = e2frac.find(*it3)->second;
							double add = it1->first * it2->first *frac;
							//if (addindex ==1){
							//	cout << add << " "<< it1->first << " "<< frac << " "<< p1 << " "<< p2 << " X -2\n";
							//}
							add = -2*add;
							//if (index  == 12 && addindex   == 0) cout << "addind for overlap "<< add << "\n";
							gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
						}
					}
				}
			}
			map<string, int> tmp1;
			map<string, int> tmp2;
			tmp1.insert(make_pair(p1, index));
			tmp2.insert(make_pair(p2, index));
			popnames2index.insert(make_pair(p1, tmp2));
			popnames2index.insert(make_pair(p2, tmp1));
			index++;
		}
	}
	set_branches_ls_f2_precompute();
	current_llik = llik();
}

void GraphState2::print_X(){
	int p = X_current->size2;
	int n = (current_npops* (current_npops-1))/2;
	for (int i = 0 ; i < n ; i++){
		cout << gsl_vector_get(y_current, i) << " " << "y"<<" ";
		for (int j = 0; j < p ; j++){
			cout << gsl_matrix_get(X_current, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n\n";
	for (map<Graph::edge_descriptor, int>::iterator it = e2index.begin(); it != e2index.end(); it++){
		cout << tree->g[source(it->first, tree->g)].index << " "<< tree->g[target(it->first, tree->g)].index << " "<< e2frac[it->first] << " "<< it->second << "\n";
	}
	cout << "\n\n";
	for (map<string, set<pair<double, set<Graph::edge_descriptor> > > >::iterator it = popname2paths.begin(); it != popname2paths.end(); it++){
		cout << it->first << "\n";
		for (set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			cout << it2->first <<"\n";
			for (set<Graph::edge_descriptor>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++){
				cout << tree->g[source(*it3, tree->g)].index << "->" << tree->g[target(*it3, tree->g)].index << ";";
			}
			cout <<"\n";
		}
	}
}
void GraphState2::set_branches_ls_f2_precompute_old(){
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = X_current->size2; // p is the number of branches lengths to be estimated

	/*for (int i = 0 ; i < n ; i++){
		cout << gsl_vector_get(y_current, i) << " " << "y"<<" ";
		for (int j = 0; j < p ; j++){
			cout << gsl_matrix_get(X_current, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n\n";
	*/
	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);

	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	//cout << "Set up\n"; cout.flush();
	gsl_multifit_linear(X_current, y_current, c, cov, &chisq, work);
	//cout << "estimated\n"; cout.flush();
	//and put in the solutions in the graph
	negsum =0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = e2index[*it];
		double frac = e2frac[*it];
		double l = gsl_vector_get(c, i);
		tree->g[*it].len = l*frac;
		if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
		if (l < 0) 	negsum += -l;

	}

	//and the corrected covariance matrix
	int index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X_current, index, k) * gsl_vector_get(c, k);
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//free memory
	gsl_multifit_linear_free(work);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}

void GraphState2::set_branches_ls_f2_precompute(){
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = X_current->size2; // p is the number of branches lengths to be estimated

	/*for (int i = 0 ; i < n ; i++){
		cout << gsl_vector_get(y_current, i) << " " << "y"<<" ";
		for (int j = 0; j < p ; j++){
			cout << gsl_matrix_get(X_current, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n\n";
	*/
	//set up the workspace


	NNLS_SOLVER nnls(n, p);


	double * A = new double[n*p]; //A holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root.

	double * b  = new double[n];  // y contains the entries of the empirical covariance matrix
	double * x = new double[p];   // x will be estimated, contains the fitted branch lengths
	double rNorm;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < p; j++){
			A[j*n+i] = gsl_matrix_get(X_current, i, j);
		}
	}
	for (int i = 0; i < n ; i++) b[i] = gsl_vector_get(y_current, i);

	//fit NNLS
	bool converged = nnls.solve(A, p, b, x, rNorm);


	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = e2index[*it];
		double frac = e2frac[*it];
		double l = x[i];
		tree->g[*it].len = l*frac;
		//if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
	}

	//and the corrected covariance matrix
	int index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += A[ k * n + index ] * x[k];
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//free memory

	delete [] A;
	delete [] b;
	delete [] x;

}

void GraphState2::update_mig(Graph::edge_descriptor e, double w){
	//cout << "entering update\n"; cout.flush();
	tree->g[e].weight = w;
	set<int> tips = e2tips[e];
	vector<string> tipnames;
	set<string> tipnames_set;
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	for (set<int>::iterator it = tips.begin(); it != tips.end(); it++){
		//cout << "getting path\n";
		Graph::vertex_descriptor v = i2v[*it];
		string name = tree->g[v].name;
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(v);


		map<string, set<pair<double, set<Graph::edge_descriptor> > > >::iterator todelete = popname2paths.find(name);
		if (todelete == popname2paths.end()){
			cerr << "ERROR: cannot find " << tree->g[v].name << "\n";
			exit(1);
		}
		popname2paths.erase(todelete);
		popname2paths.insert(make_pair(name, tmpset));
		tipnames_set.insert(name);
		tipnames.push_back(name);
	}
	//cout << "got paths\n"; cout.flush();
	map<string, int> pair2index;
	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){

			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			if (tipnames_set.find(p1) == tipnames_set.end() && tipnames_set.find(p2) == tipnames_set.end()){
				index++;
				continue;
			}
			for (int k = 0; k < X_current->size2; k++) gsl_matrix_set(X_current, index, k, 0);
			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = popname2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = popname2paths[p2];

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
					if ( tree->g[*it3].is_mig) continue;
					double frac = e2frac.find(*it3)->second;
					int addindex = e2index.find(*it3)->second;
					double add = it1->first * it1->first *frac;
					set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
					it4++;
					while (it4 != paths_1.end()){
						if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
						it4++;
					}
					gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
				}
			}

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						double frac = e2frac.find(*it3)->second;
						int addindex = e2index.find(*it3)->second;
						double add = it1->first * it1->first *frac;
						set<pair<double, set<Graph::edge_descriptor> > >::iterator it4 = it1;
						it4++;
						while (it4 != paths_2.end()){
							if (it4->second.find(*it3) != it4->second.end())	add += 2 * it1->first * it4->first * frac;
							it4++;
						}
						//if (p1 == "pop3" && p2 == "pop10") cout << index << " "<< addindex << " "<< add << "\n";
						//if (index  == 12 && addindex   == 0) cout << "addind for pop "<< p2 << " "<< add << "\n";
						gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
					}
				}

			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if ( tree->g[*it3].is_mig) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = e2index.find(*it3)->second;
							double frac = e2frac.find(*it3)->second;
							double add = it1->first * it2->first *frac;
							add = -2*add;
							//if (index  == 12 && addindex   == 0) cout << "addind for overlap "<< add << "\n";
							gsl_matrix_set(X_current, index, addindex, gsl_matrix_get(X_current, index, addindex)+add);
						}
					}
				}
			}
			pair2index.insert(make_pair(p1+p2, index));
			pair2index.insert(make_pair(p2+p1, index));
			index++;
		}
	}
	//cout << "updated_mig\n"; cout.flush();
	//cout << "here4\n"; cout.flush();
	set_branches_ls_f2_precompute();
	//current_llik = llik();

}

void GraphState2::clean_negedge(){

	vector<Graph::edge_descriptor> m = tree->get_mig_edges();
	int i = 0;
	while (i < m.size()){
		//cout << i <<" in trim\n";
		Graph::edge_descriptor e = m[i];
		Graph::vertex_descriptor v = source( e, tree->g);
		Graph::in_edge_iterator it = in_edges(v, tree->g).first;
		Graph::edge_descriptor e2 = *it;
		Graph::vertex_descriptor c = tree->get_child_node_mig(v);
		double l = tree->get_parent_node(c).second;
		if ( !tree->g[source(e2, tree->g)].is_mig && tree->g[e2].len < 0 && l > 0){
			Graph::vertex_descriptor v2 = source(e2, tree->g);
			Graph::vertex_descriptor t = target(e, tree->g);
			tree->remove_mig_edge(e);
			Graph::edge_descriptor enew = tree->add_mig_edge(v2, t);
			initialize_migupdate();
			optimize_weight_quick(enew);
			m = tree->get_mig_edges();
			i = 0;
		}

		i++;
	}

}

void GraphState2::target_pop(){
	// do global search for a single population allowing one migration

	// find the relevant vertex
	map<string, Graph::vertex_descriptor> tips = tree->get_tips_nomig(tree->root);
	string p = params->target;
	if (tips.find(p) == tips.end()){
		cerr << "ERROR: no population "<< p << "  to target" << "\n";
		exit(1);
	}
	Graph::vertex_descriptor t = tips[p];
	Graph::vertex_descriptor tp = tree->get_parent_node(t).first;
	if (tree->g[tp].is_root){
		cerr << "ERROR: cannot target root\n";
		exit(1);
	}
	int ind = tree->g[t].index;
	int indp = tree->g[tp].index;
	// make backup of tree structure
	double max = current_llik;
	double lik_bk = current_llik;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);


	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);
	//do pairwise comparisons

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal_noroot(countdata->npop);
	vector<int> all_indices;
	for (vector<Graph::vertex_descriptor>::iterator it = inorder.begin(); it != inorder.end(); it++) all_indices.push_back(tree->g[*it].index);

	for (int i = 0; i < all_indices.size(); i++){
		//cout << i  << "\n";
		int tomove = all_indices[i];
		if (tomove == ind || tomove == indp) continue;
		//cout << "moving\n"; cout.flush();
		//tree->print("test1");
		//cout << "printed\n"; cout.flush();
		rearrange(ind, tomove);
		//cout << "moved\n"; cout.flush();
		for (int j = i+1; j < all_indices.size(); j++){

			int indmig = all_indices[j];
			//cout << i << " "<< j << " "<< tomove << " "<< indmig << "\n";
			if (indmig == ind) continue;
			//cout << "adding\n"; cout.flush();
			pair<bool, Graph::edge_descriptor> added = add_mig(indmig, ind);
			//cout << "added\n"; cout.flush();
			if (!added.first) continue;
			Graph::edge_descriptor e = added.second;
			//initialize_migupdate();
			//optimize_weight_quick(e);

    		Graph::vertex_descriptor p1 = source( e, tree->g);
     		p1 = tree->get_child_node_mig(p1);
     		Graph::vertex_descriptor p2 = tree->get_parent_node( target(e, tree->g)).first;
     		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch = tree->get_child_nodes(p2);
     		if (tree->g[ ch.first ].index  == ind) p2 = ch.second;
     		else p2 = ch.first;
     		cout << tree->get_newick_format(p1) << " ---> ";
     		cout << p;
     		cout << " <---- "<< tree->get_newick_format(p2) << "\n";
     		cout << current_llik << " "<< max << "\n";
			if (current_llik > max){
				tree_bk2->copy(tree);
				gsl_matrix_memcpy( tmpfitted2, sigma_cor);
				max = current_llik;
			}
			tree->remove_mig_edge(e);
		}
	}

	tree->copy(tree_bk2);
	gsl_matrix_memcpy( sigma_cor, tmpfitted2);
	current_llik = max;
	//cout << "DONE\n"; cout.flush();
	many_local_hillclimb_wmig_all();
	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
}

void GraphState2::add_mig(string p1, string p2){
	map<string, Graph::vertex_descriptor> pop2node = tree->get_tips(tree->root);
	Graph::vertex_descriptor n1 = pop2node[p1];
	Graph::vertex_descriptor n2 = pop2node[p2];
	add_mig( tree->g[n1].index, tree->g[n2].index);
}


// Start of functions added by EKM
void GraphState2::mlno_print_edge_w_params(Graph::edge_descriptor e) {
	Graph::vertex_descriptor sv, tv;
	sv = source(e, tree->g); 
	tv = target(e, tree->g);

	// Source vertex info
	cout << " edge " << tree->g[sv].index;
	if (tree->g[sv].is_mig) cout << " (is_mig f=" << tree->g[sv].mig_frac << ")";
	if (in_degree(sv, tree->g) > 1) cout << " (mig target)";

	// Target vertex info
	cout << " -> " << tree->g[tv].index;
	if (tree->g[tv].is_tip) cout << "=" << tree->g[tv].name;
	if (tree->g[tv].is_mig) cout << " (is_mig f=" << tree->g[tv].mig_frac << ")";
	if (in_degree(tv, tree->g) > 1) cout << " (mig target)";

	// Edge info
	cout << " : l = " << tree->g[e].len << ", w = " << tree->g[e].weight;
	if (tree->g[e].is_mig) cout << " (is_mig)";	

	cout << "\n";
	cout.flush();
}


void GraphState2::mlno_print_graph_w_params() {
	// Define variables
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	queue<Graph::vertex_descriptor> vqueue;
	Graph::vertex_descriptor v, tv;	

	cout << " graph with root " << tree->g[tree->root].index << "\n";

	// Print tree edges in preorder traversal
	cout << "+ Base tree:\n";
	vqueue.push(tree->root);
	while(vqueue.size() > 0) {
		v = vqueue.front();
		out_eiter = out_edges(v, tree->g);
		for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
			if (!tree->g[*out_ei].is_mig) {
				mlno_print_edge_w_params(*out_ei);
				tv = target(*out_ei, tree->g);
				vqueue.push(tv);
			}
		}
		vqueue.pop();
	}

	// Print migration edges
	cout << "+ Migration edges:\n";
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (tree->g[*ei].is_mig)
			mlno_print_edge_w_params(*ei);
	}

	cout.flush();
}


void GraphState2::mlno_print_vertex(Graph::vertex_descriptor v) {
	cout << " vertex ";
	cout << tree->g[v].index;
	if (tree->g[v].is_tip) cout << " (tip " << tree->g[v].name << ")";
	if (tree->g[v].is_root) cout << " (root)";
	if (in_degree(v, tree->g) > 1) cout << " (mig target)";
	if (tree->g[v].is_mig) cout << " (is_mig)";
	cout << "\n";
	cout.flush();
}


void GraphState2::mlno_print_edge(Graph::edge_descriptor e) {
	Graph::vertex_descriptor sv, tv;
	sv = source(e, tree->g); 
	tv = target(e, tree->g);
	cout << " edge " << tree->g[sv].index;
	if (in_degree(sv, tree->g) > 1) cout << " (mig target)";
	if (tree->g[sv].is_mig) cout << " (is_mig)";
	cout << " -> " << tree->g[tv].index;
	if (tree->g[tv].is_tip) cout << " (tip " << tree->g[tv].name << ")";
	if (in_degree(tv, tree->g) > 1) cout << " (mig target)";
	if (tree->g[e].is_mig) cout << " (is_mig)";
	cout  << "\n";
	cout.flush();
}


void GraphState2::mlno_print_graph() {
	// Define variables
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	queue<Graph::vertex_descriptor> vqueue;
	Graph::vertex_descriptor v, tv;	

	cout << " graph with root " << tree->g[tree->root].index << "\n";

	// Print tree edges in preorder traversal
	cout << "+ Base tree:\n";
	vqueue.push(tree->root);
	while(vqueue.size() > 0) {
		v = vqueue.front();
		out_eiter = out_edges(v, tree->g);
		for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
			if (!tree->g[*out_ei].is_mig) {
				mlno_print_edge(*out_ei);
				tv = target(*out_ei, tree->g);
				vqueue.push(tv);
			}
		}
		vqueue.pop();
	}

	// Print migration edges
	cout << "+ Migration edges:\n";
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (tree->g[*ei].is_mig)
			mlno_print_edge(*ei);
	}

	cout.flush();
}


void GraphState2::mlno_print_edges(vector<Graph::edge_descriptor> &edges) {
	Graph::edge_descriptor e;
	Graph::vertex_descriptor sv, tv;
	cout << " edges = [";
	for (int i = 0; i < edges.size(); i++) {
		sv = source(edges[i], tree->g);
		tv = target(edges[i], tree->g);
		cout << " " << tree->g[sv].index << "->" << tree->g[tv].index;
	}
	cout << " ]\n";
}


void GraphState2::mlno_fit_graph(bool do_clear) {
	/*
	 * Fits branch lengths and migration weights to graph
	 *
	 * Added by EKM in January 2021
	 */
 
	if (tree->isbinary) {
		cout << "ERROR in mlno_fit_graph: Graph is not in the correct format!\n";
		exit(1);
	}

	//mlno_print_graph_w_params();

	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	Graph::vertex_descriptor sv;

	// Re-set graph parameters
	/*if (do_clear) {
		eiter = edges(tree->g);
		for (ei = eiter.first; ei != eiter.second; ++ei) {
			sv = source(*ei, tree->g);
			if (tree->g[*ei].is_mig) {
				tree->g[sv].mig_frac = 0.5;
				tree->g[*ei].weight = 0.0;
				tree->g[*ei].len = 0.0;
			} else {
				tree->g[sv].mig_frac = 0.0;
				tree->g[*ei].weight = 1.0;
				tree->g[*ei].len = 1.0;
			}
		}
	}*/

	if (params->f2) set_branches_ls_f2();
	else set_branches_ls();
	//cout << "Done initial branch length optimization\n"; mlno_print_graph_w_params();

	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (tree->g[*ei].is_mig) {
			tree->g[*ei].weight = 0.1;  // or 0.05
			if (params->f2) {
				initialize_migupdate();
				optimize_weight_quick(*ei);
			} else {
				optimize_weight(*ei);
			}
		}
	}
	//cout << "Done inital weight optimization\n"; mlno_print_graph_w_params();

	if (params->f2) {
		set_branches_ls_f2();
		initialize_migupdate();
		optimize_weights_quick();
	} else {
		set_branches_ls();
		optimize_weights();
	}
	//cout << "Done final paramerter optimization\n"; mlno_print_graph_w_params();

	current_llik = llik(); // Set log-likelihood

	// IMPORTANT for avoiding segfaults!
	e2index.clear();
	e2frac.clear();
	e2tips.clear();
	popname2paths.clear();
	popnames2index.clear();
}


void GraphState2::mlno_fit_graph_mlbt() {
	pair<bool, double> outpair;
	bool is_treebased;

	outpair = mlno_find_best_base_tree();
	is_treebased = outpair.first;
	if (!is_treebased) {
		cerr << "Only call function mlno_fit_graph_mlbt if graph is tree-based!!\n";
		exit(1);
	}

	mlno_fit_graph();
	current_llik = llik(); 
}


bool GraphState2::mlno_fit_graph_mlno() {
	bool is_reoriented = mlno_find_best_orientation();
	mlno_fit_graph();
	current_llik = llik();
	return is_reoriented;
}


pair<bool, double> GraphState2::mlno_find_best_base_tree() {
	/*
 	 * The likelihood of the graph seems to depend on which
 	 * edges are selected to be in the base tree and which
 	 * edges are selected to be migrations.
 	 * 
 	 * This implements an exhaustive search for the best scoring
 	 * base tree -- it should be replaced with something smarter.
	 *
	 * Recall that TreeMix requires the graph to be tree-based.
	 * The requirements are as follows.
 	 *
 	 * The root
 	 * + has indegree 0 and outdegree 2.
 	 * + Both (outgoing) edges incident must be 'tree edges'.
	 *
	 * The leaves
	 * + Leaf vertices have indegree >= 1 and outdegree 0.
	 * + Exactly one of the incoming edges must be a 'tree edge'.
	 * 
	 * The internal vertices
	 * + Internal vertices have indegree >= 1 and outdegree 2.
 	 * + Exactly one of the incoming edges must be a 'tree edge'.
 	 * + At least one of the outgoing edges must be a 'tree edge'.
	 *
 	 * Consider a single migration target vertex. We need to pick
 	 * exactly one of incoming edge to be a 'tree edge'. So if
 	 * there is one migration target vertex in the graph with indegree
 	 * 3, there are 3 potential base trees, although not all of them
 	 * will necessarily result in a tree-based graph. For example, 
 	 * suppose the migration target vertex is incident to the root,
 	 * then the edge incident to the root *must* be selected to be
 	 * in the base tree. If there are n migration targets, then there
 	 * are |indegree(m1)| * |indegree(m2)| * ... * | indegree(mn)|
 	 * potential base trees to evaluate.
	 * 
	 * NOTE: This function finds a topology. The following functions can
	 * be called to update the graph state.
	 *     mlno_fit_graph();
	 *     current_llik = llik();
	 * See also implemented in mlno_fit_graph_mlbt();
	 *
	 * Added by EKM in January 2021
	 */

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	list<vector<Graph::edge_descriptor> > combos;
	vector<vector<Graph::edge_descriptor> > evecs;
	vector<Graph::edge_descriptor> combo, keep_combo;
	Graph::edge_descriptor e;

	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	Graph::vertex_descriptor sv, tv, tmp_sv, tmp_tv;

	pair<bool, double> toreturn;
	double keep_llik;
	int total;
	bool is_treebased;
	bool found_base_tree = false;

	vector<double> keepllikvec;

	//cout << "Entering find_best_base_tree\n";

	keep_llik = -1 * numeric_limits<double>::max();

	if (tree->isbinary) {
		//Transform graph from binary to TreeMix format
		mlno_binary_to_treemix_graph();
		//cout << "\nTransformed binary to TreeMix"; mlno_print_graph_w_params();
	} else {
		// Graph in correct format - remove 'is_mig' labels
		eiter = edges(tree->g);
		for (ei = eiter.first; ei != eiter.second; ++ei) {
			sv = source(*ei, tree->g);
			tree->g[sv].is_mig = false;
			tree->g[*ei].is_mig = false;
		}
    }

    //cout << "\nStarting"; mlno_print_graph_w_params();

	// Find migration edges by migration target vertex
	//cout << "Finding migration edges (i.e. edges incident to migration target)\n";
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (in_degree(*vi, tree->g) > 1) {
			vector<Graph::edge_descriptor> new_evec;
			in_eiter = in_edges(*vi, tree->g);
			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
				new_evec.push_back(*in_ei);
			}
			//mlno_print_edges(new_evec);
			evecs.push_back(new_evec);
		}
	}

	// Find combinations for selecting one tree edge 
	// i.e. setting one is_mig = false for one incomine edge
	for (int j = 0; j < evecs[0].size(); j++) {
		vector<Graph::edge_descriptor> new_combo;
		new_combo.push_back(evecs.at(0).at(j));
		combos.push_back(new_combo);
	}
	for (int i = 1; i < evecs.size(); i++) {
		total = combos.size();
		for (int k = 0; k < total; k++) {
			vector<Graph::edge_descriptor> tmp_combo = combos.front();
			for (int j = 0; j < evecs[i].size(); j++) {
				vector<Graph::edge_descriptor> new_combo;
				new_combo.resize(tmp_combo.size());
				copy(tmp_combo.begin(), tmp_combo.end(), new_combo.begin());
				new_combo.push_back(evecs.at(i).at(j));
				combos.push_back(new_combo);
			}
			combos.pop_front();
		}
	}

	//cout << "Finding edge combos\n";
	//for (std::list<vector<Graph::edge_descriptor> >::iterator it=combos.begin(); it != combos.end(); ++it) {
	//	combo = *it;
	//	mlno_print_edges(combo);
	//}

	for (std::list<vector<Graph::edge_descriptor> >::iterator it=combos.begin(); it != combos.end(); ++it) {
		combo = *it;
		//cout << "Looking at"; mlno_print_edges(combo); cout.flush();

		// Label migration edges
		for (int j = 0; j < combo.size(); j++) {
			e = combo.at(j);
			in_eiter = in_edges(target(e, tree->g), tree->g);
			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
				// Remove mig info for all vertices
				sv = source(*in_ei, tree->g);
				tree->g[sv].is_mig = false;
				tree->g[sv].mig_frac = 0.0;

				// Label edge as migration
				if (*in_ei == e) {
					tree->g[*in_ei].is_mig = false;
					//cout << "Found tree"; mlno_print_edge(*in_ei);
				} else {
					tree->g[*in_ei].is_mig = true;
					//cout << "Found migration"; mlno_print_edge(*in_ei);
				}
			}
		}

		eiter = edges(tree->g);
		for (ei = eiter.first; ei != eiter.second; ++ei) {
			if (tree->g[*ei].is_mig) {
				// Set vertices to be migration sources
				sv = source(*ei, tree->g);
				tree->g[sv].is_mig = true;
				tree->g[sv].mig_frac = 0.5;

				tree->g[*ei].weight = 0.0;
				tree->g[*ei].len = 0.0;
			} else {
				tree->g[*ei].weight = 1.0;
				tree->g[*ei].len = 1.0;
			}
		}

		//cout << "Labeled"; mlno_print_graph_w_params();
		is_treebased = mlno_is_labeled_treebased();
		if (!is_treebased) continue;
		found_base_tree = true;

		// Fit branch lengths
		mlno_fit_graph(false);
		keepllikvec.push_back(current_llik);

		// Compute log-likelihood score
		if ((current_llik - keep_llik) > 1e-8) {
			keep_combo = combo;
			keep_llik = current_llik;
		}
	}

	// Label migration edges and the source vertices
	for (int j = 0; j < keep_combo.size(); j++) {
		e = keep_combo[j];
		tv = target(e, tree->g);
		in_eiter = in_edges(tv, tree->g);
		for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
			sv = source(*in_ei, tree->g);
			if (*in_ei == e) {
				tree->g[sv].is_mig = false;
				tree->g[*in_ei].is_mig = false;
			} else {
				tree->g[sv].is_mig = true;
				tree->g[*in_ei].is_mig = true;
			}
		}
	}

	if (params->doscore) {
		cout << "\nFound " << keepllikvec.size() << " base trees with scores: ";
		if (keepllikvec.size() > 0) { 
			cout << keepllikvec[0];
			for (int i = 1; i < keepllikvec.size(); i++) cout << "," << keepllikvec[i];
		}
	    cout << "\n";
	} else {
		if (found_base_tree) cout << " is treebased"; cout.flush();
		//else cout << " is NOT treebased"; cout.flush();
	}

	toreturn.first = found_base_tree;
	toreturn.second = keep_llik;

	//cout << setprecision(12) << "\nBest base tree has llik of " << keep_llik << "\n";

	return toreturn;
}


bool GraphState2::mlno_is_labeled_treebased() {
	/*
 	 * Checks if graph is tree-based.
 	 *
 	 * TreeMix requires the graph to be tree-based, that is,
 	 * the edges of the graph have a tree-based labeling.
	 * 
 	 * The label attribute is called 'is_mig'
 	 * + 0 means the edge is in the base tree
 	 * + 1 means the edge is 'gene flow' or 'migration'
 	 *
 	 * The requirements of a tree-based labeling are now given.
 	 *
 	 * The root has
 	 * + has indegree 0 and outdegree 2.
 	 * + both (outgoing) edges incident labeled 0.
	 *
	 * Each leaf vertex has
	 * + indegree >= 1 and outdegree 0.
	 * + exactly one of its incoming edges labeled 0.
	 * 
	 * Each internal vertex has
	 * + indegree >= 1 and outdegree 2.
 	 * + exactly one of its incoming edges labeled 0.
 	 * + at least one of its outgoing edges labeled 0.
 	 *
	 * Added by EKM in January 2021
	 */

	if (tree->isbinary) {
		cout << "ERROR in mlno_is_labeled_treebased: Graph is not in the proper state!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	Graph::vertex_descriptor tv;
	int count_tree_node;

	// Check edges incident to the root are part of the base tree.
	out_eiter = out_edges(tree->root, tree->g);
	for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
		if (tree->g[*out_ei].is_mig) return false;
	}

	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		//cout << "Processing "; mlno_print_vertex(*vi);

		// Check exactly one incoming edge is part of the base tree.
		if (!tree->g[*vi].is_root) {
			count_tree_node = 0;
			in_eiter = in_edges(*vi, tree->g);
			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
				if (!tree->g[*in_ei].is_mig) count_tree_node++;
			}
			if (count_tree_node != 1) return false;
		}

		// Check at least one of the two outgoing edge is part of the base tree.
		if (!tree->g[*vi].is_tip) {
			count_tree_node = 0;
			out_eiter = out_edges(*vi, tree->g);
			for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
				if (!tree->g[*out_ei].is_mig) {
					count_tree_node++;
					break;
				}
			}
			if (count_tree_node < 1) return false;
		}
	}

	return true;
}


bool GraphState2::mlno_find_best_orientation() {
	/*
 	 * Run exhaustive search for maximum likelihood network orientation
 	 *
 	 * NOTE: This function finds a topology. The following functions can
	 * be called to update the graph state.
	 *     mlno_fit_graph();
	 *     current_llik = llik();
	 * See also implemented in mlno_fit_graph_mlno();
 	 *
 	 * Added by EKM in January 2021
 	 */ 

	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;	
	vector<pair<int, int> > root_einds;
	vector<set<int> > admixture_vind_combos;  // Change to unordered set

	pair<bool, double> outpair;
	double tmp_llik, max_llik;
	int nmig, total, k;
	bool is_treebased, is_orientable, is_reoriented;

	// Store current graph and its likelihood
	cout << setprecision(12) << "Entering mlno_doit with llik " << current_llik << "\n"; cout.flush();

	// Refit graph and compute its likelihood
	mlno_fit_graph();

	// Store current graph and its likelihood
	max_llik = current_llik;
	tree_maxllik->copy(tree);
	is_reoriented = false;

	//cout << "Starting"; mlno_print_graph_w_params();

	// Transform graph from TreeMix format to binary format
	nmig = mlno_treemix_to_binary_graph();
	//cout << "Transformed TreeMix to binary"; mlno_print_graph();

	if (nmig == 0) {
		cout << "ERROR in mlno_doit: No migration events found!\n";
		exit(1);
	}
	
	// Specify orientations by selecting admixture vertices and root edges
	// TO DO: Remove use of indicies?? Can probably use admixture vertices if
	// using tree_mlno?
	mlno_get_root_einds(root_einds);
	mlno_get_admixture_vind_combos(nmig, admixture_vind_combos);

	// Evaluate likelihood of all valid orientations
	total = root_einds.size() * admixture_vind_combos.size();
	k = 1;
	for (int i = 0; i < root_einds.size(); i++) {
		for (int j = 0; j < admixture_vind_combos.size(); j++) {
			cout << "\n";
			cout << "Orientation " << k << " / " << total; 
			cerr << "Orientation " << k << " / " << total << "\n";

			tree->copy(tree_binary);
			k++;

			// Apply reorientation algorithm
			is_orientable = mlno_reorient_huber2019(root_einds[i], admixture_vind_combos[j]);
			if (!is_orientable) continue;
			cout << " is valid"; cout.flush();
			//cout << "\nReoriented"; mlno_print_graph();

			// Fit branch lengths and migration weights and compute likelihood
			outpair = mlno_find_best_base_tree();
			is_treebased = outpair.first;
			if (!is_treebased) continue;
			tmp_llik = outpair.second;
			//cout << "\nScored"; mlno_print_graph();

			// Compare log-likelihood score
			cout << setprecision(12) << " has llik " << tmp_llik;  cout.flush();
			if ((tmp_llik - max_llik) > 1e-8) {
				cout << " which is higher!";  cout.flush();
				is_reoriented = true;
				tree_maxllik->copy(tree);
				max_llik = tmp_llik;
			}
		}
	}

	tree->copy(tree_maxllik);

	cout << setprecision(12) << "\nBest orientation has llik of " << current_llik << "\n";
	

	return is_reoriented;
}


bool GraphState2::mlno_reroot_at_outgroup() {
	/*
 	 * Reroot graph at so that the root is incident to the outgroup.
 	 *
 	 * To do this, we reorient the graph so that it is rooted at the
 	 * outgroup, but we DO NOT change the migration target vertices.
 	 * This is only possible if the outgroup is NOT a migration
 	 * target.
 	 *
 	 * Added by EKM in January 2021
 	 */ 

	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	set<int> admixture_vinds;  // Change to unordered set
	pair<int, int> root_eind;
	pair<bool, int> outpair;
	double tmp_llik, bk_llik;
	int nmig;
	bool is_treebased, is_orientable;

	cout << "Checking position of the root\n";

	// Store current graph and its log-likelihood
	tree_bk->copy(tree);
	bk_llik = current_llik;

	//cout << "Starting"; mlno_print_graph_w_params();

	// Transform graph from TreeMix format to binary format
	nmig = mlno_treemix_to_binary_graph();
	//cout << "Transformed TreeMix to binary"; mlno_print_graph();

	// Find edge incident to root
	root_eind = mlno_get_root_eind_to_outgroup();

	if (root_eind.first == tree->g[tree->root].index) {
		tree->copy(tree_bk);
		current_llik = bk_llik;
		cout << "...fine.\n";
		return true;
	}
	if (nmig == 0) {
		place_root(params->root);
		current_llik = bk_llik;
		cout << "...placed root at outgroup.\n";
		return true;
	}

	// Find admixture target vertices
	// TODO: Make this into a function
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (in_degree(*vi, tree->g) == 2)
			admixture_vinds.insert(tree->g[*vi].index);
	}

	// Apply reorientation algorithm
	//cout << "Orientation "; cout.flush();
	is_orientable = mlno_reorient_huber2019(root_eind, admixture_vinds);
	if (!is_orientable) {
		cout << "\n...unable to place root at outgroup.\n";
		tree->copy(tree_bk);
		current_llik = bk_llik;
		return false;
	}
	//cout << "Reoriented "; mlno_print_graph_w_params();

	// Find best base tree
	outpair = mlno_find_best_base_tree();
	is_treebased = outpair.first;
	if (!is_treebased) {
		cout << "\n...unable to place root at outgroup.\n";
		tree->copy(tree_bk);
		current_llik = bk_llik;
		return false;
	}

	cout << "\n...placed root at outgroup.\n";

	// Fit graph
	mlno_fit_graph();
	//cout << "Best base tree "; mlno_print_graph_w_params();

	if ((bk_llik - current_llik) > 1e-8) {
		cout << "WARNING: Re-rooting graph at outgroup changed llik from "
		     << bk_llik << " to " << current_llik << "!\n";
	}

	return true;
}


int GraphState2::mlno_check_if_nonbinary() {
	/*
 	 * Checks that network does not violate conditions for algorithms
 	 * added by EKM to work. These conditions are that the tree is
 	 * binary and that leaf nodes have in-degree 1.
 	 *
 	 * Added by EKM in January 2021
 	 */

	// Define variables
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	int deg, indeg, nmig;

	nmig = 0;

	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		indeg = in_degree(*vi, tree->g);
		deg = indeg + out_degree(*vi, tree->g);

		if (tree->g[*vi].is_root) {
			if (deg != 2) return -1;
		} else if (tree->g[*vi].is_tip) {
			if (deg != 1) return -1;
		} else {
			if (deg != 3) return -1;
		}

		if (indeg == 2) nmig++;
	}

	return nmig;
}


int GraphState2::mlno_treemix_to_binary_graph() {
	/*
 	 * The data structure used by TreeMix is a contracted version
 	 * of a binary phylogenetic network, meaning that all edges
 	 * whose target vertex has outdegree 1 are contracted!
 	 * 
 	 * This can result in leaf vertices with indegree > 1. This is not
 	 * allowed by many algorithms, including the algorithm from Huber
 	 * et al., 2019. Therefore, we refine network prior to reorienting. 
 	 *
 	 * If a binary network has 'stacked' admixture nodes, then the
 	 * contracted versino of this network will have a vertex has indegree
 	 * > 2, so multiple refinements are possible. Our current approach
 	 * is to refine these vertices arbitrarily. There could be issues due
 	 * to the arbitary branching order introduced. However, we expect that
 	 * orientation is most helpful during the *first* edge addition for the
 	 * grouping of stacked vertices.
 	 * 
 	 * To come back to the main point, this function copies 'tree' to
 	 * 'tree_binary' and then updates 'tree_binary' so that its
 	 * topology is binary.
 	 *
 	 * Added by EKM in January 2021
 	 */

	if (tree->isbinary) {
		cout << "ERROR in mlno_treemix_to_binary_graph: Graph is not in the correct format!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	vector<Graph::edge_descriptor> evec;
	Graph::edge_descriptor e;

	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
        vector<Graph::vertex_descriptor> vvec;
	Graph::vertex_descriptor v, sv, sv1, sv2, tv, newv;

	int nmig, is_nonbinary;
	bool is_incident_to_root;

	// IMPORTANT for avoiding segfaults
	e2index.clear();
	e2frac.clear();
	e2tips.clear();
	popname2paths.clear();
	popnames2index.clear();

	//cout << "Preparing tree for MLNO algorithms\n";
	tree_binary->copy(tree);

	// Count number of migration edges
	nmig = 0;
	eiter = edges(tree_binary->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (tree_binary->g[*ei].is_mig) nmig++;
	}

	// Find migration target (i.e. admixture vertices)
	viter = vertices(tree_binary->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (in_degree(*vi, tree->g) > 1)
			vvec.push_back(*vi);
	}

	// Resolve migration targets (i.e. admixture vertices)
	// so they out-degree 1 (instead of 0 or 2)
	for (int i = 0; i < vvec.size(); i++) {
		v = vvec[i];
		//cout << "Resolving migration target"; mlno_print_vertex(v);

		// Create new admixture node
		newv = add_vertex(tree_binary->g);
		tree_binary->g[newv].name = "NA";
		tree_binary->g[newv].is_root = false;
		tree_binary->g[newv].is_tip = false;
		tree_binary->g[newv].index = tree_binary->indexcounter;
		tree_binary->indexcounter++;
		//cout << "  Created new"; mlno_print_vertex(newv);

		// Attach admixture edges to new admixture node
		evec.clear();
		in_eiter = in_edges(v, tree_binary->g);
		for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
			sv = source(*in_ei, tree_binary->g);
			e = add_edge(sv, newv, tree_binary->g).first;
			//cout << "  Adding"; mlno_print_edge(e);
			evec.push_back(*in_ei);
		}

		// Remove admixture edges from old admixture node
		for (int j = 0; j < evec.size(); j++) {
			//cout << "  Removing"; mlno_print_edge(evec[j]);
			remove_edge(evec[j], tree_binary->g);
		}

		// Attach old admixture node to new admixture node
		e = add_edge(newv, v, tree_binary->g).first;
		//cout << "  Adding"; mlno_print_edge(e);
	}

	// Resolve migration target (i.e. admixture) vertices so they
	// in-degree 2 (instead of >2)
	for (int i = 0; i < vvec.size(); i++) {
		in_eiter = in_edges(vvec[i], tree_binary->g);
		v = source(*in_eiter.first, tree->g);

		if (in_degree(v, tree_binary->g) > 2) {
			//cout << "Again resolving migration target"; mlno_print_vertex(v);

			cout << "WARNING: More than one migration edges are entering the same target!\n";

			// Get incoming edges
			queue<Graph::edge_descriptor> equeue;
			is_incident_to_root = false;
			in_eiter = in_edges(v, tree_binary->g);
			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
				if (source(*in_ei, tree_binary->g) == tree_binary->root) 
					is_incident_to_root = true;
				else equeue.push(*in_ei);
			}
			if (!is_incident_to_root) equeue.pop();

			while (equeue.size() > 1) {
				e = equeue.front();
				sv1 = source(e, tree_binary->g);
				//cout << "  Removing"; mlno_print_edge(e);
				remove_edge(e, tree->g);
				equeue.pop();

				e = equeue.front();
				sv2 = source(e, tree_binary->g);
				//cout << "  Removing"; mlno_print_edge(e);
				remove_edge(e, tree->g);
				equeue.pop();

				newv = add_vertex(tree_binary->g);
				tree_binary->g[newv].name = "NA";
				tree_binary->g[newv].is_root = false;
				tree_binary->g[newv].is_tip = false;
				tree_binary->g[newv].index = tree_binary->indexcounter;
				tree_binary->indexcounter++;
				//cout << "  Created new"; mlno_print_vertex(newv);

				e = add_edge(sv1, newv, tree_binary->g).first;
				//cout << "  Adding"; mlno_print_edge(e);
				e = add_edge(sv2, newv, tree_binary->g).first;
				//cout << "  Adding"; mlno_print_edge(e);
				e = add_edge(newv, v, tree_binary->g).first;
				//cout << "  Adding"; mlno_print_edge(e);

				equeue.push(e);
			}
		}
	}

	// Process edges
	eiter = edges(tree_binary->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		tree_binary->g[*ei].weight = 1;
		tree_binary->g[*ei].len = 0;
		tree_binary->g[*ei].is_labeled = false;
		tree_binary->g[*ei].is_mig = false;
		tree_binary->g[*ei].is_oriented = false;
	}

	// Process vertices
	viter = vertices(tree_binary->g); 
	for (vi = viter.first; vi != viter.second; ++vi) {
		tree_binary->g[*vi].desired_indegree = 1;
		tree_binary->g[*vi].current_indegree = 0;
		tree_binary->g[*vi].is_mig = false;
		tree_binary->g[*vi].mig_frac = 0;  // Set to 0.5 if is_mig is true
		tree_binary->g[*vi].height = 0;    // How is this variable used?
		tree_binary->g[*vi].rev = false;
	}

	tree_binary->g[tree->root].desired_indegree = 0;
	tree_binary->g[tree->root].current_indegree = 0;
	tree_binary->isbinary = true;

	tree->copy(tree_binary);

	// Check graph is now binary
	is_nonbinary = mlno_check_if_nonbinary();
	if (is_nonbinary != nmig) {
		cout << "ERROR in mlno_treemix_to_binary_graph!\n";
		exit(1);
	}

	return nmig;
}

void GraphState2::mlno_binary_to_treemix_graph() {
	/* The data structure used by TreeMix is a contracted version
 	 * of a binary phylogenetic network, meaning that all edges
 	 * whose target vertex has outdegree 1 are contracted!
 	 *
 	 * This function takes checks if 'tree' is binary, and if so,
 	 * it contracts edges so that the resulting network is suitable
 	 * for TreeMix. 
 	 *
 	 * We loose information about the tree-based labeling, i.e. the 'is_mig'
 	 * labels.
	 *
	 * Added by EKM in January 2021
 	 */

	if (!tree->isbinary) {
		cout << "ERROR in mlno_binary_to_treemix_graph: Graph is not in the correct format!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	vector<Graph::edge_descriptor> evec;
	Graph::edge_descriptor e;

	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	vector<Graph::vertex_descriptor> vvec;
	Graph::vertex_descriptor sv, tv;

	// Suppress migration target vertices (note: migration
	// target vertices can be stacked so it's important to do
	// this carefully)
	viter = vertices(tree->g);
        for (vi = viter.first; vi != viter.second; ++vi) {
		if ((in_degree(*vi, tree->g) > 1) && (out_degree(*vi, tree->g) == 1)) {
			//cout << "  Processing migration "; mlno_print_vertex(*vi);
			vvec.push_back(*vi);  // Remove outside of loop

			// Process outgoing edges
			out_eiter = out_edges(*vi, tree->g);
			evec.push_back(*out_eiter.first);
			tv = target(*out_eiter.first, tree->g);

			// Process incoming edges
			in_eiter = in_edges(*vi, tree->g);
			for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
				sv = source(*in_ei, tree->g);
				e = add_edge(sv, tv, tree->g).first;
				tree->g[e].weight = 1;
				tree->g[e].len = 0;
				tree->g[e].is_labeled = false;
				tree->g[e].is_mig = false;
				tree->g[e].is_oriented = false;
				//cout << "    Adding "; mlno_print_edge(e);
				evec.push_back(*in_ei);  // Remove outside of loop
			}

			// Remove incoming and outgoing edges
			for (int i = 0; i < evec.size(); i++) {
				//cout << "    Removing"; mlno_print_edge(evec[i]);
				remove_edge(evec[i], tree->g);
			}
			evec.clear();
		}
	}

	// Remove detached vertices from graph
	for (int i = 0; i < vvec.size(); i++) {
		//cout << "  Removing migration "; mlno_print_vertex(vvec[i]);
		remove_vertex(vvec[i], tree->g);
	}

	tree->isbinary = false;
}


bool GraphState2::mlno_reorient_huber2019(pair<int, int> &root_eind, set<int> &admixture_vinds) {
	/*
 	 * Reorient graph using algorithm 1 in Huber et al. 2019
 	 *
 	 * Added by EKM in January 2021
 	 */  

	if (!tree->isbinary) {
		cout << "ERROR in mlno_reorient_huber2019: Graph is not in the proper state!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	vector<Graph::vertex_descriptor> vvec;
	queue<Graph::vertex_descriptor> vqueue;
	Graph::vertex_descriptor v, sv, tv;

	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	vector<Graph::edge_descriptor> evec;
	Graph::edge_descriptor e, new_root_edge;

	set<int>::const_iterator it;

	// Print message
	cout << " root edge = " << root_eind.first << " -> " << root_eind.second
	     << " admixture vertices = {";
	for (it = admixture_vinds.begin(); it != admixture_vinds.end(); ++it)
		cout << " " << *it;
	//cout << "\n";
	cout << " }";
	cout.flush();

	// Process admixture vertices
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (admixture_vinds.count(tree->g[*vi].index) > 0) {
    			tree->g[*vi].desired_indegree = 2;
		}
	}

	// Process root vertex
	mlno_move_root(root_eind);

	out_eiter = out_edges(tree->root, tree->g);
	for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
		//cout << "    Orientating "; mlno_print_edge(*out_ei);
		tree->g[*out_ei].is_oriented = true;

		tv = target(*out_ei, tree->g);
		tree->g[tv].current_indegree += 1;

		if (tree->g[tv].current_indegree == tree->g[tv].desired_indegree) {
			//cout << "    Adding to queue "; mlno_print_vertex(tv);
			vqueue.push(tv);
		}	
	}

	if (vqueue.size() == 0) return false;

	// Process remaining vertices
	while (vqueue.size() > 0) {
		v = vqueue.front();

		//cout << "  Processing "; mlno_print_vertex(v);

		// For each unoriented edge
		out_eiter = out_edges(v, tree->g);
        	for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
			if (!tree->g[*out_ei].is_oriented) {
				tv = target(*out_ei, tree->g);
				//cout << "    Orientating "; mlno_print_edge(*out_ei);

				// Edge already oriented away from root
				tree->g[*out_ei].is_oriented = true;

				// Check if edge's target vertex can  be added to queue
				tree->g[tv].current_indegree += 1;
				if (tree->g[tv].current_indegree == tree->g[tv].desired_indegree) {
					//cout << "    Adding to queue "; mlno_print_vertex(tv);
					vqueue.push(tv);
				}
			}
        	}

		in_eiter = in_edges(v, tree->g);
        	for (in_ei = in_eiter.first; in_ei != in_eiter.second; ++in_ei) {
			if (!tree->g[*in_ei].is_oriented) {
				sv = source(*in_ei, tree->g);
				//cout << "    Orientating "; mlno_print_edge(*in_ei);

				// Reverse edge so that it is oriented away from root
				e = add_edge(v, sv, tree->g).first;
				tree->g[e].weight = 1;
				tree->g[e].len = 0;
				tree->g[e].is_labeled = false;
				tree->g[e].is_mig = false;
				tree->g[e].is_oriented = true;
				evec.push_back(*in_ei);  // Remove outside of loop

				// Check if edge's target vertex can be added to queue
				tree->g[sv].current_indegree += 1;
				if (tree->g[sv].current_indegree == tree->g[sv].desired_indegree) {
					//cout << "    Adding to queue "; mlno_print_vertex(sv);
					vqueue.push(sv);
				}
			}
        	}
		for (int i = 0; i < evec.size(); i++) remove_edge(evec[i], tree->g);
		evec.clear();

		vqueue.pop();
	}

	// Check that all edges have been oriented
	eiter = edges(tree->g);
        for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (!tree->g[*ei].is_oriented) return false;
	}

	return true;
}


void GraphState2::mlno_move_root(pair<int, int> &root_eind) {
	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	vector<Graph::edge_descriptor> evec;
	Graph::edge_descriptor e, root_edge;
	Graph::vertex_descriptor sv, tv, v1, v2;
	bool found_root_edge;

	// Process new root edge
	found_root_edge = false;
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		sv = source(*ei, tree->g);
		tv = target(*ei, tree->g);
		if ((tree->g[sv].index == root_eind.first) && (tree->g[tv].index == root_eind.second)) {
			//cout << "Do not need to move root\n";
			found_root_edge = true;
			root_edge = *ei;
			break;
		}
	}

	if (!found_root_edge) {
		cout << "ERROR: New root edge "
		     << root_eind.first << " -> " << root_eind.second
		     << " does not exist in graph!\n";
		exit(1);
	}

	if (tree->root == source(root_edge, tree->g)) return;

	// Suppress root vertex in graph
	out_eiter = out_edges(tree->root, tree->g);
	v1 = target(*out_eiter.first, tree->g);
	out_eiter.first++;
	v2 = target(*out_eiter.first, tree->g);
	e = add_edge(v1, v2, tree->g).first;
	evec.push_back(e);
	remove_edge(tree->root, v1, tree->g);
	remove_edge(tree->root, v2, tree->g);

	// Subdivide root edge with root vertex
	sv = source(root_edge, tree->g);
	tv  = target(root_edge, tree->g);
	e = add_edge(tree->root, sv, tree->g).first;
	evec.push_back(e);
	e = add_edge(tree->root, tv, tree->g).first;
	evec.push_back(e);
	remove_edge(sv, tv, tree->g);

	for (int i = 0; i < evec.size(); i++) {
		e = evec[i];
		tree->g[e].weight = 1;
		tree->g[e].len = 0;
		tree->g[e].is_labeled = false;
		tree->g[e].is_mig = false;
		tree->g[e].is_oriented = false;
	}
}


pair<int, int> GraphState2::mlno_get_root_eind_to_outgroup() {
	/*
 	 * Get index of edge incident to outgroup
 	 *
 	 * Added by EKM in January 2021
 	 */  
	if (!params->set_root) {
		cerr << "ERROR in mlno_get_root_eind_to_outgroup: Not outgroup specified by user!\n";
		exit(1);
	}

	// Define variables
        pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
        Graph::edge_iterator ei;
        Graph::vertex_descriptor sv, tv;
        pair<int, int> eind;
        bool found_out_group;
	vector<pair<int, int> > root_eind;

	// Find outgroup
	found_out_group = false;
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		tv = target(*ei, tree->g);
		if (params->root.compare(tree->g[tv].name) == 0) {
			found_out_group = true;
			sv = source(*ei, tree->g);
			break;
		}
	}

	if (!found_out_group) {
		cerr << "ERROR in mlno_get_root_eind_to_outgroup: Outgroup does not exist!\n";
		exit(1);
	}

	// Store indices
	eind = make_pair(tree->g[sv].index, tree->g[tv].index);
	return eind;
}


void GraphState2::mlno_get_root_einds(vector<pair<int, int> > &root_einds) {
	/* 
 	 * Get indices of edges for re-orienting the graph
 	 *
 	 * Added by EKM in January 2021
 	 */ 

	if (!tree->isbinary) {
		cerr << "ERROR in mlno_get_root_einds: Graph is not in the correct format!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	vector<Graph::edge_descriptor> evec;
	Graph::vertex_descriptor sv, tv;
	pair<int, int> eind;
	bool found_root_edge;

	// If rooted at an outgroup
	if (params->set_root) {
		eind = mlno_get_root_eind_to_outgroup();
		root_einds.push_back(eind);
		return;
	}

	// Otherwise get all edges,
	// minus one of the edges incident to the root
	found_root_edge = false;
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		sv = source(*ei, tree->g);
		tv = target(*ei, tree->g);
		eind = make_pair(tree->g[sv].index, tree->g[tv].index);

		if (tree->g[sv].is_root) {
			if (!found_root_edge) {
				found_root_edge = true;
				root_einds.push_back(eind);
			}
		} else {
			root_einds.push_back(eind);
		}
	}
}


void get_combos_util(vector<vector<int> > &combos, vector<int> &tmp, int n, int left, int k) {
	/*
 	 * Taken from https://www.geeksforgeeks.org/make-combinations-size-k/
 	 *
 	 * Added by EKM in January 2021
 	 */
  
	// Pushing this vector to a vector of vector 
	if (k == 0) {
        	combos.push_back(tmp);
        	return;
    	}

	// i iterates from left to n. First time 
	// left will be 1 
	for (int i = left; i <= n; ++i) {
		tmp.push_back(i);
		get_combos_util(combos, tmp, n, i + 1, k - 1);

		// Popping out last inserted element 
		// from the vector 
		tmp.pop_back();
	}
}


void GraphState2::mlno_get_admixture_vind_combos(int &nmig, vector<set<int> > &admixture_vind_combos) {
	/* 
 	 * Get sets of admixture vertices for reorientating graph
 	 *
 	 * NOTE: Here we use admixture is the *target* vertex of a migration edge.
 	 *
 	 * Added by EKM in January 2021
 	 */ 

	if (!tree->isbinary) {
		cerr << "ERROR in mlno_get_admixture_vind_combos: Graph is not in the correct format!\n";
		exit(1);
	}

	// Define variables
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	vector<vector<int> > int_combos;
	vector<int> vinds, tmp;
	set<int> vind_set;

	// Get indices of allowed vertices
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (!tree->g[*vi].is_tip && !tree->g[*vi].is_root)
			vind_set.insert(tree->g[*vi].index);
	}
	vinds.resize(vind_set.size());
	copy(vind_set.begin(), vind_set.end(), vinds.begin());

	// Find (n choose k) combinations of integers 1, 2, ..., n
	int n = vinds.size();
	int k = nmig;
	get_combos_util(int_combos, tmp, n, 1, k);

	// Store as admixture vertex sets
	for (int i = 0; i < int_combos.size(); i++) {
		vind_set.clear();
		for (int j = 0; j < int_combos[i].size(); j++)
			vind_set.insert(vinds[int_combos.at(i).at(j) - 1]);
		admixture_vind_combos.push_back(vind_set); 
	}
}


void GraphState2::mlno_get_vertices_below(Graph::vertex_descriptor v, set<Graph::vertex_descriptor> &vset) {
	/*
 	 * Put vertex and all vertices below it into a set
 	 *
 	 * Added by EKM in January 2021
 	 */

	// Define variables
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	queue<Graph::vertex_descriptor> vqueue;
	Graph::vertex_descriptor u, tu;

	// Add edges at and below to set
	vqueue.push(v);
	while (vqueue.size() > 0) {
		u = vqueue.front();
		out_eiter = out_edges(u, tree->g);
		for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
			tu = target(*out_ei, tree->g);
			vset.insert(tu);
			vqueue.push(tu);
		}
		vqueue.pop();
	} 
}


pair<bool, pair<int, int> > GraphState2::mlno_add_mig_to_base_tree_exhaustive() {
	/*
	 * Tries all legal ways of adding a migration edge to the base tree
	 *
	 * Added by EKM in January 2021
	 */ 

	// Define variables
	Graph::edge_descriptor e;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	vector<pair<Graph::vertex_descriptor, Graph::vertex_descriptor> > vpairs;
	pair<Graph::vertex_descriptor, Graph::vertex_descriptor> max_vpair;
	vector<Graph::vertex_descriptor> vvec;
	set<Graph::vertex_descriptor> vset; 
	Graph::vertex_descriptor v1, v2;
	pair<bool, pair<int, int> > toreturn;
	double max_llik;

	// Store current likelihood
	max_llik = current_llik;
	toreturn.first = false;

	// Get all non-root vertex indices
	// TO DO: Make this into a function
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (*vi != tree->root)
			vvec.push_back(*vi);
	}

	// Find legal edges to add
	// TO DO: Make this into a function
	for (int i = 0; i < vvec.size(); i++) {
		v1 = vvec[i];
		mlno_get_vertices_below(v1, vset);
		for (int j = 0; j < vvec.size(); j++) {
			v2 = vvec[j];
			if (tree->g[v2].is_mig) continue;
			if (vset.find(v2) == vset.end())
				vpairs.push_back(make_pair(v2, v1));
		}
		vset.clear();
	}

	// Evaluate the likelihood of all legal edges
	for (int i = 0; i < vpairs.size(); i++) {
		v1 = vpairs[i].first;
		v2 = vpairs[i].second;

		//cout << "Evaluating " << tree->g[v1].index << " " << tree->g[v2].index << "\n"; cout.flush();

		e = tree->add_mig_edge(v1, v2);
		if (params->f2) {
			initialize_migupdate();
			optimize_weight_quick(e);
		} else {
			optimize_weight(e);
		}

		if (current_llik > max_llik) {
			max_llik = current_llik;
			max_vpair.first = v1;
			max_vpair.second = v2;
			toreturn.second.first = tree->g[v1].index;
			toreturn.second.second = tree->g[v2].index;
			toreturn.first = true;
		}

		tree->remove_mig_edge(e);
		tree->indexcounter--;
	}

	if (toreturn.first == true) {
		e = tree->add_mig_edge(max_vpair.first, max_vpair.second);
		if (params->f2) {
			initialize_migupdate();
			optimize_weight_quick(e);
		} else {
			optimize_weight(e);
		}

		// Better safe than sorry
		mlno_compute_sigma_cor();
		current_llik = llik();
	}

	return toreturn;
}


void GraphState2::mlno_compute_sigma_cor() {
	/* Computes the expected value of covariances or f2-statistics, storing
         * the result in 'sigma_cor', for current graph state 'tree', *without*
         * refitting the graph.
	 *
         * Useful whenever 'sigma_cor' matrix does not match current graph state,
         * for example when a graph is read from an input file.
         *
         * Added by EKM in March 2021
         */
	if (params->f2) mlno_compute_sigma_cor_f2();
	else mlno_compute_sigma_cor_cov();
}


void GraphState2::mlno_compute_sigma_cor_cov() {
	/*
	 * Computes the expected value of covariances, storing the result in
         * 'sigma_cor', for the current graph state 'tree', *without* refitting
         * the graph.
         *
	 * Added by EKM in March 2021 - based on set_branches_ls() function
         */

	map<Graph::edge_descriptor, int> edge2index;
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;
	map<Graph::edge_descriptor, double> edge2frac;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			edge2index.insert(make_pair(*it, index));
			edge2frac.insert(make_pair(*it, 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		edge2index[*it] = joint_index;
		edge2frac.insert(make_pair(*it, 0.5));
	}
	index++;

	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2 + current_npops; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated

	//set up the workspace
	// solve Ax = b  with x >=0
	// A = weighted branches
	// x = branch lengths (to be solved)
	// b = observed f_2 matrix

	double * A = new double[n*p]; //A holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                              //   is in the path to the root.

	double * b  = new double[n];  // y contains the entries of the empirical covariance matrix
	double * x = new double[p];   // x will be estimated, contains the fitted branch lengths
	double rNorm;

	set_branch_coefs_nnls(A, b, n, p, &edge2index, &edge2frac);

	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		//double l = x[i];
		//tree->g[*it].len = l*frac;
		x[i] = tree->g[*it].len / frac;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += A[ k * n + index ] * x[k];
			}
			//cout << i << " " << j << " "<< pred << "\n";
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	//free memory
	delete [] A;
	delete [] b;
	delete [] x;
}


void GraphState2::mlno_compute_sigma_cor_f2() {
	/*
	 * Computes the expected value of f2-statistics, storing the result in
	 * 'sigma_cor', for the current graph state 'tree', *without* refitting
	 * the graph.
	 *
	 * Added by EKM in March 2021 - based on set_branches_ls_f2() function
	 */

	map<Graph::edge_descriptor, int> edge2index;
	set<Graph::edge_descriptor> root_adj = tree->get_root_adj_edge(); //get the ones next to the root
	vector<Graph::edge_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;
	map<Graph::edge_descriptor, double> edge2frac;
	for (Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (root_adj.find(*it) == root_adj.end() && tree->g[*it].is_mig == false) {
			i_nodes2.push_back( *it );
			edge2index.insert(make_pair(*it, index));
			edge2frac.insert(make_pair(*it, 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::edge_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) {
		edge2index[*it] = joint_index;
		edge2frac.insert(make_pair(*it, 0.5));
	}
	index++;

	//initialize the workspace
	int n = (current_npops * (current_npops-1))/2; //n is the total number of entries in the covariance matrix
	int p = index; // p is the number of branches lengths to be estimated

	//set up the workspace
	// solve Ax = b  with x >=0
	// A = weighted branches
	// x = branch lengths (to be solved)
	// b = observed f_2 matrix

	double * A = new double[n*p]; //A holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                              //   is in the path to the root.

	double * b  = new double[n];  // y contains the entries of the empirical covariance matrix
	double * x = new double[p];   // x will be estimated, contains the fitted branch lengths
	double rNorm;

	set_branch_coefs_f2_nnls(A, b, n, p, &edge2index, &edge2frac);

	//and put in the solutions in the graph
	negsum = 0;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		//double l = x[i];
		//tree->g[*it].len = l*frac;
		x[i] = tree->g[*it].len / frac;
		//if (tree->g[source(*it, tree->g)].is_mig  && tree->g[target(*it, tree->g)].is_mig ) continue;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i+1; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += A[ k * n + index ] * x[k];
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
			gsl_matrix_set(sigma_cor, j, i, pred);
			index++;
		}

	}

	// free memory
	delete [] A;
	delete [] b;
	delete [] x;
}


void GraphState2::print_treeout(string stem) {
	/*
	 * Added by EKM in March 2021 - Moved from TreeMix main function
	 */
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	Graph::vertex_descriptor p1, p2;
	pair<double, double> se;
	double test, pval;

	string outfile = stem + ".treeout.gz";
	ogzstream treeout(outfile.c_str());
	
	if (params->sample_size_correct == false) treeout << get_trimmed_newick() << "\n";
	else treeout << tree->get_newick_format() << "\n";

	eiter = edges(tree->g);
    ei = eiter.first;
    params->smooth_lik = false;

    while (ei != eiter.second){
    	if (tree->g[*ei].is_mig){
			double w = tree->g[*ei].weight;

			treeout << tree->g[*ei].weight << " ";

			if (params->calc_se){
				p1 = source(*ei, tree->g);
				p1 = tree->get_child_node_mig(p1);
				p2 = target(*ei, tree->g);

				cout << tree->get_newick_format(p1) << " ";
				cout << tree->get_newick_format(p2) << "\n";
				cout.flush();

				params->neg_penalty = 0;
				se = calculate_se(*ei);
				treeout << se.first << " " << se.second << " ";

				test = se.first / se.second;
				pval = 1 - gsl_cdf_gaussian_P(test, 1);
				if (pval < DBL_MIN){
					pval = DBL_MIN;
					treeout << "<" << pval << " ";
				} else treeout << pval << " ";
			}
			else treeout << "NA NA NA ";

     		tree->g[*ei].weight = w;
     		if (params->f2) set_branches_ls_f2();
     		else set_branches_ls();

     		p1 = source(*ei, tree->g);
     		p1 = tree->get_child_node_mig(p1);
     		p2 = target(*ei, tree->g);
     		treeout << tree->get_newick_format(p1) << " ";
     		treeout << tree->get_newick_format(p2) << "\n";
    	}
		ei++;
    }
}


bool GraphState2::mlno_apply_c1_treebased(Graph::edge_descriptor e, bool is_mig) {
	/*
 	 * Target vertex of edge e has in-degree 2, so in order for e to
 	 * be labeled 'is_mig', we need to apply condition C1 from
 	 * Francies & Steel, 2015. Recall that C1 says that if one of the
	 * if one incoming arc is labeled 'tree', then label the other
	 * incoming arc 'migration' and vice versa.
	 *
	 * Added by EKM in January 2021
	 */

	//cout << "    Checking C1 for label " << is_mig << " and "; mlno_print_edge(e);

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;	
	Graph::edge_descriptor ne;

	// Get relevant neighboring edge
	in_eiter = in_edges(target(e, tree->g), tree->g);
	if (e == *in_eiter.first) in_eiter.first++;
	ne = *in_eiter.first;

	// Process label for neighboring edge
	if (tree->g[ne].is_labeled) {
		if (tree->g[e].is_mig == tree->g[ne].is_mig)
			return false;
	} else {
		mlno_label_and_propagate_treebased(ne, !is_mig);
	}

	return true;
}


bool GraphState2::mlno_apply_c2_treebased(Graph::edge_descriptor e, bool is_mig) {
	/*
 	 * Source vertex of edge e has out-degree 2, so in order for e to
 	 * be labeled 'is_mig', we need to apply condition C2 from
 	 * Francies & Steel, 2015. Recall that C2 says that if one of the
 	 * outgoing arcs has label 'migration', then label the other
 	 * outgoing arc 'tree'.
 	 *
 	 * Added by EKM in January 2021
	 */

	//cout << "    Checking C2 for label " << is_mig << " and "; mlno_print_edge(e);

	// Define variables
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::edge_descriptor ne;

	// Get relevant neighboring edge
	out_eiter = out_edges(source(e, tree->g), tree->g);
	if (e == *out_eiter.first) out_eiter.first++;
	ne = *out_eiter.first;

	// Process labeling
	if (tree->g[ne].is_labeled) {
		// Validate C2
		if (tree->g[ne].is_mig && tree->g[e].is_mig)
			return false;
	} else {
		// Propagate label if necessay
		if (is_mig) {
			if (!mlno_label_and_propagate_treebased(ne, false))
				return false;
		}
	}

	return true;
}


bool GraphState2::mlno_label_and_propagate_treebased(Graph::edge_descriptor e, bool is_mig) {
	/*
 	 * Assign label 'is_mig' to edge e and propagate the label to the
	 * surronding edges using conditions C1 and C2 from Francies &
	 * Steel, 2015.
	 *
	 * Added by EKM in January 2021 
	 */
	if (tree->g[e].is_labeled) {
		if (tree->g[e].is_mig != is_mig) return false;
		return true;
	}

	//cout << "  Labeling "; mlno_print_edge(e);

	// Label edge
	tree->g[e].is_labeled = true;	
	tree->g[e].is_mig = is_mig;

	// Propagate labels using conditions C1 and C2
	if (in_degree(target(e, tree->g), tree->g) == 2)
		return mlno_apply_c1_treebased(e, is_mig);
	else if (out_degree(source(e, tree->g), tree->g) == 2)
		return mlno_apply_c2_treebased(e, is_mig);
}


bool GraphState2::mlno_label_topology_treebased() {
	/*
 	 * Labels edges of graph topology as either 'base tree edges' or
 	 * 'gene flow edges' so that the graph is tree-based. If this
 	 * is possible, 'true' is returned; otherwise, 'false' is
 	 * returned.
 	 *
 	 * This is important because TreeMix can only fit branch lengths
 	 * if the graph is tree-based. NOTE that if there are multiple
 	 * migration targets intersecting the same tree edge, we suppress
 	 * migration target vertices, and save this information as the
 	 * height (I think). In this representation all non-root vertices
 	 * have exactly one incoming tree edge and can have any number of
 	 * incoming migration edges.
 	 *
 	 * We use the algorithm from Francis & Steel, 2015:
 	 * https://doi.org/10.1093/sysbio/syv037
	 *
	 * Added by EKM in January 2021
	 */

	if (!tree->isbinary) {
		cerr << "ERROR in mlno_label_topology_treebased: Graph is not in the proper state!\n";
		exit(1);
	}

	//cout << "Checking if tree-based\n";

	// Define variables
	pair<Graph::edge_iterator, Graph::edge_iterator> eiter;
	Graph::edge_iterator ei;
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	Graph::edge_descriptor e;

	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	Graph::vertex_descriptor sv, tv;

	// Label arcs that MUST be tree edges i.e. not migration edges
	eiter = edges(tree->g);
	for (ei = eiter.first; ei != eiter.second; ++ei) {
		sv = source(*ei, tree->g);
		tv = target(*ei, tree->g);
		if ((tree->g[tv].is_tip) ||
		    (out_degree(sv, tree->g) == 1) ||
		    (tree->g[sv].is_root)) {
			if (!mlno_label_and_propagate_treebased(*ei, false))
				return false;
		}
	}

	// Ensure edges incident to vertices with in-degree 2 are labeled
	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		if (in_degree(*vi, tree->g) == 2) {
			in_eiter = in_edges(*vi, tree->g);
			e = *in_eiter.first;
			if (!tree->g[e].is_labeled) {
				if (!mlno_label_and_propagate_treebased(e, false))
					return false;
			}
		}
	}

	// Label migration source vertices
	eiter = edges(tree->g);
        for (ei = eiter.first; ei != eiter.second; ++ei) {
		if (tree->g[*ei].is_mig)
			tree->g[source(*ei, tree->g)].is_mig = true;
	}

	return true;
}


bool GraphState2::mlno_label_topology_treechild() {
	/*
 	 * Labels edges of graph topology as either 'tree edges' or
 	 * 'migration edges' so that the graph is tree-child. If this
 	 * is possible, 'true' is returned; otherwise, 'false' is
 	 * returned. 
 	 *
	 * Added by EKM in January 2021
	 */

	if (!tree->isbinary) {
		cerr << "ERROR in mlno_label_topology_treechild: Graph is not in the proper state!\n";
		exit(1);
	}

	//cout << "Checking if tree-child\n"; cout.flush();

	// Define variables
	pair<Graph::in_edge_iterator, Graph::in_edge_iterator> in_eiter;
	Graph::in_edge_iterator in_ei;
	pair<Graph::out_edge_iterator, Graph::out_edge_iterator> out_eiter;
	Graph::out_edge_iterator out_ei;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> viter;
	Graph::vertex_iterator vi;
	Graph::vertex_descriptor tv;
	bool found_tree_node;

	viter = vertices(tree->g);
	for (vi = viter.first; vi != viter.second; ++vi) {
		//cout << "Processing "; mlno_print_vertex(*vi);

		if (!tree->g[*vi].is_tip) {
			found_tree_node = false;
			out_eiter = out_edges(*vi, tree->g);
			for (out_ei = out_eiter.first; out_ei != out_eiter.second; ++out_ei) {
				tv = target(*out_ei, tree->g);
				if (in_degree(tv, tree->g) == 1) {
					found_tree_node = true;
					break;
				}
			}
			if (!found_tree_node) return false;
		}

		// Label migration edges and source vertices
		if (in_degree(*vi, tree->g) == 2) {
			in_eiter = in_edges(*vi, tree->g);
			tree->g[*in_eiter.first].is_mig = true;
			tree->g[source(*in_eiter.first, tree->g)].is_mig = true;
		}
	}

	return true;
}
// End of functions added by EKM
