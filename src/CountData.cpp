/*
 * CountData.cpp
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */
#include "CountData.h"



CountData::CountData(string infile, PhyloPop_params* p){
	params = p;
	if (p->micro) {
		read_micro_data(infile);
		//cout << "read\n"; cout.flush();
	}
	else read_counts(infile);
	if (p->restrict_pop) npop = p->pops2use;
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop, npop);
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	cov_var2 = gsl_matrix_alloc(npop, npop);
	U = gsl_matrix_alloc(npop-1, npop);
	scatter_prime = gsl_matrix_alloc(npop-1, npop-1);
	nblock = nsnp/ p->window_size;
	ncomp = (npop* (npop-1))/2+ npop;
	//cout << nwind << " "<< ncomp << "\n";
	//cov_samp = gsl_matrix_alloc(nblock, ncomp);
	//gsl_matrix_set_zero(cov_samp);
	cov_cov = gsl_matrix_alloc(ncomp, ncomp);
	if (p->micro) set_alfreqs_micro();
	else set_alfreqs();
	scale_alfreqs();
	if (p->read_hzy) set_hzy_fromfile(p->hzyfile);
	//set_scatter();
	//process_scatter();
	if (p->f2) {
		set_cov_f2();
		////if (p->cor_f2){
		//	correct_f2s(p->f2_corpop, p->migfracs, p->f2_mixdist);
		//}
	}
	else set_cov();
	///set_cov2();
	//set_ne();
	//set_ne2();
	//set_ncomp_ef();
	//cout << "Effective number of SNPs: "<< ne << "\n";
	//process_cov();
}


CountData::CountData(string infile){
	cerr << "ERROR: Do not use this constructor for CountData!\n"; exit(1);
	read_counts(infile);
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop, npop);
	set_alfreqs();
	//scale_alfreqs(which);
	set_scatter();
	process_scatter();
}


CountData::CountData(CountData * c, vector<string> names, gsl_matrix* model, PhyloPop_params* p, gsl_rng *r){
	// set covariance matrix to a random Wishart with covariance from a model, copy over the standard errors
	params= p;
	npop = names.size();
	nsnp = c->nsnp;
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	ne = c->ne;
	for (int i = 0; i < names.size(); i++){
		pop2id.insert(make_pair(names[i], i));
	}
	set_cov_ran(model, r);
	//for(int i = 0; i < names.size(); i++){
	//	for (int j = 0; j < names.size(); j++){
	//		gsl_matrix_set(cov_var, i, j, c->get_cov_var(names[i], names[j] ));
	//	}
	//}
}

void CountData::set_cov_ran(gsl_matrix* model,gsl_rng* r){

	// 1. Take SVD of model
	gsl_matrix * U2 = gsl_matrix_alloc(npop-1,npop);
	gsl_matrix * A = gsl_matrix_alloc(npop,npop);
	gsl_matrix * VT = gsl_matrix_alloc(npop,npop);
	gsl_matrix * model_prime = gsl_matrix_alloc(npop-1,npop-1);
	gsl_vector * S = gsl_vector_alloc(npop);
	gsl_vector * work = gsl_vector_alloc(npop);
	gsl_matrix_memcpy( A, model);
	gsl_matrix_set_zero(cov);
	gsl_matrix_set_zero(cov_var);
	gsl_linalg_SV_decomp(A, VT, S, work);


	// Now copy the first npop-1 eigenvectors to U

	for (int i = 0; i < npop-1; i++){
		for(int j = 0; j < npop; j++){
			gsl_matrix_set(U2, i, j, gsl_matrix_get(A, i, j));
		}
	}


	// multiply U model U^T
	gsl_matrix * US = gsl_matrix_alloc(npop-1, npop);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U2, model, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U2, 0.0, model_prime);
	int nblock = nsnp/ params->window_size;

	//initialize blocks
	vector< vector< vector<double> > > allsamp;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		allsamp.push_back(tmp1);
	}

	//sample nblock Wisharts
	for (int i = 0; i < nblock; i++){
		gsl_matrix * cov_prime = gsl_matrix_alloc(npop-1,npop-1);
		gsl_matrix * tmpcov = gsl_matrix_alloc(npop, npop);
		//get random wishart
		//const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
		rwishart(r, npop-1, ne , model_prime, cov_prime );
		//multiply U^T cov_prime U
		gsl_matrix * UTC = gsl_matrix_alloc(npop, npop-1);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U2, cov_prime, 0.0, UTC);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UTC, U2, 0.0, tmpcov);

		//cout << gsl_matrix_get(tmpcov, 1, 3) / ne << "\n";
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				gsl_matrix_set(cov, i, j, gsl_matrix_get(cov, i, j)+gsl_matrix_get(tmpcov, i, j)/ (ne/nblock));
				gsl_matrix_set(cov, j, i, gsl_matrix_get(cov, j, i)+gsl_matrix_get(tmpcov, j, i)/ (ne/nblock));
				allsamp[i][j].push_back(gsl_matrix_get(tmpcov, i, j)/ (ne/nblock));
			}
		}

		gsl_matrix_free(cov_prime);
		gsl_matrix_free(tmpcov);
	}
	//cout << gsl_matrix_get(model, 1, 3) << " model\n";
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double mean =  gsl_matrix_get(cov, i, j)/ (double) nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);
			double sum = 0;
			vector<double> all_covs = allsamp[i][j];
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= (*it-mean)*(*it-mean);
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /(double) nblock;
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);

		}
	}
	//print_cov("rancov.gz");
	//print_cov_var("rancovvar.gz");
	gsl_matrix_free(U2);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);
}
string CountData::get_pops(){
	string toreturn = "(";
	map<string , int>::iterator it = pop2id.begin();
	map<string , int>::iterator it2 = pop2id.end();
	it2--;
	toreturn+= it->first +":0.1";
	it++;
	while (it != it2){
		toreturn+=",("+it->first+":0.1";
		it++;
	}
	toreturn+=","+it->first +":0.1";
	for (int i = 0; i < npop; i++)	toreturn+= "):0.1";


	toreturn = toreturn.substr(0, toreturn.size()-4);
	toreturn+= ";";
	return toreturn;
}


vector<string> CountData::list_pops(){
	vector<string> toreturn;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++) {
		if ( params->restrict_pop == true ){
			if (it->second < params->pops2use ) toreturn.push_back( it->first);
		}
		else toreturn.push_back( it->first);
	}
	return toreturn;
}


void CountData::read_scatter(string infile){
	npop =0;
	gsl_matrix_free(scatter);
	vector<vector<double> > tmpcov;
	ifstream in(infile.c_str()); //only gzipped files
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
             buf.clear();
             stringstream ss(st);
             line.clear();
             while (ss>> buf){
                     line.push_back(buf);
             }
             vector<double> tmp;
             for (vector<string>::iterator it = line.begin(); it != line.end(); it++) tmp.push_back(atof(it->c_str()));
             tmpcov.push_back(tmp);
    }
    npop = tmpcov.size();
    scatter = gsl_matrix_alloc(npop, npop);
    for (int i = 0; i < npop; i++){
    	for (int j = 0; j< npop; j++) gsl_matrix_set(scatter, i, j, tmpcov[i][j]);
    }

}



void CountData::read_alfreqs(string infile){
	npop =0;
	gsl_matrix_free(alfreqs);
	gsl_matrix_free(scatter);
	pop2id.clear();
	vector<vector<double> > tmpalfreqs;
	igzstream in(infile.c_str()); //only gzipped files
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
             buf.clear();
             stringstream ss(st);
             line.clear();
             while (ss>> buf){
                     line.push_back(buf);
             }
             vector<double> tmp;
             for (vector<string>::iterator it = line.begin(); it != line.end(); it++) tmp.push_back(atof(it->c_str()));
             tmpalfreqs.push_back(tmp);
    }
    nsnp = tmpalfreqs.size();
    npop = tmpalfreqs[0].size();
    scatter = gsl_matrix_alloc(npop, npop);
    alfreqs = gsl_matrix_alloc(nsnp, npop);
    for (int i = 0; i < nsnp; i++){
    	for (int j = 0; j< npop; j++) gsl_matrix_set(alfreqs, i, j, tmpalfreqs[i][j]);
    }
    for(int i = 0; i < npop; i++){
    	stringstream ss;
    	ss << "pop";
    	ss << i;
    	string name = ss.str();
    	pop2id.insert(make_pair(name, i));
    }
    scale_alfreqs();
    set_scatter();

    //print_scatter("testout_scatter.gz");
	//cout << "here\n";
    set_cov();
    process_cov();

}

void CountData::read_counts(string infile){
    allele_counts.clear();
    pop2id.clear();
    id2pop.clear();
    npop = 0;
    nsnp = 0;
    rss.clear();
    chr.clear();
    pos.clear();
    a1.clear();
    a2.clear();
    string ext = infile.substr(infile.size()-3, 3);
    if (ext != ".gz"){
    	std::cerr << "ERROR: " << infile << " is not gzipped (only .gz files accepted)\n";
    	exit(1);
    }
	igzstream in(infile.c_str()); //only gzipped files
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;

    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << infile << "\n";
            exit(1);
    }

    /*
     * header contains population names
     */
    getline(in, st);
    stringstream ss(st);
    line.clear();
    while (ss>> buf){
    	line.push_back(buf);
     }
    /*
     * make map from header, number populations according to order
     */
    int start = 0;
    if (params->snpinfo) start = 5;
    for(int i = start; i < line.size(); i++) {
    	pop2id.insert(make_pair(line[i], i-start));
    	id2pop.insert(make_pair(i-start, line[i]));
    	npop ++;
    }
    int headsize = line.size();
    /*
     * read counts, store in allele_counts
     */
    while(getline(in, st)){
            buf.clear();
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            vector<pair<int, int> > topush;

            if (params->snpinfo){
            	rss.push_back(line[0]);
            	chr.push_back(line[1]);
            	pos.push_back(line[2]);
            	a1.push_back(line[3]);
            	a2.push_back(line[4]);
            }
            if (line.size() != headsize){
            	cerr << "ERROR: Line "<< nsnp <<" has "<< line.size() << " entries. Header has "<< headsize <<"\n";
            	exit(1);
            }
            for ( int i = start; i < line.size(); i++){
            	//cout <<  line[i] << "\n";
                typedef boost::tokenizer<boost::char_separator<char> >
                tokenizer;
                boost::char_separator<char> sep(",");
                tokenizer tokens(line[i], sep);
                vector<int> tmpcounts;
                for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
                        int tmp = atoi(tok_iter->c_str());
                        tmpcounts.push_back(tmp);
                }
                if (tmpcounts.size() != 2){
                	std::cerr << "ERROR: "<< line[i] << " does not have two alleles (expecting SNP data)\n";
                	exit(1);
                }
                topush.push_back(make_pair(tmpcounts[0], tmpcounts[1]));
            }
            allele_counts.push_back(topush);
            nsnp++;
    }
}

void CountData::read_micro_data(string infile){
    micro_lens.clear();
    pop2id.clear();
    id2pop.clear();
    npop = 0;
    nsnp = 0;
    rss.clear();
    chr.clear();
    pos.clear();
    a1.clear();
    a2.clear();
    string ext = infile.substr(infile.size()-3, 3);
    if (ext != ".gz"){
    	std::cerr << "ERROR:" << infile << " is not gzipped (only .gz files accepted)\n";
    	exit(1);
    }
	igzstream in(infile.c_str()); //only gzipped files
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;

    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << infile << "\n";
            exit(1);
    }

    /*
     * header contains population names
     */
    getline(in, st);
    stringstream ss(st);
    line.clear();
    while (ss>> buf){
    	line.push_back(buf);
     }
    /*
     * make map from header, number populations according to order
     */
    int start = 0;
    if (params->snpinfo) start = 5;
    for(int i = start; i < line.size(); i++) {
    	pop2id.insert(make_pair(line[i], i-start));
    	id2pop.insert(make_pair(i-start, line[i]));
    	npop ++;
    }
    int headsize = line.size();
    /*
     * read counts, store in allele_counts
     */
    while(getline(in, st)){
            buf.clear();
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            vector< vector<float> > topush;

            if (params->snpinfo){
            	rss.push_back(line[0]);
            	chr.push_back(line[1]);
            	pos.push_back(line[2]);
            	a1.push_back(line[3]);
            	a2.push_back(line[4]);
            }
            if (line.size() != headsize){
            	cerr << "ERROR: Line "<< nsnp <<" has "<< line.size() << " entries. Header has "<< headsize <<"\n";
            	exit(1);
            }
            for ( int i = start; i < line.size(); i++){
            	//cout <<  line[i] << "\n";
                typedef boost::tokenizer<boost::char_separator<char> >
                tokenizer;
                boost::char_separator<char> sep(",");
                tokenizer tokens(line[i], sep);
                vector<float> tmpcounts;
                for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
                        int tmp = atof(tok_iter->c_str());
                        tmpcounts.push_back(tmp);
                }
                if (tmpcounts.size() != 3){
                	std::cerr << "ERROR: "<< line[i] << " does not have three entries [expecting micosat data]\n";
                	exit(1);
                }
                vector<float> tmplen;
                tmplen.push_back(tmpcounts[0]); tmplen.push_back(tmpcounts[1]); tmplen.push_back(tmpcounts[2]);
                topush.push_back(tmplen);
            }
            micro_lens.push_back(topush);
            nsnp++;
    }
}

pair< vector<string>, vector<double> > CountData::get_freqs(int i){
	pair<vector<string>, vector<double> > toreturn;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		toreturn.first.push_back(it->first);
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j);
		toreturn.second.push_back(f);
	}
	return toreturn;


}


pair< vector<string>, vector<double> > CountData::get_centered_freqs(int i){
	pair<vector<string>, vector<double> > toreturn;
	double mean = 0;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j);
		if (!isnan(f)) mean+= f;
		//toreturn.second.push_back(f);
	}
	mean = mean/ (double) npop;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		toreturn.first.push_back(it->first);
		int j = it->second;
		double f = gsl_matrix_get(alfreqs, i, j) - mean;
		toreturn.second.push_back(f);
	}
	return toreturn;
}

void CountData::set_alfreqs_micro(){
	mean_ninds.clear();
	mean_var.clear();
	id2nsnp.clear();
	for (int i = 0; i < npop; i++){
		mean_ninds.insert(make_pair(i, 0.0));
		mean_var.insert(make_pair(i, 0.0));
		id2nsnp.insert(make_pair(i, 0.0));
	}
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			double m = micro_lens[i][j].at(0);
			double v = micro_lens[i][j].at(1);
			double n = micro_lens[i][j].at(2);
			if ( n < 1){
                bool display_warn = false;
                bool info_msg = false;
                if((params->num_warnings) > 0) {
                    display_warn = true;
                    --(params->num_warnings);
                    if((params->num_warnings) == 0) {
                        info_msg = true;
                    }
                } else if((params->num_warnings) == -1) {
                    display_warn = true;
                }

                if(display_warn) {
                    cerr << "WARNING: no alleles at locus "<< i << " population "<< j <<"\n";
                    if(info_msg) {
                        cerr << "WARNING: Suppressing further warnings\n";    
                    }
                }
                
				float h = 0.0;
				m = 0/h;
				//cout << m << "0 n \n"; cout.flush();
				gsl_matrix_set(alfreqs, i, j, m);
				continue;
			}
			//cout << i << " "<< j << " "<< m << "\n";
			gsl_matrix_set(alfreqs, i, j, m);
			mean_ninds[j] += n;
			mean_var[j] += v;
			id2nsnp[j]++;
		}
	}
	for (int i = 0; i < npop; i++){
		mean_ninds[i] = mean_ninds[i]/ id2nsnp[i];
		mean_var[i] = mean_var[i]/ id2nsnp[i];
		//cout << id2pop[i] << " "<< mean_hzy[i] << "\n";
	}
}


void CountData::set_alfreqs(){
	mean_ninds.clear();
	mean_hzy.clear();
	id2nsnp.clear();

	for (int i = 0; i < npop; i++){
		mean_ninds.insert(make_pair(i, 0.0));
		mean_hzy.insert(make_pair(i, 0.0));
		id2nsnp.insert(make_pair(i, 0));
	}
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			int c1 = allele_counts[i][j].first;
			int c2 = allele_counts[i][j].second;
			double f = (double) c1 / ( (double) c1 + (double) c2 );
			if ( c1+c2 < 1) {
                bool display_warn = false;
                bool info_msg = false;
                if((params->num_warnings) > 0) {
				    display_warn = true;
                    --(params->num_warnings);
                    if((params->num_warnings) == 0) {
                        info_msg = true;
                    }
                } else if((params->num_warnings) == -1) {
                    display_warn = true;
                }

                if(display_warn) {
                    cerr << "WARNING: no counts at SNP "<< i << " population "<< j <<"\n";
                    if(info_msg) {
                        cerr << "WARNING: Suppressing further warnings\n";    
                    }
                }
				gsl_matrix_set(alfreqs, i, j, f);
				continue;
			}
			gsl_matrix_set(alfreqs, i, j, f);
			mean_ninds[j] += ((double) c1+ (double) c2)/2.0;
			double tmp2 = (double) c2 / ((double) c1+ (double) c2 - 1.0);
			double tmphzy = 2* f * tmp2;
			if (c1+c2 < 2){
				tmphzy = 2*f*(1-f);
			}
			//if (id2pop[j] == "San"){
			//	cout << i << " "<< c1 << " "<< c2 << " "<< f << " "<< tmp2 << " "<< tmphzy << " "<< mean_hzy[j]<< "\n";
		//	}
			mean_hzy[j] += tmphzy; //2*f*(1-f);
			id2nsnp[j]++;
		}
	}
	for (int i = 0; i < npop; i++){
		mean_ninds[i] = mean_ninds[i]/ id2nsnp[i];
		mean_hzy[i] = mean_hzy[i]/ id2nsnp[i];
		//cout << id2pop[i] << " "<< mean_hzy[i] << "\n";
	}
}


void CountData::scale_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		double total = 0;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			double scaled;
			if (params->alfreq_scaling == 1) scaled = asin(sqrt(f));
			else scaled = f;
			total = total+scaled;
			gsl_matrix_set(alfreqs, i, j, scaled);
		}

		double m = total/ (double) npop;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			if (params->micro) gsl_matrix_set(alfreqs, i, j, f-m);
			else{
				double f = gsl_matrix_get(alfreqs, i, j);
				if (params->alfreq_scaling == 3) {
					gsl_matrix_set(alfreqs, i, j, (f-m)/sqrt(m *(1-m)) );
					if (m < 1e-8) gsl_matrix_set(alfreqs, i, j, 0);
				}
				else if (params->alfreq_scaling == 4) gsl_matrix_set(alfreqs, i, j, f);
				else gsl_matrix_set(alfreqs, i, j, f-m);
			}
		}
	}
}

void CountData::set_scatter(){
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;
			for (int k = 0; k < nsnp; k++){
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c += toadd;
			}
			gsl_matrix_set(scatter, i, j, c);
			gsl_matrix_set(scatter, j, i, c);
		}
	}
}

void CountData::set_cov(){
	/*
	 * Calculate covariance matrix in blocks on SNPs
	 *  cov[i,j] = mean( cov[i,j]_k ) over all k blocks
	 */
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_var);
	//cout << npop << "\n";
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	//cov_var2 = gsl_matrix_alloc(npop, npop);
	//initialize block estimation of covariance matrix
	vector<vector<vector<double> > > cov_block;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		cov_block.push_back(tmp1);
	}
	vector<string> popnames = list_pops();
	cov_samp.clear();
	for (int i = 0; i < popnames.size(); i++){
		map<string, vector<double> > tmp;
		for (int j = 0; j < popnames.size(); j++){
			vector<double> tmp2;
			tmp.insert(make_pair(popnames[j], tmp2));
		}
		cov_samp.insert(make_pair(popnames[i], tmp));
	}
	//if trimming covariances, get the amount to trim
	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = pop2id.begin(); it!= pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double mean_n = mean_ninds.find(id)->second;
		double t;
		if (!params->micro) {
			double meanhzy = mean_hzy.find(id)->second;
			if (params->print_hzy) cout << pop << " "<< meanhzy << "\n";
			t = meanhzy / (4.0* mean_n);
		}
		else{
			double mvar = mean_var.find(id)->second;
			t = mvar/ mean_n;
			//cout << pop << " "<< mvar << " "<< mean_n << "\n";
		}
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		//cout << pop  << " "<< t << "\n";
		if (!params->sample_size_correct){
			trim.insert(make_pair(pop, 0));
		}
		else trim.insert(make_pair(pop, t));
	}
	if (!params->sample_size_correct) sumtrim = 0;
	//calculate the covariance matrix in each block
	cout << "Estimating covariance matrix in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		int index = 0;
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				double c = 0;
				for (int n = k*params->window_size; n < (k+1)*params->window_size; n++){
					if (isnan(gsl_matrix_get(alfreqs, n, i)) || isnan(gsl_matrix_get(alfreqs, n, j))) continue;
					double toadd = gsl_matrix_get(alfreqs, n, i) * gsl_matrix_get(alfreqs, n, j);
					//cout << toadd << "\n"; cout.flush();
					c+= toadd;
				}
				double cov = c/ ((double) params->window_size);
				//cout << k << " "<< index << " "<< cov << "\n";
				string p1 = id2pop[i];
				string p2 = id2pop[j];

				double bias1 = trim[p1];
				double bias2 = trim[p2];
				cov = cov + bias1/ (double) npop + bias2/ (double) npop - sumtrim/( (double) npop* (double) npop);
				if (i ==j) cov = cov - trim[p1];

				cov_samp[p1][p2].push_back(cov);
				if (p1 != p2) cov_samp[p2][p1].push_back(cov);
				//gsl_matrix_set(cov_samp, k, index, cov);
				//cout << k << " "<< index << " "<< gsl_matrix_get(cov_samp, k, index) << "\n";
				cov_block[i][j].push_back(cov);
				index++;
			}
		}
	}
	//ofstream tout("test");
	//calculate the mean, standard error of covariance estimates
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			vector<double> all_covs = cov_block[i][j];
			double sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) {
				if (!isnan(*it)) sum+= *it;
			}
			double mean = sum/nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);

			// and standard error
			sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) {
				if (!isnan(*it)) sum+= (*it-mean)*(*it-mean);
			}
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /sqrt((double) (nblock-1) * (double)nblock);
			//cout << i << " "<< j << " "<< sd << " "<< c << "\n";
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);
		}
	}

}


void CountData::set_cov_f2(){
	/*
	 * Calculate matrix of f_2 statistics in blocks
	 *
	 */
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_var);

	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);

	//initialize block estimation of covariance matrix
	vector<vector<vector<double> > > cov_block;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		cov_block.push_back(tmp1);
	}
	vector<string> popnames = list_pops();
	cov_samp.clear();
	for (int i = 0; i < popnames.size(); i++){
		map<string, vector<double> > tmp;
		for (int j = 0; j < popnames.size(); j++){
			vector<double> tmp2;
			tmp.insert(make_pair(popnames[j], tmp2));
		}
		cov_samp.insert(make_pair(popnames[i], tmp));
	}
	//if trimming covariances, get the amount to trim
	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = pop2id.begin(); it!= pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = mean_hzy.find(id)->second;
		double mean_n = mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		trim.insert(make_pair(pop, t));
	}

	//calculate the covariance matrix in each block
	cout << "Estimating f_2 matrix in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		int index = 0;
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				double c = 0;
				string p1 = id2pop[i];
				string p2 = id2pop[j];
				for (int n = k*params->window_size; n < (k+1)*params->window_size; n++){
					if (isnan(gsl_matrix_get(alfreqs, n, i))) continue;
					if (isnan(gsl_matrix_get(alfreqs, n, j))) continue;
					double toadd = (gsl_matrix_get(alfreqs, n, i) - gsl_matrix_get(alfreqs, n, j));
					toadd = toadd*toadd;
					c+= toadd;
				}
				double cov = c/ (double) params->window_size;
				if (params->sample_size_correct && p1 != p2){

					double bias1 = trim[p1];
					double bias2 = trim[p2];

					cov = cov - bias1 -bias2;
				}



				cov_samp[p1][p2].push_back(cov);
				if (p1 != p2) cov_samp[p2][p1].push_back(cov);

				cov_block[i][j].push_back(cov);
				index++;
			}
		}
	}
	//ofstream tout("test");
	//calculate the mean, standard error of covariance estimates
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			vector<double> all_covs = cov_block[i][j];
			double sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) {
				if (isnan(*it)) continue;
				sum+= *it;
			}
			double mean = sum/nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);

			// and standard error
			sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) {
				if (isnan(*it)) continue;
				sum+= (*it-mean)*(*it-mean);
			}
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /(double) nblock;
			//cout << i << " "<< j << " "<< sd << " "<< c << "\n";
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);
		}
	}
	//get the covariance in the estimates of the covariance matrix
	gsl_matrix_set_zero(cov_cov);
}


void CountData::set_cov2(){
	/*
	 * Get SE of covariance matrix without blocks
	 */
	gsl_matrix_free(cov_var2);
	cov_var2 = gsl_matrix_alloc(npop, npop);
	gsl_matrix *tmpcov = gsl_matrix_alloc(npop, npop);
	//calculate the covariance matrix

	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;
			for (int k = 0; k < nsnp; k++){
				if (isnan(gsl_matrix_get(alfreqs, k, i))) continue;
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c+= toadd;
			}
			double cov = c/(double) nsnp;
			//cout << i << " "<< j << " "<< cov << "\n";
			gsl_matrix_set(tmpcov, i, j, cov);
			gsl_matrix_set(tmpcov, j, i, cov);
		}
	}

	//calculate the SE
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double mean = gsl_matrix_get(tmpcov, i, j);
			double sum = 0;
			for (int k = 0; k < nsnp; k++){
				double s = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);

				sum+= (s-mean)*(s-mean);
			}

			//double sd = sqrt(sum/ (double) nsnp);
			double c = sqrt(sum)/ (double) nsnp;
			//if (i == 0 && j ==2){
			//	cout << "sum "<< sum << " "<< nsnp << " "<< c << "\n";
			//}
			gsl_matrix_set(cov_var2, i, j, c);
			gsl_matrix_set(cov_var2, j, i, c);
		}
	}
	gsl_matrix_free(tmpcov);
}

void CountData::print_scatter(string outfile){
	ogzstream out(outfile.c_str());
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		out << it->first;
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(scatter, it->second, it2->second);
		out << "\n";
	}

}

void CountData::print_cov_cov(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			out << id2pop[i]<< "."<< id2pop[j] << " ";
		}
	}
	out <<  "\n";
	int index1 = 0;
	for (int l = 0; l < npop; l++){
		for (int m = l; m< npop; m++){
			out << id2pop[l]<< "."<< id2pop[m] << " ";
			int index =0;
			for (int i = 0; i < npop; i++){
				for (int j = i; j < npop; j++){
					//cout << index1 << " " << index << "\n";
					out << gsl_matrix_get(cov_cov, index1, index)<< " ";
					index++;
				}
			}
			index1++;
			out << "\n";
		}
	}
}



void CountData::print_cov_samp(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> popnames = list_pops();
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			out << popnames[i]<< "."<< popnames[j] << " ";
		}
	}
	out <<  "\n";
	for (int k = 0; k < nblock; k++){
		for (int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				string p1 = popnames[i];
				string p2 = popnames[j];
				out << cov_samp[p1][p2].at(k) << " ";
			}
		}

		out << "\n";
	}

}


void CountData::print_alfreqs(string outfile){
	ogzstream out(outfile.c_str());
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for (int i = 0; i < nsnp ; i++){
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(alfreqs, i, it2->second);
		out << "\n";
	}
}

void CountData::print_fst(string outfile){
	gsl_matrix* fst = gsl_matrix_alloc(npop, npop);
	ogzstream out(outfile.c_str());
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double total = 0;
			for (int k = 0; k < nsnp; k++){
				double f1 = gsl_matrix_get(alfreqs, k, i);
				double f2 = gsl_matrix_get(alfreqs, k, j);
				double f_hat = (f1+f2)/2;
				double between = 2*f_hat*(1-f_hat);
				double within = f1*(1-f1)+ f2*(1-f2);
				double fst_ijk;
				if (between < 1e-8) fst_ijk = 0;
				else fst_ijk = (between-within)/between;
				//cout << fst_ijk
				total+= fst_ijk;
			}
			total = total / (double) nsnp;
			gsl_matrix_set(fst, i, j, total);
			gsl_matrix_set(fst, j, i, total);
		}
	}
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		out << it->first;
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(fst, it->second, it2->second);
		out << "\n";
	}
}

void CountData::process_scatter(){
	// project scatter matrix into (npop-1) space
	// get matrix for transformation
	// get determinant and gamma

	scatter_det = 0;
	scatter_gamma = 0;
	int n = nsnp-1;
	size_t pop = npop;
	int s;

	//first do SVD on scatter matrix
	gsl_matrix * A = gsl_matrix_alloc(pop,pop);
	gsl_matrix * VT = gsl_matrix_alloc(pop,pop);
	gsl_vector * S = gsl_vector_alloc(pop);
	gsl_vector * work = gsl_vector_alloc(pop);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);

	// Now copy the first npop=1 eigenvectors to U

	for (int i = 0; i < npop-1; i++){
		for(int j = 0; j < npop; j++){
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}

	// And transform scatter into m-1 space
	// S' = U S U^T

	gsl_matrix * US = gsl_matrix_alloc(pop-1, pop);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, scatter, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U, 0.0, scatter_prime);

	// Now do LU decomposition and get determinant of scatter_prime
	gsl_matrix_free(A);
	A = gsl_matrix_alloc(npop-1, npop-1);
	gsl_matrix_memcpy( A, scatter_prime );
	gsl_permutation * p = gsl_permutation_alloc(pop-1);
	gsl_linalg_LU_decomp( A, p, &s );
	scatter_det = gsl_linalg_LU_lndet( A );


	//get the log sum of the gammas
	//cout << npop << " "<< n << "\n";
	scatter_gamma = ( (double) (pop-1) * ( (double)  npop-2.0) /4.0) * log (M_PI);
	//cout << scatter_gamma << " sg1\n";
	for (int i = 1; i <= pop-1; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);
	//cout << scatter_gamma << " sg2\n";
	cout << "scatter_gamma "<< scatter_gamma << "\n";
	cout << "scatter_det "<< scatter_det << "\n";
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_matrix_free(US);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_permutation_free(p);

}

void CountData::process_cov(){
	cov_var = gsl_matrix_alloc(npop, npop);
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;

			double tmp_cov = gsl_matrix_get(cov, i, j);
			for (int k = 0; k < nsnp; k++){
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c += (toadd-tmp_cov) * (toadd-tmp_cov);
			}
			gsl_matrix_set(cov_var, i, j, c/nsnp);
			gsl_matrix_set(cov_var, j, i, c/nsnp);
		}
	}
	//double tmp = 0;
	//int todivide =0;
	//for (int i = 0; i < npop; i++){
	//	for (int j = i; j < npop; j++)	{
	//		//cout << i << " "<< j << " "<< gsl_matrix_get(var_matrix, i, j)<< "\n";
	//		tmp+= gsl_matrix_get(var_matrix, i, j);
	//		todivide++;
	//	}
	//}
	//tmp = tmp/ (double) todivide;
	//cout << "tmp "<< tmp << "\n";
	//cov_var = tmp;
}

double CountData::get_cov(string pop1, string pop2){

	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(cov, p1, p2);
	return toreturn;
}


double CountData::get_scatter(string pop1, string pop2){

	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(scatter, p1, p2);
	return toreturn;
}

double CountData::get_cov_var(string pop1, string pop2){
	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "ERROR: No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(cov_var, p1, p2);
	return toreturn;
}

void CountData::print_cov(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::print_cov_var(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov_var, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::print_cov_var2(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov_var2, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::set_ne(){
	ne = 0;
	double tmpne = 0;
	int pairs = 0;
	for (int i = 0; i < npop; i++){
		for(int j = i ; j < npop; j++){
			double tmp =gsl_matrix_get(cov_var2, i, j)/gsl_matrix_get(cov_var, i, j);
			tmp = tmp*tmp;
			//if (i == 0&& j == 2){
			//	cout << "tmp "<< tmp <<"\n";
			//}
			tmpne+= tmp;
			pairs++;
		}
	}
	//cout << tmpne << " "<< pairs << "\n";
	tmpne = tmpne/ (double) pairs;
	tmpne = tmpne* nsnp;
	ne = int(tmpne);
}


void CountData::set_ne2(){
	ne2 = 0;
	int p = npop;
	gsl_matrix * A = gsl_matrix_alloc(npop, npop);
	gsl_matrix * VT = gsl_matrix_alloc(p,p);
	gsl_vector * S = gsl_vector_alloc(p);
	gsl_vector * work = gsl_vector_alloc(p);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);

	double s  = 0;
	double s2 = 0;
	//for (int i = 0; i < p; i++){
	//	double eig = gsl_vector_get(S, i);
	//	cout << eig << "\n";
	//}
	for (int i = 0; i < p-1; i++){
		double eig = gsl_vector_get(S, i);
		eig = eig*eig;
		s+= eig;
		s2+= eig*eig;
		//cout << eig << " "<< eig*eig << "\n";
	}
	//cout << s<< " "<< s2 << "\n";
	double num = ( (double) p+1.0)* s*s;
	double denom = ( ((double) p-1.0)* s2 ) - (s*s);
	//cout << "num denom " << num << " "<< denom << "\n";
	double tmp = num/denom;
	ne2 = int(tmp);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);

}

int rwishart(gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
  /* Wishart distribution random number generator */
  /*
   *    n        gives the dimension of the random matrix
   *    dof      degrees of freedom
   *    scale    scale matrix of dimension n x n
   *    result   output variable with a single random matrix Wishart distribution generation
   */
  int k,l;
  gsl_matrix *work = gsl_matrix_calloc(n,n);

  for(k=0; k<n; k++){
    gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
    for(l=0; l<k; l++){
      gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
    }
  }
  gsl_matrix_memcpy(result,scale);
  gsl_linalg_cholesky_decomp(result);
  gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
  gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);
  gsl_matrix_free(work);
  return 0;
}

string CountData::get_pop_in_index(int index){
	string toreturn;
	map<string , int>::iterator it = pop2id.begin();
	while (it != pop2id.end()){
		if (it->second == index) return it->first;
		it++;
	}
	if (it == pop2id.end()) {
		cerr << "ERROR: Trying to get index "<< index << " in CountData, none found\n";
		exit(1);
	}
	return toreturn;
}

/*
void CountData::set_ncomp_ef(){
	ncomp_ef = 0;
	gsl_matrix * S = gsl_matrix_alloc(ncomp, ncomp);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, cov_samp, cov_samp, 0.0, S);

	gsl_matrix * A = gsl_matrix_alloc(ncomp,ncomp);
	gsl_matrix * VT = gsl_matrix_alloc(ncomp,ncomp);
	gsl_vector * Sv = gsl_vector_alloc(ncomp);
	gsl_vector * work = gsl_vector_alloc(ncomp);
	gsl_matrix_memcpy( A, S );
	gsl_linalg_SV_decomp(A, VT, Sv, work);
	//for (int i = 0; i < ncomp ; i++){
	//	cout << gsl_vector_get(Sv, i) << "\n";
	//}
	while (ncomp_ef < ncomp && (gsl_vector_get(Sv, ncomp_ef)  > 1e-10)) ncomp_ef++; //this is the number of eigenvectors to use
	//cout << "Effective number of comparisons "<< ncomp_ef << "\n";
	gsl_matrix_free(S);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(Sv);
	gsl_vector_free(work);
}
*/

pair <double, double> CountData::f4ratio(string pop1, string pop2, string pop3, string pop4, string pop5){
	pair<double, double> toreturn;
	double mean, se;
	vector<double> ratios_block;

	int index1, index2, index3, index4, index5;
	if (pop2id.find(pop1) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop1 << "\n";
		exit(1);
	}
	index1 = pop2id[pop1];
	if (pop2id.find(pop2) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop2 << "\n";
		exit(1);
	}
	index2 = pop2id[pop2];
	if (pop2id.find(pop3) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop3 << "\n";
		exit(1);
	}
	index3 = pop2id[pop3];
	if (pop2id.find(pop4) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop4 << "\n";
		exit(1);
	}
	index4 = pop2id[pop4];

	if (pop2id.find(pop5) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop5 << "\n";
		exit(1);
	}
	index5 = pop2id[pop5];

	for (int i = 0; i < nblock ; i++){
		double c1 = 0;
		double c2 = 0;
		int tmp_nsnp = 0;
		for (int n = 0; n < nsnp; n++){

			if (isnan(gsl_matrix_get(alfreqs, n, index1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index2))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index3))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index4))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index5))) continue;
			//cout << gsl_matrix_get(alfreqs, n, index1) <<"\n";
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;

			double toadd = (gsl_matrix_get(alfreqs, n, index1) - gsl_matrix_get(alfreqs, n, index2))*(gsl_matrix_get(alfreqs, n, index3) - gsl_matrix_get(alfreqs, n, index4) );
			c1 += toadd;

			toadd = (gsl_matrix_get(alfreqs, n, index1) - gsl_matrix_get(alfreqs, n, index2))*(gsl_matrix_get(alfreqs, n, index5) - gsl_matrix_get(alfreqs, n, index4) );
			c2 += toadd;

			tmp_nsnp++;
		}
		double cov1 = c1/ (double) tmp_nsnp;
		double cov2 = c2/ (double) tmp_nsnp;

		double ratio = cov2/cov1;
		ratios_block.push_back(ratio);

	}

	// and standard error
	double sum = 0;

	for (vector<double>::iterator it = ratios_block.begin(); it != ratios_block.end(); it++) sum  += *it;
	mean = sum/ (double) nblock;


	sum = 0;
	for (vector<double>::iterator it = ratios_block.begin(); it != ratios_block.end(); it++) {
		sum+= (*it-mean)*(*it-mean);
	}
	se = ( (double) nblock- 1.0) / (double) nblock  * sum;
	se= sqrt(se);


	toreturn = make_pair(mean, se);

	return toreturn;
}


pair <double, double> CountData::f4ratio_sub(string pop1, string pop2, string pop3, string pop4, string pop5){
	pair<double, double> toreturn;
	double mean, se;
	vector<double> ratios_block;

	int index1, index2, index3, index4, index5, index6;
	if (pop2id.find(pop1) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop1 << "\n";
		exit(1);
	}
	index1 = pop2id[pop1];
	if (pop2id.find(pop2) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop2 << "\n";
		exit(1);
	}
	index2 = pop2id[pop2];
	if (pop2id.find(pop3) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop3 << "\n";
		exit(1);
	}
	index3 = pop2id[pop3];
	if (pop2id.find(pop4) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop4 << "\n";
		exit(1);
	}
	index4 = pop2id[pop4];

	if (pop2id.find(pop5) == pop2id.end()){
		cout << "ERROR: Cannot find population "<< pop5 << "\n";
		exit(1);
	}
	index5 = pop2id[pop5];

	//if (pop2id.find(pop6) == pop2id.end()){
	//	cout << "ERROR: Cannot find population "<< pop6 << "\n";
	//	exit(1);
	//}
	//index6 = pop2id[pop6];

	for (int i = 0; i < nblock ; i++){
		//double f4_x = 0;
		//double f4_y = 0;
		//double f4_2 = 0;
		//double f4_3 = 0;
		double f4_num =0;
		double f4_denom = 0;
		int tmp_nsnp = 0;
		for (int n = 0; n < nsnp; n++){

			if (isnan(gsl_matrix_get(alfreqs, n, index1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index2))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index3))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index4))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, index5))) continue;
			//if (isnan(gsl_matrix_get(alfreqs, n, index6))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;

			double f1 = gsl_matrix_get(alfreqs, n, index1);
			double f2 = gsl_matrix_get(alfreqs, n, index2);
			double f3 = gsl_matrix_get(alfreqs, n, index3);
			double f4 = gsl_matrix_get(alfreqs, n, index4);
			double f5 = gsl_matrix_get(alfreqs, n, index5);
			//double f6 = gsl_matrix_get(alfreqs, n, index6);
			//cout << f1 << " "<< f2 << " " << f3 << " "<< f4 << " "<< f5 <<  " "<< f6 << "\n";

			double toadd = (f2 - f4)*(f5 - f3);
			//f4_x += toadd;
			f4_num += toadd;
			//cout <<"add1 "<< toadd << " ";
			toadd = (f1 - f4)*(f2 -f3);
			f4_denom += toadd;
			//toadd = (f1 - f5)*(f2 -f3 );
			//f4_2 += toadd;

			//toadd = (f1 - f5)*(f3 -f4 );
			//f4_3 += toadd;
			//cout << "add2 "<< toadd << "\n";

			tmp_nsnp++;
		}
		//f4_x = f4_x / (double) tmp_nsnp;
		//f4_y = f4_y / (double) tmp_nsnp;
		//f4_2 = f4_2 / (double) tmp_nsnp;
		//f4_3 = f4_3 / (double) tmp_nsnp;

		f4_num =f4_num / (double) tmp_nsnp;
		f4_denom = f4_denom/ (double) tmp_nsnp;
		//cout << f4_num << " "<< f4_denom << "\n";
		double ratio= (f4_num/f4_denom);
		ratios_block.push_back(ratio);

	}

	// and standard error
	double sum = 0;

	for (vector<double>::iterator it = ratios_block.begin(); it != ratios_block.end(); it++) sum  += *it;
	mean = sum/ (double) nblock;


	sum = 0;
	for (vector<double>::iterator it = ratios_block.begin(); it != ratios_block.end(); it++) {
		sum+= (*it-mean)*(*it-mean);
	}
	se = ( (double) nblock- 1.0) / (double) nblock  * sum;
	se= sqrt(se);


	toreturn = make_pair(mean, se);

	return toreturn;
}


set<pair<string, pair<double, double> > > CountData::calculate_f4(int i0, int i1, int i2, int i3){
	set<pair<string, pair<double, double> > > toreturn;
	double mean1, se1, mean2, se2, mean3, se3;
	vector<double> f4_1;
	vector<double> f4_2;
	vector<double> f4_3;

	vector<double> f4_block_1;
	vector<double> f4_block_2;
	vector<double> f4_block_3;
	//calculate the covariance matrix in each block
	cout << "Estimating f_4 in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	mean1 = 0;
	mean2 = 0;
	mean3 =0;
	int total_nsnp = 0;
	for (int i = 0; i < nsnp; i++){
		if (isnan(gsl_matrix_get(alfreqs, i, i0))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i1))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i2))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i3))) continue;
		double toadd = (gsl_matrix_get(alfreqs, i, i1) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i2) );
		f4_1.push_back(toadd);
		mean1 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i1) );
		f4_2.push_back(toadd);
		mean2 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i1) );
		f4_3.push_back(toadd);
		mean3 += toadd;
		total_nsnp ++;
	}

	cout << "total_nsnp "<< total_nsnp << " nsnp "<< nsnp << "\n";
	mean1 = mean1 / (double) total_nsnp;
	mean2 = mean2 / (double) total_nsnp;
	mean3 = mean3 / (double) total_nsnp;
	for (int i = 0; i < nblock ; i++){
		double c1 = 0;
		double c2 = 0;
		double c3 = 0;
		int tmp_nsnp = 0;
		for (int n = 0; n < nsnp; n++){

			if (isnan(gsl_matrix_get(alfreqs, n, i0))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i2))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i3))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;
			double toadd = (gsl_matrix_get(alfreqs, n, i1) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i2) );
			c1+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i1) );
			c2+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i1) );
			c3+= toadd;
			tmp_nsnp++;
		}
		double cov1 = c1/ (double) tmp_nsnp;
		double cov2 = c2/ (double) tmp_nsnp;
		double cov3 = c3/ (double) tmp_nsnp;
		//cout << i <<  " "<< cov1 << " "<< cov2 << " "<< cov3 << "\n";
		f4_block_1.push_back(cov1);
		f4_block_2.push_back(cov2);
		f4_block_3.push_back(cov3);

	}

	// and standard error
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	for (vector<double>::iterator it = f4_block_1.begin(); it != f4_block_1.end(); it++) sum1  += *it;
	mean1 = sum1/ (double) nblock;

	for (vector<double>::iterator it = f4_block_2.begin(); it != f4_block_2.end(); it++) sum2  += *it;
	mean2 = sum2/ (double) nblock;

	for (vector<double>::iterator it = f4_block_3.begin(); it != f4_block_3.end(); it++) sum3  += *it;
	mean3 = sum3/ (double) nblock;

	sum1 = 0;
	sum2 = 0;
	sum3 = 0;
	for (vector<double>::iterator it = f4_block_1.begin(); it != f4_block_1.end(); it++) {
		//cout << mean1 << " "<< *it << "\n";
		sum1+= (*it-mean1)*(*it-mean1);
	}
	se1 = ( (double) nblock- 1.0) / (double) nblock  * sum1;
	se1= sqrt(se1);

	for (vector<double>::iterator it = f4_block_2.begin(); it != f4_block_2.end(); it++) {
		//cout << mean2 << " "<< *it << "\n";
		sum2+= (*it-mean2)*(*it-mean2);
	}
	se2 = ( (double) nblock- 1.0) / (double) nblock  * sum2;
	se2 = sqrt(se2);

	for (vector<double>::iterator it = f4_block_3.begin(); it != f4_block_3.end(); it++) sum3+= (*it-mean3)*(*it-mean3);
	se3 = ( (double) nblock- 1.0) / (double) nblock  * sum3;
	se3 = sqrt(se3);

	pair<double, double> tmp = make_pair(mean1, se1);
	pair<string, pair<double, double> > tmp2 = make_pair( id2pop[i0] +","+id2pop[i1]+";"+id2pop[i2]+","+id2pop[i3], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean2, se2);
	tmp2 = make_pair( id2pop[i0] +","+id2pop[i2]+";"+id2pop[i1]+","+id2pop[i3], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean3, se3);
	tmp2 = make_pair( id2pop[i0] +","+id2pop[i3]+";"+id2pop[i1]+","+id2pop[i2], tmp );
	toreturn.insert(tmp2);


	return toreturn;
}


set<pair<string, pair<double, double> > > CountData::calculate_f3(int i0, int i1, int i2){
	set<pair<string, pair<double, double> > > toreturn;
	double mean1, se1, mean2, se2, mean3, se3;
	vector<double> f3_1;
	vector<double> f3_2;
	vector<double> f3_3;

	double meanh0 = mean_hzy.find(i0)->second;
	double meanh1 = mean_hzy.find(i1)->second;
	double meanh2 = mean_hzy.find(i2)->second;
	double mean_n0 = mean_ninds.find(i0)->second;
	double mean_n1 = mean_ninds.find(i1)->second;
	double mean_n2 = mean_ninds.find(i2)->second;
	double t0 = meanh0 / (4.0* mean_n0);
	double t1 = meanh1 / (4.0* mean_n1);
	double t2 = meanh2 / (4.0* mean_n2);


	vector<double> f3_block_1;
	vector<double> f3_block_2;
	vector<double> f3_block_3;
	//calculate the covariance matrix in each block
	cout << "Estimating f_3 in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	mean1 = 0;
	mean2 = 0;
	mean3 =0;
	int total_nsnp = 0;
	for (int i = 0; i < nsnp; i++){
		if (isnan(gsl_matrix_get(alfreqs, i, i0))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i1))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i2))) continue;
		double toadd = (gsl_matrix_get(alfreqs, i, i0) - gsl_matrix_get(alfreqs, i, i1))*(gsl_matrix_get(alfreqs, i, i0) - gsl_matrix_get(alfreqs, i, i2) );
		f3_1.push_back(toadd);
		mean1 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i1) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i1) - gsl_matrix_get(alfreqs, i, i2) );
		f3_2.push_back(toadd);
		mean2 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i1) );
		f3_3.push_back(toadd);
		mean3 += toadd;
		total_nsnp ++;
	}

	cout << "total_nsnp "<< total_nsnp << " nsnp "<< nsnp << "\n";
	mean1 = mean1 / (double) total_nsnp;
	mean2 = mean2 / (double) total_nsnp;
	mean3 = mean3 / (double) total_nsnp;

	mean1 = mean1- t0;
	mean2 = mean2 - t1;
	mean3 = mean3 - t2;

	for (int i = 0; i < nblock ; i++){
		double c1 = 0;
		double c2 = 0;
		double c3 = 0;
		int tmp_nsnp = 0;
		for (int n = 0; n < nsnp; n++){

			if (isnan(gsl_matrix_get(alfreqs, n, i0))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i2))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;
			double toadd = (gsl_matrix_get(alfreqs, n, i0) - gsl_matrix_get(alfreqs, n, i1))*(gsl_matrix_get(alfreqs, n, i0) - gsl_matrix_get(alfreqs, n, i2) );
			c1+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i1) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i1) - gsl_matrix_get(alfreqs, n, i2) );
			c2+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i1) );
			c3+= toadd;
			tmp_nsnp++;
		}
		double cov1 = c1/ (double) tmp_nsnp;
		double cov2 = c2/ (double) tmp_nsnp;
		double cov3 = c3/ (double) tmp_nsnp;

		cov1 = cov1 - t0;
		cov2 = cov2 - t1;
		cov3 = cov3 - t2;

		f3_block_1.push_back(cov1);
		f3_block_2.push_back(cov2);
		f3_block_3.push_back(cov3);

	}

	// and standard error
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	for (vector<double>::iterator it = f3_block_1.begin(); it != f3_block_1.end(); it++) sum1  += *it;
	mean1 = sum1/ (double) nblock;

	for (vector<double>::iterator it = f3_block_2.begin(); it != f3_block_2.end(); it++) sum2  += *it;
	mean2 = sum2/ (double) nblock;

	for (vector<double>::iterator it = f3_block_3.begin(); it != f3_block_3.end(); it++) sum3  += *it;
	mean3 = sum3/ (double) nblock;

	sum1 = 0;
	sum2 = 0;
	sum3 = 0;
	for (vector<double>::iterator it = f3_block_1.begin(); it != f3_block_1.end(); it++) {
		//cout << mean1 << " "<< *it << "\n";
		sum1+= (*it-mean1)*(*it-mean1);
	}
	se1 = ( (double) nblock- 1.0) / (double) nblock  * sum1;
	se1= sqrt(se1);

	for (vector<double>::iterator it = f3_block_2.begin(); it != f3_block_2.end(); it++) {
		//cout << mean2 << " "<< *it << "\n";
		sum2+= (*it-mean2)*(*it-mean2);
	}
	se2 = ( (double) nblock- 1.0) / (double) nblock  * sum2;
	se2 = sqrt(se2);

	for (vector<double>::iterator it = f3_block_3.begin(); it != f3_block_3.end(); it++) sum3+= (*it-mean3)*(*it-mean3);
	se3 = ( (double) nblock- 1.0) / (double) nblock  * sum3;
	se3 = sqrt(se3);

	pair<double, double> tmp = make_pair(mean1, se1);
	pair<string, pair<double, double> > tmp2 = make_pair( id2pop[i0]+";"+id2pop[i1]+","+id2pop[i2], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean2, se2);
	tmp2 = make_pair( id2pop[i1] +";"+id2pop[i0]+","+id2pop[i2], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean3, se3);
	tmp2 = make_pair( id2pop[i2] +";"+id2pop[i0]+","+id2pop[i1], tmp );
	toreturn.insert(tmp2);


	return toreturn;
}

double CountData::calculate_f2(int p1, int p2){
	double toreturn = 0;
	for (int i = 0; i < nsnp; i++){
		double n1 = mean_ninds[p1];
		double n2 = mean_ninds[p2];
		double f1 = gsl_matrix_get(alfreqs, i, p1);
		double f2 = gsl_matrix_get(alfreqs, i, p2);

		double diff = f1-f2;
		double hz1 = f1*(1-f1);
		double hz2 = f2*(1-f2);
		double toadd = diff*diff - hz1/(2*n1) - hz2/ (2*n2);
		toreturn += toadd;
	}
	toreturn = toreturn/ (double) nsnp;
	return toreturn;
}

pair<double, double> CountData::calculate_drift(int p1, int p2){
	pair<double, double> toreturn;
	string pop1 = id2pop[p1];
	double meanhzy1 = mean_hzy.find(p1)->second;
	double mean_n1 = mean_ninds.find(p1)->second;
	double t1 = meanhzy1 / (4.0* mean_n1);
	vector<double> blocks;
	for (int i = 0; i < nblock ; i++){
		int tmp_nsnp = 0;
		double tmpnum = 0;
		double tmpdenom = 0;
		for (int n = 0; n < nsnp; n++){
			if (isnan(gsl_matrix_get(alfreqs, n, p1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, p2))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;
			double f1 = gsl_matrix_get(alfreqs, n, p1);
			double f2 = gsl_matrix_get(alfreqs, n, p2);
			tmpnum += f1*(1-f1);
			tmpdenom += f1*(1-f2);
			tmp_nsnp++;
		}
		tmpnum = tmpnum/ (double) tmp_nsnp;
		tmpdenom = tmpdenom/ (double) tmp_nsnp;
		tmpnum += t1;
		//cout << pop1 << " "<< meanhzy1 << " "<< mean_n1 << " "<< t1 << " "<< tmpnum << " "<< tmpdenom << "\n";
		double tmpc = tmpnum/tmpdenom;
		blocks.push_back(tmpc);
	}

	/*
	double num =0;
	double denom = 0;
	for (int i = 0; i < nsnp; i++){
		if (isnan(gsl_matrix_get(alfreqs, i, p1))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, p2))) continue;
		double f1 = gsl_matrix_get(alfreqs, i, p1);
		double f2 = gsl_matrix_get(alfreqs, i, p2);
		num += f1*(1-f1);
		denom += f1*(1-f2);
		//cout << f1 << " "<< f2 << "\n";
	}/
	num /= (double) nsnp;
	denom /= (double) nsnp;
	num += t1;
	//cout << num << " "<< denom << " "<< t1 << "\n";
	toreturn = num/denom;
	*/

	double mean = 0;
	for (vector<double>::iterator it = blocks.begin(); it != blocks.end(); it++) mean += *it;
	mean = mean/ (double) nblock;
	toreturn.first = mean;

	double se_sum = 0;
	for (vector<double>::iterator it = blocks.begin(); it != blocks.end(); it++) {
		se_sum+= (*it-mean)*(*it-mean);
	}

	double se1 = ( (double) nblock- 1.0) / (double) nblock * se_sum;
	se1= sqrt(se1);
	toreturn.second = se1;
	return toreturn;
}


pair<double, double> CountData::calculate_mean(int p1){
	pair<double, double> toreturn;
	string pop1 = id2pop[p1];
	double meanhzy1 = mean_hzy.find(p1)->second;
	double mean_n1 = mean_ninds.find(p1)->second;
	//double t1 = meanhzy1 / (4.0* mean_n1);
	vector<double> blocks;
	for (int i = 0; i < nblock ; i++){
		int tmp_nsnp = 0;
		double tmpnum = 0;
		double tmpdenom = 0;
		for (int n = 0; n < nsnp; n++){
			if (isnan(gsl_matrix_get(alfreqs, n, p1))) continue;
			//if (isnan(gsl_matrix_get(alfreqs, n, p2))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;
			double f1 = gsl_matrix_get(alfreqs, n, p1);
			//double f2 = gsl_matrix_get(alfreqs, n, p2);
			tmpnum += f1;
			tmp_nsnp++;
		}
		tmpnum = tmpnum/ (double) tmp_nsnp;
		//tmpdenom = tmpdenom/ (double) tmp_nsnp;
		//tmpnum += t1;
		//cout << pop1 << " "<< meanhzy1 << " "<< mean_n1 << " "<< t1 << " "<< tmpnum << " "<< tmpdenom << "\n";
		//double tmpc = tmpnum/tmpdenom;
		blocks.push_back(tmpnum);
	}


	double mean = 0;
	for (vector<double>::iterator it = blocks.begin(); it != blocks.end(); it++) mean += *it;
	mean = mean/ (double) nblock;
	toreturn.first = mean;

	double se_sum = 0;
	for (vector<double>::iterator it = blocks.begin(); it != blocks.end(); it++) {
		se_sum+= (*it-mean)*(*it-mean);
	}

	double se1 = ( (double) nblock- 1.0) / (double) nblock * se_sum;
	se1= sqrt(se1);
	toreturn.second = se1;
	return toreturn;
}

void CountData::set_cov_jackknife(int which){
	gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = 0;
			for (int k = 0; k < nblock; k++){
				if (k == which) continue;
				m+= cov_samp[p1][p2].at(k);
			}
			m = m/ (double) (nblock-1);
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}
}


void CountData::set_cov_bootstrap(gsl_rng *r){
	gsl_matrix_set_zero(cov);
	gsl_matrix_set_zero(cov_var);
	map<string, map<string, vector<double> > > samples;

	//initialize
	for (int i = 0; i < npop; i++){
		string p1 = id2pop[i];
		map<string, vector<double> > tmp;
		for (int j = i; j < npop; j++){
			string p2 = id2pop[j];
			vector<double> tmp2;
			tmp.insert(make_pair(p2, tmp2 ));
		}
		samples.insert(make_pair(p1, tmp));
	}

	//sample
	for (int i = 0; i < nblock; i++){
		int rint = gsl_rng_uniform_int(r, nblock);
		for (int j = 0; j < npop; j++){
			for (int k = j; k < npop; k++){
				string p1 = id2pop[j];
				string p2 = id2pop[k];
				samples[p1][p2].push_back( cov_samp[p1][p2].at(rint) );
			}
		}
	}

	// get mean
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = 0;
			for (int k = 0; k < nblock; k++){
				if (!isnan(samples[p1][p2].at(k))) m+= samples[p1][p2].at(k);
			}
			m = m / (double) nblock;
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}

	// and s.e.
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = gsl_matrix_get(cov, i, j);
			double s = 0;
			for (int k = 0; k < nblock; k++){

				double samp = samples[p1][p2].at(k);
				if (!isnan(samp)) s += (samp-m)*(samp-m);


			}
			double c = sqrt(s) /sqrt((double) (nblock-1) * (double)nblock);
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);
		}
	}
}


void CountData::set_cov_fromsamp(int which){
	gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = cov_samp[p1][p2].at(which);
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}

}


void CountData::set_cov_singlesnp(int which){
	gsl_matrix_set_zero(cov);
	double m = 0;
	for (int i = 0; i < npop; i++) {
		if (!isnan(gsl_matrix_get(alfreqs, which, i))) m+= gsl_matrix_get(alfreqs, which, i);
		cout << id2pop[i] << " "<< gsl_matrix_get(alfreqs, which, i) << "\n";
	}
	m = m / (double) npop;
	cout << "mean "<< m << "\n";
	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = pop2id.begin(); it!= pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = mean_hzy.find(id)->second;
		double mean_n = mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double f1 = gsl_matrix_get(alfreqs, which, i) -m;
			double f2 = gsl_matrix_get(alfreqs, which, j) -m ;
			double t1 = trim.find(p1)->second;
			double t2 = trim.find(p2)->second;
			double c = f1*f2;
			if (i == j) c = c - t1;
			c = c + t1 / (double) npop + t2/ (double) npop;
			c = c - sumtrim / (double) (npop*npop);

			//double m = cov_samp[p1][p2].at(which);
			gsl_matrix_set(cov, i, j, f1*f2);
			gsl_matrix_set(cov, j, i, f1*f2);
		}
	}

}

void CountData::set_hzy_fromfile(string infile){

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
             buf.clear();
             stringstream ss(st);
             line.clear();
             while (ss>> buf){
                     line.push_back(buf);
             }
             string pop = line[0];
             double hzy = atof(line[1].c_str());
             if (pop2id.find(pop) == pop2id.end()){
            	 cerr << "ERROR: Cannot find population "<< pop << " when reading hzy\n";
            	 exit(1);
             }
             int popid = pop2id[pop];
             mean_hzy[popid] = hzy;
             cout << "Set hzy in " << pop << " to "<< hzy << "\n";
    }
}

void CountData::correct_f2s(string mixpop, map<string, float> fracs, float d_mix){
	vector<vector<vector<double> > > corrected_blocks;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp;
		for (int j = 0; j < npop; j++){
			vector<double> tmp2;
			tmp.push_back(tmp2);
		}
		corrected_blocks.push_back(tmp);
	}
	cout << "Correcting blocked f2 stats for admixture\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		//cout <<  k << "\n";
		int index = 0;
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				double c = 0;
				string p1 = id2pop[i];
				string p2 = id2pop[j];
				//cout << p1 << " "<< p2 << "\n";
				if (p1 == p2) {

					corrected_blocks[i][j].push_back(c);
				}
				else if(p1 == mixpop && fracs.find(p2) != fracs.end()){
					c = cov_samp[p1][p2].at(k);
					c = (c-d_mix)/ ( (1-fracs[p2]) * (1- fracs[p2]));
					corrected_blocks[i][j].push_back(c);
					corrected_blocks[j][i].push_back(c);
				}
				else if (p2 == mixpop && fracs.find(p1) != fracs.end()){
					c = cov_samp[p2][p1].at(k);
					c = (c-d_mix)/ ( (1-fracs[p1])  *(1- fracs[p1]));
					corrected_blocks[i][j].push_back(c);
					corrected_blocks[j][i].push_back(c);
				}
				else if (fracs.find(p1) != fracs.end() && fracs.find(p2) != fracs.end()){
					//f2 between two admixed populations, A and B, if mixed with D
					// f2(A,B) = f2(A', B') + w_A^2 f2(A',D') + w_B^2 f2(B', D') - 2 w_A f3(A';D,B') - 2 w_B f3(B'; D, A') - 2 w_A w_B f3(D'; A',B')

					double f2_ab = cov_samp[p1][p2].at(k);
					double f2_ad = cov_samp[p1][mixpop].at(k);
					double f2_bd = cov_samp[p2][mixpop].at(k);
					double f3_adb = (f2_ad + f2_ab - f2_bd)/2;
					double f3_bda = (f2_bd + f2_ab  - f2_ad)/2;
					double f3_dab = (f2_ad+ f2_bd - f2_ab)/2;
					double wa = fracs[p1];
					double wb = fracs[p2];

					double f3_dab_cor = (f3_dab - d_mix) / ( ( 1-wb)*(1-wa) );
					double f2_ad_cor = (f2_ad - d_mix)/ ( (1 - wa)* (1-wa) );
					double f2_bd_cor = (f2_bd - d_mix) / ( (1- wb)*(1-wb) );

					double f3_adb_cor = (f3_adb + wa*wb*f3_dab_cor - wb*f3_dab_cor + wa*f2_ad_cor - wa*wa*f2_ad_cor)/ (1- wa);
					double f3_bda_cor = (f3_bda + wa*wb*f3_dab_cor - wa*f3_dab_cor + wb*f2_bd_cor - wb*wb*f2_bd_cor)/ (1- wb);

					double f2_ab_cor = f2_ab + 2*wa*wb*f3_dab_cor + 2*wb*f3_bda_cor + 2*wa*f3_adb_cor - wb*wb*f2_bd_cor - wa*wa*f2_ad_cor;
					corrected_blocks[i][j].push_back(f2_ab_cor);
					corrected_blocks[j][i].push_back(f2_ab_cor);
				}
				else if ( fracs.find(p1) != fracs.end()){
					//f2 between an admixed population A and an unadmixed population Y
					// f2(A,Y) = f2(A',Y) + w^2 f2(D', A) - 2w f3(A'; Y, D)
					double f2_ay = cov_samp[p1][p2].at(k);
					double f2_ad = cov_samp[p1][mixpop].at(k);
					double f2_yd = cov_samp[p2][mixpop].at(k);
					double f3_dya = (f2_yd+ f2_ad - f2_ay)/2;
					double f3_ady = (f2_ad + f2_ay - f2_yd)/2;
					double wa = fracs[p1];


					double f3_dya_cor = (f3_dya - d_mix)/ wa;

					double f2_ad_cor = (f2_ad - d_mix)/ ( (1 - wa)* (1-wa) );
					double f3_ady_cor = f3_ady - wa* f3_dya_cor + 2*wa*f2_ad_cor - wa*wa*f2_ad_cor;

					double f2_ay_cor = f2_ay + 2*wa * f3_ady_cor - wa*wa*f2_ad_cor;
					corrected_blocks[i][j].push_back(f2_ay_cor);
					corrected_blocks[j][i].push_back(f2_ay_cor);
				}
				else if (fracs.find(p2) != fracs.end()){

					double f2_ay = cov_samp[p1][p2].at(k);
					double f2_ad = cov_samp[p2][mixpop].at(k);
					double f2_yd = cov_samp[p1][mixpop].at(k);
					double f3_dya = (f2_yd+ f2_ad - f2_ay)/2;
					double f3_ady = (f2_ad + f2_ay - f2_yd)/2;
					double wa = fracs[p2];

					double f3_dya_cor = (f3_dya - d_mix)/ wa;

					double f2_ad_cor = (f2_ad - d_mix)/ ( (1 - wa)* (1-wa) );
					double f3_ady_cor = f3_ady - wa* f3_dya_cor + 2*wa*f2_ad_cor - wa*wa*f2_ad_cor;

					double f2_ay_cor = f2_ay + 2*wa * f3_ady_cor - wa*wa*f2_ad_cor;
					corrected_blocks[i][j].push_back(f2_ay_cor);
					corrected_blocks[j][i].push_back(f2_ay_cor);

				}
				else{
					c = cov_samp[p1][p2].at(k);
					corrected_blocks[i][j].push_back(c);
					corrected_blocks[j][i].push_back(c);
				}
			}
		}
	}
	//Now replace the entries in the f2 matrix with the corrected versions

	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double mean = 0;
			for (int k = 0; k < nblock; k++){
				mean += corrected_blocks[i][j][k];
			}
			mean = mean/ (double) nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);
		}
	}
	//And get the SEs
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			if (i == j) gsl_matrix_set(cov_var, i, j, 0);
			else{
				double c = 0;
				double mean = gsl_matrix_get(cov, i, j);
				for (int k = 0; k < npop; k++){
					double tmp = corrected_blocks[i][j][k];
					c += (tmp - mean)*(tmp- mean);
				}

				c = sqrt(c) /(double) nblock;
				gsl_matrix_set(cov_var, i, j, c);
				gsl_matrix_set(cov_var, j, i, c);
			}
		}
	}

}
