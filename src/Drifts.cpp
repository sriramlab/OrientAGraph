/*
 * Drifts.cpp
 *
 *  Created on: Jul 18, 2012
 *      Author: pickrell
 */





#include "CountData.h"
#include "PhyloPop_params.h"
#include "CmdLine.h"


void printopts(){
	cout << "drifts v. 0.1\n";
	cout << "Options:\n";
	cout << "-i [file name] input file\n";
	cout << "-k [int] number of SNPs per block for estimation of standard errors (1)\n";
	cout << "\n";
}

string infile;
int main(int argc, char *argv[])
{
	PhyloPop_params p;
    CCmdLine cmdline;
	if (cmdline.SplitLine(argc, argv) < 1){
		printopts();
		exit(1);
	}
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
     	printopts();
     	exit(1);
    }
    if (cmdline.HasSwitch("-k"))	p.window_size = atoi(cmdline.GetArgument("-k", 0).c_str());
    p.alfreq_scaling = 4;
	CountData counts(infile.c_str(), &p);

	for (int i = 0; i < counts.npop; i++){
		for (int j = 0; j < counts.npop; j++){
			if (i == j) continue;
			pair<double, double> drift = counts.calculate_drift(i, j);
			pair<double, double> m1 = counts.calculate_mean(i);
			pair<double, double> m2 = counts.calculate_mean(j);
			cout << counts.id2pop[i] << " "<< counts.id2pop[j] << " " << drift.first << " "<< drift.second << " "<< m1.first << " "<< m1.second << " "<< m2.first << " "<< m2.second << "\n";
		}
	}
	return 0;
}
