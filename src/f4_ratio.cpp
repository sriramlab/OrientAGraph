/*
 * f4_ratio.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: pickrell
 */



#include "GraphState2.h"
#include "PhyloPop_params.h"
#include "CmdLine.h"


void printopts(){
	cout << "\nf4ratio v. 0.1\n";
	cout << "Options:\n";
	cout << "-i [file name] input file\n";
	cout << "-k [int] number of SNPs per block for estimation of standard errors (1)\n";
	cout << "-d [file name] populations for denominator of f4 ratio\n";
	cout << "-n [file name] populations to estimate\n";
	cout << "\n";
}

string infile;
string denomfile;
string numfile;
int main(int argc, char *argv[])
{
	PhyloPop_params p;
	p.alfreq_scaling = 4;
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

    if (cmdline.HasSwitch("-d")) denomfile = cmdline.GetArgument("-d", 0).c_str();
    else{
     	printopts();
     	exit(1);
    }
    if (cmdline.HasSwitch("-n")) numfile = cmdline.GetArgument("-n", 0).c_str();
     else{
      	printopts();
      	exit(1);
     }

	CountData counts(infile.c_str(), &p);
	string pop1;
	string pop2;
	string pop3;
	string pop4;
	ifstream denomin(denomfile.c_str());
	ifstream numin(numfile.c_str());
	string st, buf;
	vector<string> line;

	getline(denomin, st);
	buf.clear();
	stringstream ss(st);
	line.clear();
	while (ss>> buf){
		line.push_back(buf);
	}
	pop1 = line[0];


	getline(denomin, st);
	buf.clear();
	stringstream ss1(st);
	line.clear();
	while (ss1>> buf){
		line.push_back(buf);
	}
	pop2 = line[0];

	getline(denomin, st);
	buf.clear();
	stringstream ss2(st);
	line.clear();
	while (ss2>> buf){
		line.push_back(buf);
	}
	pop3 = line[0];

	getline(denomin, st);
	buf.clear();
	stringstream ss3(st);
	line.clear();
	while (ss3>> buf){
		line.push_back(buf);
	}
	pop4 = line[0];
	//cout << pop1 << " "<< pop2 << " "<< pop3 << " "<< pop4 << "\n";
	while (getline(numin, st)){
		buf.clear();
		stringstream ss4(st);
		line.clear();
		while (ss4>> buf){
			line.push_back(buf);
		}
		string testpop = line[0];
		//cout << testpop << "\n";
		pair<double, double> tmp = counts.f4ratio(pop1, pop2, pop3, pop4, testpop);
		cout << "[["<< pop1 << "," << pop2 << "],[" << pop3 << ","<< pop4 << "]] "<< testpop <<  " "<< tmp.first << " "<< tmp.second << "\n";
	}

	return 0;
}

