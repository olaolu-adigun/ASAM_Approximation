#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <math.h>
#include <typeinfo>

using namespace std;

#include "SAMUtils.h"

/*
Ideas: 
	- Compute MAD instead of MSE ???
*/

vector<double> xin, fx;  //should use a map or valarray container instead.
/* To Work in subdirs:*/
//string md = "md ";		//make dir on unices or windows+cygwin
//string rm = "rm -r ";	//Remove dir on unices or windows+cygwin

string md = "mkdir ";		//make dir on linux
string rm = "rm -rf ";	//Remove dir on linux
string mkcmd;
string rmcmd;
string bname ("fuzzyF");	//basename for all fuzzy approx files
string bname1 ("fuzzyV");	//basename for variance


int main(int argc, char *argv[])
{
	int nRules = 5; //5; // # of Rules. Formerly = 12;
	int epochSize = 2000; // Report min error solution after epoch # of steps

	//int adaptIters = epochSize*30;
        int adaptIters = epochSize*10;
	int defaultPrec = cout.precision();
	string iname, line;
	string filename = "InputFxn.dat";
	string errlog = "Errors.dat";

	//enum Fitfxn {gauss, cauchy, tanhyp, laplace, triangle, sinc};
	double data1=0, data2=0;
	vector<double> errors(6, 100);
	int fxnpts, k;
	
	// Set input information
        double min_x = 0.00, max_x = 1.00;
	int N = 1000;
	double step_size = (max_x - min_x) / N;
	cout << step_size << endl;
	vector<double> xin(N,1);
	vector<double> fx(N,1);
	
	ofstream errfile(errlog.data(), ios::out);
	ofstream pdfile(filename.data(), ios::out);

	if ( pdfile.fail() || errfile.fail() ){
		std::cout<<"file i/o error.\n";
		system("PAUSE");
		return EXIT_SUCCESS; 
	}
	cout<<"File Opening done \n";

	// Add Noise;
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,1.0);
	

	//Need to change precision of I/O pipes here...
	pdfile.precision(9);

	// Generate the samples
	for (int i = 0; i <= N; i++){
		xin[i] = min_x + (i*step_size);
		fx[i] = sin(xin[i]); 
		pdfile << xin[i] << "\t" << fx[i] << endl;
	}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

	cout << xin.size() << " " << fx.size() << endl;
	fxnpts = fx.size();
	InitializeAll(nRules, (int) (0.5*fxnpts), (int) fxnpts);
	InitializeFxn(xin, fx);
	std::cout << endl << "Number of Testing Points: "<< ::des.size() <<"\t" << ::NUMDES;
	std::cout << endl << "Number of Training Samples: "<< ::sample.size() <<"\t" << ::NUMSAM ;
	std::cout << endl << "Number of Rules: "<< ::NUMPAT ;
	std::cout << endl << "Range of x: " << xin.front() << "<-->" << xin.back() << endl;
	std::cout << endl << "Epoch Size: "<< epochSize << endl;

	//(Make) Dirs for each fit fxn.
	// Reset Record by removing dirs.
	errfile << "Iter# ";
	for(int t = 0; t < 6; t++){ // t<5 omits Sinc SAM
		rmcmd = rm + name[t];
		mkcmd = md + name[t];	
		system(rmcmd.data());  //Reset record for new runs.
		system(mkcmd.data()); 
		errfile << "\t    " << name[t];
	}
	errfile << endl; cout << "\n \n" << endl;

	k = 0; 
	vector<double>::const_iterator i;
	double minerr; int loc; bool minQ;
	ASAMsInitialize();
    	while( /*vecmin(errors) > 0.0001*/  k < adaptIters ) { //Error Criterion or Iteration Limit
		ASAMsLearn();
		if (k%epochSize == 0){
			errors = ASAMsApprox();
			minerr = vecmin(errors);

			WriteEpoch(bname, k );

			errfile << k;
			minQ = false; loc = 0;
			for ( i = errors.begin(); i < errors.end() ; i++ ){ //Log MSEs & Locate Min.
				errfile << "\t" << *i ;
				if ( (!minQ) && (*i != minerr) ) loc++;
				else minQ = true;
			}
			errfile << endl;
			cout << "iter# " << k << ": Min. Error = " << minerr 
				<< " using " << name[loc] << " fit function." << endl;			
		}
		k++;
	}
	
	cout << endl;
	//fxnfile.close();
	errfile.close();
	pdfile.close();	
	//system("PAUSE");
    return EXIT_SUCCESS;
}
