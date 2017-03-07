// 2DASAM-v2 main function
// Linear samples version
// 17-Apr-2012
// Osonde Osoba

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Eigen/Core"

using namespace std;
using namespace Eigen;

#include "ASAM.h"
#include "GaussianASAM.h"
#include "SincASAM.h"
#include "CauchyASAM.h"

int main(int argc, char *argv[]){

	int numRules = 39;
	int epochSize = 500; int adaptIters = epochSize*20;
	matrix fxn;
	VectorXd xin, yin, fxyin;
	const int fnums = 3; int temp;

	string line=" ";
	istringstream buf(line);

	system("cd");

	string rd, mk;
	string mkcmd = "mkdir -p ";
	string rmcmd = "rm -r ";

	string fnames[3] = {"FxnGen2D.dat", "OutFxn.dat", "Errors.dat"};

        std::fstream fxnio[3];
        fxnio[0].open(fnames[0].data(), ios::in);
        fxnio[1].open(fnames[1].data(), ios::out);
        fxnio[2].open(fnames[2].data(), ios::out);


	fxnio[1].precision(9);fxnio[2].precision(9);

	if ( fxnio[0].fail() || fxnio[1].fail() || fxnio[2].fail()){
		cout<<"file i/o error.\n";
		//system("PAUSE");  This system command does not exist in linux. You can achieve 
		//		    same fanctionality with the linux command read 
		return EXIT_FAILURE; 
	}

	// Determine dynamic matrix size
	int nrows=0;
	while(getline(fxnio[0], line)) nrows++;	

	fxn.resize(nrows, 3);
	fxn.fill(0);
	cout << "\n# of Input Data Points= " << fxn.rows() << endl;

	// Read in function values
	fxnio[0].clear(); fxnio[0].seekg(0);
	int r=0, c=0;
	double tmp; 
	while(getline(fxnio[0], line)){
		buf.clear(); buf.str(line);
		while(buf >> tmp){
			fxn(r,c)=tmp; 
			fxnio[1] << tmp << "\t";
			c++;
		}
		fxnio[1] << endl;
		r++; c = 0;
	}	

////////////////////////////////////////////////////////////////////////////

	// Set input information
	int N = 1000;
	double min_x = 0.00, max_x = 6.414;	
	double step_size_x = (max_x - min_x) / N;

	double min_y = 0.00, max_y = 6.414;	
	double step_size_y = (max_y - min_y) / N;
	
	vector<double> Xin(N,1);
	vector<double> Yin(N,1);
	vector<double> Fxyin(N,1);

	///////// DEFINE PARAMETER FOR ADDED NOISE
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
    	std::mt19937 gen(rd());
	// Define the parameter for Uniform noise
	double a = 0.0 , b = 0.1;
	uniform_real_distribution<double> uniform(a,b);
	// Define the parameter for the Cauchy noise
	double m = 0.0, d = 0.001;
	cauchy_distribution<double> cauchy(m, d);
	// Define the paramter for Gaussian noise
	double mu = 0.0, sigma = 0.01;
	normal_distribution<double> gaussian(mu, sigma);
	
	int noise_status = 0;  // Noise status is 0 or 1



///////////////////////////////////////////////////////////////////////////

	xin = (fxn.col(0)); vector<double> xv = eigvec2stdvec(xin) ;	
	yin = (fxn.col(1)); vector<double> yv = eigvec2stdvec(yin) ;	
	fxyin = (fxn.col(2)); vector<double>fxyv = eigvec2stdvec(fxyin) ;	

	//...ASAM(xvals, yvals, fxyvals, _numpat, _numsam, _numdes)
	//ASAM* sam =  new CauchyASAM( xv, yv, fxyv, numRules, (int) (0.5*xin.size()), xin.size() );
	//delete sam;

	GaussianASAM gsam( xv, yv, fxyv, numRules, (int) (0.5*xin.size()), xin.size() ); 
	SincASAM ssam( xv, yv, fxyv , numRules, (int) (0.5*xin.size()), xin.size() );
	CauchyASAM csam( xv, yv, fxyv , numRules, (int) (0.5*xin.size()), xin.size() );

	//(Make) Dirs for each fit fxn.
	// Reset Record by removing dirs.
	string name[] = {"Gauss", "Sinc", "Cauchy"};
	fxnio[2] << "Iter# ";
	for(int t = 0; t < 3 ; t++){
		rd = rmcmd + name[t];
		mk = mkcmd + name[t];
		system(rd.data());  //Reset record for new runs.
		system(mk.data()); 
		fxnio[2] << "\t    " << name[t];
	}
	fxnio[2] << endl;
	cout << "Epoch Size: "<< epochSize << endl;

	int k = 0; vector<double> errors;
	vector<double>::const_iterator i;
	double minerr; int loc; bool minQ;
	do {
		gsam.Learn(); 
		ssam.Learn();
		csam.Learn();
		if (k% epochSize == 0){	
			errors.clear();
			errors.push_back(gsam.Approx()); errors.push_back(ssam.Approx()); errors.push_back(csam.Approx());
			gsam.WriteEpoch(k); ssam.WriteEpoch(k); csam.WriteEpoch(k);

			minerr = *(std::min_element(errors.begin(), errors.end()));			
			fxnio[2] << k ;
			minQ = false; loc = 0;
			for ( i = errors.begin(); i < errors.end() ; i++ ){ //Log MSEs & Locate Min.
				fxnio[2] << "\t" << *i ;
				if ( (!minQ) && (*i != minerr) ) loc++;
				else minQ = true;
			}
			fxnio[2]<< endl;
			cout << "iter# " << k << ": Min. Error = " << minerr 
				<< " using " << name[loc] << " fit function." << endl;
		}
		k++;
	} while( k <= adaptIters || minerr < 1e-7 );

	for(temp = 0; temp < fnums ; temp++) fxnio[temp].close();	
	std::cout << "\nDone. Final min MSE = " << errors[loc] << endl;
	//system("PAUSE");  This system command does not exist in linux. You can achieve same fanctionality with the linux command read
	return EXIT_SUCCESS;
}
