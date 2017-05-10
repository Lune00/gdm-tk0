#include <math.h>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>

//g++ -lm -lfftw3  fftwexample.cpp && ./a.out
using namespace std;

unsigned int getcolumnsnumber(const char * fname)
{
	string line;
	ifstream is(fname);
	getline(is,line);
	is.close();

	std::size_t pos = line.find(' ');
	std::size_t initpos = 0 ;
	unsigned int nc = 1 ;
	do
	{
		initpos = pos + 1 ;
		pos = line.find(' ',initpos);
		nc++;
	}while(pos != std::string::npos);

	return nc ;
}

int main(){


	system("mkdir -p FFT");
	for(unsigned int j = 4 ; j < 8 ; j++)
	{
		char inputfile[100];
		char outputfile[100];
		sprintf(inputfile,"%.2d/Analyse/stress.txt",j);
		sprintf(outputfile,"FFT/stressFT%.2d.txt",j);

		double value;

		int nc = getcolumnsnumber(inputfile);
		double donnees [nc];
		std::vector<double> values;
		ifstream is(inputfile);
		if(is.is_open())
		{
			while(is)
			{
				for(unsigned int i = 0 ; i < nc; i++)
				{
					is >> value;
					donnees[i]=value;
				}

				if(is.eof()) break;
				double mu = donnees[7]/donnees[8];
				values.push_back(mu);
			}
			is.close();
		}

		cout<<values.size()<<endl;

		//On recupere la taille
		int n = values.size();
		//int n = 128 ;
		//On rentre les valeurs dans un tableau 
		fftw_complex in[n];
		fftw_complex out[n];
		fftw_plan p;
		//On fait un test d'abord
		for(unsigned int i = 0 ; i < n ; i++)
		{
			in[i][0] = values[i];// cos(3 * 2*M_PI*i/n);
			in[i][1] = 0. ;
		}
		p = fftw_plan_dft_1d(n,in,out,FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);

		ofstream outfft(outputfile);
		for(unsigned int i = 0 ; i < n/2+1 ; i++)
		{
			outfft<<i<<" "<<out[i][0]<<" "<<out[i][1]<<" "<<in[i][0]<<" "<<in[i][1]<<endl;
		}
		outfft.close();

		fftw_destroy_plan(p);
	}
	return 0;
}
