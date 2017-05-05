#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;


unsigned int getcolumnsnumber(const string line)
{
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

class Dataset{
	private:
		std::vector<double> values_;
		double mean_;
		double variance_;
		double min_, max_;
	public:
		Dataset();
		~Dataset();
		void add(double);
		void write(ostream&);
		void printvalues();
		void read(const char *, unsigned int);
		void readstress(const char * );
		void extractValues();
		void normalize(double);
		double mean() const {return mean_;}
		double min() const {return min_;}
		double max() const {return max_;}
		double variance() const {return variance_;}

};

Dataset::Dataset(){
	mean_ =  0. ;
	min_ = 0. ;
	max_ = 0. ;
	variance_ = 0. ;
	values_.clear();
}

Dataset::~Dataset(){
	values_.clear();
}

void Dataset::add(double value)
{
	values_.push_back(value);
}

void Dataset::extractValues()
{
	if(values_.size()==0) return;

	unsigned int N = values_.size();
	double min = values_[0];
	double max = values_[0];
	double mean = 0. ;
	double meancarre = 0. ;

	for(unsigned int i = 0 ; i != values_.size(); i ++)
	{
		(values_[i] < min) ? (min_ = values_[i]) : min_ = min ; 
		(values_[i] > max) ? (max_ = values_[i]) : max_ = max ; 
		mean += values_[i];
		meancarre += values_[i] * values_[i];
	}

	mean_ = mean / (double) N ;
	variance_ = meancarre / (double) N - mean_ * mean_ ;
}

void Dataset::read(const char * fname, unsigned int col )
{
	//Detect the max number of columns et check si le nombre demande est bon
	double value;
	string line;
	ifstream is(fname);
	getline(is,line);
	int nc = getcolumnsnumber(line);

	if(col > nc ) { cerr<<"Mauvais choix de colonne ! "<<endl; return ; }

	double donnees [ nc ] ;
	//Lire directement les donnees dans le tableau col
	is.close();
	is.open(fname);
	values_.clear();

	ofstream test("reecriture.txt");
	if(is.is_open()){

		while(is)
		{
			for(unsigned int i = 0 ; i < nc; i++)
			{
				is >> value;
				donnees[i]=value;
			}
			
			if(is.eof()) break;
			values_.push_back(donnees[col]);

		}

		is.close();
	}
	test.close();
	cout<<"Donnees extraites du fichier "<<fname<<", colonne "<<col<<" sur "<<nc<<", "<<values_.size()<<" valeurs."<<endl;
}

void Dataset::readstress(const char * fname)
{
	//Detect the max number of columns et check si le nombre demande est bon
	double value;
	string line;
	ifstream is(fname);
	getline(is,line);
	int nc = getcolumnsnumber(line);

	double donnees [ nc ] ;
	//Lire directement les donnees dans le tableau col
	is.close();
	is.open(fname);
	values_.clear();

	if(is.is_open()){

		while(is)
		{
			for(unsigned int i = 0 ; i < nc; i++)
			{
				is >> value;
				donnees[i]=value;
			}
			
			if(is.eof()) break;
			values_.push_back(donnees[8]/donnees[9]);

		}

		is.close();
	}
}
void Dataset::printvalues()
{
	for(std::vector<double>::iterator it = values_.begin() ; it != values_.end();it++)
	{
		cout<<*it<<endl;
	}
}


class Pdf{
	private:

		unsigned int Nbins_;
		double min_, max_;
		double * histogram_ ;
		double dx_;
	public:
		Pdf();
		~Pdf();
		double Nbins() const {return Nbins_;}
};


int main(){
	Dataset stress;
	Dataset hauteur;
	char output[100] = "stress.tx";
	for(unsigned int i = 4 ; i < 5 ; i++)
	{
		char fname[100];
		char fname2[100];
		//sprintf(fname,"%.2dbidon.txt",i);
		sprintf(fname,"%.2d/Analyse/stress.txt",i);
		sprintf(fname2,"%.2d/Analyse/system.txt",i);
		stress.readstress(fname);
		stress.printvalues();
	}
	cerr<<"Resultats enregistres dans le fichier "<<output<<endl;
	return 0;
}








