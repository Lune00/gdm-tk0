#ifndef _dataSet_hpp
#define _dataSet_hpp


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "math.h"
#include <algorithm>
#include "pointSet.hpp"

using namespace std;

class pointSet;

class DataSet
{

private:

	string name_;
	vector <double> set_; // ... a voir
	double max_,min_,mean_,variance_/*,meandev_*/;

public:

	DataSet();
	DataSet(const char * name);
	DataSet(string name);
	~DataSet();

	void add(double);
	void write(ostream&, unsigned int n);
	void write(char * );
	void write(const char *);

	void read(istream&);
	void addRead(istream&);

	void read(const char * fname); 


// Extract some informations
	double Mean();
	double Min();
	double Max();
	double Variance();
	void   extractValues();
//particular mean
	//double Mean(unsigned int);

// Modify the set
	void Normalize(double ref);
	void Logarithmize();
	void Sort();
	void DecreasingSort();

// Build new set
	DataSet * PDF(unsigned int );
	DataSet * RuningAverage(DataSet&,unsigned int nb);
	DataSet * Derivative(DataSet&,unsigned int,unsigned int);
	DataSet * Period(unsigned int nb);

	pointSet Rich_PDF(unsigned int);//Number of class of same number of event
	pointSet Rich_PDF2(unsigned int );//Number of event per class of same number of events
	pointSet slidingPdf(unsigned int );//sliding window
	pointSet kernelPdf(unsigned int,double);
	void upFilter( double );//filtre les valeurs superieures a l'entree
	void downFilter( double );//filtre les valeurs inferieure a l'entree

	double linearCorrelation( DataSet *);

	void meanCentering();
	
	void exportDat( );

	double  set( unsigned int i )  const { return set_[i];}
	vector<double>  set( )   { return set_;}

	double & max()       { return max_; }
	double   max() const { return max_; }
	double & min()       { return min_; }
	double   min() const { return min_; }
	double & mean()       { return mean_; }
	double   mean() const { return mean_; }
	double & variance()       { return variance_; }
	double   variance() const { return variance_; }

	/*double & meandev()		{return meandev_;}
	double	 meandev() const {return meandev_;}*/
	unsigned int setSize() const { return set_.size();}
	void setClear()  { set_.clear();}

	struct DecreasingOrder
	{
		bool operator () (const double & a,const double & b) const
			{return (a > b) ;}
	};

	double rectangular( double u)//u compris en -.5 et .5
	{
		if( fabs(u)<=1 ) return .5;
		else return 0;
	}
	double triangular( double u )// u compris entre -1 et 1
	{
		if( fabs(u)<= 1. ) return (1. - fabs(u));
		else return 0;
	}

	DataSet::DataSet operator + (DataSet a)
	{
		DataSet temp;

		temp.set()=this->set();

		for(unsigned int i=0;i<a.setSize();++i)
			temp.add(a.set(i));

		return temp;
	}

};

#endif // _dataSet_hpp
