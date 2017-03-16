#ifndef _pointSet_hpp
#define _pointSet_hpp

#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "math.h"
#include <algorithm>
#include "point.hpp"

using namespace std;


class pointSet
{

private:

	
	string name_;
	vector <dtk::point > pset_;
	double xmean_,xvar_;
	double ymean_,yvar_;
	double xmin_,xmax_;
	double ymin_,ymax_;
	double xsum_,ysum_;
	bool periodic_;

public:

	pointSet( int Nbin ): pset_(Nbin),periodic_(false) { 	}
	
	pointSet():periodic_(false) { };
	
	void add(double, double);
	void read(istream &);
	void addRead( istream &);
	
	double & px( unsigned int i )       { return this->pset_[i].x(); }
	double   px( unsigned int i ) const { return this->pset_[i].x(); }
	double & py( unsigned int i )       { return this->pset_[i].y(); }
	double   py( unsigned int i ) const { return this->pset_[i].y(); }
	
	vector <dtk::point>   pset() const { return this->pset_;}
	vector <dtk::point> & pset()       { return this->pset_;}
	
	double & ysum(  )       { return this->ysum_; }
	double   ysum(  ) const { return this->ysum_; }
	
	double & ymean(  )       { return this->ymean_; }
	double   ymean(  ) const { return this->ymean_; }
	
	double & xmean(  )       { return this->xmean_; }
	double   xmean(  ) const { return this->xmean_; }
	
	bool & periodic(  )       { return this->periodic_; }
	bool   periodic(  ) const { return this->periodic_; }
	
	void upXFilter( double );//filtre les valeurs superieures a l'entree
	void downXFilter( double );//filtre les valeurs inferieure a l'entree
	
	
	void increasingXSort();
	void increasingYSort();
	void decreasingXSort();
	void decreasingYSort();
	
	void Mean();
	void Variance();
	void extractValues();
	void xNormalise( double );
	void yNormalise( double );
	void xoffSet( double );
	void yoffSet( double );
	void Max();
	void Min();
	void yCumul();
	void yAbs();
	
	void circularise(double,int);// duplique les donnes en graphe polaire
	void closeCircular();//duplique le premier point apres le dernier pour fermer la courbe
	
	
	pointSet histoBin( unsigned int );//Nbin
	pointSet histoNumberBin( unsigned int );
	pointSet histoBinVar( unsigned int);
	pointSet mobileMean(unsigned int , unsigned int );
	pointSet mobileMean2(unsigned int , unsigned int );
	pointSet slidingHisto( unsigned int , double );
	
	double integral(void);
	

	//void histoEvents( unsigned int );//Nevents per Bin
	//Si le y est nul, on calcule le pdf des x;
	
	void  write (const char * file);
	void  write (std::string file);

	dtk::point   pset( unsigned int i )  const { return this->pset_[i];}

	unsigned int psetSize() const { return this->pset_.size();}
	
	void psetClear()  { this->pset_.clear();}
	
	struct DecreasingXOrder
	{
		bool operator () (const dtk::point & a1, const dtk::point & a2) const 
	      { 
	        return (a1.x() > a2.x());
	      }

	      bool operator () (const dtk::point * a1, const dtk::point * a2) const 
	      { 
	        return (a1->x() > a2->x());
	      }
	};
	struct DecreasingYOrder
	{
		bool operator () (const dtk::point & a1, const dtk::point & a2) const 
	      { 
	        return (a1.y() > a2.y());
	      }

	      bool operator () (const dtk::point * a1, const dtk::point * a2) const 
	      { 
	        return (a1->y() > a2->y());
	      }
	};
	struct IncreasingXOrder
	{
		bool operator () (const dtk::point & a1, const dtk::point & a2) const 
	      { 
	        return (a1.x() < a2.x());
	      }

	      bool operator () (const dtk::point * a1, const dtk::point * a2) const 
	      { 
	        return (a1->x() < a2->x());
	      }
	};
	struct IncreasingYOrder
	{
		bool operator () (const dtk::point & a1, const dtk::point & a2) const 
	      { 
	        return (a1.y() < a2.y());
	      }

	      bool operator () (const dtk::point * a1, const dtk::point * a2) const 
	      { 
	        return (a1->y() < a2->y());
	      }
	};
	
	unsigned int Weight( double u )
	{
		return ( (unsigned int) floor ( 10*u) +1 );
	}
	
	pointSet operator + (pointSet a)
	{
		pointSet temp;

		temp.pset()=this->pset();

		for(unsigned int i=0;i<a.psetSize();++i)
		temp.add(a.px(i),a.py(i));

		return temp;
	}

};

#endif // _pointSet_hpp
