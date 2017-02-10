#include "grid.hpp"

//Place disk or polygon on a square grid with ratio=height/width
void squareGrid(Sample & spl, double ratio)
{
	spl.radiusExtrema(0);
	unsigned int nBod = spl.lbody().size();
	double rmax = spl.rmax();
	double lsquare = sqrt((double) (nBod)/ratio);
	//double ratio1 = 1./ratio;
	unsigned int nw = (unsigned int) floor(lsquare);
	cout<<nw<<endl;
	
	unsigned int jh = 0;
	unsigned int i = 0;
	double offset=0;
	double dr =0;
	while ( i < nBod)
	{
		unsigned int j=0;
		while(j< nw && i< nBod)
		{
			offset =  (((double)rand()/(double)RAND_MAX)-1.0)* .1*rmax;
			dr = (((double)rand()/(double)RAND_MAX))*.15*rmax;
			//cout<<offset<<endl;
			spl.body(i)->x()=(2.2+offset)*j*rmax;
			spl.body(i)->y()=(2.2+offset)*jh*rmax;
			//dynamic_cast<disk*>(spl.body(i))->R()=rmax+dr;
			++j;
			++i;
		}
		++jh;
	}
	
}

