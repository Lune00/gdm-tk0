#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "biaxial_A.hpp"
#include "sample.hpp"
#include "dof.hpp"

					
//-----------------------------------------------------modif le 14/04/09

int main (int argc, char * argv[])
{
		Sample * theSpl = new Sample();
		string token;
		vector<dof *>ldof;
		if (argv[1]!=0)
		{
			ifstream input(argv[1]);
			
//------------------------------------Contribution de Charles
			if( ! input)
		{
			cerr << "Echec overture du fichier spl " << argv[1] << endl;
			
		}
		else
		{		
				input>>token;
				if( token=="Sample{")
				{
					ldof = theSpl->read(input);
				}
				else
				{
					cout<<"Echec ......."<<endl;
					exit(0);
				}
				
//------------------------------------------------------------			
				dof * B;
			
				for (unsigned int i=0 ; i<ldof.size() ; ++i)
				{
					B = ldof[i];
					//cout<<i<<" "<<B->lctrlBodies().size()<<endl;
				
					double dx = (B->ctrlBody(0)->x()) - (B->ctrlBody(1)->x());
					double dy = (B->ctrlBody(0)->y()) - (B->ctrlBody(1)->y());
					double d= dx * dx + dy * dy;
					d=sqrt(d);
				
					disk * D=dynamic_cast<disk *>(B->ctrlBody(0));
				
					double r = D->R();
			
					double aire = 0;
					double eta = d/(2*r);
					double alpha = 2*asin((sqrt(3)/2.)*(sqrt(1.-eta*eta)-eta/sqrt(3)));

			
					if ( eta>=sqrt(3)/2 && eta<=1 ) //vide interne considéré comme une phase solide
					{
						aire =  3*M_PI*r*r-6*r*r*(acos(eta)-eta*sqrt(1-eta*eta)) + r*r*(eta*eta*sqrt(3) - M_PI/2 +3*(acos(eta)-eta*sqrt(1-eta*eta)));
					}
	
					else if (eta >=0 && eta< .5*sqrt(3)) // 
					{
						aire = 3*M_PI*r*r-6*r*r*(acos(eta)-eta*sqrt(1-eta*eta)) + (3*sqrt(3)*r*r*0.5)*(sqrt(1-eta*eta)-eta/sqrt(3)) * (0.5-cos(alpha/2)) +3*r*r*alpha*0.5; 
					}

					//cout<<i<<" "<<B->Area()<<endl;
	
					B ->Area()=aire;
			
					//cout<<i<<"apres "<<B->Area()<<endl;
					
					
//-----------------------------------------------------------------------------------------------------------------le 16/04/09
				}
			}
	}
	ofstream datafile(argv[1]);
	if(!datafile)
	{
		cerr << "@save_data, cannot open file " << argv[1] << endl;
	}

	datafile << "Sample{" << endl;
	
	theSpl->write(datafile);
	datafile << "}" << endl << endl;

		return 0;
}
