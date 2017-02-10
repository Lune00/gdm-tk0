#include "speedProfile.hpp"


// For disk only
void speedProfile( vector < heightProbe* > & lprb,vector <double> & Xprofile,vector <double> & Yprofile, Sample& spl)
{
	unsigned int Nb = spl.lbody().size();
	unsigned int Nprb = lprb.size();
	vector <unsigned int> Nbod( lprb.size(),0);
		
	cout<<"Nprb = "<<lprb.size()<<endl;
	
	
	for (unsigned int i=0;i<Nb; ++i)
	{
		unsigned int j=0;
		while( j< Nprb )
		{
			if (  lprb[j]->containEntireBody(spl.body(i))  )
				{
					++Nbod[j];
					Xprofile[j]+=spl.body(i)->vx();
					Yprofile[j]+=spl.body(i)->vy();
					//cout<<j<<" ";
					break;
				}
				
			if ( lprb[j]->intersection( spl.body(i)) || lprb[j]->containCenter( spl.body(i)))
				{
					++Nbod[j];
					Xprofile[j]+=spl.body(i)->vx();
					Yprofile[j]+=spl.body(i)->vy();
					j++;
					//cout<<j<<" ";
				}
			else
				{
				++j;
				//cout<<j<<" ";
				}
				
			
		}
		
		//cout<<endl;
		//getchar();

	}
	
	for (unsigned int i=0;i<Nprb; ++i)
	{
		//cout<<profile[i];
		Xprofile[i]/=(double) (Nbod[i]);
		Yprofile[i]/=(double) (Nbod[i]);
		//cout<<" "<<profile[i]<<endl;

	}
	
}

