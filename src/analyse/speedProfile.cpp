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

void zProfile( vector < heightProbe* > & lprb,vector <double> & Zprofile, Sample& spl)
{
	unsigned int Nb = spl.lbody().size();
	unsigned int Nprb = lprb.size();
	vector <unsigned int> Nbod( lprb.size(),0);


	for (unsigned int i=0;i<Nb; ++i)
	{
		unsigned int j=0;
		while( j< Nprb )
		{
			if (  lprb[j]->containEntireBody(spl.body(i))  )
			{
				++Nbod[j];
				Zprofile[j]+=spl.body(i)->z();
				break;
			}

			if ( lprb[j]->intersection( spl.body(i)) || lprb[j]->containCenter( spl.body(i)))
			{
				++Nbod[j];
				Zprofile[j]+=spl.body(i)->z();
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
		if( Nbod[i] != 0 ) Zprofile[i]/=(double) (Nbod[i]);
		//cout<<" "<<profile[i]<<endl;

	}

}

void rotProfile( vector < heightProbe* > & lprb,vector <double> & ROT, Sample& spl)
{
	unsigned int Nb = spl.lbody().size();
	unsigned int Nprb = lprb.size();
	vector <unsigned int> Nbod( lprb.size(),0);

	cout<<"Nprb = "<<lprb.size()<<endl;

	for (unsigned int i=0;i<Nb; ++i)
	{


		if( spl.body(i)->bodyDof()!=NULL ) continue;

		unsigned int j=0;
		while( j< Nprb )
		{
			if (  lprb[j]->containEntireBody(spl.body(i))  )
			{
				++Nbod[j];
				ROT[j]+=(spl.body(i)->vrot());
				//cout<<j<<" ";
				break;
			}

			if ( lprb[j]->intersection( spl.body(i)) || lprb[j]->containCenter( spl.body(i)))
			{
				++Nbod[j];
				ROT[j]+=(spl.body(i)->vrot());
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
		//cout<<Nbod[i]<<" "<<lprb[i]->h2()<<endl;

		if(Nbod[i]==0) {
			ROT[i]=0.;

		}

		else{
			ROT[i]/=(double) (Nbod[i]);
		}

		//cout<<" "<<profile[i]<<endl;

	}


}

//Profils de fluctuations vitesses pour chaque composante x et y , temperature dans les deux directions de l'espace
void TemperatureProfile( vector < heightProbe* > & lprb,vector <double> & Xprofile,vector <double> & Yprofile, vector <double> & XYprofile, Sample& spl, System* sys_)
{
	unsigned int Nb = spl.lbody().size();
	unsigned int Nprb = lprb.size();
	vector <unsigned int> Nbod( lprb.size(),0);

	cout<<"Nprb = "<<lprb.size()<<endl;

	vector <double> Vxprofile(lprb.size(),0.);
	vector <double> Vyprofile(lprb.size(),0.);

	//Construction du profil de vitesse moyen:
	speedProfile( lprb, Vxprofile, Vyprofile, spl) ;

	for (unsigned int i=0;i<Nb; ++i)
	{
		unsigned int j=0;

		while( j< Nprb )
		{
			if (  lprb[j]->containEntireBody(spl.body(i))  )
			{
				double vxi=spl.body(i)->vx();
				double vyi=spl.body(i)->vy();
				double VX=Vxprofile[j];
				double VY=0.;//Vyprofile[j];

				++Nbod[j];
				Xprofile[j]+=(vxi-VX)*(vxi-VX);
				Yprofile[j]+=(vyi-VY)*(vyi-VY);
				XYprofile[j]+=(vxi-VX)*(vyi-VY);
				break;
			}

			if ( lprb[j]->intersection( spl.body(i)) || lprb[j]->containCenter( spl.body(i)))
			{
				double vxi=spl.body(i)->vx();
				double vyi=spl.body(i)->vy();
				double VX=Vxprofile[j];
				double VY=0.;//Vyprofile[j];
				++Nbod[j];
				Xprofile[j]+=(vxi-VX)*(vxi-VX);
				Yprofile[j]+=(vyi-VY)*(vyi-VY);
				XYprofile[j]+=(vxi-VX)*(vyi-VY);


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
		//cout<<Nbod[i]<<" "<<lprb[i]->h2()<<endl;

		if(Nbod[i]==0) {
			Xprofile[i]=0.;
			Yprofile[i]=0.;
			XYprofile[i]=0.;

		}

		else{
			Xprofile[i]/=(double) (Nbod[i]);
			Yprofile[i]/=(double) (Nbod[i]);
			XYprofile[i]/=(double) (Nbod[i]);
		}

		//cout<<" "<<profile[i]<<endl;

	}

}

