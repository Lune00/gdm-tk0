void shearP_CD_A::fractalDimension( ) 
{
	cout<<" fractal Dimension : "<<endl;
	double ll= sys_->spl()->xmin() + 1.*sys_->spl()->rmax();
	double rl= sys_->spl()->xmax() - 1.*sys_->spl()->rmax();
	double bl= totalProbe_.h1() + sys_->spl()->rmax();
	double tl= totalProbe_.h2() - sys_->spl()->rmax();

	double sizeprobe;
	double rmin=sys_->spl()->rmin();
	double sf;
	long int seed = 1234324321;
	pointSet solid;
	DataSet solidTemp;

	removeRattlers();
	//return;
	circularProbe test(.5*(rl+ll),.5*(tl+bl),.5*(rl-ll) );
	sf=solidFraction( test, *(sys_->spl()), *(sys_->nwk()) );
	cout<< test.x()<<" "<<test.y()<<" "<<test.R()<<endl;
	cout<<" SF macro rond "<<sf<<endl;
//	return;

	for ( unsigned int i=1; i< 2000;i+=10)
	{
		sizeprobe= rmin*(double) (i);

		unsigned int Nmes=0;
		double x,y;
		while( Nmes < 2000)
		{
			x=ll+sizeprobe+ ran1(&seed)*( (rl-ll) - 2*sizeprobe);
			y=bl+sizeprobe+ ran1(&seed)*( (tl-bl) - 2*sizeprobe);

			circularProbe temp( x , y , sizeprobe);

			if( temp.x()+temp.R()>rl) {/*cout<<"over x"<<endl; cout<<Nmes<<endl;*/Nmes++;continue; }
			if( temp.x()-temp.R()<ll) {/*cout<<"under x"<<endl;cout<<Nmes<<endl;*/Nmes++;continue; }
			if( temp.y()+temp.R()>tl) {/*cout<<"over y"<<endl; cout<<Nmes<<endl;*/Nmes++;continue; }
			if( temp.y()-temp.R()<bl) {/*cout<<"under x"<<endl;cout<<Nmes<<endl;*/Nmes++;continue; }

			sf= solidFraction(temp, *(sys_->spl()), *(sys_->nwk()) );
		// cout<<ix<<" "<<iy<<" "<<sf<<endl;
			solidTemp.add(sf);
			Nmes++;
		}
		solidTemp.extractValues();
		if( solidTemp.mean() != 0. )
		{
			cout<<"i= "<<i<<" "<<solidTemp.mean()<<" "<<sizeprobe/rmin<<" "<<test.R()/rmin<<endl;

			solid.add(sizeprobe,solidTemp.mean() );
		}

		solidTemp.setClear();
	}

	solid.write("fractal.txt");

	SF();
}

void shearP_CD_A::BboxCounting( )
{
	cout<<" fractal Dimension ( of bodies ) : BOX COUNTING "<<endl;
	double ll= sys_->spl()->xmin() + 1.*sys_->spl()->rmax();
	double rl= sys_->spl()->xmax() - 1.*sys_->spl()->rmax();
	double bl= totalProbe_.h1() ;//;+ sys_->spl()->rmax();
	double tl= totalProbe_.h2() ;//- sys_->spl()->rmax();

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int Nb= sys_->spl()->lbody().size();

	double h= min(rl-ll,tl-bl);//size of original  box
	double amp,amp2;//size of current circuler box
	unsigned int Ndiv=2;// Number of division of initial lenght
	unsigned int Nbox;//Number of intersected boxes
	unsigned int Npoint=0;
	pointSet frac;

	vector <body2d *> inside;
	Ndiv=(unsigned int ) floor( h/(1.1*rmax));

	amp=h/(double) (Ndiv);
	amp2=.5*amp;
	bool add=true;

//	removeRattlers();

	while( amp > .1* rmin && Npoint < 30)
	{
		Sample splTemp;
		//cout<<"Ndiv= "<<Ndiv<<" amp= "<<amp<<endl;
		Nbox=0;
		cout<<" Ndiv= "<<Ndiv<<" carre= "<<Ndiv*Ndiv<<endl<<"	 amp/rmax = "<<amp/rmax<<" amp/rmin="<<amp/rmin<<flush;
		for( unsigned int i = 0; i< Ndiv ; ++i)
		{
			heightProbe dirtyDetect( bl +(.5+i)*amp - amp2, bl+(.5+i)*amp + amp2 );
			//cout<<" dirty prb "<<dirtyDetect.h1()<<" "<<dirtyDetect.h2()<<endl;
			add=true;
			inside.clear();


			for( unsigned int k=0; k < Nb; ++k )
			{
				if ( sys_->spl()->body(k)->sizeVerlet()>= amp)
					if( dirtyDetect.intersection( sys_->spl()->body(k)) || dirtyDetect.containCenter( sys_->spl()->body(k)))

					inside.push_back( sys_->spl()->body(k));
			}
		//	cout<<"-----size "<<inside.size()<<endl;


			for( unsigned int j=0; j< Ndiv; ++j)
			{
				if(add && j>500 || i>500) add=false;

				circularProbe temp( ll+ (.5+ j) *amp, bl+(.5+i)*amp, amp2) ;
			//	cout<<"     circular = "<<temp.x()<<" "<<temp.y()<<endl;

				for( unsigned int k=0; k < inside.size(); ++k )
				{
					if(  /*temp.containCenter( inside[k]) || temp.probeInsideBody(inside[k]) ||*/ temp.intersection( inside[k])  )
					{
						//cout<<i<<" "<<j<<" "<<k<<endl;

						if (add) splTemp.addBody( new disk( temp.x(), temp.y(), temp.R()) );
						Nbox++;
						break;
					}
				}

			}

		}
		//cout<<" export reconstruction "<<endl;
		ofstream out("spltemp.his",ios::out);
		out<<"Sample{"<<endl;
		splTemp.write( out);
		out<<"}"<<endl;

		cout<<" Nbox/Ncell "<<(double)(Nbox)/(double)(Ndiv*Ndiv)<<" Nbox= "<<Nbox<<" Npoint= "<<Npoint<<endl;
		//if( Ndiv>1 )
		if( Nbox != 0.) frac.add(  log10(amp), log10((double) (Nbox) ));
		Ndiv+=1000;Npoint++;
		amp=h/(double) (Ndiv);
		amp2=.5*amp;
		frac.write("FBCspl.txt");
	}

}

void shearP_CD_A::NboxCounting( )
{
	cout<<" fractal Dimension ( of contact network ) : BOX COUNTING "<<endl;
	double ll= sys_->spl()->xmin() + 1.*sys_->spl()->rmax();
	double rl= sys_->spl()->xmax() - 1.*sys_->spl()->rmax();
	double bl= totalProbe_.h1() ;//;+ sys_->spl()->rmax();
	double tl= totalProbe_.h2() ;//- sys_->spl()->rmax();

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int Nc= sys_->nwk()->clist().size();

	double h= min(rl-ll,tl-bl);//size of original  box
	double amp,amp2;//size of current circuler box
	unsigned int Ndiv=2;// Number of division of initial lenght
	unsigned int Nbox;//Number of intersected boxes
	unsigned int Npoint=0;
	pointSet frac;

	vector <inter2d *> inside;
	Ndiv=(unsigned int ) floor( h/(1.1*rmax));

	amp=h/(double) (Ndiv);
	amp2=.5*amp;
	bool add=true;

//	removeRattlers();

	while( amp > .1* rmin && Npoint < 300)
	{
		Sample splTemp;
		//cout<<"Ndiv= "<<Ndiv<<" amp= "<<amp<<endl;
		Nbox=0;
		cout<<" Ndiv= "<<Ndiv<<" carre= "<<Ndiv*Ndiv<<endl<<"	 amp/rmax = "<<amp/rmax<<" amp/rmin="<<amp/rmin<<flush;
		for( unsigned int i = 0; i< Ndiv ; ++i)
		{
			heightProbe dirtyDetect( bl +(.5+i)*amp - amp2, bl+(.5+i)*amp + amp2 );
			//cout<<" dirty prb "<<dirtyDetect.h1()<<" "<<dirtyDetect.h2()<<endl;
			add=true;
			inside.clear();


			for( unsigned int k=0; k < Nc; ++k )
			{
				if(  dirtyDetect.contain( sys_->nwk()->inter( sys_->nwk()->clist(k))))

					inside.push_back( sys_->nwk()->inter(sys_->nwk()->clist(k)));
			}
		//	cout<<"-----size "<<inside.size()<<endl;


			for( unsigned int j=0; j< Ndiv; ++j)
			{
					//if(add && j>500 || i>500) add=false;

				rectangularProbe temp( ll+ (.5+ j) *amp, bl+(.5+i)*amp, amp2,amp2) ;
			//	cout<<"     circular = "<<temp.x()<<" "<<temp.y()<<endl;

				for( unsigned int k=0; k < inside.size(); ++k )
				{
					if(  /*temp.containCenter( inside[k]) || temp.probeInsideBody(inside[k]) ||*/ temp.contain( inside[k])  )
					{
						//cout<<i<<" "<<j<<" "<<k<<endl;

						if (add) splTemp.addBody( new disk( temp.x(), temp.y(), temp.hh()) );
						Nbox++;
						break;
					}
				}

			}

		}
		//cout<<" export reconstruction "<<endl;
		ofstream out("spltemp.his",ios::out);
		out<<"Sample{"<<endl;
		splTemp.write( out);
		out<<"}"<<endl;

		cout<<" Nbox/Ncell "<<(double)(Nbox)/(double)(Ndiv*Ndiv)<<" Nbox= "<<Nbox<<" Npoint= "<<Npoint<<endl;
		//if( Ndiv>1 )
		if( Nbox != 0.) frac.add(  log10(amp), log10((double) (Nbox) ));
		Ndiv+=10;Npoint++;
		amp=h/(double) (Ndiv);
		amp2=.5*amp;
		frac.write("FBCnwk.txt");
	}

}


void shearP_CD_A::BmultifractalBC( )
{
	cout<<" multi fractal Dimension ( bodies with forces max ): BOX COUNTING "<<endl;

	unsigned int Nbinf=3;//Nombre de classe de force
	cout<<" 	Nombre d'intervalles de forces : "<<Nbinf<<endl;

	vector < double > fmax = forcesMaxCorrelation();
	double fmx=0,fmn=100000;
	double fmoy=0;
	unsigned int n=0;
	for ( unsigned int i=0; i< fmax.size(); ++i)
	{
		fmx= max( fmx,fmax[i]);

		if( fmax[i]!=0.)
		{
			fmn= min( fmn,fmax[i]);
			fmoy+=fmax[i];
			n++;
		}
	}
	fmoy/=(double) (n);
	cout<<" fmax = "<<fmx<<" fmn = "<<fmn<<endl;
	double ampf=(fmx-fmn)/(double) (Nbinf);

	double ll= sys_->spl()->xmin() + 1.*sys_->spl()->rmax();
	double rl= sys_->spl()->xmax() - 1.*sys_->spl()->rmax();
	double bl= totalProbe_.h1() ;//;+ sys_->spl()->rmax();
	double tl= totalProbe_.h2() ;//- sys_->spl()->rmax();

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int Nb= sys_->spl()->lbody().size();

	double h= min(rl-ll,tl-bl);//size of original  box
	double amp,amp2;//size of current circuler box
	unsigned int Ndiv;// Number of division of initial lenght
	unsigned int Nbox;//Number of intersected boxes
	unsigned int Npoint;
	pointSet frac;

	vector <body2d *> inside;

	bool add;

	//	removeRattlers();

	for ( unsigned int N=0; N<Nbinf;++N)
	{
		double flim=fmn+ (N+1.)*ampf;
		cout<<" 	--------- flim = "<<flim<<endl;
		Ndiv=(unsigned int ) floor( h/(1.1*rmax));
		amp=h/(double) (Ndiv);
		amp2=.5*amp;
		add=true;
		Npoint=0;

		char dir[50],fich[50];
		sprintf(dir,"./%d_%d",N+1,Nbinf);		
		mkdir (dir, S_IRUSR | S_IWUSR | S_IXUSR |
			S_IRGRP |           S_IXGRP |
			S_IROTH |           S_IXOTH  );

		while( amp > .1* rmin && Npoint < 30)
		{
			Sample splTemp;
			//cout<<"Ndiv= "<<Ndiv<<" amp= "<<amp<<endl;
			Nbox=0;
			cout<<" Ndiv= "<<Ndiv<<" carre= "<<Ndiv*Ndiv<<endl<<"	 amp/rmax = "<<amp/rmax<<" amp/rmin="<<amp/rmin<<flush;
			for( unsigned int i = 0; i< Ndiv ; ++i)
			{
				heightProbe dirtyDetect( bl +(.5+i)*amp - amp2, bl+(.5+i)*amp + amp2 );
				//cout<<" dirty prb "<<dirtyDetect.h1()<<" "<<dirtyDetect.h2()<<endl;
				add=true;
				inside.clear();


				for( unsigned int k=0; k < Nb; ++k )
				{
					if( fmax[k] <= flim)
						if ( sys_->spl()->body(k)->sizeVerlet()>= amp)
						if( dirtyDetect.intersection( sys_->spl()->body(k)) || dirtyDetect.containCenter( sys_->spl()->body(k)))

						inside.push_back( sys_->spl()->body(k));
				}
			//	cout<<"-----size "<<inside.size()<<endl;


				for( unsigned int j=0; j< Ndiv; ++j)
				{
					if(add && j>500 || i>500) add=false;

					circularProbe temp( ll+ (.5+ j) *amp, bl+(.5+i)*amp, amp2) ;
				//	cout<<"     circular = "<<temp.x()<<" "<<temp.y()<<endl;

					for( unsigned int k=0; k < inside.size(); ++k )
					{
						if(  /*temp.containCenter( inside[k]) || temp.probeInsideBody(inside[k]) ||*/ temp.intersection( inside[k])  )
						{
							//cout<<i<<" "<<j<<" "<<k<<endl;

							if (add) splTemp.addBody( new disk( temp.x(), temp.y(), temp.R()) );
							Nbox++;
							break;
						}
					}

				}

			}
			//cout<<" export reconstruction "<<endl;


			sprintf(fich,"%s/spltemp.his",dir);
			ofstream out(fich,ios::out);
			out<<"Sample{"<<endl;
			splTemp.write( out);
			out<<"}"<<endl;

			cout<<" Nbox/Ncell "<<(double)(Nbox)/(double)(Ndiv*Ndiv)<<" Nbox= "<<Nbox<<" Npoint= "<<Npoint<<endl;
			//if( Ndiv>1 )
			if( Nbox != 0.) frac.add(  log10(amp), log10((double) (Nbox) ));
			Ndiv+=1000;Npoint++;
			amp=h/(double) (Ndiv);
			amp2=.5*amp;
			sprintf(fich,"%s/FBCspl.txt",dir);
			frac.write(fich);
		}
		frac.psetClear();
	}
}

void shearP_CD_A::NmultifractalBC( )
{
	cout<<" multi fractal Dimension ( contact with forces max ): BOX COUNTING "<<endl;

	unsigned int Nbinf=3;//Nombre de classe de force
	cout<<" 	Nombre d'intervalles de forces : "<<Nbinf<<endl;

	vector < double > fmax = forcesMaxCorrelation();
	double fmx=0,fmn=100000;
	double fmoy=0;
	unsigned int n=0;
	for ( unsigned int i=0; i< fmax.size(); ++i)
	{
		fmx= max( fmx,fmax[i]);

		if( fmax[i]!=0.)
		{
			fmn= min( fmn,fmax[i]);
			fmoy+=fmax[i];
			n++;
		}
	}
	fmoy/=(double) (n);
	cout<<" fmax = "<<fmx<<" fmn = "<<fmn<<endl;
	double ampf=(fmx-fmn)/(double) (Nbinf);

	double ll= sys_->spl()->xmin() + 1.*sys_->spl()->rmax();
	double rl= sys_->spl()->xmax() - 1.*sys_->spl()->rmax();
	double bl= totalProbe_.h1() ;//;+ sys_->spl()->rmax();
	double tl= totalProbe_.h2() ;//- sys_->spl()->rmax();

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int Nb= sys_->nwk()->clist().size();

	double h= min(rl-ll,tl-bl);//size of original  box
	double amp,amp2;//size of current circuler box
	unsigned int Ndiv;// Number of division of initial lenght
	unsigned int Nbox;//Number of intersected boxes
	unsigned int Npoint;
	pointSet frac;

	vector <inter2d *> inside;

	bool add;

	//	removeRattlers();

	for ( unsigned int N=0; N<Nbinf;++N)
	{
		double flim=fmn+ (N+1.)*ampf;
		cout<<" 	--------- flim = "<<flim<<endl;
		Ndiv=(unsigned int ) floor( h/(1.1*rmax));
		amp=h/(double) (Ndiv);
		amp2=.5*amp;
		add=true;
		Npoint=0;

		char dir[50],fich[50];
		sprintf(dir,"./%d_%d",N+1,Nbinf);		
		mkdir (dir, S_IRUSR | S_IWUSR | S_IXUSR |
			S_IRGRP |           S_IXGRP |
			S_IROTH |           S_IXOTH  );

		while( amp > .1* rmin && Npoint < 30)
		{
			Sample splTemp;
			//cout<<"Ndiv= "<<Ndiv<<" amp= "<<amp<<endl;
			Nbox=0;
			cout<<" Ndiv= "<<Ndiv<<" carre= "<<Ndiv*Ndiv<<endl<<"	 amp/rmax = "<<amp/rmax<<" amp/rmin="<<amp/rmin<<flush;
			for( unsigned int i = 0; i< Ndiv ; ++i)
			{
				heightProbe dirtyDetect( bl +(.5+i)*amp - amp2, bl+(.5+i)*amp + amp2 );
				//cout<<" dirty prb "<<dirtyDetect.h1()<<" "<<dirtyDetect.h2()<<endl;
				add=true;
				inside.clear();


				for( unsigned int k=0; k < Nb; ++k )
				{
					if( fmax[k] <= flim)
						if(  dirtyDetect.contain( sys_->nwk()->inter(sys_->nwk()->clist(k))))

						inside.push_back( sys_->nwk()->inter(sys_->nwk()->clist(k)));
				}
			//	cout<<"-----size "<<inside.size()<<endl;


				for( unsigned int j=0; j< Ndiv; ++j)
				{
				//	if(add && j>500 || i>500) add=false;

					rectangularProbe temp( ll+ (.5+ j) *amp, bl+(.5+i)*amp, amp2,amp2) ;
				//	cout<<"     circular = "<<temp.x()<<" "<<temp.y()<<endl;

					for( unsigned int k=0; k < inside.size(); ++k )
					{
						if(  temp.contain( inside[k])  )
						{
							//cout<<i<<" "<<j<<" "<<k<<endl;

							if (add) splTemp.addBody( new disk( temp.x(), temp.y(), temp.hh()) );
							Nbox++;
							break;
						}
					}

				}

			}
			//cout<<" export reconstruction "<<endl;


			sprintf(fich,"%s/spltemp.his",dir);
			ofstream out(fich,ios::out);
			out<<"Sample{"<<endl;
			splTemp.write( out);
			out<<"}"<<endl;

			cout<<" Nbox/Ncell "<<(double)(Nbox)/(double)(Ndiv*Ndiv)<<" Nbox= "<<Nbox<<" Npoint= "<<Npoint<<endl;
			//if( Ndiv>1 )
			if( Nbox != 0.) frac.add(  log10(amp), log10((double) (Nbox) ));
			Ndiv+=10;Npoint++;
			amp=h/(double) (Ndiv);
			amp2=.5*amp;
			sprintf(fich,"%s/FBCnwk.txt",dir);
			frac.write(fich);
		}
		frac.psetClear();
	}
}



void shearP_CD_A::fractalDensity()
{
	cout<<" fractal Density : "<<endl;

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	double ampr= (rmax-rmin)/1000.;

	unsigned int Ninside;
	double densite,rlim;
	pointSet dr;


	rlim=rmin-ampr;
	while( rlim < rmax)
	{	
		Ninside=0;
		for ( unsigned int i=0; i< sys_->spl()->lbody().size();++i)
		{
			if( totalProbe_.containCenter( sys_->spl()->body(i) ) && sys_->spl()->body(i)->sizeVerlet() > rlim)		
				Ninside++;
		}
		densite= (double) (Ninside)/totalProbe_.area();
		dr.add( rlim, densite);
		rlim+=ampr;
	}
	dr.write("fractalDenity.txt");


}
