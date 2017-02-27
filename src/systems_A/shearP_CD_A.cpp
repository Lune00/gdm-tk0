#include "shearP_CD_A.hpp"

void shearP_CD_A::plugRef()
{
	partref=sys_->spl()->body(partrefId);
	
}



void shearP_CD_A::read_parameters(istream & is)
{
	string token,temp;
	is>>token;
	double trash;
	NbinPT=NbinFN=NbinFT=Nbingranulo=NbinFA=0;
	while(is)
	{
		//cout<<token<<endl;
		if     (token== "/*")
		{
			is>>temp;
			while( temp != "*/")
			{
				cout<<" @shearP_CD_A :Ingored parameter: "<<temp<<endl;
				is >> temp;
			}
		}
		else if (token == "ZoomSample") is >> zoom_;
		else if(token== "Sample") displaySample=true;
		else if(token== "Force") displayForce=true;
		else if(token== "Solidfraction") calcsf=true;
		else if(token== "Fabric") calcFabric=true;
		else if(token== "Fractal") calcfracdim=true;
		else if(token== "Granulo") granulo=true;
		else if(token== "Inout")
		{
			calcinout=true;
			is >> Ninout;
		}
	//	else if(token== "Filtergap") filter=true;
		
		else if(token== "FA")
		{
			 calcforcesA=true;
			is>>NbinFA;
			is>>temp;
			if(temp=="plot")
			{
				plotFA=true;
			}
			else if( temp=="noplot")
			{
				plotFA=false;
			}
			else
			{
				 cout<<" @shearP_CD_A :Forcesanisotropy : missing plot command : "<<endl;
			}
		}
		else if(token== "FLA") 
		{
			 calcforcesAL=true;
			is>>NbinFA;
			is>>temp;
			if(temp=="plot")
			{
				plotFA=true;
			}
			else if( temp=="noplot")
			{
				plotFA=false;
			}
			else
			{
				 cout<<" @shearP_CD_A :Forcesanisotropy : missing plot command : "<<endl;
			}
		}
		else if(token== "FCA")
		{
			 calcforcesAC=true;
			is>>NbinFA;
			is>>temp;
			if(temp=="plot")
			{
				plotFA=true;
			}
			else if( temp=="noplot")
			{
				plotFA=false;
			}
			else
			{
				 cout<<" @shearP_CD_A :Forcesanisotropy : missing plot command : "<<endl;
			}
		}
		
		else if(token== "Export") exportDistribution=true;
		else if(token== "Granulopdf") 
		{
			calcgranulopdf=true;
			is>>Nq_;
		}
		else if(token=="Z") calcz=true;
		else if(token=="ZP") calczp=true;
		else if(token=="Zgranulo") 
		{
			calczg=true;
			if (Nbingranulo !=0 )
			{
				cout<<" @shearP_CD_A :Zgranulo : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ shearP_CD_A :Zgranulo : Nbingranulo undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Granulostress") 
		{
			calcgranulostress=true;
			if (Nbingranulo !=0 )
			{
				cout<<" @shearP_CD_A :granulostress : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ shearP_CD_A :granulostress:  Nbingranulo undefined "<<endl;
				exit(0);
			}
			is>>perG;
			if ( perG == 0 )
			{
				cout<<" @ shearP_CD_A :granulostress:  perG undefined "<<endl;
				exit(0);
			}
			is>>wG;
			if ( wG == 0 )
			{
				cout<<" @ shearP_CD_A :granulostress:  wG undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="PDFFN") 
		{
			calcfn=true;
			is >> NbinFN;

			if ( NbinFN == 0 )
			{
				cout<<" @ shearP_CD_A : NbinFN undefined "<<endl;
				exit(0);
			}
			is >> perF;
			if ( perF == 0 )
			{
				cout<<" @ shearP_CD_A : perF undefined "<<endl;
				exit(0);
			}
			is >> wF;
			if ( wF == 0 )
			{
				cout<<" @ shearP_CD_A : wF undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="PDFL") 
		{
			calcl=true;
			is >> NbinL;

			if ( NbinL == 0 )
			{
				cout<<" @ shearP_CD_A : NbinL undefined "<<endl;
				exit(0);
			}
			is >> perL;
			if ( perL == 0 )
			{
				cout<<" @ shearP_CD_A : perL undefined "<<endl;
				exit(0);
			}
			is >> wL;
			if ( wL == 0 )
			{
				cout<<" @ shearP_CD_A : wL undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Ptheta") 
		{
			calcPtheta=true;
			is >> NbinPT;

			if ( NbinPT == 0 )
			{
				cout<<" @ shearP_CD_A : NbinPT undefined "<<endl;
				exit(0);
			}
			is >> mobperiod;
			if ( mobperiod == 0 )
			{
				cout<<" @ shearP_CD_A : mobperiod undefined "<<endl;
				exit(0);
			}
			is >> mobwidth;
			if ( mobwidth == 0 )
			{
				cout<<" @ shearP_CD_A : mobwidth undefined "<<endl;
				exit(0);
			}

		}

        //else
		else if( token=="norm") normpdf=true;
		else if( token=="Gap") calcgap=true;
      
		else if( token=="Def") calcdef=true;
		else if(token=="PDFFT")
		{
			calcft=true;
			is >> NbinFT;
			is >> perF;
			is >> wF;
			if ( NbinFT == 0 )
			{
				cout<<" @ shearP_CD_A : NbinFT undefined "<<endl;
				exit(0);
			}
		} 
		else if(token=="Sprofile") 
		{
			calcSprofile=true;
			is >> Nprb_;
			if ( Nprb_ == 0 )
			{
				cout<<" @ shearP_CD_A : Nprb_ undefined "<<endl;
				exit(0);
			}
		} 
		else if(token=="SFprofile") 
		{
			calcSFprofile=true;
			is >> Nprb_;
			if ( Nprb_ == 0 )
			{
				cout<<" @ shearP_CD_A : Nprb_ undefined "<<endl;
				exit(0);
			}
		} 

		else if(token=="Globalstress") calcglobalstress=true;
		else if(token=="PartialLengthStress")
		{
			 calcpartialLengthstress=true;
				is >> Npointps;
		}
		else if(token=="PartialNormalForceStress")
		{
			 calcpartiaNormalForcestress=true;
				is >> Npointps;
		}
		else if(token=="Forcescorrelation") 
		{
			calcFC=true;
			is >> Nbincor_;
		}
		else if(token=="Sumforce") sumforce_defined=true;
		else if(token=="GranuloBin") 
		{
			if (Nbingranulo !=0 )
			{
				cout<<" @shearP_CD_A :granulobin : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ shearP_CD_A :granulobin : Nbingranulo undefined "<<endl;
				exit(0);
			}
		}
		else if(token=="Nanalyze") is>>Nanalyze();
		else if(token=="Remove")	removeR=true;
		else if( token=="RemoveVisu") removeVisu=true;
		else if( token=="Visu") visu=true;

		else if(token=="Grow") 
		{
			growR=true;
			is>>incR;
		}
		else if(token=="RFD")
		{
			calcrfd=true;
			is >> Np;
			is >> Nd;

		}
		else if(token=="}") break;
		//else if(token=="include_interrupt") incl_int=true;
		else cout<<" @ ShearP_CD_A : parametre de commande inconnu :"<<token<<endl;

		is>>token;
	}	

}


void shearP_CD_A::setFolder(std::string dir)
{
	princDir_=dir;
	cout<<princDir_<<endl;
}	


void shearP_CD_A::initAnalyse( ) 
{
    
	cout<<"-----oooo---shearP_CD_A::init()----oooo----"<<endl;

	//Particule de reference pour le calcul des deformations
	partref=sys_->ldof(1)->leftBody();
	partrefId=partref->id();
	cout<<partref->x()<<" "<<partref->y()<<endl;

	w0 = sys_->spl()->boundWidth();
//	h0 = partref->y() - origindef->y();

	X00 = partref->x();
	Y00 = sys_->ldof(0)->leftBody()->y();

	X0 = partref->x();
	Y0 = partref->y();

	strain.xx()=0.;
	strain.xy()=0.;
	strain.yx()=0.;
	strain.yy()=0.;

	//defXY=defYY=0.;


	//height probe:
	

	totalProbe_.width() = sys_->spl()->boundWidth();
	totalProbe_.h1()= (sys_->ldof(0)->lowerBody())->ymin()+sys_->spl()->rmax();//sys_->spl()->body(0)->y()+sys_->spl()->rmax();
	totalProbe_.h2()= (sys_->ldof(1)->lowerBody())->ymin()-sys_->spl()->rmax();//-sys_->spl()->rmax();


	//Rectangular probe:
	
	//On la met au centre de l'Žchantillon:
	totalProbe_R.x()=(sys_->spl()->rightBoundary() - sys_->spl()->leftBoundary())*0.5;
	totalProbe_R.y()=((sys_->ldof(1)->lowerBody())->ymin()-(sys_->ldof(0)->lowerBody())->ymin())*0.5;
	//On dŽfinit les dimensions probe:
	totalProbe_R.hh()=(sys_->ldof(1)->lowerBody())->ymin() -2*sys_->spl()->rmax()-totalProbe_R.y();
	totalProbe_R.hl()=sys_->spl()->rightBoundary()-3*sys_->spl()->rmax()-totalProbe_R.x();


	double p_xmin,p_xmax,p_ymin,p_ymax,rmax,xc,yc;

	xc=totalProbe_R.x();
	yc=totalProbe_R.y();
	p_xmin=totalProbe_R.x()-totalProbe_R.hl();
	p_xmax=totalProbe_R.x()+totalProbe_R.hl();
	p_ymin=totalProbe_R.y()-totalProbe_R.hh();
	p_ymax=totalProbe_R.y()+totalProbe_R.hh();
	rmax=sys_->spl()->rmax();
	//Ecriture des dimensions de la probe:
	
	cout<<"Probe position :"<<totalProbe_R.x()<<" "<<totalProbe_R.y()<<endl;
	cout<<string(100,'*')<<endl;
	
	Vinit_=totalProbe_.area();
	
	system("mkdir -p Analyse");
	system("mkdir -p Analyse/Anisotropies");

	ofstream strain("Analyse/strain.txt",ios::out);
	strain.close();	
	
	if(calcinout)
	{
	ofstream inout("Analyse/inout_time.txt",ios::out);
	inout.close();
	}

	if ( calcforcesA ) 
	{
		ofstream FA_out("Analyse/Anisotropies/FAnisotropy.txt",ios::out); 
		//FA_out<<"# time epsq ac an at al .5*(ac+an+at)  q/p"<<endl;
		FA_out.close();
		ofstream DPM_out("Analyse/Anisotropies/F_DPM.txt",ios::out); FA_out.close();
				// ofstream FA_out("FA.txt",ios::out); FA_out.close();
	}
	if ( calcforcesAC ) 
	{
	ofstream FA_out("Analyse/Anisotropies/FCAnisotropy.txt",ios::out);
	FA_out.close();

	ofstream DPM_out("Analyse/Anisotropies/FC_DPM.txt",ios::out);
	DPM_out.close();
	}
	if ( calcforcesAL) 
	{
	ofstream FA_out("Analyse/Anisotropies/FLAnisotropy.txt",ios::out);
	FA_out.close();

	ofstream DPM_out("Analyse/Anisotropies/FL_DPM.txt",ios::out);
	DPM_out.close();
	}
	
	if( calcfn || calcft || calcl)
	{
		system("mkdir Analyse/pdf");
	}
	
	if( calcglobalstress)
	{
		ofstream GS_out("Analyse/stress.txt",ios::out);
		GS_out.close();
		
	}

	if ( calcSprofile ) 
	{
		ofstream XS_out("Analyse/SProfile.txt",ios::out); XS_out.close();
				// ofstream FA_out("FA.txt",ios::out); FA_out.close();
	}

	if (calcz)
	{
		system("mkdir Analyse/Connect");
	//	ofstream FA_out("Connect/Zp.txt",ios::out); FA_out.close();
		if (calczp) {ofstream zp_out("Analyse/Connect/Zpt.txt",ios::out); zp_out.close();}
		ofstream scalZ_out("Analyse/Connect/scalarZ.txt",ios::out); scalZ_out.close();
		if( granulo )
		{
			ofstream zg_out("Analyse/Connect/Granulo.txt",ios::out);zg_out.close();
		}

	}

	if( calcPtheta)
	{
		system("mkdir Analyse/angDistrib");
		
	}
	
	if( calcFC)
	{
		system("mkdir Analyse/forceCorr");
		ofstream lc("Analyse/forceCorr/linearCor.txt",ios::out);lc.close();
		
	}


		//cout<<"h1= "<<totalProbe_.h1()<<" h2 = "<<totalProbe_.h2()<<endl;
		//cout<<"Compacite initiale = "<<sf()<<endl;		
		//if( calcfracdim ) removeR=true;
		
	if (displaySample || displayForce)
	{
		system("mkdir -p Analyse/PS1");
	}
		
}

void shearP_CD_A::analyse( double t, unsigned int nsi, unsigned int nsf )
{
	cout<<"------oooo><ooo-----shearP_CD_A::analyze() t= "<<t<<" ---------"<<endl;
	
	char fname[100];
	
	
	time=t;
    
    
    followparticles();

	if (calcz)  Z(calczp,Nbingranulo,granulo);
	if (calcPtheta) Ptheta( NbinPT, mobperiod, mobwidth);
	
	if( removeR) removeRattlers();
	if( growR) growRattlers( sys_->spl()->rmin()*incR);

	if (calcglobalstress) globalStress();
	if (calcpartialLengthstress) partialLengthStress(Npointps);
	if (calcpartiaNormalForcestress) partialNormalForceStress(Npointps);
	
	if (calcfn) pdfforce(true,NbinFN,normpdf,perF,wF);
	if (calcft) pdfforce(false,NbinFT,normpdf,perF,wF);	
	if (calcl)  pdflength(NbinL,normpdf,perL,wL);
	
	if (calcinout) normalForceInOut(Ninout);

    
	if (calcsf) SF();
	
	if (calcFabric)  A();
	if (calcforcesA)	forces_A(NbinFA,plotFA);
	if (calcforcesAL)	forces_AL(NbinFA,plotFA);
	if (calcforcesAC)	forces_AC(NbinFA,plotFA);
	if (calcSprofile ) 	 profiles(calcSprofile,calcSFprofile);
	if (calcgranulostress) granuloStress3( Nbingranulo);//,perG,wG);
    
	if (calcgap) Gap();
	if (calcdef) def();
	
	if (displaySample || displayForce)
	{
		sprintf(fname,"Analyse/PS1/particles%.4i.ps",Nanalyze());
		writePS(fname);
	}

	cout<<"		Nombre d'interactions = "<<sys_->nwk()->linter().size()<<endl;
	cout<<"		Nombre de contacts   = "<<sys_->nwk()->clist().size()<<endl;
	
	
	Nanalyze() ++;


}


void shearP_CD_A::Ptheta(unsigned int Nbin , unsigned int period, unsigned int width) 
{
	cout<<"	Ptheta  : "<<endl;
	vector <double> PT(Nbin,0);
	vector <double> FN(Nbin,0);
	vector <double> FT(Nbin,0);
	vector <double> L(Nbin,0);
	vector <unsigned int > NcN(Nbin,0);
	vector <unsigned int > NcT(Nbin,0);

	double fnmoy=0;
	double ftmoy=0;
	double lmoy =0,l;
	double angN,dx,dy;
	double amp = M_PI/(double) (Nbin);
	unsigned int NctotN=0,rangN,NctotT=0;
	inter2d * interc;
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));

		if( interc->fn() !=0. && totalProbe_.contain(interc))
		{
			angN= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angN = M_PI - angN;

			rangN=(unsigned int) ( floor( angN/amp));
			if (rangN==Nbin) rangN--;

			dx=interc->Vbranchx();
			dy=interc->Vbranchy();
			l=sqrt(dx*dx+dy*dy);

			NcN[ rangN ]++;
			FN [ rangN ]+=interc->fn();
			L  [ rangN ]+=l;
			NctotN++;
			
			fnmoy+=interc->fn();
			lmoy +=l;

			if( interc->ft() != 0. )
			{
				FT [ rangN ]+=interc->ft();
				NcT[ rangN ]++;
				ftmoy += interc->ft();
				NctotT++;
			}
		}


	}
	cout<<" Nctot = "<<NctotN<<endl;
	fnmoy/=(double) (NctotN);
	ftmoy/=(double) (NctotT);
	lmoy /=(double) (NctotN);
	
	fnmoy_=fnmoy;
	lmoy_=lmoy;

	pointSet fn,ft,ld,p;

	for ( unsigned int i =0;i<Nbin;++i)
	{
		p.add( (.5+i)*amp, (double) (NcN[i])/(double) (NctotN)*(double)(Nbin)/M_PI );//Distrib angulaire de prob
		fn.add( (.5+i)*amp, FN[i]/(double) (NcN[i]) );//Distrib ang de l'intensite de forces normales
		ld.add( (.5+i)*amp, L[i]/(double) (NcN[i]) );//Distrib ang des longueurs
		
		if( NctotT!=0)
			ft.add( (.5+i)*amp, FT[i]/(double) (NcT[i]));

	}

	
	if( NctotT!=0)
	{
		ft.circularise(M_PI,width);
		pointSet ftm =ft.mobileMean( period,width);
		ftm.closeCircular();
		ftm.write("Analyse/angDistrib/FTtheta_mob.txt");
		ftm.yNormalise(fnmoy);
		ftm.write("Analyse/angDistrib/FTtheta_MOD.txt");
	}   

	
	p.circularise(M_PI,width);
	fn.circularise(M_PI,width);
	ld.circularise(M_PI,width);
	
	pointSet pm  = p.mobileMean(  period,width);
	pointSet fnm = fn.mobileMean( period,width);
	pointSet ldm = ld.mobileMean( period,width);
	pm.closeCircular();
	fnm.closeCircular();
	ldm.closeCircular();
	
	pm.write("Analyse/angDistrib/Ptheta_mob.txt");
	fnm.write("Analyse/angDistrib/FNtheta_mob.txt");
	ldm.write("Analyse/angDistrib/Ltheta_mob.txt");

	pointSet  pt = PthetaInProbe( totalProbe_, *(sys_)->spl(),*(sys_)->nwk() ,Nbin);
	pt.write("Analyse/angDistrib/Ptheta_mob2.txt");
	
	fnm.yNormalise(fnmoy);
	ldm.yNormalise(lmoy);
	
	fnm.write("Analyse/angDistrib/FNtheta_MOD.txt");
	ldm.write("Analyse/angDistrib/Ltheta_MOD.txt");
	
	vector <double> Ls(Nbin,0);
	vector <unsigned int > NcNs(Nbin,0);
	vector <double> FNs(Nbin,0);
	vector <double> Lw(Nbin,0);
	vector <unsigned int > NcNw(Nbin,0);
	vector <double> FNw(Nbin,0);
	
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));

		if( interc->fn() !=0. && totalProbe_.contain(interc))
		{
			angN= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angN = M_PI - angN;

			rangN=(unsigned int) ( floor( angN/amp));
			if (rangN==Nbin) rangN--;

			dx=interc->Vbranchx();
			dy=interc->Vbranchy();
			l=sqrt(dx*dx+dy*dy);
			
			if( interc->fn() >= fnmoy)
			{
			NcNs[ rangN ]++;
			FNs [ rangN ]+=interc->fn();
			Ls  [ rangN ]+=l;
			}
			else
			{
				NcNw[ rangN ]++;
				FNw [ rangN ]+=interc->fn();
				Lw  [ rangN ]+=l;
			}
		}
	}
	
	pointSet fns,fnw,ls,lw;

	for ( unsigned int i =0;i<Nbin;++i)
	{
		fns.add( (.5+i)*amp, FNs[i]/(double) (NcNs[i]) );//Distrib ang de l'intensite de forces normales fortes
		fnw.add( (.5+i)*amp, FNw[i]/(double) (NcNw[i]) );//Distrib ang de l'intensite de forces normales fortes

		ls.add( (.5+i)*amp, Ls[i]/(double) (NcNs[i]) );//Distrib ang des longueurs
		lw.add( (.5+i)*amp, Lw[i]/(double) (NcNw[i]) );//Distrib ang des longueurs
	}	
	
	(fns.mobileMean( period,width)).write("Analyse/angDistrib/fn_theta_strong.txt");
	(fnw.mobileMean( period,width)).write("Analyse/angDistrib/fn_theta_weak.txt");
	
	(ls.mobileMean( period,width)).write("Analyse/angDistrib/l_theta_strong.txt");
	(lw.mobileMean( period,width)).write("Analyse/angDistrib/l_theta_weak.txt");
	
	
}


void shearP_CD_A::Gap()
{

	if(sys_->nwk()->linter().empty()) return;
	unsigned int i,Ngap=0;
	unsigned int Ni = sys_->nwk()->linter().size();
	double d;
	double temp;
	temp=-100.;
	gapmoy_=0.;
	ofstream gapout("Analyse/gap.txt",ios::out);
    
	for ( i=0;i< Ni;++i)
	{
		d=sys_->nwk()->inter(i)->Overlap(); // Overlap = -1 si pas d'overlap
        //cout<<"d="<<d<<" "<<i<<endl;
        // Si overlap de deux particules:
		if ( d > 0. )
		{
			gapout<<d/sys_->spl()->rmin()<<endl;
			gapmoy_+= d;
			Ngap++;
			
			if ( d > temp)
			{	
				ngap1_=sys_->nwk()->inter(i)->first();
                
				ngap2_=sys_->nwk()->inter(i)->second();
				temp=d;
				gaprmin = d/sys_->spl()->rmin();
				gaprmax = d/sys_->spl()->rmax();
			}
		}
	}
	gapout.close();
}


void shearP_CD_A::def()
{
    if(sys_->nwk()->linter().empty()) return;
    
    if ( fabs(partref->x() - X0)>.5* sys_->spl()->boundWidth())
    {
        cout<<"*****"<<endl;
        cout<<".Period completed."<<endl;
        cout<<"*****"<<endl;
        X00+=sys_->spl()->boundWidth();
    }
    
    //On calcule la dŽformation totale a partir d'une refŽrence initiale: epsxy
    //Volume variation
    dvov_=(totalProbe_.area()-Vinit_)/Vinit_;

    
    //Strain Tensor : en petite dŽformation, calcul par incrŽment (faux)
    
    double defx = ( partref->x() - X0)/sys_->spl()->boundWidth();
    strain.xx() += defx ;//0.;//( partref->x() - X0)/w0;//(partref->x()-X00); A test
    strain.xy() += ( partref->x() - X0)/(partref->y()-Y00);
    strain.yx() += 0.;//( partref->y() - Y0)/(w0);
    strain.yy() += ( partref->y() - Y0)/(partref->y()-Y00);
    
    ofstream str("Analyse/strain.txt",ios::app);
    
    X0=partref->x();
    Y0=partref->y();
    
    strain.eigenValues();
    
    epsp=strain.l1()+strain.l2();
    epsq=max(strain.l1(),strain.l2())-min(strain.l1(),strain.l2());
    epsxy=-(partref->x() - X00)/(partref->y()-Y00);
    cout<<" epsp = "<<epsp<<" epsq = "<<epsq<<" espxy = "<<epsxy<<endl;
    
    str<<time<<" "<<epsxy<<" "<<epsp<<" "<<epsq<<" "<<endl;
    
    str.close();
    
}


void shearP_CD_A::SF()
{
	sf()=solidFraction(totalProbe_,*sys_->spl(), *sys_->nwk());
	cout<<".Solid Fraction = "<<sf()<<endl;
}

void shearP_CD_A::profiles(bool Speed, bool Solfrac)
{
	unsigned int Nprobe = Nprb_;//floor( .5*( totalProbe_.h2() - totalProbe_.h1() ) / sys_->spl()->rmax() );
	double ampProbe=( totalProbe_.h2() - totalProbe_.h1() ) / (double) (Nprobe);
    
	cout<<".Profile number of probes = "<<Nprobe<<endl;
	cout<<".Profile amplitude of probes/<d> = "<< ampProbe/(2. * sys_->spl()->rmoy())<<endl;
    
	
	vector < heightProbe *> lprobe(Nprobe);


	for ( unsigned int i=0;i<Nprobe;++i)
	{
		lprobe[i]=new heightProbe( totalProbe_.h1() + (double) (i) * ampProbe,
			totalProbe_.h1() + (double) (i+1) * ampProbe);
	}

	if(Solfrac)
	{
        cout<<".Solid fraction Profile "<<flush;
		vector <double> SF(Nprobe,0.);

		for ( unsigned int i=0;i<Nprobe;++i)
		{
			SF[i]=solidFraction(*lprobe[i],*sys_->spl(),*sys_->nwk());
		}

		ofstream SFprofile ( "SFProfile.txt" , ios::out|ios::app );
		if ( ! SFprofile ) cout<<"erreur creation de SFprofile"<<endl;
		for ( unsigned int i=0;i<Nprobe;++i)
		{
			SFprofile<<time<<" "<<lprobe[i]->halfHeight()<<" "<<SF[i]<<endl;
		}
		SFprofile.close();
	}

	if( Speed)
	{
		cout<<".Speed Profile "<<endl;
		vector <double> XS(Nprobe,0.);
		vector <double> YS(Nprobe,0.);
	//cout<<" taille lpr "<<lprobe.size()<<" taille XS "<<XS.size()<<endl;

		speedProfile( lprobe  , XS , YS , *sys_->spl()  );


		ofstream Sprofile ( "Analyse/SProfile.txt" , ios::out|ios::app );
		if ( ! Sprofile ) cout<<"erreur creation de Sprofile"<<endl;


		ofstream Sinst ( "Analyse/SPinst.txt" , ios::out );
		if ( ! Sprofile ) cout<<"erreur creation de SPint"<<endl;


		for ( unsigned int i=0;i<Nprobe;++i)
		{
			Sprofile<<time<<" "<<lprobe[i]->halfHeight()<<" "<<XS[i]<<" "<<YS[i]<<endl;
			Sinst<<time<<" "<<lprobe[i]->halfHeight()<<" "<<XS[i]<<" "<<YS[i]<<endl;
		}

		Sprofile.close();
		Sinst.close();

	}
}



void shearP_CD_A::normalForceInOut(int Ninout)
{
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No InOut "<<endl;return ;}
	
	unsigned int Nc = sys_->nwk()->clist().size();
	
	unsigned int i=0;
	inter2d * interc, *interc2;
	
	DataSet fnd;
	DataSet fnind;
	DataSet fnoutd;
	DataSet Nind;
	DataSet Noutd;
	DataSet fninoutnull;
	
	double fn,fnin,fnout;
	unsigned int Nin,Nout,Noutnull=0;
	double fproj;
	
	fnin=fnout=0.;
	
	ofstream out("Analyse/Connect/inout_dist.txt",ios::out);
	
	
	for ( i=0;i< Nc;++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));
		if (totalProbe_.contain(interc))
		{
			double nx = interc->nx();
			double ny = interc->ny();
			
			Nin=Nout=0;
			fn = interc->fn();
			fnin=fnout=0;
			
			//if(interc->first()->id()==86)
			//cout<<" interaction entre "<<interc->first()->id()<<" et "<<interc->second()->id()<<endl;
			
			
			for( unsigned int j=0;j<Nc;++j)
			{
				if(j==i) continue;
				
				interc2 = sys_->nwk()->inter(sys_->nwk()->clist(j));
									
					if( (interc)->first() == (interc2)->first() )
					{
						fproj = interc2->fx() * nx + interc2->fy() * ny;
						
						if( fproj<0.) 
						{
							fnin -= fproj;
							Nin++;
						}
						else
						{
							fnout += fproj;
							Nout++;
						}
					
						//cout<<"		1"<<endl;
					}
					else if( (interc)->first() == (interc2)->second() )
					{
						// - pour A-R de j vers i
						fproj = - (interc2->fx() * nx +  interc2->fy() * ny);
						
						if( fproj<0.) 
						{
							fnin -= fproj;
							Nin++;
						}
						else
						{
							fnout += fproj;
							Nout++;
						}
						
						//cout<<"		2"<<endl;
						
					}
					else if( (interc)->second() == (interc2)->first() )
					{
						fproj =  interc2->fx() * nx +  interc2->fy() * ny;
						
						if( fproj>0.) 
						{
							fnin += fproj;
							Nin++;
						}
						else
						{
							fnout -= fproj;
							Nout++;
						}
					
						//cout<<"		3"<<endl;
						
					}
					else if( (interc)->second() == (interc2)->second() )
					{
						fproj =  -(interc2->fx() * nx +  interc2->fy() * ny);
						
						if( fproj>0.) 
						{
							fnin += fproj;
							Nin++;
						}
						else
						{
							fnout -= fproj;
							Nout++;
						}
						
						//cout<<"		4"<<endl;
						
					}
					//else cout<<" probleme gars....  i = "<<i<<" j = "<<j<<endl;
			
			}//fin pour j
			
			//cout<< " verif  "<<fn<<" "<<.5*(fnin-fnout)<<" "<<Nin+Nout<<endl;
			//ecriture dans fichier
			out<< interc->fn() << " " << fnin <<" "<<fnout<<" "<<Nin<<" "<<Nout<<endl;
			
			if( Nout !=0)
			{
			fnd.add(fn);
			fnind.add(fnin);
			fnoutd.add(fnout);
			Nind.add(Nin);
			Noutd.add(Nout);
			}
			else
			{
				fninoutnull.add(fnin);
				Noutnull++;
			}
		}
		
	}
	out.close();
	
	fnd.   extractValues();
	
	fnind .Normalize(fnd.mean());
	fnoutd.Normalize(fnd.mean());
	fnd   .Normalize(fnd.mean());
	fninoutnull. Normalize(fnd.mean());
	
	fnind. extractValues();
	fnoutd.extractValues();
	Nind.  extractValues();
	Noutd. extractValues();
	fninoutnull.extractValues();
	
	fnind .DecreasingSort();
	fnoutd.DecreasingSort();
	fnd   .DecreasingSort();
	fninoutnull   .DecreasingSort();
	
	pointSet pdfin  = fnind.Rich_PDF(Ninout);
	pointSet pdfout = fnoutd.Rich_PDF(Ninout);
	pointSet pdffn  = fnd.Rich_PDF(Ninout);
	pointSet pdfinoutnull = fninoutnull.Rich_PDF(Ninout);
	
	pdfin.write ("Analyse/Connect/pdfin.txt");
	pdfout.write("Analyse/Connect/pdfout.txt");
	pdffn.write ("Analyse/Connect/pdffn.txt");
	pdfinoutnull.write ("Analyse/Connect/pdfinoutnull.txt");
	
	
	ofstream out2("Analyse/inout_time.txt",ios::app);
	out2<<epsxy<<" "<< fnd.mean() << " " << fnind.mean() <<" "<<fnoutd.mean()<<" "<<Nind.mean()<<" "<<Noutd.mean()<<" "<<(double) Noutnull/(double)(fnd.setSize())<<endl;
	out2.close();	
		
}

void shearP_CD_A::A()
{
	an_=at_=al_=0;
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No Fabric "<<endl;return ;}
	
	cout<<"	Fabric A : ";
	gdm::Tensor2x2 * F = FabricInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk() );
	if (F != NULL ) 
	{
		F->eigenValues();
		a()=2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );
		da_=F->majorDirection();
		oa_=a()*cos( 2.*da_);
	}
		//F->adresse();
	delete F;

	cout<<" direction " << da_/M_PI*180.<<" ac= "<<a()<<endl;
}

void shearP_CD_A::forces_A(int Nbin,bool plot )
{
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	cout<<"-----> Uncorrelated anisotropies (classical) : "<<endl;
	
	double an,at,al,dn,dt,dl;
	an=at=al=dn=dt=dl=0;
	cout<<"	fn_A : ";
	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 10);
	/*//BON--gdm::Tensor2x2 * Fn = fnlAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 10);
	gdm::Tensor2x2 * L = lengthAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	*/

	if (Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		//Fn->print();
		an = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
		//	cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dn=Fn->majorDirection();
		an-=a();
		cout<<"     direction "<<dn/M_PI*180.<<" an= "<<an<<endl;
	}
	
	
	cout<<"	ft_A : ";
	gdm::Tensor2x2 * Ft = ftAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin );
	
	if( Ft != NULL )
	{
		*Ft= *Ft + *Fn;
		

		if (Ft->xx() != 0 && Ft->yy() !=0) 
		{
			Ft->eigenValues();
	//Ft->print();
			at=2.*( max(Ft->l2(),Ft->l1()) - min(Ft->l2(),Ft->l1()) )/(Fn->l1() + Fn->l2() );
	//cout<<Ft->l1()<<" "<<Ft->l2()<<endl;
			dt=Ft->majorDirection();
			//at()-=an()+a();
			//avec fnl->
			at-=a()+an;
			cout<<"     direction "<<dt/M_PI*180.<<" at= "<<at<<endl;
		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		at=0;
	}
	
	cout<<"	l_A : ";
	
	gdm::Tensor2x2 * L = lengthAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	//gdm::Tensor2x2 * F = FabricInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk() );
	
	if (L->xx() != 0. && L->yy() !=0.) 
	{
		//*L=*L+*F;
		L->eigenValues();
		al = 2.*( max(L->l2(),L->l1()) - min(L->l2(),L->l1()) )/(L->l1() + L->l2() );
		//	cout<<" max "<<max(L->l2(),L->l1())<<" "<<min(L->l2(),L->l1())<<endl;
		dl=L->majorDirection();
		al=al-a();
		cout<<"      direction "<<dl/M_PI*180.<<" al= "<<al<<endl;
		//L->print();
	}
	
	
	ofstream F_approx("Analyse/angDistrib/FA_approx.txt",ios::out);
	double amp = M_PI/(double) (Nbin);
	
	for (int i =0;i<2*Nbin+1;++i)
	{
		F_approx<<(.5+i)*amp<<" "
			<<1./M_PI*( 1.+ a() *cos( 2.*((.5+i)*amp - da_ )))<<" "
			<<fnmoy_*  (1. + an*cos( 2.*((.5+i)*amp - dn )))<<" "
			<<fnmoy_*  ( -1.*at*sin( 2.*((.5+i)*amp - dt)))<<" "
			<<lmoy_* (1.+ al*cos( 2.*((.5+i)*amp - dl )))<<" "
			<<fnmoy_<<" "<<lmoy_<<endl;
	}
	F_approx.close();
	
	ofstream F_mod("Analyse/angDistrib/FA_mod.txt",ios::out);
	
	for (int i =0;i<2*Nbin+1;++i)
	{
		F_mod<<(.5+i)*amp<<" "
			<< (1. + an*cos( 2.*((.5+i)*amp - dn )))<<" "
			<< ( -1.*at*sin( 2.*((.5+i)*amp - dt)))<<" "
			<< (1.+ al*cos( 2.*((.5+i)*amp - dl )))<<" "
			<<endl;
	}
	F_mod.close();
	
		ofstream FA_out("Analyse/Anisotropies/FAnisotropy.txt",ios::app);
		FA_out<<time<<" "<<epsxy<<" "<<a()<<" "<<an<<" "<<at<<" "<<al<<" "<<.5*( a()+an+at+al)<<" "<<qop_<<endl;
		FA_out.close();

		ofstream DPM_out("Analyse/Anisotropies/F_DPM.txt",ios::app);
		DPM_out<<time<<" "<<epsxy<<" "<<da_<<" "<<dn<<" "<<dt<<" "<<dl<<" "<<ds_<<endl;
		DPM_out.close();
	
	cout<<"		.5*(a+an+at+al) = "<<.5*(a()+an+at+al)<<endl;
	cout<<"		            q/p = "<<qop_<<endl;

	delete Fn;
	delete Ft;
	delete L;
	

	

}

void shearP_CD_A::forces_AL(int Nbin,bool plot)
{
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	cout<<"-----> Coupled Length-Forces anisotropies : "<<endl;

	double an,at,al,dn,dt;
	an=at=al=dn=dt=0;
	
	cout<<"	fn_A : ";
	gdm::Tensor2x2 * Fn = fnlAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 10);
	/*//BON--gdm::Tensor2x2 * Fn = fnlAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 10);
	gdm::Tensor2x2 * L = lengthAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	*/

	if (Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		//Fn->print();
		an = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
		//	cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dn=Fn->majorDirection();
		an-=a();
		cout<<"     direction "<<dn/M_PI*180.<<" an= "<<an<<endl;
	}
	
	
	cout<<"	ft_A : ";
	gdm::Tensor2x2 * Ft = ftlAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin );
	
	if( Ft != NULL )
	{
		*Ft= *Ft + *Fn;
		

		if (Ft->xx() != 0 && Ft->yy() !=0) 
		{
			Ft->eigenValues();
	//Ft->print();
			at=2.*( max(Ft->l2(),Ft->l1()) - min(Ft->l2(),Ft->l1()) )/(Fn->l1() + Fn->l2() );
	//cout<<Ft->l1()<<" "<<Ft->l2()<<endl;
			dt=Ft->majorDirection();
			//at()-=an()+a();
			//avec fnl->
			at-=a()+an;
			cout<<"     direction "<<dt/M_PI*180.<<" at= "<<at<<endl;
		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		at=0;
	}
	
	ofstream F_approx("Analyse/angDistrib/FLA_approx.txt",ios::out);
	double amp = M_PI/(double) (Nbin);
	double flnmoy = Fn->l1()+Fn->l2();
	for (int i =0;i<Nbin;++i)
	{
		F_approx<<(.5+i)*amp<<" "
			<<1./M_PI*( 1.+ a() *cos( 2.*((.5+i)*amp - da_ )))<<" "
			<<flnmoy*  (1. + an*cos( 2.*((.5+i)*amp - dn )))<<" "
			<<flnmoy*  ( -1.*at*sin( 2.*((.5+i)*amp - dt)))<<" "<<endl;
	}
	F_approx.close();
	
		ofstream FA_out("Analyse/Anisotropies/FLAnisotropy.txt",ios::app);
		FA_out<<time<<" "<<epsxy<<" "<<a()<<" "<<an<<" "<<at<<" "<<.5*( a()+an+at)<<" "<<qop_<<endl;
		FA_out.close();

		ofstream DPM_out("Analyse/Anisotropies/FL_DPM.txt",ios::app);
		DPM_out<<time<<" "<<epsxy<<" "<<da_<<" "<<dn<<" "<<dt<<" "<<ds_<<endl;
		DPM_out.close();
	
	cout<<"		.5*(a+anl+atl) = "<<.5*(a()+an+at)<<endl;
	cout<<"		           q/p = "<<qop_<<endl;
	
	delete Fn;
	delete Ft;

}

void shearP_CD_A::forces_AC(int Nbin ,bool plot)
{
	cout<<"-----> Correlated anisotropies with C(theta) : "<<endl;
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	double an,at,al,dn,dt,ac,dc,dl;
	an=at=al=ac=dn=dt=dl=dc=0;
	
	cout<<"	l_A : ";
	
	gdm::Tensor2x2 * L = lengthAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin);
//	gdm::Tensor2x2 * F = FabricInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk() );
	
	if (L->xx() != 0. && L->yy() !=0.) 
	{
		//*L=*L+*F;
		L->eigenValues();
		al = 2.*( max(L->l2(),L->l1()) - min(L->l2(),L->l1()) )/(L->l1() + L->l2() );
		//	cout<<" max "<<max(L->l2(),L->l1())<<" "<<min(L->l2(),L->l1())<<endl;
		dl=L->majorDirection();
		al-=a();
		cout<<"      direction "<<dl/M_PI*180.<<" al= "<<al<<endl;
		//L->print();
	}
	
	cout<<"	C_A : ";
	gdm::Tensor2x2 * Cn = CorrelationAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin);
	if (Cn->xx() != 0. && Cn->yy() !=0.) 
	{
		//*Cn=*Cn+*F;
		Cn->eigenValues();
		//Cn->print();
		ac = 2.*( max(Cn->l2(),Cn->l1()) - min(Cn->l2(),Cn->l1()) )/(Cn->l1() + Cn->l2() );
		//cout<<" ac = "<<ac<<endl;
		//	cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dc=Cn->majorDirection();
		ac-=a();
		cout<<"     direction "<<dc/M_PI*180.<<" ac= "<<ac<<endl;
	}
	
	
	
	cout<<"	fn_A : ";
	gdm::Tensor2x2 * Fn = fnlcAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), Nbin);

	if (Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		//Fn->print();
		an = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
		//	cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dn=Fn->majorDirection();
		an-=a()+al+ac;
		cout<<"     direction "<<dn/M_PI*180.<<" an= "<<an<<endl;
	}
	
	
	cout<<"	ft_A : ";
	gdm::Tensor2x2 * Ft = ftlcAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin );
	
	if( Ft != NULL )
	{
		*Ft= *Ft + *Fn;
		

		if (Ft->xx() != 0 && Ft->yy() !=0) 
		{
			Ft->eigenValues();
	//Ft->print();
			at=2.*( max(Ft->l2(),Ft->l1()) - min(Ft->l2(),Ft->l1()) )/(Fn->l1() + Fn->l2() );
	//cout<<Ft->l1()<<" "<<Ft->l2()<<endl;
			dt=Ft->majorDirection();
			//at()-=an()+a();
			//avec fnl->
			at-=a()+an+al+ac;
			cout<<"     direction "<<dt/M_PI*180.<<" at= "<<at<<endl;
		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		at=0;
	}
	
	
	double cmoy=1./(M_PI*M_PI) *(Cn->l1()+Cn->l2());
	cout<<"cmoy = "<<cmoy<<endl;
	
	
	ofstream F_approx("Analyse/angDistrib/FCA_approx.txt",ios::out);
	double amp = M_PI/(double) (Nbin);
	
	for (int i =0;i<Nbin;++i)
	{
		F_approx<<(.5+i)*amp<<" "
			<<1./M_PI*( 1.+ a() *cos( 2.*((.5+i)*amp - da_ )))<<" "
			<<fnmoy_*  (1. + an*cos( 2.*((.5+i)*amp - dn )))<<" "
			<<fnmoy_*  ( -1.*at*sin( 2.*((.5+i)*amp - dt)))<<" "
			<<lmoy_* (1.+ al*cos( 2.*((.5+i)*amp - dl )))<<" "
			<<cmoy* (1.+ ac*cos( 2.*((.5+i)*amp - dc )))<<" "<<endl;
	}
	F_approx.close();
	
	pointSet Ct = CorrelationthetaInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk() ,Nbin);
	Ct.write("Analyse/angDistrib/Ctheta.txt");
	
		ofstream FA_out("Analyse/Anisotropies/FCAnisotropy.txt",ios::app);
		FA_out<<time<<" "<<epsxy<<" "<<a()<<" "<<an<<" "<<at<<" "<<al<<" "<<ac<<" "<<.5*( a()+an+at+al)<<" "<<qop_<<endl;
		FA_out.close();

		ofstream DPM_out("Analyse/Anisotropies/FC_DPM.txt",ios::app);
		DPM_out<<time<<" "<<epsxy<<" "<<da_<<" "<<dn<<" "<<dt<<" "<<dl<<" "<<dc<<" "<<ds_<<endl;
		DPM_out.close();
	
	cout<<"		.5*(a+an+at+al+ac) = "<<.5*(a()+an+at+al+ac)<<endl;
	cout<<"		               q/p = "<<qop_<<endl;
	
	delete Cn;
	delete Fn;
	delete Ft;
	delete L;

}

//Rajouter calcul tenseur de contraintes avec les composantes propres de l'Žcoulement (xy,yy)
void shearP_CD_A::globalStress()
{
	cout<<"	Globalstress : ";
	gdm::Tensor2x2 * S = StressInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk());

	if (S != NULL ) 
	{
		S->eigenValues();
			//S->print();
		s1=S->l1();
		s2=S->l2();
		
		pressure_=.5*(s1+s2);
		q_ = .5*( max(s1,s2)-min(s1,s2));
		qop_ = q_/pressure_;
		pressure_defined=true;
		ds_=S->majorDirection();
		cout<<" direction = "<<ds_/M_PI*180.<<" q/p = "<<qop_;
		ofstream GS_out("Analyse/stress.txt",ios::app);
		GS_out<<time<<" "<<max(s1,s2)<<" "<<min(s1,s2)<<" "<<pressure_<<" "<<q_<<" "<<qop_<<" "<<ds_<<endl;
		GS_out.close();

	}

	delete S;
	cout<<endl;
}
void shearP_CD_A::partialLengthStress(unsigned int Npoint)
{
	cout<<"	Partialstress in length: ";
	gdm::Tensor2x2 * S = PartialLengthStressInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Npoint,lmoy_,sys_->spl()->rmax(),pressure_);	
	
	delete S;
	cout<<endl;
}
void shearP_CD_A::partialNormalForceStress(unsigned int Npoint)
{
	cout<<"	Partialstress in normal force: ";
	cout<<" Pressure "<<pressure_<<endl;
	gdm::Tensor2x2 * S = PartialNormalForceStressInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Npoint,fnmoy_,pressure_);	
	
	delete S;
	cout<<endl;
}

unsigned int shearP_CD_A::Z(bool exportZP, unsigned int Nbin, bool cvd)
{

	unsigned int i,ni;
	unsigned int Nb = sys_->spl()->lbody().size();
	unsigned int Nc = sys_->nwk()->clist().size();
	unsigned int Nz = 0; //Number of particles without dof with contacts ( more to one)
	unsigned int Nf = 0; //Number of "floatting" particles ( less than two contacts )

//	sys_->spl()->radiusExtrema();
	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();

	unsigned int Ncwf = 0;//Number of contacts with force 
	unsigned int Ncf0 =0; // Number of contats in probe with null force
	vector < unsigned int> Ncpp(Nb,0);//Number of contact per particles
	vector < bool > Activate(Nb,false);
	vector< unsigned int > Nw(Nb,0); 

	unsigned int NcppMax;
	unsigned int NcZ=0;
	unsigned int id1,id2;
	double fmoy=0;
	cout<<"	Z : "<<flush;
	NcppMax=0;
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		//if( totalProbe_.contain( sys_->nwk()->inter(ni)))
		//{
		id1=sys_->nwk()->inter(ni)->first()->id();
		id2=sys_->nwk()->inter(ni)->second()->id();


		if (sys_->nwk()->inter(ni)->fn() == 0.)
		{
			++Ncf0;
				//cout<<" force nulle "<<ni<<endl;
		}
		else
		{
			++Ncwf;
			fmoy+=sys_->nwk()->inter(ni)->fn();
			Ncpp[id1]+=1;
			Ncpp[id2]+=1;
			NcppMax=max(Ncpp[id1],NcppMax);
			NcppMax=max(Ncpp[id2],NcppMax);
				//cout<<sys_->nwk()->inter(ni)->type()<<" ";
				//cout<<sys_->nwk()->inter(ni)->first()->id()<<" "<<sys_->nwk()->inter(ni)->second()->id()<<endl;
		}
		//}

	}
		//double Pf0= (double) (Ncf0)/(double) (Ncf0 + Ncwf);
	fmoy/=(double) (Ncwf);
	unsigned int Ncw=0;
//	cout<<" fmoy = "<<fmoy<< " ";
	//cout<<" Nforce nulle = "<<Ncf0<<" Ncwf = "<<Ncwf<<endl;
	if (Ncwf==0) 
		{ cout<<" no forces "<<endl;return 0;}
	ofstream scalZ("Analyse/Connect/scalarZ.txt",ios::app);
	scalZ<<time<<" "<<Ncwf<<" "<<Ncf0<<" "<<flush;

	//Partial connectivity
	unsigned int N=0;
	if ( exportZP)
	{
		cout<<"  / ZP "<<flush;
		vector <unsigned int> Npnc(NcppMax+1,0);//Number of particles with n contacts
		for(i=0;i<Nb;++i)
		{ 
			if ( sys_->spl()->body(i)->bodyDof()==NULL)
			{
				Npnc[Ncpp[i]]+=1;
				N++;
				//if( Ncpp[i]==2) cout<<i<<endl;
			}
		}
		//cout<<" N = "<<N<<endl;
		//A prevoir: la numerotation successive des fichier zp
		ofstream out("Analyse/Connect/ZPt.txt",ios::app);
		ofstream out1("Analyse/Connect/ZP.txt",ios::out);
		out<<time<<" ";

		vector <double> Zp(NcppMax+1);
		for(i=0;i<NcppMax;++i)
		{ 
			Zp[i]=(double) (Npnc[i])/(double) (N);

			out<<Zp[i]<<" ";
			out1<<i<<" "<<Zp[i]<<endl;
			//cout<<i<<" "<<Zp[i]<<endl;

		}
		out<<endl;
		out.close();
		out1.close();

		Npnc.clear();
		Zp.clear();
	}
	Ncpp.clear();

	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		//if( totalProbe_.contain( sys_->nwk()->inter(ni)))
		//{
		id1=sys_->nwk()->inter(ni)->first()->id();
		id2=sys_->nwk()->inter(ni)->second()->id();

		if ( Ncpp[id1]>1)
		{
			if( sys_->nwk()->inter(ni)->fn() < fmoy)
			{
				Nw[id1]+=1;
				Ncw++;
			}
		}

		if ( Ncpp[id2]>1)
		{
			if( sys_->nwk()->inter(ni)->fn() < fmoy)
				Nw[id2]+=1;
		}

		//}

	}
	//cout<<" Fraction of weak contact : "<<(double) (Ncw)/(double)(Ncwf);
	pointSet fracw;
	for(i=0;i<Nb;++i)
	{
		if( Ncpp[i]>1 && sys_->spl()->body(i)->bodyDof()==NULL )
			fracw.add( sys_->spl()->body(i)->sizeVerlet(), (double) (Nw[i])/(double)(Ncpp[i]));
	}
	
	fracw.decreasingXSort();
	//pointSet weak= fracw.histoBinVar(200);
	pointSet weak= fracw.slidingHisto(50,.10);
	if ( exportDistribution )
	{	
		fracw.xoffSet(-rmin);
		fracw.xNormalise( rmax-rmin);
		fracw.write("Analyse/Connect/fracweak_brut.txt");
	}
	fracw.psetClear();
	//pointSet mean= weak.mobileMean(1,2);
	weak.xoffSet(-rmin);
	weak.xNormalise( rmax-rmin);
	weak.write("Analyse/Connect/fracweak.txt");


	Nz=0;
	NcppMax=0;
	NcZ=0;
	double r;
	pointSet gcu,gc,vu,nu;
	double masstu=0,masst=0;

	rattlersB.clear();
	rattlersVisu.clear();
	rattlersVisu.resize( Nb);

	for(i=0;i<Nb;++i)
	{ 
		if (sys_->spl()->body(i)->bodyDof()==NULL)
		{
			r=sys_->spl()->body(i)->sizeVerlet();

			gc.add(r , r*r);
			masst +=r*r;

			if (Ncpp[i] > 1    )    
			{	
				NcppMax=max(Ncpp[i],NcppMax);
				Nz++;
				NcZ += Ncpp[i];
				gcu.add( r, r*r);
				nu.add( r, 1.);
					//vu.add( r, r*r);
				masstu +=r*r;
				rattlersVisu[i]=false;
			}
			else if( Ncpp[i]<=1 )
			{
				++Nf;
				nu.add(r, 0.);
				rattlersB.push_back(sys_->spl()->body(i) );
				rattlersVisu[i]=true;
				//cout<<" nf i= "<<i<<endl;
			}

		}
	}

	ratlers_=(double)(Nf)/(double)(Nf+Nz);
	z()=(double) (NcZ)/(double) (Nz);
	cout<<endl<<"	 Nz = "<<Nz<<" Ncz = "<<NcZ<<" Nf = "<<Nf<<endl;
	cout<<" 	z= "<<z()<<" ratlers_= "<<ratlers_<<endl;
	
	scalZ<<Nz<<" "<<NcZ<<" "<<Nf<<" "<<ratlers_<<" "<<z()<<endl;

	if (cvd )
	{

		//cout<<" mt "<<masst<<" mtu "<<masstu<<endl;
		gcu.increasingXSort();
		gc.increasingXSort();

		pointSet histvolutil = gcu.histoNumberBin(Nbin);
		pointSet histvoltot  = gc .histoNumberBin(Nbin);


		pointSet temp;

		for( unsigned int i=0; i< Nbin; ++i)
		{
			//cout<<histvoltot.py(i)<<" "<<histvolutil.py(i)<<endl;
			temp.add( histvoltot.px(i), histvolutil.py(i)/histvoltot.py(i));
		}
		temp.write("Analyse/Connect/fvgu.txt");

		gcu.yNormalise( masstu);
		gcu.yCumul();
		gcu.write("Analyse/Connect/granutilZ.txt");
		gcu.xoffSet(-rmin);
		gcu.xNormalise( rmax-rmin);
		gcu.write("Analyse/Connect/granutilZnorm.txt");
		
		ofstream zg_out("Analyse/Connect/Granulo.txt",ios::app);
		unsigned int i=0;
		double d1,d3,d6;
		d1=d3=d6=0;
		while ( i<gcu.psetSize() )
		{
			if( gcu.py(i) >.1 && d1 ==0)
				d1 = gcu.px(i);
			if( gcu.py(i) >.3 && d3 ==0)
				d3 = gcu.px(i);
			if( gcu.py(i) >.6 && d6 ==0)
			{
				d6 = gcu.px(i);
				break;
			}
			++i;
		}
		zg_out<<time<<" "<<d1<<" "<<d3<<" "<<d6<<" ";
		
		gc.yNormalise( masst);
		gc.yCumul();
		gc.write("Analyse/Connect/grantot.txt");
		
		gc.xoffSet(-rmin);
		gc.xNormalise( rmax-rmin);
		gc.write("Analyse/Connect/grantotnorm.txt");
		
		i=0;
		d1=d3=d6=0;
		
		while ( i< gc.psetSize() )
		{
			if( gc.py(i) >.1 && d1 ==0)
				d1 = gc.px(i);
			if( gc.py(i) >.3 && d3 ==0)
				d3 = gc.px(i);
			if( gc.py(i) >.6 && d6 ==0)
			{
				d6 = gc.px(i);
				break;
			}
			++i;
		}
		zg_out<<d1<<" "<<d3<<" "<<d6<<endl;
		
		nu.decreasingXSort();
		
		pointSet histonu = nu.histoBinVar(Nbin);
		
		if ( exportDistribution )
		{	
			nu.write("Analyse/Connect/gu_brut.txt");
		}
		
		histonu.write("Analyse/Connect/fgu.txt");//fraction granulo utile
		histonu.xoffSet(-rmin);
		histonu.xNormalise( rmax-rmin);
		histonu.write("Analyse/Connect/fgu_adim.txt");//fraction granulo utile
		
		nu.xoffSet(-rmin);
		nu.xNormalise( rmax-rmin);
		if ( exportDistribution )
		{	
			nu.write("Analyse/Connect/gu_adim_brut.txt");
		}
		
		nu.psetClear();
		//histonu.psetClear();
		histvolutil.psetClear();
		histvoltot.psetClear();
		
		
	}

	//cout<<" pf01 = "<<Pf0<<" pf02 "<<(double) (Ncf0)/(double) (Ncf0 + Ncwf- Nf)<<endl;
	

	//Connectivity as a function of granulometry
	if( calczg)
	{
		cout<<" / Zg "<<flush;
		vector < DataSet*> DZ(Nbin,NULL);
		vector < DataSet*> R(Nbin,NULL); 
		vector < unsigned int > Nrattlers(Nbin,0); 
		vector < unsigned int > Ncont(Nbin,0); 
		
		for(  i=0;i<Nbin;++i)
		{
			DZ[i]=new DataSet;
			R[i]=new DataSet;
		}
		vector <DataSet*> ZD(NcppMax+1,NULL);
		for(  i=0;i<NcppMax+1;++i)
		{
			ZD[i]=new DataSet;
		}

		ofstream data("Analyse/Connect/Zid.txt",ios::out);

		double rmax=sys_->spl()->rmax();
		double rmin=sys_->spl()->rmin();
		unsigned int jc;
		double amp=(rmax-rmin )/(double) (Nbin);
		pointSet PDZ;

		for(i=0;i<Nb;++i)
		{ 
			if ( sys_->spl()->body(i)->bodyDof()==NULL )
			{	
				r=sys_->spl()->body(i)->sizeVerlet();
				jc = (unsigned int) floor ( (r-rmin)/amp);

				if ( fabs(r - rmax)< 1e-10 )  jc=Nbin-1;
				if ( fabs(r - rmin)< 1e-10 )  jc=0;
				Ncont[jc]++;
				if( Ncpp[i]>1)
				{

					//cout<<i<<" "<<r<<endl;
					data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<Ncpp[i] <<endl;

					ZD[Ncpp[i]]->add( r );

						//r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
					DZ[jc]->add( Ncpp[i]);
					R [jc]->add( r );
					PDZ.add( r,Ncpp[i]);

				}
				else
				{
					Nrattlers[jc]++;
				}
			}

		}
		data.close();
		
		ofstream data3("Analyse/Connect/RattlersR.txt",ios::out);
		for(  i=0;i<Nbin;++i)
		{
			data3<<((i+.5)*amp)/(rmax-rmin)<<" "<<(i+.5)*amp+rmin<<" "<< (double) (Nrattlers[i])/(double)(Ncont[i])<<endl;
		}
		data3.close();
		cout<<" ok "<<endl;


		bool increase=false;
		if( increase)
		{
			PDZ.  increasingXSort();
		}
		else
		{
			PDZ.  decreasingXSort();
		}

		pointSet sli = PDZ.slidingHisto(50,.05);
		//PDZ.write("connect/testZ.txt");
		sli.xoffSet(-rmin);
		sli.xNormalise( rmax-rmin);
		sli.write("Analyse/Connect/sliDZ.txt");

		ofstream output ( "Analyse/Connect/DZ2.txt" , ios::out );
		if ( ! output ) cout<<"erreur creation de DZ2"<<endl;


		DataSet dPZ,dR;
		unsigned int N=0,Nev=0,Nepc;

		Nepc=PDZ.psetSize()/Nbin;
		//cout<<"Nepc = "<<Nepc<<" "<<flush;

		while( N < PDZ.psetSize())
		{

			if( Nev < Nepc)
			{
				dPZ.add( PDZ.py(N));
				dR.add( PDZ.px(N));			
					//cout<<RP.px(N)<<" "<<RP.py(N)<<endl;
				N++;
				Nev++;
			}
			else
			{
					//cout<<P.setSize()<<endl;
				dPZ.extractValues();
				dR.extractValues();

				output<<dR.mean()<<" "<<(dR.mean()-rmin)/(rmax-rmin)<<" "<< dR.variance()<<" "
					<<dPZ.mean()<<" "<< dPZ.variance()<<endl;

				dPZ.setClear();
				dR.setClear();
				Nev=0;
			}
		}
		output.close();

		ofstream data1("Analyse/Connect/ZD.txt",ios::out);
		for( i=0;i<NcppMax;++i)
		{
			ZD[i]->extractValues();
			if( ZD[i]->mean() != 0. )
				data1<<i<<" "<<ZD[i]->mean()<<" "<<ZD[i]->variance()<<endl;
		}
		data1.close();

		ofstream data2("Analyse/Connect/DZ.txt",ios::out);
		for(i=0;i<Nbin;++i)
		{
			DZ[i]->extractValues();
			R[i]->extractValues();
			data2<<(R[i]->mean()-rmin)/(rmax-rmin)<<" "<<R[i]->mean()<<" "<<DZ[i]->mean()<<" "<<R[i]->variance()<<" "<<DZ[i]->variance()<<endl;

		}
		data2.close();


	}


	cout<<endl;
	return 1;

}

void shearP_CD_A::granuloSpeed(unsigned int Nbin)
{
	cout<<" GranuloSpeed : ";
	unsigned int Nb = sys_->spl()->lbody().size();
	vector < DataSet*> S(Nbin,NULL);
	vector < DataSet*> R(Nbin,NULL);
	double v,r;
	double amp = (sys_->spl()->rmax()-sys_->spl()->rmin() )/(double) (Nbin);
	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int jc;


	for( unsigned int i=0;i<Nbin;++i)
	{
		S[i]=new DataSet;
		R[i]=new DataSet;
	}
	for( unsigned int i=0; i < Nb;++i)
	{
		if(totalProbe_.containCenter(sys_->spl()->body(i)))
		{
			r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
			jc = (unsigned int) floor ( (r-rmin)/amp);

			if ( fabs(r - rmax)< 1e-10 )  jc=Nbin-1;
			if ( fabs(r - rmin)< 1e-10 )  jc=0;

			v=sqrt( pow(sys_->spl()->body(i)->vx(),2.)+pow(sys_->spl()->body(i)->vy(),2.) );			
			if(v != 0.)
			{
				S[jc]->add(v);
				R[jc]->add(r);
			}
		}
	}

	ofstream output ( "granulospeed.txt" , ios::out | ios::trunc);
	if ( ! output ) cout<<"erreur creation de granulospeed"<<endl;
	double vdof=sqrt( pow( partref->vx(),2) + pow( partref->vy(),2));
	Vmoy_=0.;
	for( unsigned int i=0;i<Nbin;++i)
	{
		if( S[i]->setSize() != 0)
		{
			S[i]->extractValues();
			R[i]->extractValues();
			output<<R[i]->mean()<<" "<<S[i]->min()<<" "<<S[i]->mean()<<" "<<S[i]->max()<<" "<<vdof<<endl;
			Vmoy_+=S[i]->mean();
		}
	}
	Vmoy_/=(double) (Nbin)*vdof;
	output.close();
	cout<<" OK "<<endl;

}

unsigned int shearP_CD_A::granuloStress(unsigned int Nbinstress )
{
	cout<<"	Granulostress : ";
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr );

	for( unsigned int i=0; i < Nb;++i)
	{
		//rajouter une multi par le volume specifique =1/compacitŽ
		bodstr[i].eigenValues();
		if ( bodstr[i].l1() !=0. ) noIMT=false;
	}
	if (noIMT) 
	{
		cout<<" All IMT are Null "<<endl;
		return 0;
	}
	//Granulometry dependence analysis

	double amp = (sys_->spl()->rmax()-sys_->spl()->rmin() )/(double) (Nbinstress);
	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	unsigned int jc;
	double r;
	DataSet rglob;
	vector < DataSet*> P(Nbinstress,NULL);
	vector < DataSet*> Q(Nbinstress,NULL); 
	vector < DataSet*> R(Nbinstress,NULL);
	vector < DataSet*> QP(Nbinstress,NULL);
	vector < DataSet*> DPM(Nbinstress,NULL); 

	for( unsigned int i=0;i<Nbinstress;++i)
	{
		P[i]=new DataSet;
		Q[i]=new DataSet;
		R[i]=new DataSet;
		QP[i]=new DataSet;
		DPM[i]=new DataSet;
	}

	double qtemp,ptemp,phase,dpm;
	//	double halfPI=M_PI*.5;
	ofstream data("qpid.txt",ios::out);
	for( unsigned int i=0; i < Nb;++i)
	{
		if(totalProbe_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
		{
			r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
			//cout<<" r = "<<(r-rmin)/amp<<endl;
			jc = (unsigned int) floor ( (r-rmin)/amp);
			//cout<<jc<<endl;
			//cout<<" r-rmax = "<<endl;
			if ( fabs(r - sys_->spl()->rmax())< 1e-10 )  jc=Nbinstress-1;
			if ( fabs(r - sys_->spl()->rmin())< 1e-10 )  jc=0;

			//cout<<jc<<" "<<r - sys_->spl()->rmax()<<endl;
			phase=1.;//cos( 2.*(bodstr[i].majorDirection()-halfPI));
			//cout<<bodstr[i].majorDirection()/M_PI*180.<<endl;
			qtemp=max(bodstr[i].l1(),bodstr[i].l2()) - min(bodstr[i].l1(),bodstr[i].l2());
			ptemp=bodstr[i].l1() + bodstr[i].l2();
			dpm=bodstr[i].majorDirection();
			//cout<<qtemp<<" "<<ptemp<<endl;

			data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<
				ptemp<<" "<<qtemp <<" "<<endl;
			if(ptemp>0) cout<<ptemp<<endl;

			if( qtemp!=0.)
			{			
				P[jc]->add( .5 * ptemp *phase);
				Q[jc]->add( .5 * qtemp *phase);
				R[jc]->add( r );
				QP[jc]->add( qtemp/ptemp * phase);
				DPM[jc]->add( dpm );
				rglob.add(r);
			}
			else cout<<" qtemp nul"<<endl;
		}
	}
	data.close();


	ofstream output ( "stressgranulomoy.txt" , ios::out | ios::trunc);
	if ( ! output ) cout<<"erreur creation de stressgranulomoy"<<endl;

	else
	{
		output.setf( ios::scientific);
		for (unsigned int i=0;i<Nbinstress;++i)
		{
			P[i]->extractValues();
			Q[i]->extractValues();
			QP[i]->extractValues();
			R[i]->extractValues();
			DPM[i]->extractValues();

			output<<R[i]->mean()<<" "<<(R[i]->mean()-rmin)/(rmax-rmin)<<" "<< R[i]->variance()<<" "
				<<Q[i]->mean()<<" "<< Q[i]->variance()<<" "
				<<P[i]->mean()<<" "<< P[i]->variance()<<" "
				<<DPM[i]->mean()<<" "<< DPM[i]->variance()<<" "
				<<QP[i]->mean()<<" "<< QP[i]->variance()<<" "<<endl;

		}
	}
	output.close();

	double mass=0;
	DataSet granulo;
	rglob.Sort();
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		mass+=rglob.set(i)*rglob.set(i);
	}
	ofstream gran("granutilS.txt",ios::out);
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		granulo.add( rglob.set(i)*rglob.set(i)/mass);
		gran<<rglob.set(i)<<" "<<granulo.set(i)<<endl;
	}
	gran.close();



	for( unsigned int i=0;i<Nbinstress;++i)
	{
		delete P[i];
		P.clear();
		delete Q[i];
		Q.clear();
		delete R[i];
		R.clear();
		delete DPM[i];
		DPM.clear();
		delete QP[i];
		QP.clear();
	}

	return 1;
}

unsigned int shearP_CD_A::granuloStress2(unsigned int Nc ,unsigned int period, unsigned int width)
{
	cout<<"	Granulostress2 : "<<flush;
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr );
	cout<<" imt ok"<<endl;
	for( unsigned int i=0; i < Nb;++i)
	{
	//rajouter une multi par le volume specifique =1/compacitŽ
		bodstr[i].eigenValues();
		if ( bodstr[i].l1() !=0. ) noIMT=false;
	}
	if (noIMT) 
	{
		cout<<" All IMT are Null "<<endl;
		return 0;
	}
//Granulometry dependence analysis

	pointSet RP,RQ,RQP,RDPM;

	double qtemp,ptemp,phase,dpm;
	double r,rmin,rmax;
	PGTM_.clear();
	for( unsigned int i=0; i < Nb;++i)
	{
		if(totalProbe_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
		{
			r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();

			phase=1.;//cos( 2.*(bodstr[i].majorDirection()-halfPI));
		//cout<<bodstr[i].majorDirection()/M_PI*180.<<endl;
			qtemp= -.5*(max(bodstr[i].l1(),bodstr[i].l2()) - min(bodstr[i].l1(),bodstr[i].l2()));
			ptemp= -.5*(bodstr[i].l1() + bodstr[i].l2());
			dpm=bodstr[i].majorDirection();

			if( pressure_defined)
			{
				if( ptemp > .8* pressure_)
					PGTM_.push_back(true);
				else PGTM_.push_back(false);
			}

			if( qtemp!=0.)
			{			
				RP.add(r,ptemp);
				RQ.add(r,qtemp);
				RQP.add(r,qtemp/ptemp);
				RDPM.add(r,dpm);
				//cout<<dpm<<endl;
			}
			else cout<<" qtemp nul"<<endl;
		}
		else if( pressure_defined) PGTM_.push_back(false);
	}

	bool increase=false;
	if( increase)
	{
		RP.  increasingXSort();
		RQ.  increasingXSort();
		RQP. increasingXSort();
		RDPM.increasingXSort();

		rmin = RP.px(0);
		rmax = RP.px(RP.psetSize()-1);
	}
	else
	{
		RP.  decreasingXSort();
		RQ.  decreasingXSort();
		RQP. decreasingXSort();
		RDPM.decreasingXSort();

		rmax = RP.px(0);
		rmin = RP.px(RP.psetSize()-1);
	}
	ofstream output ( "stressgranulomoy.txt" , ios::out | ios::trunc);
	if ( ! output ) cout<<"erreur creation de stressgranulomoy"<<endl;


	DataSet P,R,Q,QP,DPM,rglob;
	pointSet P2,P3,Q2,Q3,QP2,QP3,DPM2,DPM3;
	unsigned int N=0,Nev=0,Nepc;

	Nepc=RP.psetSize()/Nc;
	//cout<<"Nepc = "<<Nepc<<" ";

	while( N < RP.psetSize())
	{

		if( Nev < Nepc)
		{
			P.add( RP.py(N));
			Q.add( RQ.py(N));
			QP.add( RQP.py(N));
			DPM.add( RDPM.py(N));
			R.add( RP.px(N));
			rglob.add(RP.px(N));			
			//cout<<RP.px(N)<<" "<<RP.py(N)<<endl;
			N++;
			Nev++;
		}
		else
		{
			//cout<<P.setSize()<<endl;
			P.extractValues();
			Q.extractValues();
			QP.extractValues();
			DPM.extractValues();
			R.extractValues();
			P2.add( R.mean(),P.mean() );
			P3.add( (R.mean()-rmin)/(rmax-rmin), P.mean() );
			Q2.add( R.mean(),Q.mean() );
			Q3.add( (R.mean()-rmin)/(rmax-rmin), Q.mean() );
			QP2.add( R.mean(),QP.mean() );
			QP3.add( (R.mean()-rmin)/(rmax-rmin), QP.mean() );
			DPM2.add( R.mean(),DPM.mean() );
			DPM3.add( (R.mean()-rmin)/(rmax-rmin), DPM.mean() );

			output<<R.mean()<<" "<<(R.mean()-rmin)/(rmax-rmin)<<" "<< R.variance()<<" "
				<<Q.mean()<<" "<< Q.variance()<<" "
				<<P.mean()<<" "<< P.variance()<<" "
				<<DPM.mean()<<" "<< DPM.variance()<<" "
				<<QP.mean()<<" "<< QP.variance()<<" "<<endl;

			P.setClear();
			Q.setClear();
			QP.setClear();
			DPM.setClear();
			R.setClear();
			Nev=0;
		}
	}
	pointSet P2m = P2.mobileMean( period, width );
	pointSet P3m = P3.mobileMean( period, width );
	P2m.write("pg_mob.txt");
	P3m.write("pga_mob.txt");
	pointSet Q2m = Q2.mobileMean( period, width );
	pointSet Q3m = Q3.mobileMean( period, width );
	Q2m.write("qg_mob.txt");
	Q3m.write("qga_mob.txt");
	pointSet QP2m = QP2.mobileMean( period, width );
	pointSet QP3m = QP3.mobileMean( period, width );
	QP2m.write("qpg_mob.txt");
	QP3m.write("qpga_mob.txt");
	pointSet DPM2m = DPM2.mobileMean( period, width );
	pointSet DPM3m = DPM3.mobileMean( period, width );
	DPM2m.write("DPMg_mob.txt");
	DPM3m.write("DPMga_mob.txt");


	double mass=0;
	DataSet granulo;
	rglob.Sort();
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		mass+=rglob.set(i)*rglob.set(i);
	}
	ofstream gran("granutilS.txt",ios::out);
	double massc=0;
//	pointSet gr;
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		massc+=rglob.set(i)*rglob.set(i);
		granulo.add( massc/mass);
	//	gr.add( rglob.set(i), massc/mass);

		gran<<rglob.set(i)<<" "<<granulo.set(i)<<endl;
	}
	gran.close();

	return 1;
}

unsigned int shearP_CD_A::granuloStress3(unsigned int Nbinstress )
{
	system("mkdir Analyse/granuloS");
	cout<<"	Granulostress : ";
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr,sf_ );

	for( unsigned int i=0; i < Nb;++i)
	{
		//rajouter une multi par le volume specifique =1/compacitŽ
		bodstr[i].eigenValues();
		if ( bodstr[i].l1() !=0. ) noIMT=false;
	}
	if (noIMT) 
	{
		cout<<" All IMT are Null "<<endl;
		return 0;
	}
	//Granulometry dependence analysis

	//double amp = (sys_->spl()->rmax()-sys_->spl()->rmin() )/(double) (Nbinstress);
	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();
	//unsigned int jc;
	double r;
	DataSet rglob;

	pointSet P;
	pointSet Q; 
	pointSet QP;
	pointSet DPM; 

	double qtemp,ptemp,phase,dpm;
	//	double halfPI=M_PI*.5;
	ofstream data("Analyse/granuloS/qpid.txt",ios::out);
	for( unsigned int i=0; i < Nb;++i)
	{
		if(totalProbe_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
		{
			r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
			//cout<<" r = "<<(r-rmin)/amp<<endl;
			//cout<<jc<<endl;
			//cout<<" r-rmax = "<<endl;

			//cout<<jc<<" "<<r - sys_->spl()->rmax()<<endl;
			phase=1.;//cos( 2.*(bodstr[i].majorDirection()-halfPI));
			//cout<<bodstr[i].majorDirection()/M_PI*180.<<endl;
			qtemp=.5*(max(bodstr[i].l1(),bodstr[i].l2()) - min(bodstr[i].l1(),bodstr[i].l2()));
			ptemp=.5*(bodstr[i].l1() + bodstr[i].l2());
			dpm=bodstr[i].majorDirection();
			//cout<<qtemp<<" "<<ptemp<<endl;

			data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<(r-rmin)/(rmax-rmin)
			<<" "<<-ptemp<<" "<<qtemp <<" "<<endl;
			
			if(ptemp>0) cout<<ptemp<<endl;

			if( qtemp!=0.)
			{			
				P.add(r, - ptemp *phase);
				Q.add(r,  qtemp *phase);
				QP.add(r, qtemp/ptemp * phase);
				DPM.add(r, dpm );
				rglob.add(r);
			}
			else cout<<" qtemp nul"<<endl;
		}
	}
	data.close();

	int width=1;
	P  .increasingXSort();
	pointSet ps= (P.histoBinVar(Nbinstress)).mobileMean(1,width);//P.slidingHisto(Nbinstress,width);
	
	ps.xoffSet(-rmin);
	ps.xNormalise( rmax-rmin);
	ps.write("Analyse/granuloS/PG_adim.txt");
	
	pointSet ps2= P.slidingHisto(Nbinstress,width);
	ps2.xoffSet(-rmin);
	ps2.xNormalise( rmax-rmin);
	ps2.yNormalise( pressure_ );
	ps2.write("Analyse/granuloS/PG_xyadim.txt");
	
	pointSet ps3= P.histoBinVar(Nbinstress);
	ps3.xoffSet(-rmin);
	ps3.xNormalise( rmax-rmin);
	ps3.yNormalise( pressure_ );
	ps3.write("Analyse/granuloS/PG_xyadimVar.txt");
	
	pointSet ps4= ps2.mobileMean(1,3);
	ps4.write("Analyse/granuloS/PG_xyadimVarMob.txt");

	Q .decreasingXSort();
	pointSet qs= Q.slidingHisto(Nbinstress,width);
	qs.xoffSet(-rmin);
	qs.xNormalise( rmax-rmin);
	qs.write("Analyse/granuloS/QG_adim.txt");
	qs.yNormalise(pressure_);
	qs.write("Analyse/granuloS/Qpressure_xyadim.txt");

	/*QP .decreasingXSort();
	pointSet pqs= QP.slidingHisto(Nbinstress,width);
	pqs.xoffSet(-rmin);
	pqs.xNormalise( rmax-rmin);
	pqs.write("Analyse/granuloS/PQG_adim.txt");
	pqs.yNormalise(qop_);
	pqs.write("Analyse/granuloS/QOPG_xyadim.txt");*/
	

	DPM .decreasingXSort();
	pointSet dpms= DPM.slidingHisto(Nbinstress,width);
	dpms.xoffSet(-rmin);
	dpms.xNormalise( rmax-rmin);
	dpms.write("Analyse/granuloS/DPMG_adim.txt");


	double mass=0;
	DataSet granulo;
	rglob.Sort();
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		mass+=rglob.set(i)*rglob.set(i);
	}
	ofstream gran("Analyse/granuloS/granutilS.txt",ios::out);
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		granulo.add( rglob.set(i)*rglob.set(i)/mass);
		gran<<rglob.set(i)<<" "<<granulo.set(i)<<endl;
	}
	gran.close();

	cout<<" OK"<<endl;
	return 1;
}

int shearP_CD_A::pdfforce(bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	
	
	char  fichier[100];
	char  fichierbrut[100];
	if( fn )
	{
		sprintf(fichier,"Analyse/pdf/pdffn.txt");
		sprintf( fichierbrut,"Analyse/pdf/fn_brut.txt");
	}
	else
	{
		sprintf(fichier,"Analyse/pdf/pdfft.txt");
		sprintf(fichierbrut,"Analyse/pdf/ft_brut.txt");
	}

	/*ofstream out (fichier,ios::out);
	if ( ! out ) cout<<"@ forces : file error "<<fichier<<endl;

	out.close();*/
	
	//cout<<" debut pdf"<<endl;
	cout<<"	PDF  : " ;
	DataSet f;
	for( unsigned int i=0;i< sys_->nwk()->linter().size();++i)
	{
		if( totalProbe_.contain( sys_->nwk()->inter(i) ))
		{
			if( fn )
			{
				if( sys_->nwk()->inter(i)->fn() != 0.) 
				{
					f.add(sys_->nwk()->inter(i)->fn());
			//mean+=sys_->nwk()->inter(i)->fn();
				}
			}
			else
				if( sys_->nwk()->inter(i)->ft() != 0.) 
			{
				f.add(fabs(sys_->nwk()->inter(i)->ft()));
			//mean+=sys_->nwk()->inter(i)->ft();
			}
		}
	}
	if (f.setSize() <100)
		{cout<<" no enough forces (min 100 events ) "<<f.setSize()<<endl; return 0;}

	f.extractValues();
	if (exportDistribution) f.write(fichierbrut);
	
	if( normalize ) 
		f.Normalize( f.mean() );
	cout<<" fmoy = "<<f.mean()<<" ";
	//f.write("forces.txt");
	

	//f.downFilter( .001);
	f.DecreasingSort();
	
		
	pointSet pdf0 = f.Rich_PDF( nbin );//nombre de classe 
	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffnrich.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfftrich.txt");
	pdf0.write(fichier);
	
	
	pointSet pdf = f.kernelPdf( nbin,.05);
	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffnker.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfftker.txt");
	pdf.write(fichier);
	
	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffn_kermob.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft_kermob.txt");
	pointSet pdfm = pdf.mobileMean(period, width);
	pdfm.write(fichier);


	pointSet pdfs = (f.slidingPdf( nbin) ).mobileMean(period,width);
	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffnslid.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfftslid.txt");
	pdfs.write(fichier);


	return 1;
}

int shearP_CD_A::pdflength(int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	cout<<"	PDF L : "<<flush ;
	
	const char * fichier,*fichierbrut;
	
	fichier="Analyse/pdf/pdfL.txt";
	fichierbrut = "Analyse/pdf/L_brut.txt";
	

	/*ofstream out (fichier,ios::out);
	if ( ! out ) cout<<"@ forces : file error "<<fichier<<endl;

	out.close();*/
	
	//cout<<" debut pdf"<<endl;
	double rmoy=0;
	double rmin,rmax;
	rmin = sys_->spl()->rmin();
	rmax = sys_->spl()->rmax();
	rmoy = sys_->spl()->rmoy();

	
	inter2d * interc;
	DataSet l,lmoy,lmax,lmin;
	double temp;
	for( unsigned int i=0;i< sys_->nwk()->clist().size();++i)
	{
		interc=sys_->nwk()->inter(sys_->nwk()->clist(i));
				
		
		if( totalProbe_.contain( interc ))
		{
			//disk only
			temp=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
			l.add(temp);
			lmoy.add(temp);
			lmin.add(temp);
			lmax.add(temp);
			
		}
	}
	if (l.setSize() <100)
		{cout<<" no enough length (min 100 events ) "<<l.setSize()<<endl; return 0;}
	cout<<" ok "<<endl;

	l.extractValues();
	
	lmoy.extractValues();
	lmoy.DecreasingSort();
	lmin.extractValues();
	lmin.DecreasingSort();
	lmax.extractValues();
	lmax.DecreasingSort();
	
	if (exportDistribution) l.write(fichierbrut);
	
	lmoy.Normalize(rmoy);
	pointSet pdfs2 = (lmoy.slidingPdf( nbin) ).mobileMean(period,width);
	fichier="Analyse/pdf/pdfLslid_normRmoy.txt";
	pdfs2.write(fichier);
	
	lmin.Normalize(rmin);
	pointSet pdfs3 = (lmin.slidingPdf( nbin) ).mobileMean(period,width);
	fichier="Analyse/pdf/pdfLslid_normRmin.txt";
	pdfs3.write(fichier);
	
	lmax.Normalize(rmax);
	pointSet pdfs4 = (lmax.slidingPdf( nbin) ).mobileMean(period,width);
	fichier="Analyse/pdf/pdfLslid_normRmax.txt";
	pdfs4.write(fichier);
	
	
	if( normalize ) 
		l.Normalize( l.mean() );
	cout<<" lmoy = "<<l.mean()<<" ";
	//f.write("forces.txt");
	

	//f.downFilter( .001);
	l.DecreasingSort();
	
	pointSet pdf0 = l.Rich_PDF( nbin );//nombre de classe 
	
		fichier="Analyse/pdf/pdfLrich.txt";
	
	pdf0.write(fichier);
	
	
	pointSet pdf = l.kernelPdf( nbin,.05);
		fichier="Analyse/pdf/pdfLker.txt";
	pdf.write(fichier);
	
		fichier="Analyse/pdf/pdfL_kermob.txt";
	
	pointSet pdfm = pdf.mobileMean(period, width);
	pdfm.write(fichier);


	pointSet pdfs = (l.slidingPdf( nbin) ).mobileMean(period,width);
		fichier="Analyse/pdf/pdfLslid.txt";
	pdfs.write(fichier);

	cout<<"	OK "<<endl;
	return 1;
}


void shearP_CD_A::removeBody(body2d * rm)
{
	//Cette fonction est a complete, elle n'enleve pas les contacts 
	//dans la clist (dur avec les entier et le changement de taille de nwk)
	//--->il faudrait passer clist en pointeur d'interaction.
	vector< body2d*>::iterator itspl;
	vector< inter2d*>::iterator itnwk;
	//vector< unsigned int >::iterator itclist;	
	//unsigned int iclist=0;
	for(itnwk = sys_->nwk()->linter().begin() ; itnwk != sys_->nwk()->linter().end() ; )
	{
		if( (*itnwk)->first() == rm  || (*itnwk)->second()==rm )
		{
			itnwk = sys_->nwk()->linter().erase( itnwk);
		}
		else
		{	
			++itnwk;
		}
	}

	/*for(itclist = sys_->nwk()->clist().begin() ; itclist != sys_->nwk()->clist().end() ; )
	{
		if( (*itclist) == iclist )
		{
			itclist = sys_->nwk()->clist().erase( itclist);
		}
		else
		{	
			++itnwk;

		}
		}*/

			for(itspl = sys_->spl()->lbody().begin() ; itspl != sys_->spl()->lbody().end() ; )
		{
			if( (*itspl)==rm)
				itspl=sys_->spl()->lbody().erase( itspl);
			else
				itspl++;
		}
}

	void shearP_CD_A::removeRattlers()
	{
		cout<<" +++++++ WARNING : removing body : ";
		cout<<" taille echant : "<<sys_->spl()->lbody().size()<<" corps"<<endl;
		cout<<" taille reseau : "<<sys_->nwk()->linter().size()<<" inter"<<endl;

		for( unsigned int i=0; i< rattlersB.size();++i)
		{
			cout<<rattlersB[i]->id()<<" "<<flush;
			removeBody( rattlersB[i]);
		}
		cout<<" *********  SUPPRESSION  **********"<<endl;
		cout<<" taille echant : "<<sys_->spl()->lbody().size()<<" corps"<<endl;
		cout<<" taille reseau : "<<sys_->nwk()->linter().size()<<" inter"<<endl;
		sys_->spl()->updateBands();
		removeR=false;
	}

	void shearP_CD_A::reduceRattlers()
	{
		disk* temp;
		cout<<" +++++++ WARNING : reducing body : ";
		for( unsigned int i=0; i< rattlersB.size();++i)
		{
			cout<<rattlersB[i]->id()<<" ";
			temp= dynamic_cast< disk*>( rattlersB[i] );
			temp->R()=0.;
		//cout<<ran1(&seed)<<endl;
		//cout<<" Rayon nul pour le corps "<<rattlersB[i]->id()<<endl;
		}
		cout<<endl;

	}

	void shearP_CD_A::growRattlers( double dr)
	{
		disk* temp;
		cout<<" +++++++ WARNING : growing body : ";
		for( unsigned int i=0; i< rattlersB.size();++i)
		{
			cout<<rattlersB[i]->id()<<" "<<flush;
			temp= dynamic_cast< disk*>( rattlersB[i] );
			temp->R() += dr;
			temp->Fill(2650);

		}
		growR=false;
		sys_->spl()->radiusExtrema(0);
	}

void shearP_CD_A::writePS( const char * fname)
{
	//cout<<" WARNING : writePS running only with disk shaped body"<<endl;
	ofstream ps(fname);
	
	//cout << endl << "---\t " << sys_->spl()->xmin() << " " << sys_->spl()->xmax() << " " << sys_->spl()->boundWidth() << endl << endl; 
		
	//sys_->spl()->updateBoundaries();
	sys_->spl()->updateBands();
	sys_->spl()->radiusExtrema(1);
	
	double R = sys_->spl()->rmax();
	
	double Xmin = sys_->spl()->xmin() + R;
	double Ymin = sys_->spl()->ymin() + R;
	
	double xmin_ = sys_->spl()->xmin() - 5.*R;
	double ymin_ = sys_->spl()->ymin() - 5.*R;
	double xmax_ = sys_->spl()->xmax() + 5.*R;
	double ymax_ = sys_->spl()->ymax() + 5.*R;
	
	double zoom = zoom_;
	double x_offset = fabs(Xmin*zoom);
	double y_offset = fabs(Ymin*zoom);
	
	unsigned int id1,id2;
	
	//!------ Header for ps file	
	ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	//ps<<"%%BoundingBox:"<<" "<<"-30 -20 515 550"<<endl;
	ps<<"%%BoundingBox:"<<" "<<x_offset + (xmin_ - 2.*sys_->spl()->bandWidth())*zoom<<" "<<y_offset + (ymin_ - 2.*sys_->spl()->bandWidth())*zoom<<" "
	<<x_offset + (xmax_ + 2.*sys_->spl()->bandWidth())*zoom<<" "<<y_offset + (ymax_ + 2.*sys_->spl()->bandWidth())*zoom<<endl;
	ps<<"%%Pages: 1"<<endl;
	ps<<"0.1 setlinewidth 0. setgray "<<endl;
	
	double x_A,y_A,x_B,y_B,x_C,y_C,x_D,y_D;
	x_A = sys_->spl()->xmax() + R;//   sys_->spl()->body(3)->x()+R;
	y_A = sys_->spl()->ymin() - R;//   sys_->spl()->body(0)->y()-R;
	x_B = sys_->spl()->xmin() - R;//   sys_->spl()->body(2)->x()-R;
	y_B = sys_->spl()->ymin() - R;//   sys_->spl()->body(0)->y()-R;
	x_C = sys_->spl()->xmin() - R;//   sys_->spl()->body(2)->x()-R;
	y_C = sys_->spl()->ymax() + R;//   sys_->spl()->body(1)->y()+R;
	x_D = sys_->spl()->xmax() + R;//   sys_->spl()->body(3)->x()+R;
	y_D = sys_->spl()->ymax() + R;//   sys_->spl()->body(1)->y()+R;
	
	/*
	 //     ========1=========
	 //     =                =
 	 //     =                =
	 //     2                3
	 //     =                =
	 //     =                =
	 //     ========0=========
	 */
	
	
	
	//!------ Draw bounds of the sample		
	/*ps<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" moveto "<<x_B*zoom + x_offset<<" "<<y_B*zoom + y_offset<<" "<<" lineto "<<x_C*zoom + x_offset<<" "<<y_C*zoom + y_offset<<" lineto "<<x_D*zoom + x_offset<<" "
	 <<y_D*zoom + y_offset<<" lineto "<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" "<<" lineto stroke"<<endl;
	 */
	if (sys_->spl()->body(0)->type() ==2)
		ps<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" moveto "
		<<sys_->spl()->body(0)->xmin()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmin()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymax()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymax()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" lineto stroke"<<endl;
	
	//!------ Draw sample
	if(displaySample)
	{
		DataSet V;
		
		for(unsigned i=0 ; i<sys_->spl()->lbody().size() ; ++i)
		{
			//if (sys_->spl()->body(i)->type() ==0)
			{
				V.add(sqrt(sys_->spl()->body(i)->vx()*sys_->spl()->body(i)->vx()+sys_->spl()->body(i)->vy()*sys_->spl()->body(i)->vy()));
			}
		}
		
		V.extractValues();
		
		//!------ Draw Particles
		ps<<"/coul_disk_grey_ {1 setlinecap 0.75 0.75 0.75 setrgbcolor} def"<<endl;
		ps<<"/coul_disk_black {1 setlinecap 0.25 0.25 0.25 setrgbcolor} def"<<endl;
		ps<<"/coul_disk_green {1 setlinecap 0.00 0.76 0.00 setrgbcolor} def"<<endl;
		
		for (unsigned int j=0 ; j<sys_->spl()->lbody().size() ; ++j)
		{
			if (sys_->spl()->body(j)->type() ==0)
			{
				if (!displayForce)
				{
					if (sys_->spl()->body(j)->bodyDof()==NULL)
					{
						if (sys_->spl()->lbody().size()<= 5000)
						{
							ps	<<"/cirg {1 " 
							<<(sqrt(sys_->spl()->body(j)->vx()*sys_->spl()->body(j)->vx()+sys_->spl()->body(j)->vy()*sys_->spl()->body(j)->vy())-V.min())/(V.max()+V.min()) 
							<< " 0 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;//make circle 
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" cirg"<<endl;
						}
						else
						{
							ps	<< "/coul_disk {1 setlinecap 1 " 
							<<(sqrt(sys_->spl()->body(j)->vx()*sys_->spl()->body(j)->vx()+sys_->spl()->body(j)->vy()*sys_->spl()->body(j)->vy())-V.min())/(V.max()+V.min()) 
							<< " 0 setrgbcolor} def"<<endl;
							ps	<< "newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<< " 0.0 360.0 arc closepath 0.0 coul_disk stroke" << endl;
						}
						
					}
					else
					{
						if (sys_->spl()->lbody().size()<= 5000)
						{
							ps	<<"/cirg {0.25 0.25 0.25 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" cirg"<<endl;
						}
						else
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_black stroke" << endl;
					}
				}
				else
				{
					if (sys_->spl()->body(j)->bodyDof()==NULL)
					{
						if (sys_->spl()->lbody().size()<= 5000)
						{
							ps	<<"/cirg {0.75 0.75 0.75 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" cirg"<<endl;
						}
						else
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_grey_ stroke" << endl;
					}
					else
					{
						if (sys_->spl()->lbody().size()<= 5000)
						{
							ps	<<"/cirg {0.25 0.25 0.25 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;	
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" cirg"<<endl;
						}
						else
							ps	<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "
							<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_black stroke" << endl;
					}
				}
			}
		}
		
        
		for (unsigned int j=0 ; j<sys_->spl()->leftband().size() ; ++j)
		{
			if (sys_->spl()->body(sys_->spl()->leftband(j))->type() ==0)
			{
				if (sys_->spl()->body(sys_->spl()->leftband(j))->bodyDof()==NULL)
					ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->leftband(j))->x() + sys_->spl()->boundWidth())*zoom<<" "
					<<y_offset + sys_->spl()->body(sys_->spl()->leftband(j))->y()*zoom<<" "<<sys_->spl()->body(sys_->spl()->leftband(j))->sizeVerlet()*zoom<<" " 
					<<"0.0 360.0 arc closepath 0.0 coul_disk_green stroke" << endl;
				else
				{
					if (sys_->spl()->lbody().size()<= 5000)
						ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->leftband(j))->x() + sys_->spl()->boundWidth())*zoom<<" "
						<<y_offset + sys_->spl()->body(sys_->spl()->leftband(j))->y()*zoom<<" "
						<<sys_->spl()->body(sys_->spl()->leftband(j))->sizeVerlet()*zoom<<" cirg"<<endl;
					else
						ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->leftband(j))->x() + sys_->spl()->boundWidth())*zoom<<" "
						<<y_offset + sys_->spl()->body(sys_->spl()->leftband(j))->y()*zoom<<" "
						<<sys_->spl()->body(sys_->spl()->leftband(j))->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_black stroke" << endl;
				}
			}
		}
		
		for (unsigned int j=0 ; j<sys_->spl()->rightband().size() ; ++j)
		{
			if (sys_->spl()->body(sys_->spl()->rightband(j))->type() ==0)
			{
				if (sys_->spl()->body(sys_->spl()->rightband(j))->bodyDof()==NULL)
					ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->rightband(j))->x() - sys_->spl()->boundWidth())*zoom<<" "
					<<y_offset + sys_->spl()->body(sys_->spl()->rightband(j))->y()*zoom<<" "
					<<sys_->spl()->body(sys_->spl()->rightband(j))->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_green stroke" << endl;
				else
				{
					if (sys_->spl()->lbody().size()<= 5000)
						ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->rightband(j))->x() - sys_->spl()->boundWidth())*zoom<<" "
						<<y_offset + sys_->spl()->body(sys_->spl()->rightband(j))->y()*zoom<<" "
						<<sys_->spl()->body(sys_->spl()->rightband(j))->sizeVerlet()*zoom<<" cirg"<<endl;
					else
						ps	<<"newpath "<<x_offset + (sys_->spl()->body(sys_->spl()->rightband(j))->x() - sys_->spl()->boundWidth())*zoom<<" "
						<<y_offset + sys_->spl()->body(sys_->spl()->rightband(j))->y()*zoom<<" "
						<<sys_->spl()->body(sys_->spl()->rightband(j))->sizeVerlet()*zoom<<" 0.0 360.0 arc closepath 0.0 coul_disk_black stroke" << endl;
				}
			}
		}
	}
	
	//!------ Draw contact forces
	if(displayForce && !sys_->nwk()->clist().empty())
	{
		DataSet NEGATIVE_FN;
		DataSet POSITIVE_FN;
		
		sys_->spl()->radiusExtrema(2);
		
		for(unsigned i=0 ; i<sys_->nwk()->clist().size() ; ++i)
		{
			id1=sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id();
			id2=sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id();
			
			if(sys_->nwk()->inter(i)->type()==0 && 
			   sys_->spl()->body(id1)->bodyDof()==NULL &&
			   sys_->spl()->body(id2)->bodyDof()==NULL)
			{
				if(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()<0)
					NEGATIVE_FN.add(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
				else
					POSITIVE_FN.add(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
			}
		}
		
		ps<<"/coul_ligth_blue {1 setlinecap 0 1 1 setrgbcolor} def"<<endl;
		ps<<"/coul_blue {1 setlinecap 0 0 1 setrgbcolor} def"<<endl;
		ps<<"/coul_blanc {1 setlinecap 1 1 1 setrgbcolor} def"<<endl;
		
		ps<<0.0*(sys_->nwk()->inter(sys_->nwk()->clist(0))->fn())<<" setlinewidth coul_blanc"<<endl;
		ps<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->x()*zoom + x_offset<<" "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->y()*zoom + y_offset<<" "<<"moveto "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->second()->x()*zoom + x_offset<<" "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->second()->y()*zoom + y_offset<<" lineto stroke"<<endl;
		
		
		NEGATIVE_FN.extractValues();
		POSITIVE_FN.extractValues();
		
		for(unsigned i=0 ; i<sys_->nwk()->clist().size() ; ++i)
		{
			id1=sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id();
			id2=sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id();
			
			if(sys_->nwk()->inter(sys_->nwk()->clist(i))->type()==0 && 
			   sys_->spl()->body(id1)->bodyDof()==NULL &&
			   sys_->spl()->body(id2)->bodyDof()==NULL)
			{
				/*
				 ps<<0.01*(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn())<<" setlinewidth coul_blue"<<endl;
				 ps<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()*zoom + x_offset<<" "
				 <<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()*zoom + y_offset<<" "<<"moveto "
				 <<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()*zoom + x_offset<<" "
				 <<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()*zoom + y_offset<<" lineto stroke"<<endl;
				 */
				
				if(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn() < 0.0)
				{
					/*
					ps<<2.*(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()-NEGATIVE_FN.min())/(NEGATIVE_FN.min()+NEGATIVE_FN.max())<<" setlinewidth coul_ligth_blue"<<endl;
					ps<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()*zoom + x_offset<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()*zoom + y_offset<<" "<<"moveto"<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()*zoom + x_offset<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()*zoom + y_offset<<" "<<"lineto stroke"<<endl;
					*/
				}
				else 
				{
					ps<<10.*(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()-POSITIVE_FN.min())/(POSITIVE_FN.min()+POSITIVE_FN.max())<<" setlinewidth coul_blue"<<endl;
					ps<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()*zoom + x_offset<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()*zoom + y_offset<<" "<<"moveto"<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()*zoom + x_offset<<" "
					<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()*zoom + y_offset<<" "<<" "<<"lineto stroke"<<endl;
				}
				
			}
		}
		
		NEGATIVE_FN.setClear();
		POSITIVE_FN.setClear();
	}
	
	ps.close();
	
}


void shearP_CD_A::followparticles()
{

    cout<<".Followparticles"<<endl;
    
    unsigned int idref = 60 ;
    
    ofstream follow_("Analyse/particules.txt",ios::app);
	
	follow_.precision(10);
    
    for (unsigned int i=0; i<sys_->spl()->lbody().size();++i)
	{
        if(sys_->spl()->body(i)->id()==idref) follow_<< sys_->spl()->body(i)->sizeVerlet()<<" "<< sys_->spl()->body(i)->vx()<< " "<< sys_->spl()->body(i)->vy()<<endl;
        
    }
}
