#include "brazilian_A.hpp"

void brazilian_A::plugRef()
{
	top_    = dynamic_cast< rline*>(sys_->spl()->body(topId));
	bottom_ = dynamic_cast< rline*>(sys_->spl()->body(bottomId));
	left_   = dynamic_cast< rline*>(sys_->spl()->body(leftId));
	right_  = dynamic_cast< rline*>(sys_->spl()->body(rightId));
}

void brazilian_A::read_parameters(istream & is)
{
	string token,temp;
	is>>token;
	double trash;
	NbinPT=NbinFN=NbinFT=Nbingranulo=NbinFA=0;
	while(is)
	{
		if     (token== "/*")
		{
			is>>temp;
			while( temp != "*/")
			{
				cout<<" @shearP_CD_A :Ingored parameter: "<<temp<<endl;
				is >> temp;
			}
		}
		else if (token == "Sample")		displaySample = true;
		else if (token== "Cluster") calClust=true;
		else if(token== "Solidfraction") calcsf=true;
		else if(token== "Fabric") calcFabric=true;
		else if(token== "Fabricpolyg") {calcFabricPolyg=true;is >>NFabric;}
		else if(token== "Fractal") calcfracdim=true;
		else if(token== "Granulo") granulo=true;
		else if(token== "Forcesanisotropy") calcforcesA=true;
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
				cout<<" @brazilian_A :Zgranulo : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ brazilian_A :Zgranulo : Nbingranulo undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Granulostress") 
		{
			calcgranulostress=true;
			if (Nbingranulo !=0 )
			{
				cout<<" @brazilian_A :granulostress : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ brazilian_A :granulostress:  Nbingranulo undefined "<<endl;
				exit(0);
			}
			is>>perG;
			if ( perG == 0 )
			{
				cout<<" @ brazilian_A :granulostress:  perG undefined "<<endl;
				exit(0);
			}
			is>>wG;
			if ( wG == 0 )
			{
				cout<<" @ brazilian_A :granulostress:  wG undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="PDFFN") 
		{
			calcfn=true;
			is >> NbinFN;

			if ( NbinFN == 0 )
			{
				cout<<" @ brazilian_A : NbinFN undefined "<<endl;
				exit(0);
			}
			is >> perF;
			if ( perF == 0 )
			{
				cout<<" @ brazilian_A : perF undefined "<<endl;
				exit(0);
			}
			is >> wF;
			if ( wF == 0 )
			{
				cout<<" @ brazilian_A : wF undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Ptheta") 
		{
			calcPtheta=true;
			is >> NbinPT;

			if ( NbinPT == 0 )
			{
				cout<<" @ brazilian_A : NbinPT undefined "<<endl;
				exit(0);
			}
			is >> mobperiod;
			if ( mobperiod == 0 )
			{
				cout<<" @ brazilian_A : mobperiod undefined "<<endl;
				exit(0);
			}
			is >> mobwidth;
			if ( mobwidth == 0 )
			{
				cout<<" @ brazilian_A : mobwidth undefined "<<endl;
				exit(0);
			}
		}


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
				cout<<" @ brazilian_A : NbinFT undefined "<<endl;
				exit(0);
			}
		} 
		
	
		else if(token=="Globalstress") calcglobalstress=true;
		else if(token=="Forcescorrelation") calcFC=true;
		else if(token=="Sumforce") sumforce_defined=true;
		else if(token=="GranuloBin") 
		{
			if (Nbingranulo !=0 )
			{
				cout<<" @brazilian_A :granulobin : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ brazilian_A :granulobin : Nbingranulo undefined "<<endl;
				exit(0);
			}
		}
		else if(token=="Nanalyze") is>>Nanalyze();
		else if(token=="Remove")	removeR=true;
		else if( token=="RemoveVisu") removeVisu=true;


		else if(token=="Grow") 
		{
			growR=true;
			is>>incR;
		}
		else if(token=="Compactness") 
		{
			calcompact=true;
			is>>divisionl>>divisionh;
		}
		else if(token=="Fracture") 
			calfracture=true;
		else if(token=="}") break;
		//else if(token=="include_interrupt") incl_int=true;
		else cout<<token<<" : parametre de commande inconnu "<<endl;

		is>>token;
	}	
	
}

void brazilian_A::analyse( double t, unsigned int nsi, unsigned int nsf )
{
	cout<<"--------------brazilian_A::analyze()---------"<<endl;
	
	char fname[100];
	time=t;
	

	//Definir la probe avec les 4 premier corps du spl
	double hx=.5*(sys_->spl()->ymax()-sys_->spl()->ymin());
	double hy=.5*(sys_->spl()->xmax()-sys_->spl()->xmin());
	
	totalProbe_.x() = sys_->spl()->xmin()+ hx;
	totalProbe_.y() = sys_->spl()->ymin() + hy;
	totalProbe_.R() = 1.0*min(hx,hy);
	
	prb_.x() = sys_->spl()->xmin()+ hx;
	prb_.y() = sys_->spl()->ymin()+ hy;
	prb_.hh() = 1.*hy;
	prb_.hl() = 1.*hx;	
	
	//cout<<" aire "<<prb_.area()<< " "<<totalProbe_.area()<<endl;
	//cout<<" rmin "<<sys()->spl()->rmin()<<" rmax "<<sys()->spl()->rmax()<<endl;
	
	//def();
	//Gap();
	if (calcz)  Z(calczp,Nbingranulo,granulo);
	if( removeR) removeRattlers();
	if( growR) growRattlers( sys_->spl()->rmin()*incR);
	if (calcglobalstress) globalStress();
	if (calcfn) pdfforce(true,NbinFN,normpdf,perF,wF);
	if (calcft) pdfforce(false,NbinFT,normpdf,perF,wF);
	if (calcsf) SF();
	if (calcFabric)  A();
	if (calcforcesA)	forces_A(NbinFA);
	if (calcFabricPolyg) polygAnisotropy(NFabric);
	if (calcgranulostress) granuloStress3( Nbingranulo);//,perG,wG);

	if (calcgap) Gap();
	if (calcdef) def();

	if (calcPtheta) Ptheta( NbinPT, mobperiod, mobwidth);

	if (calcFC) forcesMaxCorrelation();
	//if (calcfracdim) BmultifractalBC();//NboxCounting();//BboxCounting();//fractalDimension( );//fractalDensity();
	if (calcgranulopdf) granulopdf(Nq_,true,NbinFN,normpdf,perF,wF);
	//writePS("Visu.ps");
	
	if (calClust) 
	{
		cluster();		
	}

	if(calcompact && (divisionh==1) && (divisionl==1))
	{
		sys_->spl()->updateBoundaries();
		
		double xmin=sys_->spl()->xmin();
		double xmax=sys_->spl()->xmax();
		double ymin=sys_->spl()->ymin();
		double ymax=sys_->spl()->ymax();
		double hight=ymax-ymin;
		double large=xmax-xmin;
		
		ofstream str("Analyse/compactness.txt",ios::app);
		str<<time<<"\t"<<compactness(hight,large)<<endl;
		str.close();
	}
	
	else if (calcompact && (nsi==nsf))
	{
		sys_->spl()->updateBoundaries();
		
		double xmin=sys_->spl()->xmin();
		double xmax=sys_->spl()->xmax();
		double ymin=sys_->spl()->ymin();
		double ymax=sys_->spl()->ymax();
		double hight=(ymax-ymin)/double(divisionh);
		double large=(xmax-xmin)/double(divisionl);
		
		ofstream str("Analyse/compactness.txt",ios::app);
		
		
		for (unsigned int i=1;i<=divisionl;i++)
			for(unsigned int j=1;j<=divisionh;j++)
			{
				double centrex=xmin+(i-0.5)*large;
				double centrey=ymin+(j-0.5)*hight;
				if (sys_->spl()->body(4)->type()==0)
					str<<time<<"\t"<<i<<"\t"<<j<<"\t"<<compactness_disk(centrex,centrey,hight,large)<<endl;
				else
					str<<time<<"\t"<<i<<"\t"<<j<<"\t"<<compactness(centrex,centrey,hight,large)<<endl;
			}
			
		str.close();
	}	
	cout<<endl;
	cout<<"		.5*(a+an+at) = "<<.5*(a()+an()+at())<<endl;
	cout<<"		         q/p = "<<(max(s1,s2)-min(s1,s2))/(s1+s2)<<endl;
	cout<<"		Nombre d'interaction = "<<sys_->nwk()->linter().size()<<endl;
	cout<<"		Nombre de contacts   = "<<sys_->nwk()->clist().size()<<endl;
	
		if( calcforcesA)
		{
			ofstream FA_out("Analyse/Anisotropy.txt",ios::app);
			FA_out<<time<<" "<<-defy<<" "<<Z_<<" "<<a()<<" "<<an()<<" "<<at()<<" "<<al()<<" "<<da_/M_PI*180.<<" "<<.5*( a()+an()+at() )<<" "<<qop_<<" "<<max(s1,s2)/sys_->grpRel()->fragmentation()->sigmac()<<endl;
			FA_out.close();

			ofstream DPM_out("Analyse/DPM.txt",ios::app);
			DPM_out<<time<<" "<<epsq<<" "<<da_<<" "<<dn_<<" "<<dt_<<" "<<dl_<<" "<<ds_<<endl;
			DPM_out.close();
	}
	
	ofstream an("Analyse/analyse.txt",ios::app);
	an.setf( ios::scientific);
	an<<time<<" "<<epsp<<" "<<epsq<<" "<<sf_<<" "<<Z_<<" "<<ratlers_<<" "<<(1-sf_)/sf_<<endl;
	
	if (displaySample)
	{
		sprintf(fname,"Analyse/PS1/particles%.4li.ps",numFile_);
		writePS(fname);	
	//	if (nsi==nsf)
	//	{
			sprintf(fname,"Analyse/PS2/particles%.4li.ps",numFile_);
			writePS2(fname);	
	//	}
		numFile_++;
	}
	
	if (calfracture) fracture();
}


void brazilian_A::initAnalyse( ) 
{ 
	cout<<"--------------brazilian_A::init()-----------"<<endl;	
	
	top_    = (dynamic_cast< brazilian*>( this->sys()))->top();
	bottom_ = (dynamic_cast< brazilian*>( this->sys()))->bottom();
	left_   = (dynamic_cast< brazilian*>( this->sys()))->left();
	right_  = (dynamic_cast< brazilian*>( this->sys()))->right();
	
	topId    = top_->id();
	bottomId = bottom_->id();
	leftId   = left_->id();
	rightId  = right_->id();
	
	cout<<" top "<<topId<<" bottom "<<bottomId<<" right "<<rightId<<" left "<<leftId<<endl;
	
	sys_->spl()->radiusExtrema(4);
	cout<<" rmin "<<sys()->spl()->rmin()<<" rmax "<<sys()->spl()->rmax()<<endl;
	
	this->sys()->spl()->updateBoundaries(); 		//Ajouter le 25/11/2011 	
	h0= sys_->spl()->ymax()-sys_->spl()->ymin();
	l0 = sys_->spl()->xmax()-sys_->spl()->xmin();
	
	defx=defy=0.;
	ngap1_=ngap2_=NULL;
	
	
	//double R = left_->R();
	double hy=.5*(sys_->spl()->ymax()-sys_->spl()->ymin());
	double hx=.5*(sys_->spl()->xmax()-sys_->spl()->xmin());
	
	totalProbe_.x() = sys_->spl()->xmin()+ hx;
	totalProbe_.y() = sys_->spl()->ymin() + hy;
	totalProbe_.R() = 1.0*min(hx,hy);
	
	prb_.x() = sys_->spl()->xmin()+ hx;
	prb_.y() = sys_->spl()->ymin()+ hy;
	prb_.hh() = 1.*hy;
	prb_.hl() = 1.*hx;	
	zoom_=min(610/hy/2.,710/hx/2.); //hinh chu nhat nam
	//zoom_=0.8*min(832/hy/2.,585/hx/2.);   //hcn dung
	
	system("mkdir Analyse");
	

	ofstream monitoring("Analyse/monitoring.txt",ios::out);
	monitoring.close();
	
	ofstream analyse("Analyse/analyse.txt",ios::out);
	analyse.close();

	ofstream strain("Analyse/strain.txt",ios::out);
	//strain<<"time	defx	defy	epsp	epsq	asin(epsp/epsq)	s1	s2	2*q	q/p	p	q"<<endl;
	strain.close();
	
	ofstream str("Analyse/compactness.txt",ios::out);
	str.close();
	// création du fichier measure
	//ofstream measure("Analyse/measure.txt",ios::out);
	//measure.close();
	
	if ( calcforcesA ) 
	{
			ofstream FA_out("Analyse/Anisotropy.txt",ios::out); 
			FA_out<<"# time epsq ac an at al .5*(ac+an+at)  q/p"<<endl;
			FA_out.close();
			ofstream DPM_out("Analyse/DPM.txt",ios::out); FA_out.close();
	}
	if( calcfn || calcft)
	{
		system("mkdir Analyse/pdf");
	}
	
	if (calcFabricPolyg)
	{
		system("mkdir Polygons");
		ofstream out("Polygons/anisotropy.txt",ios::out);out.close();
		ofstream out1("Polygons/dpm.txt",ios::out);out1.close();
		
	}
	if( calcglobalstress)
	{
		ofstream GS_out("Analyse/stress.txt",ios::out);
		//GS_out<<"time	s1	s2	max(s1,s2)-min(s1,s2)	epsq	qop_	(s1+s2)/2.	qop_*(s1+s2)/2.	ds_"<<endl;
		GS_out.close();
		
		ofstream U_compression_data_out("Analyse/U_compression.dat",ios::out);
		U_compression_data_out<<"#   P   rho   Nb_contact"<<endl;
		U_compression_data_out.close();
	}
	if( calcPtheta)
	{
		system("mkdir Analyse/angDistrib");
	}
	
	if (calcz)
	{
		system("mkdir Analyse/Connect");
		if (calczp) {ofstream zp_out("Analyse/Connect/Zpt.txt",ios::out); zp_out.close();}
		ofstream scalZ_out("Analyse/Connect/scalarZ.txt",ios::out); scalZ_out.close();
		if( granulo )
		{
			ofstream zg_out("Analyse/Connect/Granulo.txt",ios::out);zg_out.close();
		}		
	}
			
	if (calClust)
	
		{
			system("mkdir Analyse/Cluster");//creation répertoire Cluster
			system("mkdir Analyse/Cluster/pdf");
			//system("mkdir Analyse/Cluster/ptheta");
			
			
			ofstream clu_out("Analyse/Cluster/connectivity.txt",ios::out);
			clu_out <<" "<< endl;
			clu_out.close();
			
			ofstream Xcontact_out("Analyse/Cluster/Xcontact_type.txt",ios::out);
			Xcontact_out <<" # epsq ks kd kt" <<endl;
			Xcontact_out.close();
			
			/*ofstream Xdouble_out("Analyse/Cluster/Xdouble_type.txt ",ios::out);
			Xdouble_out<<"# epsq D DS"<<endl;
			Xdouble_out.close();*/
			
			//FMR Force Moyenne relative calculée avec les forces radiales
			ofstream FMR_out("Analyse/Cluster/meanRadial.txt",ios::out);
			FMR_out<<"# epsq fs fd ft"<<endl;
			FMR_out.close();
			
			ofstream FMOR_out(" Analyse/Cluster/meanOrthorad.txt",ios::out);
			FMOR_out<<" # epsq fs fd ft"<<endl;
			FMOR_out.close();
			//ofstream FMRdoub_out("Analyse/Cluster/forcemoyenne_double.txt",ios::out);
			//FMRdoub_out<<"# epsq fD fDS"<<endl;
			//FMRdoub_out.close();
			
			ofstream anisotropy_out("Analyse/Cluster/AnisotropyCluster.txt",ios::out);
			anisotropy_out<<"# epsq qoverp ac an at 0.5(ac+an+at)"<<endl;
			anisotropy_out.close();
			
			ofstream anisotropyBytype_out("Analyse/Cluster/AnisotropyClusterBytype.txt",ios::out);
			anisotropyBytype_out<<"# epsq qoverpS qoverpNS acs acns"<<endl;
			anisotropyBytype_out.close();
			
			ofstream coordinence_out("Analyse/Cluster/coordinence.txt",ios::out);
			coordinence_out<<"# epsq Z ratlers" <<endl;
			coordinence_out.close();
			
			ofstream branchLengthC_out("Analyse/Cluster/BranchC.txt",ios::out);
			branchLengthC_out<<"# Branch "<<endl;
			branchLengthC_out.close();
		}
		
	if (calcgap)
	{
		system("mkdir Analyse/Inter_geometry");
				
		ofstream branchLengthD_out("Analyse/Inter_geometry/Branch_D.txt",ios::out);
		branchLengthD_out<<"# id1 id2 vbranch" << endl;
		branchLengthD_out.close();
				
		ofstream gapOwner_out("Analyse/Inter_geometry/particle_gap.txt",ios::out);
		gapOwner_out<<" # gaprmin"<<endl;
		gapOwner_out.close();
	}
	
	if (calfracture)
	{
		ofstream fra("Analyse/fracture.txt",ios::out);
		//fra<<"#	time	Nfissure	Nfragmentation	Ncasse"<<endl;
		fra.close();
	}
			
	if (displaySample)
	{
		system("mkdir Analyse/PS1");
		system("mkdir Analyse/PS2");
	}
}

void brazilian_A::cluster()
{
	cla_ = new cluster_A(sys_);
	cla_->buildListinterdof( prb_ );
	cla_->Stress_Cluster( prb_ );
	cla_->Stress_ClusterSimple(prb_);
	cla_->Stress_ClusterNsimple(prb_);
	//cla_->Fabric_Cluster( prb_); fabrique de branche
	
	cla_->ContactFabric(prb_);
	cla_->simpleContactFabric(prb_);
	cla_->NonsimpleContactFabric(prb_);
	
	cla_->frAnisoInProbe(prb_);
	

	cla_->MeanForce();
	cla_->contact_connectivity();
	//cla_->ContactForceRatio();
	//cla_->PDFFN_Cluster();
	//cla_->PDFFT_Cluster();
	
	ofstream anisotropy("Analyse/Cluster/AnisotropyCluster.txt",ios::app);
	anisotropy<<" "<< epsq <<" "<< cla_->sigma() <<" "<< cla_->ac() <<" "<< cla_->an()<<" "<<cla_->at()<<" "<<0.5*(cla_->ac()+cla_->an()+cla_->at())<<endl;
	anisotropy.close();
	
	ofstream anisotropyBytype("Analyse/Cluster/AnisotropyClusterBytype.txt",ios::app);
	anisotropyBytype<<" "<< epsq <<" " <<cla_->sigmaSimple()<<" " <<cla_->sigmaNsimple()<<" " <<cla_->acs() <<" "<< cla_->acns() <<endl;
	anisotropyBytype.close();
	
	ofstream coordinence("Analyse/Cluster/coordinence.txt",ios::app);
	coordinence<<" "<< epsq <<" " << cla_->ZCluster() <<" " <<cla_->Xratlers()<<endl;
	coordinence.close();
	
	vector<double> Xpkc = cla_->connectivity();
	double Zclust=0;
	ofstream clu("Analyse/Cluster/connectivity.txt",ios::app);//ouverture du fichier avant la boucle
	for(unsigned int i=0;i<Xpkc.size();++i)
	{
		Zclust+= i * Xpkc[i];
		//cout<<i<<" "<<Xpkc[i]<<" "<<Zclust<<endl;
		//// clu<<i<< " "<< Zclust << endl; 

	}
	clu<<" "<<epsq<<" "<<Zclust << endl; 

	clu.close();

//------------------------------------------------------------------------Dénombrement types de contact--------------------------------------	
	
	vector<double> XType_ = cla_->contactTypes();
	ofstream Xcontact("Analyse/Cluster/Xcontact_type.txt",ios::app);
	for (unsigned int i=0 ; i<XType_.size() ; ++i)
	{
		Xcontact<<" "<<epsq<<" "<<XType_[i];
	}
	Xcontact<<endl;
	Xcontact.close();

	/*vector<double> XDouble_ = cla_->contactDouble(); // pour l'instant on ne s'occupe pas de distinguer les différents contacts double
	ofstream Xdouble("Analyse/Cluster/Xdouble_type.txt",ios::app);
	for (unsigned j=0 ; j<XDouble_.size() ; ++j)
	{
		Xdouble<<" "<<epsq<<" "<<XDouble_[j];
	}
	Xdouble<<endl;
	Xdouble.close();*/

	ofstream FMR("Analyse/Cluster/meanRadial.txt",ios::app);
	for (unsigned int j=0 ; j<XType_.size() ; ++ j)
	{
		FMR<<" "<<epsq<<" "<<cla_->MeanRadial(j); 
	}
	FMR<<endl;
	FMR.close();
	
	
	ofstream FMOR("Analyse/Cluster/meanOrthorad.txt",ios::app);
	for (unsigned int p=0 ; p<XType_.size() ; ++ p)
	{
		FMOR<<" "<<epsq<<" "<<cla_->MeanOrtho(p); 
	}
	FMOR<<endl;
	FMOR.close();

	

	
	
	
//---------------------------------------------------------------------------------on s'occupe des contacts double uniquement------------------
	
	/*DataSet meanForce_D;
	DataSet meanForce_DS;
	for (unsigned int j =0 ; j<cla_->storage_a().size() ; ++j)
	{
		meanForce_D.add(cla_->catchInterdof(j)->fresNormal());
	}
	
	for (unsigned int j=0 ; j<cla_->storage_b().size() ; ++j)
	{
		meanForce_DS.add(cla_->catchInterdof(j)->fresNormal());
	}
	
	meanForce_D.extractValues();
	meanForce_DS.extractValues();
	ofstream forcemoye("Analyse/Cluster/forcemoyenne_double.txt",ios::app);
	forcemoye<<" "<<epsq<<" "<<meanForce_D.mean()/meanForceGlobal.mean()<<" "<<meanForce_DS.mean()/meanForceGlobal.mean()<<endl;
	forcemoye.close();*/
	
	// A sauver : Xpkc et Zclust
	
	
	//vector<double> Xk = cla->rankStat();
	//for(unsigned int i=0;i<Xk.size();++i)
	//{
	//	cout<<i<<" "<<Xk[i]<<endl;
	//}
	
}

void brazilian_A::Ptheta(unsigned int Nbin , unsigned int period, unsigned int width) 
{
	cout<<"	Ptheta  : ";
	//cout<<"NbinPT = "<<Nbin<<endl;
	
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
	double amp = M_PI/(double) (Nbin);//amplitude de chaque classe

	unsigned int NctotN=0,rangN,NctotT=0;
	inter2d * interc;
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));

		if( interc->fn() !=0. && prb_.contain(interc))
		{
			angN= acos( interc->nx() ) ;//calcul de l'orientation de la normale au contact
			if(interc->ny()< 0.) angN = M_PI - angN;

			rangN=(unsigned int) ( floor( angN/amp));// choix de la classe
			if (rangN==Nbin) rangN--;

			dx=interc->Vbranchx();
			dy=interc->Vbranchy();
			l=sqrt(dx*dx+dy*dy);

			NcN[ rangN ]++;//remplissage de chaque classe avec son effectif
			FN [ rangN ]+=interc->fn();//résultante des forces normales dans chaque classe
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
	ofstream str("Analyse/angDistrib/Ptheta.txt",ios::out);
	for ( unsigned int i =0;i<Nbin;++i)
	{
		str<<(.5+i)*amp<<"\t"<<(NcN[i])<<endl;
		str<<(.5+i)*amp<<"\t"<<(double)(NcN[i])/(double)(NctotN)*Nbin/M_PI<<endl;
	}
	str.close();
	
	fnmoy/=(double) (NctotN);
	ftmoy/=(double) (NctotT);
	lmoy /=(double) (NctotN);
	//cout<<"Nctot = "<<Nctot<<endl;

	pointSet fn,ft,ld,p;

	for ( unsigned int i =0;i<Nbin;++i)
	{
		p.add( (.5+i)*amp, (double) (NcN[i])/(double) (NctotN)*(double)(Nbin)/M_PI );//Distrib angulaire de prob c'est une densité de fréquence
		fn.add( (.5+i)*amp, FN[i]/(double) (NcN[i]) );//Distrib ang de l'intensite de forces normales
		ld.add( (.5+i)*amp, L[i]/(double) (NcN[i]) );//Distrib ang de l'intensite de forces tangentielle
		if( NctotT!=0) ft.add( (.5+i)*amp, FT[i]/(double) (NcT[i]));
	} 
	//(NcN[i])/(NctotN)*(Nbin)/M_PI )

	if( NctotT!=0)
	{
		pointSet ftm =ft.mobileMean( period,width);
		ftm.yNormalise( fnmoy);
		ftm.write("ANALYSE/angDistrib/FTtheta_mob.txt");
	}   
	pointSet pm  =p.mobileMean(  period,width);
	pointSet fnm =fn.mobileMean( period,width);
	pointSet ldm = ld.mobileMean( period,width);

	fnm.yNormalise( fnmoy);
	//p.write("ANALYSE/angDistrib/Ptheta_mob.txt");
	pm.write("ANALYSE/angDistrib/Ptheta_mob.txt");
	fnm.write("ANALYSE/angDistrib/FNtheta_mob.txt");
	ldm.write("ANALYSE/angDistrib/Ltheta_mob.txt");

	ofstream F_approx("ANALYSE/angDistrib/PF_approx.txt",ios::out);
	for ( unsigned int i =0;i<Nbin;++i)
	{
		F_approx<<(.5+i)*amp<<" "
			<<1./M_PI*( 1.+ a() *cos( 2.*((.5+i)*amp - da_ )))<<" "
			<<  (1. + an()*cos( 2.*((.5+i)*amp - dn_ )))<<" "
			<<  ( -1.*at()*sin( 2.*((.5+i)*amp - dt_ )))<<endl;
			
	}
	F_approx.close();
	
	pointSet  pt = PthetaInProbe( totalProbe_, *(sys_)->spl(),*(sys_)->nwk() ,Nbin);
	pt.write("ANALYSE/angDistrib/Ptheta_mob2.txt");
	
	cout<<" ok "<<endl;
}
/*
void brazilian_A::Ptheta(unsigned int Nbin , unsigned int period, unsigned int width) 
{
	
	cout<<"	Ptheta  : ";
	//cout<<"NbinPT = "<<Nbin<<endl;


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
	double amp = M_PI/(double) (Nbin);//amplitude of each class

	unsigned int NctotN=0,rangN,NctotT=0;
	inter2d * interc;
	
	unsigned int count = 0;
	
	//ofstream Data_out("ANALYSE/ROUGH_DATA.dat",ios::out); // creates output file for rough data
	
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		count++;
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));

		if( interc->fn() !=0. && prb_.contain(interc))
		{
			angN= acos( interc->nx() ) ;//compute contact direction
			if(interc->ny()< 0.) angN = M_PI - angN;
			
			//Data_out<<angN<<endl;// writes rough data into output file

			rangN=(unsigned int) ( floor( angN/amp));// choix de la classe
			if (rangN==Nbin) rangN--;

			dx=interc->Vbranchx();
			dy=interc->Vbranchy();
			l=sqrt(dx*dx+dy*dy);

			NcN[ rangN ]++;//filling each class with number of elements
			FN [ rangN ]+=interc->fn();//normal forces resultant in class
			L  [ rangN ]+=l;//branch length in class
			NctotN++;//counts active forces that has been treated
			
			fnmoy+=interc->fn();
			lmoy +=l;

			if( interc->ft() != 0. )
			{
				FT [ rangN ]+= interc->ft();
				NcT[ rangN ]++;
				ftmoy += interc->ft();
				NctotT++;
			}
		}


	}
	fnmoy/=(double) (NctotN);
	ftmoy/=(double) (NctotT);
	lmoy /=(double) (NctotN);
	//cout<<"Nctot = "<<Nctot<<endl;
	
	cout<<"event number:"<<count<<endl;
	
	//Data_out.close();// close the output file of rough data

	pointSet fn,ft,ld,p;

	for ( unsigned int i =0;i<Nbin;++i)
	{
		p.add( (.5+i)*amp, (double) (NcN[i])/(double) (NctotN)*(double)(Nbin)/M_PI );//Angular distribution of contact normal-->this is frequency histogramm
		fn.add( (.5+i)*amp, FN[i]/(double) (NcN[i]) );//Angular distribution of normal contact forces intensity
		ld.add( (.5+i)*amp, L[i]/(double) (NcN[i]) );//Angular distribution of branch vectors length
		if( NctotT!=0)
			ft.add( (.5+i)*amp, fabs(FT[i]/(double) (NcT[i])));//Angular distribution of tangential contact forces intensity

	}

	if( NctotT!=0)
	{
		//pointSet ftm =ft.mobileMean( period,width); //make a mobile mean on rough data
		//ftm.yNormalise( fnmoy);// normalize by mean normal force
		ft.yNormalise( fnmoy);
		//ftm.write("ANALYSE/angDistrib/FTtheta_mob.txt");
		ft.write("ANALYSE/angDistrib/FTtheta_mob.txt");
	}   
	pointSet pm  =p.mobileMean(  period,width);
	pointSet fnm =fn.mobileMean( period,width);
	pointSet ldm = ld.mobileMean( period,width);

	fnm.yNormalise( fnmoy);
	
	pm.write("ANALYSE/angDistrib/Ptheta_mob.txt");
	fnm.write("ANALYSE/angDistrib/FNtheta_mob.txt");
	ldm.write("ANALYSE/angDistrib/Ltheta_mob.txt");

	
	ofstream F_approx("ANALYSE/angDistrib/PF_approx.txt",ios::out);//values of angular distribution obtained by the harmonic approximation
	F_approx.setf( ios::scientific);
	for ( unsigned int i =0;i<Nbin;++i)
	{
		F_approx<<(.5+i)*amp<<"   "
			<<1./M_PI*( 1.+ a() *cos( 2.*((.5+i)*amp - da_ )))<<"   "<<(.5+i)*amp<<"   "<<  (1. + an()*cos( 2.*((.5+i)*amp - dn_ )))<<"   "<<(.5+i)*amp<<"   "<<( -1.*at()*sin( 2.*((.5+i)*amp - dt_ )))<<endl;
			
	}
	F_approx.close();
	
	pointSet  pt = PthetaInProbe( totalProbe_, *(sys_)->spl(),*(sys_)->nwk() ,Nbin);
	pt.write("ANALYSE/angDistrib/Ptheta_mob2.txt");
	
	cout<<" ok "<<endl;
}
*/
void brazilian_A::Gap()
{

//	cout<<"gap"<<endl;
	unsigned int i,Ngap=0;
	unsigned int Ni = sys_->nwk()->linter().size();
	double d;
	double temp;
	temp=-1000.;
	gapmoy_=0.;
	
	//ofstream gapout("gap.txt",ios::out);
	for ( i=0;i< Ni;++i)
	{
		if( sys_->nwk()->inter(i)->type() != 0) continue;
		d=sys_->nwk()->inter(i)->Overlap();
		//if( sys_->nwk()->inter(i)->type()==2) cout<<d<<endl;
	//cout<<i<<" " <<d<<endl;
		if ( d > 0. )
		{
		//	gapout<<-d/sys_->spl()->rmin()<<endl;
			gapmoy_+=d;//sum gap this is a cumulated gap
			Ngap++;
	
			ngap1_ = sys_->nwk()->inter(i)->first();
			ngap2_ = sys_->nwk()->inter(i)->second();
						
			if ( d > temp)
			{	
				gaprmin = d/sys_->spl()->rmin();
				gaprmax = d/sys_->spl()->rmax();
				ngap1_ = sys_->nwk()->inter(i)->first();
				ngap2_ = sys_->nwk()->inter(i)->second();
				temp = d;
			}
		}
	}
		
	gapmoy_/=(double) (Ngap)*sys_->spl()->rmin();
	
	ofstream monitoring("Analyse/monitoring.txt",ios::app);
	monitoring<<time<<" "<<epsp<<" "<<epsq<<" "<<gaprmin<<" "<<gaprmax<<" "<<gapmoy_<< " " << endl;
	monitoring.close();
	
//gap=-1.*min(gap1,gap2);
	//cout<<"gap ok"<<endl;
		
	cout<< " " <<Ni<<endl;
	double Vx_=0;
	double Vy_=0;
	ofstream branchLengthD("Analyse/Inter_geometry/Branch_D.txt",ios::out);
	ofstream gapOwner("Analyse/Inter_geometry/particle_gap.txt",ios::out);
	for (unsigned int i=0 ; i<Ni ; ++i)
	{
		if( sys_->nwk()->inter(i)->type() != 0) continue;
		{
			Vx_ = sys_->nwk()->inter(i)->Vbranchx();
			Vy_ = sys_->nwk()->inter(i)->Vbranchy();
			branchLengthD<<" "<<sqrt(Vx_*Vx_+Vy_*Vy_)<<endl;
			gapOwner<<" "<<sys_->nwk()->inter(i)->Overlap()/sys_->spl()->rmin()<<endl;
		}
	}
	branchLengthD.close();
	gapOwner.close();
}

/*void brazilian_A::def()
{
	defx+=( right_->x()-x0)/(right_->x();
	defy+=( top_->y()-y0)/top_->y();
	
	y0 = top_ ->y();
	x0 = right_->x();

	epsp=defx+defy;
	epsq=max(defx,defy)-min(defx,defy);
	
	//cout<<" 	defx : "<<defx<<" defy : "<<defy<<" desp : "<<epsp<<" defq : "<<epsq<< endl;
	
	ofstream str("Analyse/strain.txt",ios::app);
	str<<time<<" "<<strain.xy()<<" "<<strain.yy()<<" "<<max(strain.l1(),strain.l2())<<
		" "<<min(strain.l1(),strain.l2())<<" "<<
		epsp<<" "<<epsq<<" "<<asin(epsp/epsq)<<endl;
	str.close();
	//out<<defXY+defYY<<" "<<endl;

}*/

void brazilian_A::def() 
{
	this->sys()->spl()->updateBoundaries(); 		//Ajouter le 25/11/2011 	
	double h= sys_->spl()->ymax()-sys_->spl()->ymin();
	double l = sys_->spl()->xmax()-sys_->spl()->xmin();
		
	defx=(l-l0)/l0;
	defy=(h-h0)/h0;
	
	//h0 = h;
	//l0 = l;
	inter2d * oxo;
	double force2=0.;
	double force1=0.;
	for( unsigned int i=0 ;i<sys_->nwk()->clist().size();++i)
	{
		oxo=sys_->nwk()->inter(sys_->nwk()->clist(i));
		if((oxo->first()->id()==2) || (oxo->second()->id()==2))
		{				
			for(unsigned int r=0 ; r < oxo->rang() ; ++r)
			{
				oxo->current() = r;
				force2+=oxo->fy();
			}
		}
		
		if((oxo->first()->id()==0) || (oxo->second()->id()==0))
		{				
			for(unsigned int r=0 ; r < oxo->rang() ; ++r)
			{
				oxo->current() = r;
				force1+=oxo->fy();
			}
		}
	}

	epsp=defx+defy;
	epsq=max(defx,defy)-min(defx,defy);
	ofstream str("Analyse/strain.txt",ios::app);
	str<<time<<" "<<defx<<" "<<-defy<<" "<<-defy/sys_->grpRel()->fragmentation()->epsilonc()<<" "<<epsp<<" "<<epsq<<" "<<asin(epsp/epsq)<<"	"<<s1<<" "<<s2<<" "<<max(s1,s2)/1.0e+6<<" "<<max(s1,s2)-min(s1,s2)<<" " << qop_ << " " << (s1+s2)/2. << " "<<qop_*(s1+s2)/2.<<"	"<<max(s1,s2)/sys_->grpRel()->fragmentation()->sigmac()<<"	"<<0.5*(force1-force2)/(totalProbe_.R()*M_PI*sys_->grpRel()->fragmentation()->sigmac())<<endl;
	//str<<time<<" "<<defx<<" "<<-defy<<" "<<epsp<<" "<<epsq<<" "<<asin(epsp/epsq)<<"	"<<s1<<" "<<s2<<" "<<max(s1,s2)/1.0e+6<<" "<<max(s1,s2)-min(s1,s2)<<" " << qop_ << " " << (s1+s2)/2. << " "<<qop_*(s1+s2)/2.<<"	"<<endl;

	str.close();
}

void brazilian_A::SF()
{
	//sf()=solidFraction(totalProbe_,*sys_->spl(), *sys_->nwk());
	sf()=solidFraction(prb_,*sys_->spl(), *sys_->nwk());
	
	double sf2 = solidFraction(totalProbe_,*sys_->spl(), *sys_->nwk());
	
	/*double surf=0.;
	for (unsigned int i=4;i<sys_->spl()->lbody().size();++i)
	{
		surf+=sys_->spl()->body(i)->Area();
		
	}
	double R = left_->R();
	//sf()=surf/((right_->x()-left_->x())*(top_->y()-bottom_->y()));	
	sf()=surf/((right_->x()-left_->x()-2*R)*(top_->y()-bottom_->y()-2*R));	*/
	
	cout<<"	SF : "<<sf()<<" in circle "<<sf2<<endl;
}

void brazilian_A::polygAnisotropy( unsigned int nc)
{
	cout<<"	polyg A : "<<endl;
	double as=0,das=0,oas=0;
	double ad=0,dad=0,oad=0;
	
	inter2d* interc;
	unsigned int Nc=0;
	double fnmoy=0;
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));

		if( interc->fn() !=0. && prb_.contain(interc))
		{
			Nc++;
			fnmoy+=interc->fn();
			if ( interc->rang() ==2)
			{
				interc->current()=1;
				fnmoy+=interc->fn();
				interc->current()=0;
			}
			
		/*	if( interc->rang() ==2 && interc->first()->id()==500)
			{
				cout<<" x1 = "<<interc->x()<<" y1 = "<<interc->y()<<endl;
				interc->current()=1;
				cout<<" x2 = "<<interc->x()<<" y2 = "<<interc->y()<<endl;
				interc->current()=0;
			}
		*/
		}
	}
	fnmoy/=(double) (Nc);
//	cout<<" Nombre de contact "<<Nc<<" fnmoy = "<<fnmoy<<endl;
	
	gdm::Tensor2x2 * Fs = SimpleContactFabricInProbe(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	if (Fs != NULL ) 
	{
		Fs->scalarMult(1./(double)(Nc));
		Fs->eigenValues();
		as=2.*(max(Fs->l2(),Fs->l1())- min(Fs->l2(),Fs->l1()) );//(Fs->l1()+Fs->l2());
	    das= Fs->majorDirection();
	    oas=as*cos( 2.*(das- .5*M_PI));
		cout<<" 		direction simple " << das/M_PI*180.<<" acs= "<<as<<" oacs = "<<oas<<endl;
	}
		//F->adresse();
	delete Fs;

	//cout<<" direction " << da_/M_PI*180.<<" ac= "<<a()<<endl;
	
	gdm::Tensor2x2 * Fd = DoubleContactFabricInProbe(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	if (Fd != NULL ) 
	{
		Fs->scalarMult(1./(double)(Nc));
		Fd->eigenValues();
		ad=2.*(max(Fd->l2(),Fd->l1())- min(Fd->l2(),Fd->l1()) );//(Fd->l1()+Fd->l2());
		dad= Fd->majorDirection();
		oad=ad*cos( 2.*(dad- .5*M_PI));
		cout<<" 		direction double " << dad/M_PI*180.<<" acd= "<<ad<<" oacd = "<<oad<<endl;
		cout<<" 		somme a = "<<as+ad<<" oa = "<<oad+oas<<endl;
		
	}
	delete Fd;
	ofstream out("Polygons/anisotropy.txt",ios::app);
	out<<time<<" "<<epsq<<" "<<as<<" "<<ad<<" "<<oas<<" "<<oad<<endl;
	out.close();
	
	ofstream out1("Polygons/dpm.txt",ios::app);
	out1<<time<<" "<<epsq<<" "<<das/M_PI*180.<<" "<<dad/M_PI*180.<<endl;
	out1.close();
	
	//cout<<nc<<endl;
	pointSet  pts = SimpleContactPthetaInProbe( totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),nc );
	pts.write("Polygons/SPtheta.txt");
	
	pointSet  ptd = DoubleContactPthetaInProbe( totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),nc );
	ptd.write("Polygons/DPtheta.txt");
	
}

void brazilian_A::A()
{
	an_=at_=al_=0;
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No Fabric "<<endl;return ;}
	
	cout<<"	Fabric A : ";
	gdm::Tensor2x2 * F = FabricInProbe(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	
	if (F != NULL ) 
	{
		F->eigenValues();
		a()=2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );
		da_=F->majorDirection();
		oa_=a()*cos( 2.*(da_- .5*M_PI));
		cout<<" direction " << da_/M_PI*180.<<" ac= "<<a()<<" oa= "<<oa_<<endl;		
	}
	delete F;

}

void brazilian_A::forces_A(int Nbin )
{
	
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	
	cout<<"	fn_A : ";
	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	
	if (Fn!=NULL)//(Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		an() = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
			//cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dn_=Fn->majorDirection();
		an()=an()-a();
		cout<<"     direction "<<dn_/M_PI*180.<<" an= "<<an()<<endl;
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
			at()=2.*( max(Ft->l2(),Ft->l1()) - min(Ft->l2(),Ft->l1()) )/(Fn->l1() + Fn->l2() );
	//cout<<Ft->l1()<<" "<<Ft->l2()<<endl;
			dt_=Ft->majorDirection();
			at()-=an()+a();
			cout<<"     direction "<<dt_/M_PI*180.<<" at= "<<at()<<endl;
		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		at()=0;
	}
	
	cout<<"	l_A : ";
	gdm::Tensor2x2 * L = lengthAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	if (L != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		L->eigenValues();
		al() = 2.*( max(L->l2(),L->l1()) - min(L->l2(),L->l1()) )/(L->l1() + L->l2() );
			//cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dl_=L->majorDirection();
		al()=al()-a();
		cout<<"      direction "<<dl_/M_PI*180.<<" al= "<<al()<<endl;
	}

	delete Fn;
	delete Ft;
	delete L;
	

}

void brazilian_A::globalStress()
{
	
	cout<<"	Globalstress : ";
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No Stress "<<endl;return ;}
	gdm::Tensor2x2 * S = StressInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk());	
	if (S != NULL ) 
	{
		S->eigenValues();
		//S->print();
		s1=S->l1();
		s2=S->l2();
		qop_=( max(s1,s2)-min(s1,s2) )/(s1+s2);
		pressure_=s1+s2;
		pressure_defined=true;
		ds_=S->majorDirection();
		cout<<" direction "<<ds_/M_PI*180.<<" q/p = "<<qop_;

	}
	ofstream GS_out("Analyse/stress.txt",ios::app);
	GS_out<<time<<" "<<s1<<" "<<s2<<" "<<max(s1,s2)-min(s1,s2)<<" " <<epsq<<" " << qop_ << " " << (s1+s2)/2. << " "<<qop_*(s1+s2)/2.<<" "<<ds_<<endl;
	GS_out.close();
	
	ofstream U_compression_data_app("Analyse/U_compression.dat",ios::app);
	U_compression_data_app<<"   "<<sys()->ctrl(topId).yval()<<"   "<<sf_<<"   "<<sys_->nwk()->clist().size()<<endl;
	U_compression_data_app.close();

	delete S;
	cout<<endl;
}

unsigned int brazilian_A::Z(bool exportZP, unsigned int Nbin, bool cvd)
{
	unsigned int i,ni;
	unsigned int Nb = sys_->spl()->lbody().size();//Number of particles
	unsigned int Nc = sys_->nwk()->clist().size();//Number of contact
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

//Calcul number of contacts with force et with null force, number of contact per particules, force moyence
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		//if( prb_.contain( sys_->nwk()->inter(ni)) && prb_.containCenter( sys_->nwk()->inter(ni)->first()) && prb_.containCenter( sys_->nwk()->inter(ni)->second()))
		//{
			id1=sys_->nwk()->inter(ni)->first()->id();
			id2=sys_->nwk()->inter(ni)->second()->id();

			Activate[id1]=true;
			Activate[id2]=true;

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
				//cout<<sys_->nwk()->inter(ni)->type()<<" ";
				//cout<<sys_->nwk()->inter(ni)->first()->id()<<" "<<sys_->nwk()->inter(ni)->second()->id()<<endl;
			}
	//	}

	}
		//double Pf0= (double) (Ncf0)/(double) (Ncf0 + Ncwf);
	if (Ncwf==0) 
		{ cout<<" no forces "<<endl;return 0;}
	fmoy/=(double) (Ncwf);
	unsigned int Ncw=0;
	cout<<" fmoy = "<<fmoy<< " ";
	//cout<<" Nforce nulle = "<<Ncf0<<" Ncwf = "<<Ncwf<<endl;

		
	ofstream scalZ("Analyse/Connect/scalarZ.txt",ios::app);
	scalZ<<time<<" "<<Ncwf<<" "<<Ncf0<<endl;

	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		if( prb_.contain( sys_->nwk()->inter(ni)))
		{
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

		}

	}
	cout<<" Fraction of weak contact : "<<(double) (Ncw)/(double)(Ncwf);
	pointSet fracw;
	for(i=0;i<Nb;++i)
	{
		if( Ncpp[i]>1)
			fracw.add( sys_->spl()->body(i)->sizeVerlet(), (double) (Nw[i])/(double)(Ncpp[i]));
	}
	
	fracw.decreasingXSort();
	//pointSet weak= fracw.histoBinVar(200);
	pointSet weak= fracw.slidingHisto(50,.10);
	//pointSet mean= weak.mobileMean(1,2);
	weak.xoffSet(-rmin);
	weak.xNormalise( rmax-rmin);
	weak.write("Analyse/connect/fracweak.txt");
	

	Nz=0;
	NcppMax=0;
	NcZ=0;
	double r;
	pointSet gcu,gc,nu,vu;
	double masstu=0,masst=0;

	rattlersB.clear();
	rattlersVisu.clear();
	rattlersVisu.resize(Nb);

	for(i=4;i<Nb;++i)
	{ 
		if (prb_.containCenter(sys_->spl()->body(i) ) )
		{
			r=sys_->spl()->body(i)->sizeVerlet();

			gc.add(r , r*r);
			masst +=r*r;

			if (Ncpp[i] > 1 && Activate[i]==true   )  // occurs in contact list and probe and with contact number per particles larger than 1
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
			else if( Ncpp[i]<=1 || Activate[i]==false)
			{
				++Nf;
				nu.add(r, 0.);
				rattlersB.push_back(sys_->spl()->body(i) );
				rattlersVisu[i]=true;
				//cout<<" nf i= "<<i<<endl;
			}

		}
	}
	unsigned int Ncrl =0;
	for (unsigned int k =0 ; k<sys_->nwk()->clist().size() ; ++k)
	{
		if(sys_->nwk()->inter(sys_->nwk()->clist(k))->type() == _type_dkrl && sys_->nwk()->inter(sys_->nwk()->clist(k))->fn() != 0)
		{Ncrl += 1;}
	}
	NcZ = NcZ + Ncrl;//Tong so contact cua tat ca cac hat co nhieu hon 1 contact
	ratlers_=(double)(Nf)/(double)(Nf+Nz);//Ty le so phan tu co it hon hoac bang 1 contact
	z()=(double) (NcZ)/(double) (Nz);//So contact trung binh cua moi hat (chi tinh voi cac hat co nhieu hon 1 contact)
	cout<<" z= "<<z()<<" ratlers_= "<<ratlers_<<flush;

	if (cvd )
	{

		//cout<<" mt "<<masst<<" mtu "<<masstu<<endl;
		//gcu.increasingXSort();
		//gc.increasingXSort();

		pointSet histvolutil = gcu.histoNumberBin(Nbin);
		pointSet histvoltot  = gc .histoNumberBin(Nbin);


		pointSet temp;

		for( unsigned int i=0; i< Nbin; ++i)
		{
			//cout<<histvoltot.py(i)<<" "<<histvolutil.py(i)<<endl;
			temp.add( histvoltot.px(i), histvolutil.py(i)/histvoltot.py(i));
		}
		temp.write("Analyse/connect/fvgu.txt");

		gcu.yNormalise( masstu);
		gcu.yCumul();
		gcu.write("Analyse/connect/granutilZ.txt");

		gc.yNormalise( masst);
		gc.yCumul();
		gc.write("Analyse/connect/grantot.txt");

		pointSet histonu = nu.histoBin(Nbin);
		histonu.write("Analyse/connect/fgu.txt");//fraction granulo utile
	}

	//cout<<" pf01 = "<<Pf0<<" pf02 "<<(double) (Ncf0)/(double) (Ncf0 + Ncwf- Nf)<<endl;
	//cout<<" Nz = "<<Nz<<" Ncz = "<<NcZ<<" Nf = "<<Nf<<endl;

	

	//Connectivity as a function of granulometry
	//So luong contact cho tung gia tri cua ban kinh
	if( calczg)
	{
		cout<<" / Zg ";
		vector < DataSet*> DZ(Nbin,NULL);
		vector < DataSet*> R(Nbin,NULL); 

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

		ofstream data("Analyse/connect/Zid.txt",ios::out);
		double rmax=sys_->spl()->rmax();
		double rmin=sys_->spl()->rmin();
		unsigned int jc;
		double amp=(rmax-rmin )/(double) (Nbin);
		pointSet PDZ;

		for(i=4;i<Nb;++i)
		{ 
			if ( Activate[i]==true && Ncpp[i]>1)
			{
				r=sys_->spl()->body(i)->sizeVerlet();
					//cout<<i<<" "<<r<<endl;
				data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<Ncpp[i] <<endl;

				ZD[Ncpp[i]]->add( r );

				if(prb_.containCenter(sys_->spl()->body(i)))
				{
						//r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
					jc = (unsigned int) floor ( (r-rmin)/amp);

					if ( fabs(r - rmax)< 1e-10 )  jc=Nbin-1;
					if ( fabs(r - rmin)< 1e-10 )  jc=0;

					DZ[jc]->add( Ncpp[i]);
					R [jc]->add( r );
					PDZ.add( r,Ncpp[i]);
				}
			}
		}
		data.close();
		cout<<" ok "<<flush<<endl;


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
		sli.write("Analyse/connect/sliDZ.txt");

		ofstream output ( "Analyse/connect/DZ2.txt" , ios::out );
		if ( ! output ) cout<<"erreur creation de DZ2"<<endl;


		DataSet dPZ,dR;
		unsigned int N=0,Nev=0,Nepc;

		Nepc=PDZ.psetSize()/Nbin;
		cout<<"Nepc = "<<Nepc<<" "<<flush;

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

		ofstream data1("Analyse/connect/ZD.txt",ios::out);
		for( i=0;i<NcppMax;++i)
		{
			ZD[i]->extractValues();
			if( ZD[i]->mean() != 0. )
				data1<<i<<" "<<ZD[i]->mean()<<" "<<ZD[i]->variance()<<endl;
		}
		data1.close();

		ofstream data2("Analyse/connect/DZ.txt",ios::out);
		for(i=0;i<Nbin;++i)
		{
			DZ[i]->extractValues();
			R[i]->extractValues();
			data2<<(R[i]->mean()-rmin)/(rmax-rmin)<<" "<<R[i]->mean()<<" "<<DZ[i]->mean()<<" "<<R[i]->variance()<<" "<<DZ[i]->variance()<<endl;

		}
		data2.close();


	}

	//Partial connectivity
	unsigned int N=0;
	if ( exportZP)
	{
		cout<<"  / ZP ";
		vector <unsigned int> Npnc(NcppMax+1,0);//Number of particles with n contacts
		for(i=4;i<Nb;++i)
		{ 
			if ( Activate[i]==true)
			{
				Npnc[Ncpp[i]]+=1;
				N++;
			}
		}
		//cout<<" N = "<<N<<endl;
		//A prevoir: la numerotation successive des fichier zp
		ofstream out("Analyse/connect/ZPt.txt",ios::app);
		ofstream out1("Analyse/connect/ZP.txt",ios::out);
		out<<time<<" ";

		vector <double> Zp(NcppMax+1);
		for(i=0;i<NcppMax;++i)
		{ 
			Zp[i]=(double) (Npnc[i])/(double) (N);

			out<<Zp[i]<<" ";
			out1<<i<<" "<<Zp[i]<<endl;

		}
		out<<endl;
		out.close();
		out1.close();

		Npnc.clear();
		Zp.clear();
	}
	Ncpp.clear();
	cout<<endl;
	return 1;

}

void brazilian_A::granuloSpeed(unsigned int Nbin)
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
		if(prb_.containCenter(sys_->spl()->body(i)))
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
	double vdof=0.;//sqrt( pow( partref->vx(),2) + pow( partref->vy(),2));
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

unsigned int brazilian_A::granuloStress(unsigned int Nbinstress )
{
	cout<<"	Granulostress : ";
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr );

	for( unsigned int i=0; i < Nb;++i)
	{
		//rajouter une multi par le volume specifique =1/compacité
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
		if(prb_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
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


	cout<<" OK"<<endl;
	return 1;
}

unsigned int brazilian_A::granuloStress2(unsigned int Nc ,unsigned int period, unsigned int width)
{
	cout<<"	Granulostress2 : "<<flush;
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr );
	cout<<" imt ok"<<endl;
	for( unsigned int i=0; i < Nb;++i)
	{
	//rajouter une multi par le volume specifique =1/compacité
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
		if(prb_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
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
	cout<<"Nepc = "<<Nepc<<" ";

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
//	gc.write("gc.txt");
	/*ofstream output ( "stressgranulomoy.txt" , ios::out | ios::trunc);
	if ( ! output ) cout<<"erreur creation de stressgranulomoy"<<endl;

	else
	{
		output.setf( ios::scientific);
		for (unsigned int i=0;i<Nbinstress;++i)
		{

		}
	}
	output.close();

	for( unsigned int i=0;i<Nbinstress;++i)
	{

	}
*/

	cout<<" OK"<<endl;
	return 1;
}

unsigned int brazilian_A::granuloStress3(unsigned int Nbinstress )
{
	system("mkdir granuloS");
	cout<<"	Granulostress : ";
	unsigned int Nb = sys_->spl()->lbody().size();
	bool noIMT=true;
	vector <gdm::Tensor2x2> bodstr(Nb);

	IMT_Body( *sys_->spl(), *sys_->nwk(),bodstr ,sf_);

	for( unsigned int i=0; i < Nb;++i)
	{
		//rajouter une multi par le volume specifique =1/compacité
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
	ofstream data("granuloS/qpid.txt",ios::out);
	for( unsigned int i=0; i < Nb;++i)
	{
		if(prb_.containCenter(sys_->spl()->body(i)) && bodstr[i].l1() != 0. )
		{
			r= sys_->spl()->body(i)->ymax() - sys_->spl()->body(i)->y();
			//cout<<" r = "<<(r-rmin)/amp<<endl;
			//cout<<jc<<endl;
			//cout<<" r-rmax = "<<endl;

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
				P.add(r, -.5 * ptemp *phase);
				Q.add(r, .5 * qtemp *phase);
				QP.add(r, qtemp/ptemp * phase);
				DPM.add(r, dpm );
				rglob.add(r);
			}
			else cout<<" qtemp nul"<<endl;
		}
	}
	data.close();


	P  .increasingXSort();
	pointSet ps= P.slidingHisto(Nbinstress,.05);
	ps.xoffSet(-rmin);
	ps.xNormalise( rmax-rmin);
	ps.write("granuloS/PG_adim.txt");
	
	Q .decreasingXSort();
	pointSet qs= Q.slidingHisto(Nbinstress,.05);
	qs.xoffSet(-rmin);
	qs.xNormalise( rmax-rmin);
	qs.write("granuloS/QG_adim.txt");
	
	QP .decreasingXSort();
	pointSet pqs= QP.slidingHisto(Nbinstress,.05);
	pqs.xoffSet(-rmin);
	pqs.xNormalise( rmax-rmin);
	pqs.write("granuloS/PQG_adim.txt");
	
	DPM .decreasingXSort();
	pointSet dpms= DPM.slidingHisto(Nbinstress,.05);
	dpms.xoffSet(-rmin);
	dpms.xNormalise( rmax-rmin);
	dpms.write("granuloS/DPMG_adim.txt");
		

	double mass=0;
	DataSet granulo;
	rglob.Sort();
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		mass+=rglob.set(i)*rglob.set(i);
	}
	ofstream gran("granuloS/granutilS.txt",ios::out);
	for( unsigned int i=0;i<rglob.setSize();++i)
	{
		granulo.add( rglob.set(i)*rglob.set(i)/mass);
		gran<<rglob.set(i)<<" "<<granulo.set(i)<<endl;
	}
	gran.close();

	cout<<" OK"<<endl;
	return 1;
}

int brazilian_A::pdfforce(bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	char fichier[100];
	char fichierbrut[100];
	cout<<"	PDF  : " ;
	
	if( fn )
	{
		sprintf(fichier,"Analyse/pdf/pdffn.txt");
		sprintf(fichierbrut , "Analyse/pdf/fn_brut.txt");
	}
	else
	{
		sprintf(fichier,"Analyse/pdf/pdfft.txt");
		sprintf(fichierbrut , "Analyse/pdf/ft_brut.txt");
	}
	
	DataSet f;
	double force=0.;
	inter2d * oxo;

	for( unsigned int i=0;i< sys_->nwk()->linter().size();++i)
	{
		oxo=sys_->nwk()->inter(i);
		if( prb.contain(oxo) && (oxo->rang()>0) && sys_->grpRel()->fragmentation()->fco(oxo)!=0.)
		{
			if(fn)
			{
				for (unsigned int r=0;r<oxo->rang();++r)
				{
					oxo->current()=r;
					if((sys_->grpRel()->fragmentation()->res_traction(oxo)!=0.) && (oxo->fn()!=0.))
					f.add(oxo->fn()/sys_->grpRel()->fragmentation()->res_traction(oxo));
				}
				oxo->current()=0;
			}
			else
			{
				for (unsigned int r=0;r<oxo->rang();r++)
				{
					oxo->current()=r;
					force=oxo->ft()/(sys_->grpRel()->fragmentation()->res_cisaillement(oxo)+sys_->grpRel()->fragmentation()->frict()*oxo->fn());
				}
				oxo->current()=0;

				if(force!= 0.) 
				{
					f.add(fabs(force));//stockage dans f des forces de contacts à traiter
				}
			}
		}
	}
	if (f.setSize() <100)
		{cout<<" no more forces (min 100 events ) : "<<f.setSize()<<endl; return 0;}

	f.extractValues();
		
	f.DecreasingSort();
	
	if(fn)
		sprintf(fichier,"Analyse/pdf/pdffn.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft.txt");
	
	f.write(fichier);
	
	cout<<"	OK "<<endl;
	return 1;
}

/*int brazilian_A::pdfforce(bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	char * fichier,*fichierbrut;
	cout<<"	PDF  : " ;
	
	if( fn )
	{
		fichier="Analyse/pdf/pdffn.txt";
		fichierbrut = "Analyse/pdf/fn_brut.txt";
	}
	else
	{
		fichier="Analyse/pdf/pdfft.txt";
		fichierbrut = "Analyse/pdf/ft_brut.txt";
	}
	
	DataSet f;
	double force;
	for( unsigned int i=0;i< sys_->nwk()->linter().size();++i)
	{
		if( prb.contain( sys_->nwk()->inter(i) ))
		{
			force=0.;
			if( fn )
			{
				for (unsigned int r=0;r<sys_->nwk()->inter(i)->rang();r++)
				{
					sys_->nwk()->inter(i)->current()=r;
					force+=sys_->nwk()->inter(i)->fn();
					if(fabs(force)>1.0e+20) {cout<<"force:="<<force<<"	i:="<<i<<"	fn:="<<sys_->nwk()->inter(i)->fn()<<"	first:="<<sys_->nwk()->inter(i)->first()->id()<<"	second:="<<sys_->nwk()->inter(i)->second()->id()<<"	rang:="<<sys_->nwk()->inter(i)->rang()<<"	current:="<<sys_->nwk()->inter(i)->current()<<endl;exit(0); }
				}
				sys_->nwk()->inter(i)->current()=0;

				if(force!= 0.) 
				{
					f.add(force);//stockage dans f des forces de contacts à traiter
				}			
			}
			else
			{
				for (unsigned int r=0;r<sys_->nwk()->inter(i)->rang();r++)
				{
					sys_->nwk()->inter(i)->current()=r;
					force+=sys_->nwk()->inter(i)->ft();
				}
				sys_->nwk()->inter(i)->current()=0;

				if(force!= 0.) 
				{
					f.add(fabs(force));//stockage dans f des forces de contacts à traiter
				}
			}
		}
	}
	if (f.setSize() <100)
		{cout<<" no more forces (min 100 events ) : "<<f.setSize()<<endl; return 0;}

	f.extractValues();

	if( normalize )//considération des forces normalisées 
		f.Normalize( f.mean() );
	cout<<" fmoy = "<<f.mean()<<" ";
	
	if( fn )
		fichier="Analyse/pdf/pdffnrich.txt";
	else
		fichier="Analyse/pdf/pdfftrich.txt";
		
	pointSet pdf0 = f.Rich_PDF2( f.setSize()/20 );
	pdf0.write(fichier);
	
	f.DecreasingSort();
	pointSet pdf = f.kernelPdf( nbin,.01);
	
	if(fn)
		fichier="Analyse/pdf/pdffn.txt";
	else
		fichier="Analyse/pdf/pdfft.txt";
	
	pdf.write(fichier);

	pointSet pdfs = (f.slidingPdf( nbin) ).mobileMean(period,width);

	if( fn )
		fichier="Analyse/pdf/pdffnslid.txt";
	else
		fichier="Analyse/pdf/pdfftslid.txt";
	pdfs.write(fichier);


	if( fn )
		fichier="Analyse/pdf/pdffn_mob.txt";
	else
		fichier="Analyse/pdf/pdfft_mob.txt";
	pdf.periodic()=false;
	pointSet pdfm = pdf.mobileMean2(period, width);
	pdfm.write(fichier);

	cout<<"	OK "<<endl;
	return 1;
}*/

int brazilian_A::granulopdf(unsigned int Nq,bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	//cout<<" debut pdf"<<endl;
	
	cout<<"	granuloPDF  : WARNING --- pdffn must be well parametrized and called "<<endl ;
	//determination des dx
	//number of mass quartiles Nq
	double ampdx=1./(double) (Nq);
	unsigned int N=1; 
	pointSet df,dm;
	double ftemp;
	double r;

	for( unsigned int i=0; i < sys_->spl()->lbody().size();++i)
	{
		if(prb_.containCenter(sys_->spl()->body(i))  )
		{	
			r= sys_->spl()->body(i)->sizeVerlet();
			dm.add( r, r*r);
		}
	}

	dm.increasingXSort();
	dm.yCumul();
	dm.yNormalise(dm.ysum());

	for( unsigned int i=0; i < sys_->nwk()->linter().size();++i)
	{
		if( fn ) ftemp=sys_->nwk()->inter(i)->fn();
		else     ftemp=sys_->nwk()->inter(i)->ft();

		if(ftemp!=0.)
		{
			if( prb_.containCenter(sys_->nwk()->inter(i)->first()) )
				df.add( sys_->nwk()->inter(i)->first()->sizeVerlet(), ftemp);

			if( prb_.containCenter(sys_->nwk()->inter(i)->second()) )
				df.add( sys_->nwk()->inter(i)->second()->sizeVerlet(), ftemp);
		}		
	}


	df.increasingXSort();

	DataSet f;
	unsigned int i=0,j=0;
	double dlim=0;

	char dir[50],fichier[50];
	sprintf(dir,"./granuloPDF");		
	mkdir (dir, S_IRUSR | S_IWUSR | S_IXUSR |
		S_IRGRP |           S_IXGRP |
		S_IROTH |           S_IXOTH  );

	while ( N <= Nq)
	{
		while ( dm.py(i)< N*ampdx)
		{
			++i;
		}
		dlim=dm.px(i);


		f.setClear(); 
		while( df.px(j)<= dlim && j<df.psetSize())
		{
			f.add( df.py(j));
			++j;
		}
		//cout<<" N= "<<N<<" dlim = "<<dlim<<" dm.py(i) = "<<dm.py(i)<<" Nev= "<<f.setSize()<<endl;

		if (f.setSize() ==0)
			{cout<<" no forces "<<endl; return 0;}
		f.extractValues();

		if( normalize ) 
			f.Normalize( f.mean() );


		if( fn )
			sprintf(fichier,"%s/pdffn%d.txt",dir,N);

		else
			sprintf(fichier,"%s/pdfft%d.txt",dir,N);

		pointSet pdf = f.Rich_PDF(nbin);
		pdf.write(fichier);
		
		if( fn )
			sprintf(fichier,"%s/pdffn_mob%d.txt",dir,N);

		else
			sprintf(fichier,"%s/pdfft_mob%d.txt",dir,N);

		pointSet pdfm = pdf.mobileMean(period, width);
		pdfm.write(fichier);

		++N;
	}
	cout<<"	OK "<<endl;	

	return 1;
}

void brazilian_A::forcesMaxCorrelation( )
{
	if ( sys_->nwk()->clist().empty()) return ;
	system("mkdir forceCorr");
	cout<<"	forcesCorrelations  : "<<endl ;
	DataSet f1,d1,d2;
	double fn,r1,r2;
	unsigned int id1,id2,Nbin=100;
	unsigned int Nb= sys_->spl()->lbody().size();

	sys_->spl()->radiusExtrema(4);
	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();

	vector<double> fnmax(Nb,0.);
	vector<double> fnmin(Nb,10000.);
	vector<double> ftot(Nb,0.);//Force totale apps a chaque particules
	vector<double> ftotS(Nb,0.);
	vector<double> ftotW(Nb,0.);
	vector< unsigned int > NC(Nb,0);
	vector< bool > activate(Nb,false);

	pointSet fsum,fdiff,fdtot,fmaxd;;
	inter2d * inter;

	for(unsigned int i=0;i< sys_->nwk()->clist().size();++i)
	{
		inter=sys_->nwk()->inter(sys_->nwk()->clist(i));
		if( inter->fn() != 0.) 
		{
			id1=inter->first()->id();
			id2=inter->second()->id();
			r1= inter->first()->sizeVerlet();
			r2= inter->second()->sizeVerlet();
			fn= inter->fn();

			fnmax[id1]= max( fnmax[id1],fn);
			fnmax[id2]= max( fnmax[id2],fn);
			fnmin[id1]= min( fnmin[id1],fn);
			fnmin[id2]= min( fnmin[id2],fn);
			ftot[id1]+= fn;
			ftot[id2]+= fn;
			NC[id1]+=1;
			NC[id2]+=1;
			activate[id1]=true;
			activate[id2]=true;

			fsum .add( r1+r2 , fn);
			fdiff.add( max(r1,r2)-min(r1,r2),fn);

		}
	}

	for (unsigned int i=0;i<Nb;++i)
	{
		if( activate[i] && ftot[i]!=0.)
		{
			fdtot.add( sys_->spl()->body(i)->sizeVerlet(), ftot[i]/(double) (NC[i]) );//Force totale en fonction de la taille
		}
	}

	fdtot.extractValues();
	sumforce_defined=true;

	if( sumforce_defined)
	{
		FSGTM_.clear(); 
		for (unsigned int i=0;i<Nb;++i)
		{
			if( ftot[i]> fdtot.ymean())
				FSGTM_.push_back(true);
			else 
				FSGTM_.push_back(false);
		}
	}

	fsum. decreasingXSort();
	fdiff.decreasingXSort();
	fdtot.decreasingXSort();

	fsum.extractValues();

	fnmoy_=fsum.ymean();
	cout<<" fnmoy= "<<fnmoy_<<endl;

//	pointSet diff= fdiff.histoBinVar(40);
	pointSet diff= fdiff.slidingHisto(50,.05);
	diff.write("forceCorr/lddiff.txt");

//	pointSet sum = fsum.histoBinVar(200,.05);
	pointSet sum = fsum.slidingHisto(50,.05);
	sum.yNormalise( fsum.ymean() );
	sum.xNormalise( fdtot.xmean());
	sum.write("forceCorr/ldsumm.txt");


	//fdtot.yNormalise( fdtot.ymean());
//  pointSet htot= fdtot.histoBinVar(Nbin,.05);
	pointSet htot= fdtot.slidingHisto(Nbin,.05);
//	pointSet  mean= htot.mobileMean(1,3);
	htot.write("forceCorr/ftotdSW.txt");

//	htot.yNormalise( fdtot.ymean());
//	htot.write("forceCorr/ftotdSW_yNorm.txt");

	htot.xoffSet(-rmin);
	htot.xNormalise( rmax-rmin);

	htot.write("forceCorr/ftotdSW_xyadim.txt");

	//Somme des forces faibles et fortes
	vector< unsigned int > NCS(Nb,0);
	vector< unsigned int > NCW(Nb,0);

	for( unsigned int i=0;i< sys_->nwk()->clist().size();++i)
	{
		inter=sys_->nwk()->inter(sys_->nwk()->clist(i));
		if( inter->fn() != 0.) 
		{
			id1=inter->first()->id();
			id2=inter->second()->id();
			fn= inter->fn();

			if( fn >= fnmoy_)
			{
				ftotS[id1]+= fn;
				ftotS[id2]+= fn;
				NCS[id1]+=1;
				NCS[id2]+=1;
			}
			else
			{
				ftotW[id1]+= fn;
				ftotW[id2]+= fn;
				NCW[id1]+=1;
				NCW[id2]+=1;
			}
		}
	}
	
	ofstream out("inert.txt",ios::out);
	for (unsigned int i=0;i<Nb;++i)
	{
		double v= sqrt(sys_->spl()->body(i)->vx()*sys_->spl()->body(i)->vx()+
			sys_->spl()->body(i)->vy()*sys_->spl()->body(i)->vy());
		double f= sqrt(sys_->spl()->body(i)->fx()*sys_->spl()->body(i)->fx()+
			sys_->spl()->body(i)->fy()*sys_->spl()->body(i)->fy());
		out<<i<<" "<<sys_->spl()->body(i)->sizeVerlet()<<" "<<sys_->spl()->body(i)->mass()*v<<" "<<f<<endl;
	}
	
	pointSet fdtotS,fdtotW,fs,fw;
	for (unsigned int i=0;i<Nb;++i)
	{
		if( activate[i] )
		{
			fdtotS.add( sys_->spl()->body(i)->sizeVerlet(), ftotS[i]/(double) (NC[i]) );//Force totale en fonction de la taille
			fdtotW.add( sys_->spl()->body(i)->sizeVerlet(), ftotW[i]/(double) (NC[i]) );//Force totale en fonction de la taille
			
			if( NCS[i]!=0) fs.add( sys_->spl()->body(i)->sizeVerlet(), ftotS[i]/(double) (NCS[i]));//Force totale en fonction de la taille
			if( NCW[i]!=0) fw.add( sys_->spl()->body(i)->sizeVerlet(), ftotW[i]/(double) (NCW[i]));//Force totale en fonction de la taille
			//cout<<NCS[i]<<" "<<NCW[i]<<" "<<NC[i]<<endl;
		}
		//cout<<ftot[i]<<" "<<ftotS[i]+ftotW[i]<<endl;
	}

	fdtotS.extractValues();
	fdtotW.extractValues();
	
	fs.extractValues();
	fw.extractValues();

	fdtotS.decreasingXSort();
	fdtotW.decreasingXSort();

//	pointSet S= fdtotS.histoBinVar(Nbin,.05);
	pointSet S= fdtotS.slidingHisto(Nbin,.05);
	//pointSet Smean=S.mobileMean(1,3);
	//S->yNormalise( fdtotS.ymean());
	//S->write("ftotdS.txt");
	S.write("forceCorr/ftotdS.txt");

//	pointSet W= fdtotW.histoBinVar(Nbin,.05);
	pointSet W= fdtotW.slidingHisto(Nbin,.05);
	//pointSet Wmean=W.mobileMean(1,3);
	//W->yNormalise( fdtotW.ymean());
	//W->write("ftotdW.txt");
	W.write("forceCorr/ftotdW.txt");

	S.yNormalise( fdtot.ymean());
	S.xoffSet(-rmin);
	S.xNormalise( rmax-rmin);
	S.write("forceCorr/ftotdS_xyadim.txt");

	W.yNormalise( fdtot.ymean());
	W.xoffSet(-rmin);
	W.xNormalise( rmax-rmin);
	W.write("forceCorr/ftotdW_xyadim.txt");
	
	pointSet st=fs.slidingHisto(Nbin,.05);
	pointSet we=fw.slidingHisto(Nbin,.05);
	
	st.yNormalise( fs.ymean());
	st.xoffSet(-rmin);
	st.xNormalise( rmax-rmin);
	st.write("forceCorr/fs_adim.txt");
	
	we.yNormalise( fw.ymean());
	we.xoffSet(-rmin);
	we.xNormalise( rmax-rmin);
	we.write("forceCorr/fw_adim.txt");

	//double f1d1=f1.linearCorrelation(& d1);
	//double f1d2=f1.linearCorrelation(& d2);


		//pointSet fmind(0);
/*	ofstream out("FC2.txt", ios::out);
	for( unsigned int i=0;i< Nb;++i)
		{	//if (i<50) cout<<i<<" "<<fnmax[i]<<" "<<fnmin[i]<<endl;
	fmaxd.add(0.,0.);
	if( fnmax[i]!=0.)
	{

		fmaxd.add( sys_->spl()->body(i)->sizeVerlet(),fnmax[i] );
		out<<sys_->spl()->body(i)->sizeVerlet()<<" "<<fnmax[i]<<" "<<fnmin[i]<<endl;

	}
			
	if( fnmin[i]!=10000.)
	{
			//	cout<<"min non nul"<<endl;
			//fmind.add( sys_->spl()->body(i)->sizeVerlet(),fnmin[i] );
	}
			
}
out.close();*/
pointSet hist= fmaxd.histoBin( 20);
		//fmind.histoBin( 20);
hist.write("forceCorr/FCmean.txt");

return ;
	//cout<<" r1= "<<f1d1<<" r2= "<<f1d2<<endl;
}


/*void brazilian_A::rfd( unsigned int Npoint, unsigned int Nrmean)
{
	//double rmax = sys_->spl()->rmax();
	double rmean = sys_->spl()->rmoy();
	double P  = sys_->spl()->boundWidth();
	double ll= sys_->spl()->xmin() ;
	double rl= sys_->spl()->xmax() ;
	double band = Nrmean*rmean*1.1;

	unsigned int Nbod=0;
	cout<<"rmean " <<rmean<<" P/rmean = "<<.5* P/rmean<<endl;
	cout<<" ll "<<ll<<" "<<rl<<" "<<endl;
	heightProbe lim(prb_.h1() + band, prb_.h2() - band ,P);

	vector <body2d*> inlimit;
	vector <body2d*> inlimitP;


	for (unsigned int i=0; i<sys_->spl()->lbody().size();++i)
	{
		if ( lim.containCenter( sys_->spl()->body(i) ) )
		{ 
			inlimit.push_back(sys_->spl()->body(i));
			++Nbod;
		}
	}

	for (unsigned int i=0; i<sys_->spl()->lbody().size();++i)
	{

		if( sys_->spl()->body(i)->x() < ll + band )
		{
			inlimitP.push_back(sys_->spl()->body(i)->duplicate());
			inlimitP.back()->x()+=P;
			//cout<<" per G ";
		}

		if( sys_->spl()->body(i)->x() > rl - band)
		{
			inlimitP.push_back(sys_->spl()->body(i)->duplicate());
			inlimitP.back()->x()-=P;
			//cout<<" per D ";
		}
			//cout<<endl;
	}
	//Nbod=inlimit.size();
	double dx,dy,d;
	unsigned int rang;

	double amp=Nrmean*rmean/(double) Npoint;
	//double amp=.1*rmean;
	vector<unsigned int > N(Npoint,0); 
	cout<<" NBOD = "<<Nbod<<endl;

	for( unsigned int i=0 ; i< Nbod ;++i)
	{
		double xi= inlimit[i]->x();
		double yi= inlimit[i]->y();
		unsigned int idi = inlimit[i]->id();
		//cout<<idi<<endl;
		for( unsigned int j = 0 ; j< sys_->spl()->lbody().size() ;++j)
		{	
			if( idi==j) 
			{
				continue;
			}
			dx= xi - sys_->spl()->body(j)->x();
			dy= yi - sys_->spl()->body(j)->y();
			d=sqrt(dx*dx+dy*dy);

			rang=(unsigned int) floor( d/amp );

			//cout<<"j "<<j<<" d "<<d<<" rg "<<rang<<endl;
			if( rang > Npoint)
			{
				//cout<<"rang " <<rang<<endl;
				continue;
			}	
			else
			{

				N[rang]++;
			}
		}
	//	cout<<"fin norm"<<endl;

		for( unsigned int j = 0 ; j< inlimitP.size();++j)
		{
			if( idi == inlimitP[j]->id()) continue;
			//cout<<inlimitP[j]->id()<<endl;
			dx=xi - inlimitP[j]->x();
			dy=yi - inlimitP[j]->y();
			d=sqrt(dx*dx+dy*dy);

			rang=(unsigned int) floor( d/amp );
			if( rang > Npoint)
			{
				//cout<<"rang " <<rang<<endl;
				continue;
			}	
			else
			{
				N[rang]++;
			}
		}
	}

	double meandensity= (double) (Nbod)/lim.area();
	vector<double> density(Npoint);
	double r;
	density[0]= N[0]/Nbod/( M_PI*amp*amp);
	for( unsigned int i=1;i<Npoint;++i)
	{	
		r=(i+1.)*amp;
		density[i] = (double) (N[i])/(double)( Nbod);
		density[i]/=( M_PI*(2.*r*amp+ amp*amp ));

	//	cout<<r/rmean<<" "<<N[i]<<" "<<density[i]<<" "<<meandensity<<endl;
	}

	ofstream rfd("rfd.txt",ios::out);
	for( unsigned int i=0;i<Npoint;++i)
	{	
		r=(i+1)*amp;
		rfd<<.5*r/rmean<<" "<<density[i]/meandensity<<endl;
	}


}*/

void brazilian_A::removeBody(body2d * rm)
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

void brazilian_A::removeRattlers()
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

void brazilian_A::reduceRattlers()
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

	void brazilian_A::growRattlers( double dr)
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
		sys_->spl()->radiusExtrema(4);
	}

	void brazilian_A::writePS( const char * fname)
	{
				
		ofstream ps(fname);
		/*
		this->sys()->spl()->updateBoundaries(); 		//Ajouter le 25/11/2011 	
		double height= sys_->spl()->ymax()-sys_->spl()->ymin();
		double width = sys_->spl()->xmax()-sys_->spl()->xmin();*/
		
		this->sys()->spl()->updateBoundaries(); 
		double height= sys_->spl()->ymax()-sys_->spl()->ymin();
		double width = sys_->spl()->xmax()-sys_->spl()->xmin();
		double periode = 0.;
		double xoffset=5.;
		double c;
		double s;
		double v=0;
		double largeur;
		double decal;
	
		DataSet fn,fnabs,ft,frot;
		for (unsigned int i=0;i<sys_->nwk()->clist().size();++i)
		{
			if( (sys_->nwk()->inter(sys_->nwk()->clist(i))->fn() != 0.)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()>3))
				fn.add( sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
			
			if( (sys_->nwk()->inter(sys_->nwk()->clist(i))->ft() != 0.)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()>3))
				ft.add( sys_->nwk()->inter(sys_->nwk()->clist(i))->ft());
			
			//if( (sys_->nwk()->inter(sys_->nwk()->clist(i))->frot() != 0.)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()>3))
			//	frot.add( sys_->nwk()->inter(sys_->nwk()->clist(i))->frot());
			
			if( (sys_->nwk()->inter(sys_->nwk()->clist(i))->fn() != 0.)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()>3))
				fnabs.add( fabs(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()));
		}
		fn.extractValues();
		double fmin = fn.min();
		double fmax = fn.max();
		//double fmoy = fn.mean(); 

		ft.extractValues();
		double ftmin = ft.min();
		double ftmax = ft.max();
		
		//frot.extractValues();
		//double frmin = frot.min();
		//double frmax = frot.max();
		
		fnabs.extractValues();
		double fabsmin = fnabs.min();
		double fabsmax = fnabs.max();
		
		cout<<"fnmin:= "<<fmin<<"    fnmax:="<<fmax<<endl;		
		cout<<"ftmin:= "<<ftmin<<"    ftmax:="<<ftmax<<endl;
		//cout<<"frotmin:= "<<frmin<<"    frotmax:="<<frmax<<endl;
	
		// pour les disks
		if (sys_->spl()->body(4)->type()==0)
		{
			if( sys_->spl()->isMonoPeriodic()) 
			{	
				width  += 2.*sys_->spl()->bandWidth(); 
				periode = sys_->spl()->boundWidth();
				xoffset=sys_->spl()->bandWidth();
			} 

			//Pour les rectang
			ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
			ps<<"%%BoundingBox: 0 0 842 595"<<endl;
			ps<<"0.5 setlinewidth"<<endl;
			ps<<"0.0 setgray"<<endl;
		
			ps<<"newpath"<<endl;
			ps<<xoffset<<" "<<xoffset<<" moveto"<<endl;
			ps<<zoom_*width<<" 0 rlineto"<<endl;
			ps<<"0 "<<zoom_*height<<" rlineto"<<endl;
			ps<<-zoom_*width<<" "<<0<<" rlineto"<<endl;
			ps<<"closepath"<<endl;
			ps<<"1 setlinewidth"<<endl;
			ps<<"stroke"<<endl;
			
			//Pour les disks
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C: trace un cercle de centre xy"<<endl;
			ps<<"% de rayon r et niveau de gris g  :"<<endl;
			ps<<"% x y r g C "<<endl;
			ps<<"/C{newpath 4 1 roll 0 360 arc gsave setgray fill  grestore stroke }def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C2: trace un cercle plein de centre xy"<<endl;
			ps<<"% de rayon r  :"<<endl;
			ps<<"% x y r  C "<<endl;
			ps<<"/C2{newpath  0 360 arc gsave  fill grestore stroke}def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C3: trace un cercle vide de centre xy"<<endl;
			ps<<"% de rayon r  :"<<endl;
			ps<<"% x y r  C "<<endl;
			ps<<"/C3{newpath  0 360 arc gsave grestore stroke}def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			

			for( unsigned int i=sys_->lctrl().size() ;i<sys_->spl()->lbody().size();++i)
			{
				
				ps<<xoffset+(sys_->spl()->body(i)->x()-sys_->spl()->xmin())*zoom_<<" "<<xoffset+(sys_->spl()->body(i)->y()-sys_->spl()->ymin())*zoom_<<" "
				<<sys_->spl()->body(i)->sizeVerlet()*zoom_<<" 0.85 C"<<endl;
				ps<<"/C{newpath 4 1 roll 0 360 arc gsave setgray fill  grestore stroke }def"<<endl;
				ps<<"0.5 setlinewidth"<<endl;
				ps<<"gsave 0.5 setgray fill grestore"<<endl;
				ps<<"stroke"<<endl;
			}
			
			//Pour les forces
			for( unsigned int i=0 ;i<sys_->nwk()->clist().size();++i)
			{
				if ((sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()>3)&&(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()!=0))
				{
					ps<<"1 0 0 setrgbcolor"<<endl;
					ps<<"newpath"<<endl;
					if (fmax==fmin)
						largeur = 2.1;
					else
						largeur = 0.1 + 5*(fabs(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()) - fabsmin)/(fabsmax-fabsmin);
					
					ps<<largeur<<" setlinewidth"<<endl;
					
					if ( sys_->nwk()->inter(sys_->nwk()->clist(i))->type()== _type_dkdkP) 
					{
						decal = periode;
					}
					else decal = 0.;
					
					ps<<xoffset+zoom_*(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()-sys_->spl()->xmin())<<" "
				 	  <<xoffset+zoom_*(sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()-sys_->spl()->ymin())<<" "<<" moveto"<<endl;
					ps<<xoffset+ zoom_*(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()-sys_->spl()->xmin() + decal)<<" "
					  <<xoffset+zoom_*(sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()-sys_->spl()->ymin())<<" "<<" lineto"<<endl;
					ps<<"stroke"<<endl;
				}
			}
			ps<<"0.5 setlinewidth"<<endl;	
			ps<<"stroke"<<endl;
			ps<<"showpage"<<endl;
			
		}
		
		// pour les polygones
		else
		{
			ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
			ps<<"%%BoundingBox: 0 0 720 620"<<endl;
			//ps<<"%%BoundingBox: 0 0 595 842"<<endl;
			ps<<"1.5 setlinewidth"<<endl;
			ps<<"0.0 setgray"<<endl;
					
			//Tracer le rectang
			ps<<"newpath"<<endl;
			ps<<xoffset<<" "<<xoffset<<" moveto"<<endl;
			ps<<zoom_*width<<" 0 rlineto"<<endl;
			ps<<"3 setlinewidth"<<endl;
			ps<<"stroke"<<endl;

			ps<<"newpath"<<endl;
			ps<<xoffset<<" "<<xoffset+zoom_*height<<" moveto"<<endl;
			ps<<xoffset+zoom_*width<<" "<<xoffset+zoom_*height<<" lineto"<<endl;
			ps<<"3 setlinewidth"<<endl;
			ps<<"stroke"<<endl;

			for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
			v+=sys_->spl()->body(i)->Area();
			
			//Tracer les polygones
			//ps<<"0 0 1 setrgbcolor"<<endl;
			ps<<"0.0 setgray"<<endl;
			for(unsigned int i=sys_->lctrl().size();i<sys_->spl()->lbody().size();i++)
			{
				//polyg * polygo = new polyg;
				polyg * polygo;//Attention reutiliser meme mémoire
				polygo=dynamic_cast<polyg*> (sys_->spl()->body(i));	
				
				c=cos(polygo->rot());
				s=sin(polygo->rot());
				
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(polygo->x()+polygo->Vertex(0).x()*c-polygo->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(polygo->y()+polygo->Vertex(0).x()*s+polygo->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<polygo->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(polygo->x()+polygo->Vertex(j).x()*c-polygo->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(polygo->y()+polygo->Vertex(j).x()*s+polygo->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}
 				//ps<<"gsave 0 1 1 setrgbcolor fill grestore"<<endl;
 				ps<<"gsave 0.85 setgray fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.15 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
			}
			//Tracer la ligne entre deux centre des bodies	

			inter2d * oxo;
			DataSet f1;
			double force;
			
			for( unsigned int i=0 ;i<sys_->nwk()->clist().size();++i)
			{
				force=0.;
				oxo=sys_->nwk()->inter(sys_->nwk()->clist(i));
				
				for(unsigned int r=0 ; r < oxo->rang() ; ++r)
				{
					oxo->current() = r;
					force+=oxo->fn();
				//	cout<<"rang= "<<r<<" fn="<<oxo->fn()<<"force:="<<force<<endl;
				}
				force/=oxo->rang();
				//cout<<"force:="<<force<<endl<<endl;
				f1.add(fabs(force));
			}
			
			f1.extractValues();
			double f1min=f1.min();
			double f1max=f1.max();
				
			for( unsigned int i=0 ;i<sys_->nwk()->clist().size();++i)
			{	
				force=0.;
				oxo=sys_->nwk()->inter(sys_->nwk()->clist(i));
				
				for(unsigned int r=0 ; r < oxo->rang() ; ++r)
				{
					oxo->current() = r;
					force+=oxo->fn();
				}
				
				force/=oxo->rang();
				
				if ((oxo->first()->id()>3)&&(oxo->second()->id()>3))
				{
					if (force==0)
						ps<<"0 0 1 setrgbcolor"<<endl;
					else if(force>0)
						ps<<"1 0 0 setrgbcolor"<<endl;
					else
						ps<<"0 1 0 setrgbcolor"<<endl;
						
					ps<<"newpath"<<endl;
					
 					if (f1max==f1min)
							largeur = 2.6;
					else
							largeur = 0.1 + 5*(fabs(force) - f1min)/(f1max-f1min);
						
					ps<<largeur<<" setlinewidth"<<endl;
					
					ps<<xoffset+zoom_*(oxo->first()->x()-sys_->spl()->xmin())<<" "
				 	  <<xoffset+zoom_*(oxo->first()->y()-sys_->spl()->ymin())<<" "<<" moveto"<<endl;
					ps<<xoffset+zoom_*(oxo->second()->x()-sys_->spl()->xmin())<<" "
					  <<xoffset+zoom_*(oxo->second()->y()-sys_->spl()->ymin())<<" "<<" lineto"<<endl;
					ps<<"stroke"<<endl;
					/*
					ps<<"/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def"<<endl;
					ps<<"/Times-Roman findfont 10 scalefont setfont"<<endl;
					ps<<"1 0 0 setrgbcolor"<<endl;
					ps<<"("<<sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()<<") "<<xoffset+zoom_*(0.5*sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()+0.5*sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()-sys_->spl()->xmin())<<" "<<xoffset+zoom_*(0.5*sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()+0.5*sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()-sys_->spl()->ymin())<<" cprint"<<endl;					
					*/
				}	
				
					//Tracer position des contacts
				//if (oxo->rang()==1)
				/*for(unsigned int r=0 ; r < oxo->rang() ; ++r)
				{
					oxo->current() = r;

				ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
				ps<<"% Procedure cprint: position des contacts"<<endl;
				ps<<"% (texte) x y cprint  :"<<endl;
				ps<<"/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def"<<endl;
				ps<<"/Times-Roman findfont 10 scalefont setfont"<<endl;
				ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
				ps<<"0 0 1 setrgbcolor"<<endl;
				
				ps<<"("<<round(oxo->fn())<<") "<<(1.0+0.01*r)*(xoffset+zoom_*(oxo->x()-sys_->spl()->xmin()))<<" "<<(1.0+0.01*r)*(xoffset+zoom_*(oxo->y()-sys_->spl()->ymin()))<<" cprint"<<endl;
				
				}*/
				
				/*double sum=0.;
				double xsum=0.;
				double ysum=0.;
				for(unsigned int r=0 ; r < oxo->rang() ; ++r)
				{
					oxo->current() = r;
					sum+=oxo->fn();
					xsum+=oxo->x();
					ysum+=oxo->y();
				}
				//sum/=oxo->rang();
				xsum/=oxo->rang();
				ysum/=oxo->rang();
				
				ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
				ps<<"% Procedure cprint: position des contacts"<<endl;
				ps<<"% (texte) x y cprint  :"<<endl;
				ps<<"/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def"<<endl;
				ps<<"/Times-Roman findfont 10 scalefont setfont"<<endl;
				ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
				ps<<"0 0 1 setrgbcolor"<<endl;
				
				ps<<"("<<sum<<") "<<1.0*(xoffset+zoom_*(xsum-sys_->spl()->xmin()))<<" "<<1.0*(xoffset+zoom_*(ysum-sys_->spl()->ymin()))<<" cprint"<<endl;
				*/
				
				
			}
			cout<<endl;
		}	
		
		//Tracer le text
			/*ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure cprint: texte centre"<<endl;
			ps<<"% (texte) x y cprint  :"<<endl;
			ps<<"/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def"<<endl;
			ps<<"/Times-Roman findfont 10 scalefont setfont"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"1 0 0 setrgbcolor"<<endl;
			for(unsigned int i=sys_->lctrl().size();i<sys_->spl()->lbody().size();i++)
			{
				ps<<"("<<sys_->spl()->body(i)->id()-4<<") "<<1.01*(xoffset+zoom_*(sys_->spl()->body(i)->x()-sys_->spl()->xmin()))<<" "<<1.01*(xoffset+zoom_*(sys_->spl()->body(i)->y()-sys_->spl()->ymin()))<<" cprint"<<endl;
			}*/
	ps.close();
}

//Les particules en noir pour étudier la distribution des volumes
void brazilian_A::writePS2( const char * fname)
	{
				
		ofstream ps(fname);

		this->sys()->spl()->updateBoundaries(); 		//Ajouter le 25/11/2011 	
		double height= sys_->spl()->ymax()-sys_->spl()->ymin();
		double width = sys_->spl()->xmax()-sys_->spl()->xmin();
		double xoffset=5.;
		zoom_=min((595-2*xoffset)/height,(842-2*xoffset)/width);
		double c;
		double s;
		double v=0;

		DataSet fn;
		for (unsigned int i=0;i<sys_->nwk()->clist().size();++i)
		{
			if( sys_->nwk()->inter(sys_->nwk()->clist(i))->fn() != 0.)
				fn.add( sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
		}
		fn.extractValues();
		if (sys_->spl()->body(4)->type()==0)
		{

			ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
			ps<<"%%BoundingBox: 0 0 842 595"<<endl;

			ps<<"%%Pages: 2"<<endl;
			ps<<"1 setlinecap"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure cprint: texte centre"<<endl;
			ps<<"% (texte) x y cprint  :"<<endl;
			ps<<"/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def"<<endl;
			ps<<"/Times-Roman findfont 15 scalefont setfont"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C: trace un cercle de centre xy"<<endl;
			ps<<"% de rayon r et niveau de gris g  :"<<endl;
			ps<<"% x y r g C "<<endl;
			ps<<"/C{newpath 4 1 roll 0 360 arc gsave setgray fill  grestore stroke }def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C2: trace un cercle plein de centre xy"<<endl;
			ps<<"% de rayon r  :"<<endl;
			ps<<"% x y r  C "<<endl;
			ps<<"/C2{newpath  0 360 arc gsave  fill grestore stroke}def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure C3: trace un cercle vide de centre xy"<<endl;
			ps<<"% de rayon r  :"<<endl;
			ps<<"% x y r  C "<<endl;
			ps<<"/C3{newpath  0 360 arc gsave grestore stroke}def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure Lic: ligne intercentre de largeur donnee"<<endl;
			ps<<"% largeur t, du point 0 au point 1"<<endl;
			ps<<"% x0 y0 x1 y1 t Lic"<<endl;
			ps<<"/Lic{gsave setrgbcolor 0.1 mul setlinewidth newpath moveto lineto stroke grestore}def"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"% Procedure Lis: ligne de largeur fixe"<<endl;
			ps<<"% du point 0 au point 1"<<endl;
			ps<<"% x0 y0 x1 y1 Lis"<<endl;
			ps<<"/Lis{gsave newpath moveto lineto stroke grestore}def"<<endl; 
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			//ps<<"100 300 translate 0.2 dup scale"<<endl;
			ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
			//******************

					
			//Tracer le rectang
			ps<<"newpath"<<endl;
			ps<<xoffset<<" "<<xoffset<<" moveto"<<endl;
			ps<<zoom_*width<<" 0 rlineto"<<endl;
			ps<<"0 "<<zoom_*height<<" rlineto"<<endl;
			ps<<-zoom_*width<<" "<<0<<" rlineto"<<endl;
			ps<<"closepath"<<endl;
			ps<<"1 setlinewidth"<<endl;
			ps<<"stroke"<<endl;
			
			ps<<"/C2{newpath  0 360 arc gsave  fill grestore stroke}def"<<endl;
			ps<<"0.0 setgray"<<endl;
			for( unsigned int i=sys_->lctrl().size();i<sys_->spl()->lbody().size();++i)
			ps<<xoffset+(sys_->spl()->body(i)->x()-sys_->spl()->xmin())*zoom_<<" "<<xoffset+(sys_->spl()->body(i)->y()-sys_->spl()->ymin())*zoom_<<" "
			<<sys_->spl()->body(i)->sizeVerlet()*zoom_<<" C2"<<endl;
			
			ps<<"stroke"<<endl;
			ps<<"showpage"<<endl;
			ps.close();
		}
		else
		{
			ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
			ps<<"%%BoundingBox: 0 0 842 595"<<endl;
			ps<<"0.3 setlinewidth"<<endl;
			ps<<"0.3 setgray"<<endl;
					
			//Tracer le rectang
			ps<<"newpath"<<endl;
			ps<<xoffset<<" "<<xoffset<<" moveto"<<endl;
			ps<<zoom_*width<<" 0 rlineto"<<endl;
			ps<<"0 "<<zoom_*height<<" rlineto"<<endl;
			ps<<-zoom_*width<<" "<<0<<" rlineto"<<endl;
			ps<<"closepath"<<endl;
			ps<<"1 setlinewidth"<<endl;
			ps<<"stroke"<<endl;

			for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
			v+=sys_->spl()->body(i)->Area();
			
			//Tracer les polygones
			ps<<"0.0 setgray"<<endl;
			for(unsigned int i=sys_->lctrl().size();i<sys_->spl()->lbody().size();i++)
			{
				//polyg * polygo = new polyg;
				polyg * polygo;//Attention reutiliser meme mémoire
				polygo=dynamic_cast<polyg*> (sys_->spl()->body(i));	
				
				c=cos(polygo->rot());
				s=sin(polygo->rot());
				
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(polygo->x()+polygo->Vertex(0).x()*c-polygo->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(polygo->y()+polygo->Vertex(0).x()*s+polygo->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<polygo->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(polygo->x()+polygo->Vertex(j).x()*c-polygo->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(polygo->y()+polygo->Vertex(j).x()*s+polygo->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}
 				ps<<"gsave 0.0 setgray fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.01 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
			}
		}	
}

	
//sortie VTK
void brazilian_A::display_data(char * fname)
{
	//DataSet fn;
	disk* d_;
	
	ofstream O_VTK(fname,ios::app);
	
	O_VTK << "# vtk DataFile Version 3.0" << endl << "Sortie Grains sans" << endl << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID" << endl;
	O_VTK << "POINTS" << " " << sys_->spl()->lbody().size() - 4 << " " << "float" << endl;
	
	//écriture de la position des particules
	for (unsigned long int i=0; i<sys_->spl()->lbody().size(); i++)
	{
		if (sys_->spl()->body(i)->type() == 0)
		{
			O_VTK << sys_->spl()->body(i)->x() << " " << sys_->spl()->body(i)->y() << " 0." << endl;
		}
	}
	
	O_VTK << "POINT_DATA" << " " << sys_->spl()->lbody().size() - 4 << endl;
	
	//affectation d'un vecteur contenant les rayons des particules 
	O_VTK << "VECTORS Radius float" << endl;
	for (unsigned long int i=0; i<sys_->spl()->lbody().size(); i++)
	{
		if (sys_->spl()->body(i)->type() == 0)
		{
			d_ = dynamic_cast< disk*>(sys_->spl()->body(i));
			O_VTK << "0. 0. " << d_->R() << endl;
		}
	}
	
	//affichage en fonction de la vitesse des particules
	O_VTK << "VECTORS Velocity float" << endl;
	for (unsigned long int i=0; i<sys_->spl()->lbody().size(); i++)
	{
		if (sys_->spl()->body(i)->type() == 0)
		{
			O_VTK << sys_->spl()->body(i)->vx() << " " << sys_->spl()->body(i)->vy() << " 0." << endl;
		}
	}
	
	O_VTK.close();
}

bool brazilian_A::Inter (polyg *p, polyg *q)
{
	bool Inter = false;
	double c,s;
	polyg *r=new polyg; //polygone qui est parti commun de p et q
	gdm::vertex q1, q2, vtex;
	vector<gdm::vertex> Psection; //Point section entre deux polygone
	
	r->y()=0;
	r->x()=0;
	r->rot()=0.0;
	r->Vertex().clear();
	//garder les vertex de p qui appartient q
//	cout<<"ligne 3057:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<p->Vertex().size();i++)
	{
		c=cos(p->rot());
		s=sin(p->rot());
		vtex.x()=p->x()+p->Vertex(i).x()*c-p->Vertex(i).y()*s;
		vtex.y()=p->y()+p->Vertex(i).x()*s+p->Vertex(i).y()*c;
		if (q->Contain(vtex)) r->Vertex().push_back(vtex);
	}
	//garder les vertex de q qui appartient p
//	cout<<"ligne 3067:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<q->Vertex().size();i++)
	{
		c=cos(q->rot());
		s=sin(q->rot());
		vtex.x()=q->x()+q->Vertex(i).x()*c-q->Vertex(i).y()*s;
		vtex.y()=q->y()+q->Vertex(i).x()*s+q->Vertex(i).y()*c;
		if (p->Contain(vtex)) r->Vertex().push_back(vtex);
	}
	//ajouter les points section entre les segement du p et q
//	cout<<"ligne 3077:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<q->Vertex().size();i++)
	{
		c=cos(q->rot());
		s=sin(q->rot());
		unsigned int j=(i+1)%q->Vertex().size();
		q1.x()=q->x()+q->Vertex(i).x()*c-q->Vertex(i).y()*s;
		q1.y()=q->y()+q->Vertex(i).x()*s+q->Vertex(i).y()*c;
		q2.x()=q->x()+q->Vertex(j).x()*c-q->Vertex(j).y()*s;
		q2.y()=q->y()+q->Vertex(j).x()*s+q->Vertex(j).y()*c;

		Psection=p->SectionSegment(q1,q2);
		if (!Psection.empty())
		for (unsigned int k=0; k<Psection.size();k++) r->Vertex().push_back(Psection[k]);
	}
//	cout<<"ligne 3092:= "<<r->Vertex().size()<<endl;
	
	if (r->Vertex().size() > 2) Inter =true ;
	return Inter;
}


polyg * brazilian_A::Intersection (polyg *p, polyg *q)
{
	double c,s;
	polyg *r=new polyg; //polygone qui est parti commun de p et q
	gdm::vertex q1, q2, vtex;
	vector<gdm::vertex> Psection; //Point section entre deux polygone

	r->y()=0;
	r->x()=0;
	r->rot()=0.0;
	r->Vertex().clear();
	//garder les vertex de p qui appartient q
//	cout<<"ligne 3111:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<p->Vertex().size();i++)
	{
		c=cos(p->rot());
		s=sin(p->rot());
		vtex.x()=p->x()+p->Vertex(i).x()*c-p->Vertex(i).y()*s;
		vtex.y()=p->y()+p->Vertex(i).x()*s+p->Vertex(i).y()*c;
		if (q->Contain(vtex)) r->Vertex().push_back(vtex);
	}
	//garder les vertex de q qui appartient p
//	cout<<"ligne 3121:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<q->Vertex().size();i++)
	{
		c=cos(q->rot());
		s=sin(q->rot());
		vtex.x()=q->x()+q->Vertex(i).x()*c-q->Vertex(i).y()*s;
		vtex.y()=q->y()+q->Vertex(i).x()*s+q->Vertex(i).y()*c;
		if (p->Contain(vtex)) r->Vertex().push_back(vtex);
	}
	//ajouter les points section entre les segement du p et q
//	cout<<"ligne 3131:= "<<r->Vertex().size()<<endl;
	for(unsigned int i=0;i<q->Vertex().size();i++)
	{
		c=cos(q->rot());
		s=sin(q->rot());
		unsigned int j=(i+1)%q->Vertex().size();
		q1.x()=q->x()+q->Vertex(i).x()*c-q->Vertex(i).y()*s;
		q1.y()=q->y()+q->Vertex(i).x()*s+q->Vertex(i).y()*c;
		q2.x()=q->x()+q->Vertex(j).x()*c-q->Vertex(j).y()*s;
		q2.y()=q->y()+q->Vertex(j).x()*s+q->Vertex(j).y()*c;

		Psection=p->SectionSegment(q1,q2);
		if (!Psection.empty())
		for (unsigned int k=0; k<Psection.size();k++) r->Vertex().push_back(Psection[k]);
	}
//	cout<<"ligne 3146:= "<<r->Vertex().size()<<endl;
	
	//arranger les vertex du polygone r 
	r->Arrangement();
	//Invert ordre vertex polygone if CW
	if (r->ispolysimple()<0) r->Invert();
	r->adjustCenter();
	r->Fill(1);
	return r;
}

double brazilian_A::compactness(double &h,double &l)
{
	double compact=0;
	
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	compact+=sys_->spl()->body(i)->Area();
	compact/=(h*l);
	return compact;
}


double brazilian_A::compactness(double &x, double &y, double &h, double &l)
{
	double compact=0;
	polyg *Rec=new polyg;
	polyg *poly=new polyg;
	gdm::vertex vertex_;
	//cout<<"hight:="<<h<<"   large:="<<l<<endl;
	
	Rec->x()=0;
	Rec->y()=0;
	Rec->rot()=0;
	Rec->Vertex().clear();
	vertex_.x()=x-0.5*l; vertex_.y()=y-0.5*h;
	Rec->Vertex().push_back(vertex_);
	vertex_.x()=x+0.5*l; vertex_.y()=y-0.5*h;
	Rec->Vertex().push_back(vertex_);
	vertex_.x()=x+0.5*l; vertex_.y()=y+0.5*h;
	Rec->Vertex().push_back(vertex_);
	vertex_.x()=x-0.5*l; vertex_.y()=y+0.5*h;
	Rec->Vertex().push_back(vertex_);
	
		
	Rec->adjustCenter();
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	{
		polyg  * polyg_=new polyg;
		polyg_=dynamic_cast<polyg*> (sys_->spl()->body(i));
		if (Inter(Rec,polyg_)) 
		{
			poly=Intersection(Rec,polyg_);
			if ((poly->Area()>0) ) compact+=poly->Area();
			else 
			{
				cout<<"Aire polygone negatif ! "<<poly->Area()<<endl; 
				//exit(0);
			}	
		}
	}
	compact/=(h*l);
	//cout<<"Compactness: "<<compact<<endl;
	return compact;
}

double brazilian_A::compactness_disk(double &x, double &y, double &h, double &l)
{
	double compact=0;
	double d=0;
	double R,X,Y;
	double a,b,c;
	double alpha;
	double area=0;
	
	prb.x()=x;
	prb.y()=y;
	prb.hh()=0.5*h;
	prb.hl()=0.5*l;
	
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	{
		body_=sys_->spl()->body(i);
		R=body_->sizeVerlet();
		X=body_->x();
		Y=body_->y();
		if (prb.containEntireBody(body_)) area=body_->Area();
		else if (prb.intersection(body_))
		{
			if (y-0.5*h>body_->ymin())
			{
				if(x+0.5*l<body_->xmax())
				{
					a=x+0.5*l-(X-sqrt(R*R-(y-0.5*h-Y)*(y-0.5*h-Y)));
					b=(Y+sqrt(R*R-(x+0.5*l-X)*(x+0.5*l-X)))-(y-0.5*h);
					c=sqrt(a*a+b*b);
					alpha=asin(c/(2*R));
					area=alpha*R*R-0.5*R*R*sin(2*alpha)+0.5*a*b;
				}
				else if (x-0.5*l>body_->xmin())
				{
					a=(X+sqrt(R*R-(y-0.5*h-Y)*(y-0.5*h-Y)))-(x-0.5*l);
					b=(Y+sqrt(R*R-(x-0.5*l-X)*(x-0.5*l-X)))-(y-0.5*h);
					c=sqrt(a*a+b*b);
					alpha=asin(c/(2*R));
					area=alpha*R*R-0.5*R*R*sin(2*alpha)+0.5*a*b;
				}
				else
				{
					d=y-0.5*h-body_->y();
					alpha =acos(fabs(d)/R);
					area=alpha*R*R-0.5*R*R*sin(2*alpha);
					if (d<0) area=M_PI*R*R-area;
				}
			}
			
			else if (y+0.5*h<body_->ymax())
			{
				if(x+0.5*l<body_->xmax())
				{
					a=x+0.5*l-(X-sqrt(R*R-(y+0.5*h-Y)*(y+0.5*h-Y)));
					b=y+0.5*h-(Y-sqrt(R*R-(x+0.5*l-X)*(x+0.5*l-X)));
					c=sqrt(a*a+b*b);
					alpha=asin(c/(2*R));
					area=alpha*R*R-0.5*R*R*sin(2*alpha)+0.5*a*b;
				}
				else if (x-0.5*l>body_->xmin())
				{
					a=(X+sqrt(R*R-(y+0.5*h-Y)*(y+0.5*h-Y)))-(x-0.5*l);
					b=y+0.5*h-(Y-sqrt(R*R-(x-0.5*l-X)*(x-0.5*l-X)));
					c=sqrt(a*a+b*b);
					alpha=asin(c/(2*R));
					area=alpha*R*R-0.5*R*R*sin(2*alpha)+0.5*a*b;
				}
				else
				{
					d=y+0.5*h-body_->y();
					alpha =acos(fabs(d)/R);
					area=alpha*R*R-0.5*R*R*sin(2*alpha);
					if (d>0) area=M_PI*R*R-area;
				}
			}
			
			else if (x+0.5*l<body_->xmax())
			{
				d=x+0.5*l-body_->x();
				alpha =acos(fabs(d)/R);
				area=alpha*R*R-0.5*R*R*sin(2*alpha);
				if (d>0) area=M_PI*R*R-area;
			}
			
			else if (x-0.5*l>body_->xmin()) 
			{
				d=x-0.5*l-body_->x();
				alpha =acos(fabs(d)/R);
				area=alpha*R*R-0.5*R*R*sin(2*alpha);
				if (d<0) area=M_PI*R*R-area;
			}
		}
		else area=0;
		compact+=area;
	}
	cout<<endl;
	compact/=(h*l);
	cout<<"Compactness: "<<compact<<endl;
	return compact;
}

void brazilian_A::fracture()
{
	ofstream fra("Analyse/fracture.txt",ios::app);
	fra<<time<<"\t"<<sys_->grpRel()->fragmentation()->numfissure()<<"\t"<<sys_->grpRel()->fragmentation()->ncontact()<<"\t"<<double(sys_->grpRel()->fragmentation()->numfissure())/double(sys_->grpRel()->fragmentation()->ncontact())<<"\t"<<-defy<<"\t"<<Z_<<"\t"<<a()<<"\t"<<an()<<"\t"<<at()<<"\t"<<al()<<"\t"<<da_/M_PI*180.<<"\t"<<max(s1,s2)/sys_->grpRel()->fragmentation()->sigmac()<<endl;
	fra.close();

	cout<<"Nombre fissure : "<<sys_->grpRel()->fragmentation()->numfissure()<<endl;
}		
