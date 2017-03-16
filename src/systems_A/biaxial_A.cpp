#include "biaxial_A.hpp"

void Biaxial_A::plugRef()
{
	top_    = dynamic_cast< rline*>(sys_->spl()->body(topId));
	bottom_ = dynamic_cast< rline*>(sys_->spl()->body(bottomId));
	left_   = dynamic_cast< rline*>(sys_->spl()->body(leftId));
	right_  = dynamic_cast< rline*>(sys_->spl()->body(rightId));
}

void Biaxial_A::read_parameters(istream & is)
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
        else if (token == "ContactMesh_")   ContactMesh_ = true;
		else if (token == "ZoomSample") is >> zoom_;
		else if (token == "Sample")		displaySample = true;
		else if (token== "Cluster") calClust=true;
        else if(token== "Force") displayForce=true;
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
				cout<<" @Biaxial_A :Zgranulo : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ Biaxial_A :Zgranulo : Nbingranulo undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Granulostress") 
		{
			calcgranulostress=true;
			if (Nbingranulo !=0 )
			{
				cout<<" @Biaxial_A :granulostress : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ Biaxial_A :granulostress:  Nbingranulo undefined "<<endl;
				exit(0);
			}
			is>>perG;
			if ( perG == 0 )
			{
				cout<<" @ Biaxial_A :granulostress:  perG undefined "<<endl;
				exit(0);
			}
			is>>wG;
			if ( wG == 0 )
			{
				cout<<" @ Biaxial_A :granulostress:  wG undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="PDFFN") 
		{
			calcfn=true;
			is >> NbinFN;

			if ( NbinFN == 0 )
			{
				cout<<" @ Biaxial_A : NbinFN undefined "<<endl;
				exit(0);
			}
			is >> perF;
			if ( perF == 0 )
			{
				cout<<" @ Biaxial_A : perF undefined "<<endl;
				exit(0);
			}
			is >> wF;
			if ( wF == 0 )
			{
				cout<<" @ Biaxial_A : wF undefined "<<endl;
				exit(0);
			}

		}
		else if(token=="Ptheta") 
		{
			calcPtheta=true;
			is >> NbinPT;

			if ( NbinPT == 0 )
			{
				cout<<" @ Biaxial_A : NbinPT undefined "<<endl;
				exit(0);
			}
			is >> mobperiod;
			if ( mobperiod == 0 )
			{
				cout<<" @ Biaxial_A : mobperiod undefined "<<endl;
				exit(0);
			}
			is >> mobwidth;
			if ( mobwidth == 0 )
			{
				cout<<" @ Biaxial_A : mobwidth undefined "<<endl;
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
				cout<<" @ Biaxial_A : NbinFT undefined "<<endl;
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
				cout<<" @Biaxial_A :granulobin : Nbingranulo already defined : "<<Nbingranulo<<endl;
				is>>trash;
			}
			else
				is>>Nbingranulo;

			if ( Nbingranulo == 0 )
			{
				cout<<" @ Biaxial_A :granulobin : Nbingranulo undefined "<<endl;
				exit(0);
			}
		}
		else if(token=="Nanalyze") is>>Nanalyze();
		else if(token=="RFD") {calrdf=true; is>>npoint>>nrmean;}
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
		else if(token=="Granulometrique") 
		{	
			is>>Nbingranulometrique;
			calfgranulometrique=true;
		}
		else if(token=="Heterogeneity") 
		{
			is>>Ni>>Nf;
			calheterogeneity=true;
		}
		else if(token=="Clusterdouble") 
		{
			calclusterdouble=true;
		}
		else if(token=="}") break;
		//else if(token=="include_interrupt") incl_int=true;
		else cout<<token<<" : parametre de commande inconnu "<<endl;

		is>>token;
	}	
	
}

void Biaxial_A::analyse( double t, unsigned int nsi, unsigned int nsf )
{
	cout<<"--------------Biaxial_A::analyze()---------"<<endl;
	
	char fname[100];
	time=t;
	

	//Definir la probe avec les 4 premier corps du spl
	double R = left_->R();
	double hx=.5*(right_->x() - left_->x())   -R;
	double hy=.5*(top_->y()   - bottom_->y()) -R;
	
	totalProbe_.x() = left_->x()+ hx+ R;
	totalProbe_.y() = bottom_->y() + hy+ R ;
	totalProbe_.R() = .9*min(hx,hy);
	
	prb_.x() = left_->x()+ hx+ R;
	prb_.y() = bottom_->y() + hy+ R ;
	prb_.hh() = 1.*hy;
	prb_.hl() = 1.*hx;
	
	prb2.x()	= prb_.x();
	prb2.y()	= prb_.y();
	prb2.hh() = 0.9*hy;
	prb2.hl() = 0.9*hx;
	
	unsigned int Nb=sys_->spl()->lbody().size()-sys_->lctrl().size();
	double vmoy=0.;
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();	i++)
		vmoy+=sqrt(pow(sys_->spl()->body(i)->vx(),2)+pow(sys_->spl()->body(i)->vy(),2));
	vmoy=vmoy/Nb;
	
	ofstream paroi("Analyse/paroi.txt",ios::app);
		paroi<<time<<"\t"<<-sys_->ctrl(topId).yval()/(right_->xmin()-left_->xmax())<<"\t"<<top_->vy()<<"\t"<<vmoy<<endl;
		paroi.close();
    
    //
    if (ContactMesh_)   ContactMesh();
    
    
    
	if (calcz) Z(calczp,Nbingranulo,granulo); 
	if( removeR) removeRattlers();
	if( growR) growRattlers( sys_->spl()->rmin()*incR);

	if (calcglobalstress) globalStress();

	if (calcfn) pdfforce(true,NbinFN,normpdf,perF,wF);
	if (calcft) pdfforce(false,NbinFT,normpdf,perF,wF);
	if (calcsf) SF();
	
	if (calcFabric)  {A(); A2();}
	if (calcforcesA)	{forces_A(NbinFA);forces2_A(NbinFA);}
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
	if (calclusterdouble)
	{
		clustercontactdouble();
		cluster_df();
		cluster_volume();
		cluster_shape();
	}
	
	if(calrdf)
	{
		rfd(npoint,nrmean);
	}
	if (calheterogeneity && (nsi==nsf) )
	{
		heterogeneity(Ni,Nf);
	}
	
	if (calcompact && (divisionl==0)) 
	{
		double ratio_=divisionh/10.;
		compactness(ratio_);
	}
	else if(calcompact && (divisionh==1) && (divisionl==1))
	{
		sys_->spl()->updateBoundaries();
		
		double xmin=sys_->spl()->xmin();
		double xmax=sys_->spl()->xmax();
		double ymin=sys_->spl()->ymin();
		double ymax=sys_->spl()->ymax();
		double hight=ymax-ymin;
		double large=xmax-xmin;
		if(sys_->spl()->body(4)->type()==0)
		{
			const double g=10;
			double dmoy=0;
			double mmoy=0;
			double vmax=0;
			double v;
			double Ec=0;
			unsigned int Nb;
			//Calcule vmax,Ec moyence
			Nb=sys_->spl()->lbody().size()-sys_->lctrl().size();
			for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
			{
				dmoy+=2*sys_->spl()->body(i)->sizeVerlet();
				mmoy+=sys_->spl()->body(i)->mass();
				v=sqrt(pow(sys_->spl()->body(i)->vx(),2)+pow(sys_->spl()->body(i)->vy(),2));
				vmax=(v>vmax) ? v : vmax;
				Ec+=0.5*sys_->spl()->body(i)->mass()*pow(v,2);
			}
			dmoy/=Nb;
			mmoy/=Nb;
			Ec=Ec/(Nb*mmoy*g*dmoy);
			vmax/=sqrt(g*dmoy);//vmaxcritique
			ofstream str("Analyse/compactness.txt",ios::app);
			str<<time<<"\t"<<compactness(hight,large)<<"\t"<<vmax<<"\t"<<Ec<<endl;
			str.close();		
		}
		else if ((sys_->ctrl(rightId).x()==1) && (sys_->ctrl(rightId).xval()==0) && (sys_->ctrl(topId).yval()!=0))
		{
			ofstream indice("Analyse/indice_de_vide.txt",ios::app);
			//indice<<time<<" "<<-sys_->ctrl(topId).yval()/(right_->xmin()-left_->xmax())/sys_->grpRel()->fragmentation()->sigmac()<<" "<<-sys_->ctrl(topId).yval()/(right_->xmin()-left_->xmax())<<" "<<1./compactness(hight,large)-1.<<"	"<<compactness(hight,large)<<endl;
			indice<<time<<" "<<(0.1+0.1*t/1.0E-07*0.002)<<" "<<(0.1+0.1*t/1.0E-07*0.002)/10.<<" "<<1./compactness(hight,large)-1.<<"	"<<compactness(hight,large)<<endl;
			indice.close();
		}
		else 
		{
			ofstream str("Analyse/compactness.txt",ios::app);
			str<<time<<"\t"<<compactness(hight,large)<<endl;
			str.close();
		}
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
		double centrex;
		double centrey;
		double comp=0.;
		vector<double> compact;
		vector <vector<double> > Compact;
		
		for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
		comp+=sys_->spl()->body(i)->Area();
		comp/=((ymax-ymin)*(xmax-xmin));
		cout<<"Compact:="<<comp<<endl;
	
		Compact.clear();
		for (unsigned int i=0;i<divisionl;i++)
		{
			compact.clear();
			for(unsigned int j=0;j<divisionh;j++)
			{
				centrex=xmin+(i+0.5)*large;
				centrey=ymin+(j+0.5)*hight;
				if (sys_->spl()->body(4)->type()==0)
					compact.push_back(compactness_disk(centrex,centrey,hight,large));
				else
					compact.push_back(compactness(centrex,centrey,hight,large));
			}
			Compact.push_back(compact);
		}
		
		ofstream str("Analyse/compactness.txt",ios::out);

		for (unsigned int i=0;i<divisionl;i++)
			for(unsigned int j=0;j<divisionh;j++)
				str<<(i+0.5)/double(divisionl)<<"\t"<<(j+0.5)/(divisionh)<<"\t"<<Compact[i][j]/comp<<endl;	
		str.close();
		
		ofstream strl("Analyse/compactl.txt",ios::out);
		double compactl;
		for (unsigned int i=0;i<divisionl;i++)
		{
			compactl=0.;
			for(unsigned int j=0;j<divisionh;j++) compactl+=Compact[i][j];
			compactl/=double(divisionh);
			strl<<(i+0.5)/double(divisionl)<<"\t"<<compactl/comp<<endl;	
		}
		strl.close();

		ofstream strh("Analyse/compacth.txt",ios::out);
		double compacth;
		for (unsigned int j=0;j<divisionh;j++)
		{
			compacth=0.;
			for(unsigned int i=0;i<divisionl;i++) compacth+=Compact[i][j];
			compacth/=double(divisionl);
			strh<<compacth/comp<<"\t"<<(j+0.5)/double(divisionh)<<"\t"<<compacth/comp<<endl;	
		}
		strh.close();

	}	

	cout<<endl;
	cout<<"		.5*(a+an+at+aln+alt) = "<<.5*(a()+an()+at()+aln()+alt())<<endl;
	cout<<"		         q/p = "<<(max(s1,s2)-min(s1,s2))/(s1+s2)<<endl;
	cout<<"		Nombre d'interaction = "<<sys_->nwk()->linter().size()<<endl;
	cout<<"		Nombre de contacts   = "<<sys_->nwk()->clist().size()<<endl;
	
		if( calcforcesA)
		{
			ofstream FA_out1("Analyse/Anisotropy1.txt",ios::app);
			FA_out1<<time<<" "<<epsq<<" "<<a()<<" "<<an()<<" "<<at()<<" "<<aln()<<" "<<alt()<<" "<<.5*( a()+an()+at()+aln()+alt())<<" "<<qop_<<" " << (s1+s2)/2. << " "<<qop_*(s1+s2)/2.<<endl;
			FA_out1.close();

			ofstream FA_out2("Analyse/Anisotropy2.txt",ios::app);
			FA_out2<<time<<" "<<epsq<<" "<<a2()<<" "<<an2()<<" "<<at2()<<" "<<al2()<<" "<<.5*(a2()+an2()+at2()+al2())<<" "<<qop_<<" " << (s1+s2)/2. << " "<<qop_*(s1+s2)/2.<<endl;
			FA_out2.close();

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
		if (nsi==nsf)
		{
			sprintf(fname,"Analyse/PS2/particles%.4li.ps",numFile_);
			writePS2(fname);	
		}
		numFile_++;
	}
	if (calfracture) fracture();
	if (calfgranulometrique) granulometrique(Nbingranulometrique);
    
    Nanalyze() ++;
}


void Biaxial_A::initAnalyse( ) 
{ 
	cout<<"--------------Biaxial_A::init()-----------"<<endl;	
	
	top_    = (dynamic_cast< Biaxial*>( this->sys()))->top();
	bottom_ = (dynamic_cast< Biaxial*>( this->sys()))->bottom();
	left_   = (dynamic_cast< Biaxial*>( this->sys()))->left();
	right_  = (dynamic_cast< Biaxial*>( this->sys()))->right();
	
	topId    = top_->id();
	bottomId = bottom_->id();
	leftId   = left_->id();
	rightId  = right_->id();
	
	cout<<" top "<<topId<<" bottom "<<bottomId<<" right "<<rightId<<" left "<<leftId<<endl;
	
	sys_->spl()->radiusExtrema(4);
	cout<<" rmin "<<sys()->spl()->rmin()<<" rmax "<<sys()->spl()->rmax()<<" rmoy "<<sys()->spl()->rmoy()<<endl;
	double Areamoy=0.,massmoy=0.;
	
	for (unsigned int i=0; i<4;i++) 
	cout<<"i:=	"<<i<<"	grp:=	"<<sys_->spl()->body(i)->grp()<<"	Area:=		"<<sys_->spl()->body(i)->Area()<<"	mass:=	"<<sys_->spl()->body(i)->mass()<<endl;
	
	for (unsigned int i=4; i<sys_->spl()->lbody().size();i++) 
	{
		Areamoy+=sys_->spl()->body(i)->Area();
		massmoy+=sys_->spl()->body(i)->mass();
	}
	
	Areamoy/=(sys_->spl()->lbody().size()-4.0);
	massmoy/=(sys_->spl()->lbody().size()-4.0);

	cout<<"Area moyen:	"<<Areamoy<<"	massmoy:=	"<<massmoy<<endl;
	
	y0 = top_ ->y();
	x0 = right_->x();
	l0=right_->x() - left_->x() - 2*left_->R();
	h0=top_->y()   - bottom_->y() - 2*left_->R();
	
	defx=defy=0.;
	ngap1_=ngap2_=NULL;
	 
	double R = left_->R();
	double hx=.5*(right_->x() - left_->x())   -R;
	double hy=.5*(top_->y()   - bottom_->y()) -R;
	
	totalProbe_.x() = left_->x()+ hx+ R;
	totalProbe_.y() = bottom_->y() + hy+ R ;
	totalProbe_.R() = 1.0*min(hx,hy);
	
	prb_.x() = left_->x()+ hx+ R;
	prb_.y() = bottom_->y() + hy+ R ;
	prb_.hh() = 1.*hy;
	prb_.hl() = 1.*hx;	
	//zoom_=min(585/hy/2.,832/hx/2.); //hinh chu nhat nam
	//zoom_=1.0*min(832/hy/2.,585/hx/2.);   //hcn dung
	
	system("mkdir Analyse");
	system("mkdir ContactMesh");

	ofstream monitoring("Analyse/monitoring.txt",ios::out);
	monitoring.close();
	
	ofstream analyse("Analyse/analyse.txt",ios::out);
	analyse.close();

	ofstream strain("Analyse/strain.txt",ios::out);
	strain.close();
	
	ofstream str("Analyse/compactness.txt",ios::out);
	str.close();
	
	ofstream indice("Analyse/indice_de_vide.txt",ios::out);
	indice.close();
	
	
	// création du fichier measure
	//ofstream measure("Analyse/measure.txt",ios::out);
	//measure.close();
	
	if(calrdf)
	{
		ofstream rfd("Analyse/rfd.txt",ios::out);
		rfd.close();
	}
	
	if ( calcforcesA ) 
	{
			ofstream FA_out1("Analyse/Anisotropy1.txt",ios::out); 
			FA_out1<<"time epsq ac an at aln alt .5*(ac+an+at...)  q/p"<<endl;
			FA_out1.close();
			
			ofstream FA_out2("Analyse/Anisotropy2.txt",ios::out); 
			//FA_out<<"time epsq ac an at aln alt .5*(ac+an+at...)  q/p"<<endl;
			FA_out2.close();
			
			ofstream DPM_out("Analyse/DPM.txt",ios::out); 
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
		GS_out<<"time	s1	s2	max(s1,s2)-min(s1,s2)	epsq	qop_	(s1+s2)/2.	qop_*(s1+s2)/2.	ds_"<<endl;
		GS_out.close();
		
		ofstream U_compression_data_out("Analyse/U_compression.dat",ios::out);
		U_compression_data_out<<"#   P   rho   Nb_contact"<<endl;
		U_compression_data_out.close();
		
		ofstream DC_out("Analyse/decompte.txt",ios::out);
		DC_out<<"time epsq qop qopss_	qopsv_ ass_ asv_ afnss_ afnsv_ aftss_ aftsv_ alnss_ alnsv_ altss_ altsv_ a2ss_ a2sv_ afn2ss_ afn2sv_ aft2ss_ aft2sv_ al2ss_ al2sv_"<<endl;
		DC_out.close();

	}
	if( calcPtheta)
	{
		system("mkdir Analyse/angDistrib");
	}
	
	if (calcz)
	{
		system("mkdir Analyse/Connect");
		ofstream outz("Analyse/connect/Z.txt",ios::out); outz.close();
		ofstream out("Analyse/connect/ZPs.txt",ios::out); out.close();
		ofstream out1("Analyse/connect/ZPd.txt",ios::out); out1.close();
		ofstream out3("Analyse/connect/ZP.txt",ios::out); out3.close();
		ofstream out4("Analyse/connect/ZP2.txt",ios::out); out4.close();		
		
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
			
	if (displaySample)
	{
		system("mkdir Analyse/PS1");
		system("mkdir Analyse/PS2");
	}
	
	if (calfracture)
	{
		ofstream fra("Analyse/fracture.txt",ios::out);
		fra<<"#	time	Nfissure	Nfragmentation	Ncasse"<<endl;
		fra.close();
	}
	
	if (calfgranulometrique)
	{
		ofstream gra("Analyse/granulometrique.txt",ios::out);
		gra.close();
	//	granulometrique(Nbingranulometrique);
	}
	
	if (calclusterdouble)
	{
		system("mkdir Analyse/Cluster");
		ofstream clust("Analyse/Cluster/cluster_number.txt",ios::out);
		clust.close();
		
		ofstream clust_df("Analyse/Cluster/cluster_df.txt",ios::out);
		clust_df.close();

	}
		ofstream paroi("Analyse/paroi.txt",ios::out);
		paroi<<"time	contraint	vitesse_top	vitesse_moy"<<endl;
		paroi.close();
		
}

void Biaxial_A::ContactMesh()

{
    

char fname[100];
char name_spl[50];


// Ecriture variables particules dans un autre fichier pour post traitement fortran
/////////////////////////
system("mkdir ContactMesh ");

sprintf(name_spl,"ContactMesh/spl_%04d.his",Nanalyze());

ofstream datafile__(name_spl);

	for (long int i = 0; i<sys_->spl()->lbody().size(); i++) {
        if (sys_->spl()->body(i)->type() == _type_disk && sys_->spl()->body(i)->bodyDof()==NULL)
	 {
		// sys_->spl()->body(i)->write(datafile__);
		sys_->spl()->body(i)->writeM(datafile__);
	 }
    }

datafile__.close();


sprintf(fname,"ContactMesh/contact%04d.dat",Nanalyze());
//writePS(fname);

cout<<"clist size: "<<sys_->nwk()->clist().size()<<endl;
    
ofstream Mesh_(fname,ios::out);

for(unsigned i=0 ; i<sys_->nwk()->clist().size() ; ++i)
{

		//enlever les contacts avec parois, pas bon!
       // if(sys_->nwk()->inter(sys_->nwk()->clist(i))->type()==0)
    
        {
            
        Mesh_
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->id()     << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->id()    << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()      << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()      << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()     << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()     << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->x()               << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->y()               << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->nx()               << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->ny()               << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()               << " "
        << sys_->nwk()->inter(sys_->nwk()->clist(i))->ft()               <<endl;
        }
    }

Mesh_.close();
    
}



void Biaxial_A::cluster()
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

void Biaxial_A::Ptheta(unsigned int Nbin , unsigned int period, unsigned int width) 
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
	bool active;
	unsigned int NctotN=0,rangN,NctotT=0;
	inter2d * interc;
	for( unsigned int i=0;i<sys_->nwk()->clist().size();++i)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(i));
		active=false;
		
		for ( unsigned int r=0; r<interc->rang(); ++r)
		{
			interc->current()=r;
			if (interc->fn() !=0.) active=true;
		}
		interc->current()=0;
			
		if( active && prb2.contain(interc))
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

			if ( interc->rang()==2)
			{
				interc->current()=1;
				FN[rangN]+=interc->fn();
				fnmoy+=interc->fn();
				interc->current()=0;
			}
			
			fnmoy+=interc->fn();
			lmoy +=l;

			if( interc->ft() != 0. )
			{
				FT [ rangN ]+=interc->ft();
				NcT[ rangN ]++;
				ftmoy += interc->ft();
				NctotT++;
				if ( interc->rang()==2)
				{
					interc->current()=1;
					FT[rangN]+=interc->ft();
					ftmoy+=interc->ft();
					interc->current()=0;
				}
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
	for ( unsigned int i =0;i<Nbin;++i)
	{
		p.add( (.5+i)*amp+M_PI, (double) (NcN[i])/(double) (NctotN)*(double)(Nbin)/M_PI );//Distrib angulaire de prob c'est une densité de fréquence
		fn.add( (.5+i)*amp+M_PI, FN[i]/(double) (NcN[i]) );//Distrib ang de l'intensite de forces normales
		ld.add( (.5+i)*amp+M_PI, L[i]/(double) (NcN[i]) );//Distrib ang de l'intensite de forces tangentielle
		if( NctotT!=0) ft.add( (.5+i)*amp+M_PI, FT[i]/(double) (NcT[i]));
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
	cout<<"SIZE:="<<pm.psetSize()<<endl;
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
void Biaxial_A::Ptheta(unsigned int Nbin , unsigned int period, unsigned int width) 
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

		if( interc->fn() !=0. && prb2.contain(interc))
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
void Biaxial_A::Gap()
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

/*void Biaxial_A::def()
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

void Biaxial_A::def()
{
	/*defx+=( right_->x()-x0)/l0;
	defy+=( top_->y()-y0)/h0;
	
	y0 = top_ ->y();
	x0 = right_->x();
	l0=right_->x() - left_->x() - 2*left_->R();
	h0=top_->y()   - bottom_->y() - 2*left_->R();*/
	// Dilatations et non deformations:
	defx=( right_->x()-left_->x()-2*left_->R())/l0;
	defy=( top_->y()-bottom_->y()-2*left_->R())/h0;

	epsp=defx+defy;
	epsq=max(defx,defy)-min(defx,defy);
	ofstream str("Analyse/strain.txt",ios::app);
	str<<time<<" "<<defx<<" "<<defy<<" "<<epsp<<" "<<epsq<<" "<<epsp/epsq<<" "<<asin(epsp/epsq)<<endl;
	str.close();
}

void Biaxial_A::SF()
{
	//sf()=solidFraction(totalProbe_,*sys_->spl(), *sys_->nwk());
	sf()=solidFraction(prb2,*sys_->spl(), *sys_->nwk());
	
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

void Biaxial_A::polygAnisotropy( unsigned int nc)
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

void Biaxial_A::A()
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

void Biaxial_A::A2()
{
	an2_=at2_=al_=0;
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No Fabric "<<endl;return ;}
	
	cout<<"	Fabric A2 : ";
	gdm::Tensor2x2 * F = Fabric2InProbe(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	
	if (F != NULL ) 
	{
		F->eigenValues();
		a2()=2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );
		da2_=F->majorDirection();
		cout<<" direction " << da2_/M_PI*180.<<" ac2= "<<a2()<<endl;
	}
	delete F;

}

void Biaxial_A::forces_A(int Nbin )
{
	
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	
	cout<<"	fn_A : ";
	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	
	if (Fn!=NULL)//(Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		an() = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
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
		dl_=L->majorDirection();
		al()=al()-a();
		cout<<"      direction "<<dl_/M_PI*180.<<" al= "<<al()<<endl;
	}

	cout<<"	ln_A : ";
	gdm::Tensor2x2 * Ln = lnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	if (Ln != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		Ln->eigenValues();
		aln() = 2.*( max(Ln->l2(),Ln->l1()) - min(Ln->l2(),Ln->l1()) )/(Ln->l1() + Ln->l2() );
		dln_=Ln->majorDirection();
		aln()=aln()-a();
		cout<<"      direction "<<dln_/M_PI*180.<<" aln:= "<<aln()<<endl;
	}

	cout<<"	lt_A : ";
	gdm::Tensor2x2 * Lt = ltAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	if (Lt != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		*Lt= *Lt + *Ln;
		if(Lt->xx()!=0 && Lt->yy() !=0)
		{
			Lt->eigenValues();
			alt() = 2.*( max(Lt->l2(),Lt->l1()) - min(Lt->l2(),Lt->l1()) )/(Ln->l1() + Ln->l2() );
			dlt_=Lt->majorDirection();
			alt()=alt()-aln()-a();
			cout<<"      direction "<<dlt_/M_PI*180.<<" alt:= "<<alt()<<endl;
		}
	}

	delete Fn;
	delete Ft;
	delete L;
	delete Ln;
	delete Lt;
}

void Biaxial_A::forces2_A(int Nbin )
{
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No anisotropy "<<endl;return ;}
	
	cout<<"	fn2_A : ";
	gdm::Tensor2x2 * Fn = fn2AnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);

	if (Fn!=NULL)//(Fn->xx() != 0. && Fn->yy() !=0.) 
	{
		Fn->eigenValues();
		an2() = 2.*( max(Fn->l2(),Fn->l1()) - min(Fn->l2(),Fn->l1()) )/(Fn->l1() + Fn->l2() );
		dn2_=Fn->majorDirection();
		an2()=an2()-a2();
		cout<<"     direction "<<dn2_/M_PI*180.<<" an2= "<<an2()<<endl;
	}
	
	
	cout<<"	ft2_A : ";
	gdm::Tensor2x2 * Ft = ft2AnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),Nbin );
	
	if( Ft != NULL )
	{
		*Ft= *Ft + *Fn;
		

		if (Ft->xx() != 0 && Ft->yy() !=0) 
		{
			Ft->eigenValues();
			at2()=2.*( max(Ft->l2(),Ft->l1()) - min(Ft->l2(),Ft->l1()) )/(Fn->l1() + Fn->l2() );
			dt2_=Ft->majorDirection();
			at2()-=an2()+a2();
			cout<<"     direction "<<dt2_/M_PI*180.<<" at2= "<<at2()<<endl;
		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		at2()=0;
	}
	
	cout<<"	l2_A : ";
	gdm::Tensor2x2 * L = length2AnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	if (L != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		L->eigenValues();
		al2() = 2.*( max(L->l2(),L->l1()) - min(L->l2(),L->l1()) )/(L->l1() + L->l2() );
			//cout<<" max "<<max(Fn->l2(),Fn->l1())<<" "<<min(Fn->l2(),Fn->l1())<<endl;
		dl2_=L->majorDirection();
		al2()=al2()-a2();
		cout<<"      direction "<<dl2_/M_PI*180.<<" al2= "<<al2()<<endl;
	}


	delete Fn;
	delete Ft;
	delete L;
}

void Biaxial_A::globalStress()
{
	
	cout<<"	Globalstress : ";
	if ( sys_->nwk()->clist().empty()) {cout<<" No contact = No Stress "<<endl;return ;}
	
	gdm::Tensor2x2 * S = StressInProbe(prb_, *(sys_)->spl(),*(sys_)->nwk());	
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
		cout<<" direction "<<ds_/M_PI*180.<<" q/p = "<<qop_<<endl;

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
/*
void Biaxial_A::Decompte()
{
	double ssv1, ssv2, qopsv, sss1, sss2, qopss;
	double ass,afnss,aftss,alnss,altss,a2ss,afn2ss,aft2ss,al2ss;
	double asv,afnsv,aftsv,alnsv,altsv,a2sv,afn2sv,aft2sv,al2sv;
	
	//qopsv
	gdm::Tensor2x2 * Ssv = StressInProbe_sv(prb_, *(sys_)->spl(),*(sys_)->nwk());	
	if (Ssv != NULL ) 
	{
		Ssv->eigenValues();
		ssv1=Ssv->l1();
		ssv2=Ssv->l2();
		qopsv=( max(ssv1,ssv2)-min(ssv1,ssv2) )/pressure_;
	}
	else
	{
		qopsv=0.;
	}
	//qopss
	gdm::Tensor2x2 * Sss = StressInProbe_ss(prb_, *(sys_)->spl(),*(sys_)->nwk());	
	if (Sss != NULL ) 
	{
		Sss->eigenValues();
		sss1=Sss->l1();
		sss2=Sss->l2();
		qopss=( max(sss1,sss2)-min(sss1,sss2) )/pressure_;
	}
	else
	{
		qopss=0.;
	}
	
	afnss=aftss=alnss=altss=0;
	gdm::Tensor2x2 * Fss = FabricInProbe_ss(prb_, *(sys_)->spl(),*(sys_)->nwk() );	
	if (Fss != NULL ) 
	{
		Fss->eigenValues();
		ass=2.*(max(Fss->l2(),Fss->l1())- min(Fss->l2(),Fss->l1()));
	}
	else
	{
		ass=0.;
	}
	delete Fss;

	afnsv=aftsv=alnsv=altsv=0;
	gdm::Tensor2x2 * Fsv = FabricInProbe_sv(prb_, *(sys_)->spl(),*(sys_)->nwk() );	
	if (Fsv != NULL ) 
	{
		Fsv->eigenValues();
		asv=2.*(max(Fsv->l2(),Fsv->l1())- min(Fsv->l2(),Fsv->l1()) );
	}
	else
	{
		asv=0.;
	}
	delete Fsv;

	gdm::Tensor2x2 * Fn = fnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);	
	gdm::Tensor2x2 * Fnss = fnAnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Fnsv = fnAnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);

	if (Fn!=NULL)
	{
		Fn->eigenValues();
		Fnss->eigenValues();
		Fnsv->eigenValues();
		
		afnss = 2.*( max(Fnss->l2(),Fnss->l1()) - min(Fnss->l2(),Fnss->l1()) )/(Fn->l1() + Fn->l2() );
		afnsv = 2.*( max(Fnsv->l2(),Fnsv->l1()) - min(Fnsv->l2(),Fnsv->l1()) )/(Fn->l1() + Fn->l2() );
		afnss=afnss-ass;
		afnsv=afnsv-asv;
	}

	gdm::Tensor2x2 * Ft = ftAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),30 );
	gdm::Tensor2x2 * Ftss = ftAnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),30 );
	gdm::Tensor2x2 * Ftsv = ftAnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),30 );

	if( Ft != NULL )
	{
		*Ft= *Ft + *Fn;
		*Ftss= *Ftss + *Fnss;
		*Ftsv= *Ftsv + *Fnsv;
		
		if (Ft->xx() != 0 && Ft->yy() !=0) 
		{
			Ft->eigenValues();
			Ftss->eigenValues();
			Ftsv->eigenValues();
			aftss=2.*(max(Ftss->l2(),Ftss->l1()) - min(Ftss->l2(),Ftss->l1()))/(Fn->l1() + Fn->l2());
			aftss-=afnss+ass;
			aftsv=2.*(max(Ftsv->l2(),Ftsv->l1()) - min(Ftsv->l2(),Ftsv->l1()))/(Fn->l1() + Fn->l2());
			aftsv-=afnsv+asv;

		}
	}
	else 
	{
		cout<<"     No tangential forces "<<endl;
		aftss=0;
		aftsv=0;
	}
	
	gdm::Tensor2x2 * Ln = lnAnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Lnss = lnAnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Lnsv = lnAnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);

	if (Ln != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		Ln->eigenValues();
		Lnss->eigenValues();
		Lnsv->eigenValues();
		alnss = 2.*( max(Lnss->l2(),Lnss->l1()) - min(Lnss->l2(),Lnss->l1()) )/(Ln->l1() + Ln->l2() );
		alnsv = 2.*( max(Lnsv->l2(),Lnsv->l1()) - min(Lnsv->l2(),Lnsv->l1()) )/(Ln->l1() + Ln->l2() );
		alnss=alnss-ass;
		alnsv=alnsv-asv;
	}

	gdm::Tensor2x2 * Ltss = ltAnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Ltsv = ltAnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	
	if (Ln != NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		*Ltss= *Ltss + *Lnss;
		*Ltsv= *Ltsv + *Lnsv;
		
		if(Ln->xx()!=0 && Ln->yy() !=0)
		{
			Ltss->eigenValues();
			Ltsv->eigenValues();
			altss = 2.*( max(Ltss->l2(),Ltss->l1()) - min(Ltss->l2(),Ltss->l1()) )/(Ln->l1() + Ln->l2());
			altsv = 2.*( max(Ltsv->l2(),Ltsv->l1()) - min(Ltsv->l2(),Ltsv->l1()) )/(Ln->l1() + Ln->l2());
			altss=altss-alnss-ass;
			altsv=altsv-alnsv-asv;
		}
	}

	delete Fn;
	delete Fnss;
	delete Fnsv;
	delete Ft;
	delete Ftss;
	delete Ftsv;
	delete Ln;
	delete Lnss;
	delete Lnsv;
	delete Ltss;
	delete Ltsv;

	afn2ss=afn2sv=aft2ss=aft2sv=al2ss=al2sv=0;
	a2ss=0.;
	a2sv=0.;

	gdm::Tensor2x2 * F2ss = Fabric2InProbe_ss(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	
	if (F2ss != NULL ) 
	{
		F2ss->eigenValues();
		a2ss=2.*(max(F2ss->l2(),F2ss->l1())- min(F2ss->l2(),F2ss->l1()) );
	}
	else
	{
		a2ss=0.;
	}
	delete F2ss;
	cout<<"a2ss:="<<a2ss;
	gdm::Tensor2x2 * F2sv = Fabric2InProbe_sv(prb_, *(sys_)->spl(),*(sys_)->nwk() );
	if (F2sv != NULL ) 
	{
		F2sv->eigenValues();
		a2sv=2.*(max(F2sv->l2(),F2sv->l1())- min(F2sv->l2(),F2sv->l1()) );
	}
	delete F2sv;
	cout<<"	a2sv:="<<a2sv<<endl;
	
	gdm::Tensor2x2 * Fn2 = fn2AnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Fn2ss = fn2AnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);
	gdm::Tensor2x2 * Fn2sv = fn2AnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);

	if (Fn2!=NULL)
	{
		Fn2->eigenValues();
		Fn2ss->eigenValues();
		Fn2sv->eigenValues();
		afn2ss = 2.*( max(Fn2ss->l2(),Fn2ss->l1()) - min(Fn2ss->l2(),Fn2ss->l1()) )/(Fn2->l1() + Fn2->l2());
		afn2sv = 2.*( max(Fn2sv->l2(),Fn2sv->l1()) - min(Fn2sv->l2(),Fn2sv->l1()) )/(Fn2->l1() + Fn2->l2());
		afn2ss=afn2ss-a2ss;
		afn2sv=afn2sv-a2sv;
	}
	
	
	gdm::Tensor2x2 * Ft2ss = ft2AnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),30 );
	gdm::Tensor2x2 * Ft2sv = ft2AnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(),30 );
	
	if( Fn2 != NULL )
	{
		*Ft2ss= *Ft2ss + *Fn2ss;
		*Ft2sv= *Ft2sv + *Fn2sv;

		if (Fn2->xx() != 0 && Fn2->yy() !=0) 
		{
			Ft2ss->eigenValues();
			Ft2sv->eigenValues();
			aft2ss=2.*( max(Ft2ss->l2(),Ft2ss->l1()) - min(Ft2ss->l2(),Ft2ss->l1()) )/(Fn2->l1() + Fn2->l2() );
			aft2sv=2.*( max(Ft2sv->l2(),Ft2sv->l1()) - min(Ft2sv->l2(),Ft2sv->l1()) )/(Fn2->l1() + Fn2->l2() );
			aft2ss-=afn2ss+a2ss;
			aft2sv-=afn2sv+a2sv;
		}
	}
	else 
	{
		aft2ss=0;
		aft2sv=0;
	}

	gdm::Tensor2x2 * L2 = length2AnisoInProbe(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);		
	gdm::Tensor2x2 * L2ss = length2AnisoInProbe_ss(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);	
	gdm::Tensor2x2 * L2sv = length2AnisoInProbe_sv(totalProbe_, *(sys_)->spl(),*(sys_)->nwk(), 30);	

	if (L2!= NULL)//->xx() != 0. && L->yy() !=0.) 
	{
		L2->eigenValues();
		L2ss->eigenValues();
		L2sv->eigenValues();

		al2ss = 2.*( max(L2ss->l2(),L2ss->l1()) - min(L2ss->l2(),L2ss->l1()) )/(L2->l1() + L2->l2());
		al2sv = 2.*( max(L2sv->l2(),L2sv->l1()) - min(L2sv->l2(),L2sv->l1()) )/(L2->l1() + L2->l2());
		al2ss=al2ss-a2ss;
		al2sv=al2sv-a2sv;
	}

	delete Fn2;
	delete Fn2ss;
	delete Fn2sv;
	delete Ft2ss;
	delete Ft2sv;
	delete L2;
	delete L2ss;
	delete L2sv;

	ofstream DC_out("Analyse/decompte.txt",ios::app);
  //DC_out<<"time epsq qop qopsv_	qopss_ ass_ asv_ afnss_ afnsv_ aftss_ aftsv_ alnss_ alnsv_ altss_ altsv_ a2ss_ a2sv_ afn2ss_ afn2sv_ aft2ss_ aft2sv_ al2ss_ al2sv_"<<endl;
	DC_out<<time<<"	"<<epsq<<"	"<<qop_<<"	"<<qopsv<<"	"<<qopss<<"	"<<ass<<"	"<<asv<<"	"<<afnss<<"	"<<afnsv<<"	"<<aftss<<"	"<<aftsv<<"	"<<alnss<<"	"<<alnsv<<"	"<<altss<<"	"<<altsv<<"	";
	DC_out<<a2ss<<"	"<<a2sv<<"	"<<afn2ss<<"	"<<afn2sv<<"	"<<aft2ss<<"	"<<aft2sv<<"	"<<al2ss<<"	"<<al2sv<<endl;
	DC_out.close();
}
*/

unsigned int Biaxial_A::Z(bool exportZP, unsigned int Nbin, bool cvd)
{
	if (sys_->spl()->body(4)->type()==!0)
	{
	unsigned int i,ni;
	unsigned int Nb = sys_->spl()->lbody().size();//Number of particles
	unsigned int Nc = sys_->nwk()->clist().size();//Number of contact
	unsigned int Nz = 0; //Number of particles without dof with contacts ( more to one)
	unsigned int Nf = 0; //Number of "floatting" particles ( less than two contacts )

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();

	unsigned int Ncwf = 0;//Number of contacts with force 
	unsigned int Ncf0 =0; // Number of contats in probe with null force
	vector < unsigned int> Ncpp(Nb,0);//Number of contact per particles
	vector < double > Force(Nb,0.);//Force total per particles
	vector < bool > Activate(Nb,false);
	vector< unsigned int > Nw(Nb,0); 

	unsigned int NcZ=0;
	unsigned int id1,id2;
	double fmoy=0;
	double force=0.;
	double fn1, fn2;

	vector < unsigned int> Ncpps(Nb,0);//Number of contact simple per particles
	vector < unsigned int> Ncppd(Nb,0);//Number of contact double per particles

	unsigned int NcppMaxs,NcppMaxd,NcppMax,NcppMax2;
	unsigned int NcZs=0;
	unsigned int NcZd=0;
	cout<<"Nc:="<<Nc<<endl;	
	
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		if((prb2.containCenter( sys_->nwk()->inter(ni)->first()) || prb2.containCenter( sys_->nwk()->inter(ni)->second())))
		{
			id1=sys_->nwk()->inter(ni)->first()->id();
			id2=sys_->nwk()->inter(ni)->second()->id();
			
			if(sys_->nwk()->inter(ni)->rang()==1)
			{
				if (sys_->nwk()->inter(ni)->fn()!=0.)
				{
					++Ncwf;
					fmoy+=sys_->nwk()->inter(ni)->fn();
					if(prb2.containCenter( sys_->nwk()->inter(ni)->first())) {Ncpps[id1]+=1;Ncpp[id1]+=1;Force[id1]+=sys_->nwk()->inter(ni)->fn();}		
					if(prb2.containCenter( sys_->nwk()->inter(ni)->second())) {Ncpps[id2]+=1;Ncpp[id2]+=1;Force[id2]+=sys_->nwk()->inter(ni)->fn();}
				}
				else {++Ncf0;}
			}
			else if(sys_->nwk()->inter(ni)->rang()==2)
			{
				sys_->nwk()->inter(ni)->current()=0;
				fn1=sys_->nwk()->inter(ni)->fn();
				sys_->nwk()->inter(ni)->current()=1;
				fn2=sys_->nwk()->inter(ni)->fn();
				sys_->nwk()->inter(ni)->current()=0;
				if(fn1*fn2!=0.)
				{
					++Ncwf;
					fmoy+=fn1+fn2;
					if(prb2.containCenter( sys_->nwk()->inter(ni)->first())) {Ncppd[id1]+=1;Ncpp[id1]+=1;Force[id1]+=fn1+fn2;}				
					if(prb2.containCenter( sys_->nwk()->inter(ni)->second())) {Ncppd[id2]+=1;Ncpp[id2]+=1;Force[id2]+=fn1+fn2;}
				}	
				else if ((fn1!=0.)||(fn2!=0.))
				{
					++Ncwf;
					fmoy+=fn1+fn2;
					if(prb2.containCenter( sys_->nwk()->inter(ni)->first())) {Ncpps[id1]+=1;Ncpp[id1]+=1;Force[id1]+=fn1+fn2;}
					if(prb2.containCenter( sys_->nwk()->inter(ni)->second())) {Ncpps[id2]+=1;Ncpp[id2]+=1;Force[id2]+=fn1+fn2;}
				}
				else { ++Ncf0;}
			}
		}
	}

	NcppMaxs=0;
	NcppMaxd=0;
	NcppMax=0;
	NcppMax2=0;
	Nz=0;
	NcZ=0;
	double r;
	pointSet gcu,gc,nu,vu;
	double masstu=0,masst=0;

	rattlersB.clear();
	rattlersVisu.clear();
	rattlersVisu.resize(Nb);

	for(i=4;i<Nb;++i)
	{ 
		if (prb2.containCenter(sys_->spl()->body(i) ) )
		{
			r=sys_->spl()->body(i)->sizeVerlet();
			gc.add(r , r*r);
			masst +=r*r;

			if (Ncpp[i]>1)  // occurs in contact list and probe and with contact number per particles larger than 1
			{	
				NcppMaxs=max(Ncpps[i],NcppMaxs);
				NcppMaxd=max(Ncppd[i],NcppMaxd);
				NcppMax=max(Ncpp[i],NcppMax);
				NcppMax2=max(Ncpps[i]+2*Ncppd[i],NcppMax2);

				Nz++;
				NcZs += Ncpps[i];
				NcZd += Ncppd[i];
				NcZ += Ncpp[i];
				Activate[i]=true;
				gcu.add( r, r*r);
				nu.add( r, 1.);
				masstu +=r*r;
				rattlersVisu[i]=false;
			}	
			else
			{
				++Nf;
				nu.add(r, 0.);
				rattlersB.push_back(sys_->spl()->body(i) );
				rattlersVisu[i]=true;
			}
		}
	}
	double zs,zd;
	ratlers_=(double)(Nf)/(double)(Nf+Nz);//Ty le so phan tu co it hon hoac bang 1 contact
	z()=(double) (NcZ)/(double) (Nz);//So contact trung binh cua moi hat (chi tinh voi cac hat co nhieu hon 1 contact)
	zs=(double) (NcZs)/(double) (Nz);//Nombre contact simple moyenne->Number of connectivity simple
	zd=(double) (NcZd)/(double) (Nz);//Nombre contact double moyenne->Number of connectivity double
	cout<<"Nombre contact:="<<NcZs+NcZd<<"	Nombre particule active:="<<Nz<<"	Nombre de coordination:="<<zs+zd<<"	Nombre de connectivity:="<<zs+2*zd<<" ratlers_= "<<ratlers_<<endl;

	//Partial connectivity
	double hight=2*prb_.hh();
	double large=2*prb_.hl();
	ofstream outz("Analyse/connect/Z.txt",ios::app);
	outz<<time<<" "<<epsq<<" "<<ratlers_<<" "<<compactness(hight,large)<<"	"<<zs<<"	"<<zd<<" "<<zs+zd<<" "<<zs+2*zd<<"	"<<zs/(zs+zd)<<" "<<zd/(zs+zd)<<endl;
	outz.close();
	
	//Partial connectivity
	if ( exportZP)
	{
		vector <unsigned int> Npncs(NcppMaxs+1,0);//Number of particles with n contacts simple
		vector <unsigned int> Npncd(NcppMaxd+1,0);//Number of particles with n contacts double
		vector <unsigned int> Npnc(NcppMax+1,0);//Number of particles with n contacts
		vector <unsigned int> Npnc2(NcppMax2+1,0);//Number of particles with n contacts- contact double compte two time
		
		for(i=4;i<Nb;++i)
		{ 
			if ( Activate[i]==true)
			{
				Npncs[Ncpps[i]]+=1;
				Npncd[Ncppd[i]]+=1;
				Npnc[Ncpp[i]]+=1;
				Npnc2[Ncpps[i]+2*Ncppd[i]]+=1;
			}
		}

		ofstream out("Analyse/connect/ZPs.txt",ios::app);
		vector <double> Zps(NcppMaxs+1);
		for(i=0;i<NcppMaxs+1;++i)
		{ 
			Zps[i]=(double) (Npncs[i])/(double) (Nz);

			out<<i<<"	"<<Zps[i]<<endl;
		}
		out.close();
		
		ofstream out1("Analyse/connect/ZPd.txt",ios::app);
		vector <double> Zpd(NcppMaxd+1);
		for(i=0;i<NcppMaxd+1;++i)
		{ 
			Zpd[i]=(double) (Npncd[i])/(double) (Nz);

			out1<<i<<"	"<<Zpd[i]<<endl;
		}
		out1.close();
		
		ofstream out2("Analyse/connect/ZPt.txt",ios::app);
		ofstream out3("Analyse/connect/ZP.txt",ios::app);
		out2<<time<<" ";

		vector <double> Zp(NcppMax+1);
		for(i=0;i<NcppMax+1;++i)
		{ 
			Zp[i]=(double) (Npnc[i])/(double) (Nz);

			out2<<Zp[i]<<" ";
			out3<<i<<" "<<Zp[i]<<endl;
		}
		out2<<endl;
		out2.close();
		out3.close();

		ofstream out4("Analyse/connect/ZP2.txt",ios::app);
		vector <double> Zp2(NcppMax2+1);
		for(i=0;i<NcppMax2+1;++i)
		{ 
			Zp2[i]=(double) (Npnc2[i])/(double) (Nz);
			out4<<i<<"	"<<Zp2[i]<<endl;
		}
		out4.close();

		Npncs.clear();
		Zps.clear();
		Npncd.clear();
		Zpd.clear();
		Npnc.clear();
		Zp.clear();
	}
	Ncpps.clear();
	Ncppd.clear();
	Ncpp.clear();

	if (Ncwf==0) 
		{ cout<<" no forces "<<endl;return 0;}
	fmoy/=(double) (Ncwf);
	unsigned int Ncw=0;
	cout<<" fmoy = "<<fmoy<< " ";
		
	ofstream scalZ("Analyse/Connect/scalarZ.txt",ios::app);
	scalZ<<time<<" "<<Ncwf<<" "<<Ncf0<<endl;
	
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);
		if( prb2.contain( sys_->nwk()->inter(ni)))
		{
			id1=sys_->nwk()->inter(ni)->first()->id();
			id2=sys_->nwk()->inter(ni)->second()->id();
			
			force=0.;
			for(unsigned r=0; r<sys_->nwk()->inter(ni)->rang(); ++r)
			{
				sys_->nwk()->inter(ni)->current()=r;
				force+=sys_->nwk()->inter(ni)->fn();
			}
			sys_->nwk()->inter(ni)->current()=0;
			
			if ( Ncpp[id1]>1)
			{
				if( force < fmoy)
				{
					Nw[id1]+=1;
					Ncw++;
				}
			}

			if ( Ncpp[id2]>1)
			{
				if( force < fmoy)
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
	pointSet weak= fracw.slidingHisto(50,.10);
	weak.xoffSet(-rmin);
	weak.xNormalise( rmax-rmin);
	weak.write("Analyse/connect/fracweak.txt");
	
	if (cvd )
	{
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
	
	//Calcul force sur un particule.


	//Connectivity as a function of granulometry
	if( calczg)
	{
		cout<<" / Zg ";
		vector < DataSet*> DZ(Nbin,NULL);
		vector < DataSet*> R(Nbin,NULL);
		vector < DataSet*> DF(Nbin,NULL);

		for(  i=0;i<Nbin;++i)
		{
			DZ[i]=new DataSet;
			R[i]=new DataSet;
			DF[i]=new DataSet;
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
				data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<Ncpp[i] <<endl;

				ZD[Ncpp[i]]->add( r );

				if(prb2.containCenter(sys_->spl()->body(i)))
				{
					jc = (unsigned int) floor ( (r-rmin)/amp);

					if ( fabs(r - rmax)< 1e-10 )  jc=Nbin-1;
					if ( fabs(r - rmin)< 1e-10 )  jc=0;

					DZ[jc]->add( Ncpp[i]);
					R [jc]->add( r );
					PDZ.add( r,Ncpp[i]);
					DF[jc]->add( Force[i]);
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
		double fmoy=0.;
		for(i=0;i<Nbin;++i)
		{
			DF[i]->extractValues();
			fmoy+=DF[i]->mean();
		}
		fmoy/=Nbin;
		for(i=0;i<Nbin;++i)
		{
			DZ[i]->extractValues();
			R[i]->extractValues();
			data2<<(R[i]->mean()-rmin)/(rmax-rmin)<<" "<<R[i]->mean()<<" "<<DZ[i]->mean()<<" "<<DF[i]->mean()<<" "<<DF[i]->mean()/fmoy<<" "<<R[i]->variance()<<" "<<DZ[i]->variance()<<endl;
		}
		data2.close();
	}

	cout<<endl;
	}
	else // pour les disque
	{
	unsigned int i,ni;
	unsigned int Nb = sys_->spl()->lbody().size();//Number of particles
	unsigned int Nc = sys_->nwk()->clist().size();//Number of contact
	unsigned int Nz = 0; //Number of particles without dof with contacts ( more to one)
	unsigned int Nf = 0; //Number of "floatting" particles ( less than two contacts )

	double rmin=sys_->spl()->rmin();
	double rmax=sys_->spl()->rmax();

	unsigned int Ncwf = 0;//Number of contacts with force 
	unsigned int Ncf0 =0; // Number of contats in probe with null force
	vector < unsigned int> Ncpp(Nb,0);//Number of contact per particles
	vector < double > Force(Nb,0.);//Force total per particles
	vector < bool > Activate(Nb,false);
	vector< unsigned int > Nw(Nb,0); 

	unsigned int NcZ=0;
	unsigned int id1,id2;
	double fmoy=0;
	double force=0.;

	unsigned int NcppMax;
	cout<<"Nc:="<<Nc<<endl;	
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);

		if((prb2.containCenter( sys_->nwk()->inter(ni)->first()) || prb2.containCenter( sys_->nwk()->inter(ni)->second())))
		{
			id1=sys_->nwk()->inter(ni)->first()->id();
			id2=sys_->nwk()->inter(ni)->second()->id();

			if (sys_->nwk()->inter(ni)->fn()!=0.)
			{
				++Ncwf;
				fmoy+=sys_->nwk()->inter(ni)->fn();
				if(prb2.containCenter( sys_->nwk()->inter(ni)->first())) {Ncpp[id1]+=1;Force[id1]+=sys_->nwk()->inter(ni)->fn();}		
				if(prb2.containCenter( sys_->nwk()->inter(ni)->second())) {Ncpp[id2]+=1;Force[id2]+=sys_->nwk()->inter(ni)->fn();}
			}
			else {++Ncf0;}
		}
	}

	NcppMax=0;
	Nz=0;
	double r;
	pointSet gcu,gc,nu,vu;
	double masstu=0,masst=0;

	rattlersB.clear();
	rattlersVisu.clear();
	rattlersVisu.resize(Nb);

	for(i=4;i<Nb;++i)
	{ 
		if (prb2.containCenter(sys_->spl()->body(i) ) )
		{
			r=sys_->spl()->body(i)->sizeVerlet();
			gc.add(r , r*r);
			masst +=r*r;

			if (Ncpp[i]>1)  // occurs in contact list and probe and with contact number per particles larger than 1
			{	
				NcppMax=max(Ncpp[i],NcppMax);
				Nz++;
				NcZ += Ncpp[i];
				Activate[i]=true;
				gcu.add( r, r*r);
				nu.add( r, 1.);
				masstu +=r*r;
				rattlersVisu[i]=false;
			}	
			else
			{
				++Nf;
				nu.add(r, 0.);
				rattlersB.push_back(sys_->spl()->body(i) );
				rattlersVisu[i]=true;
			}
		}
	}
	ratlers_=(double)(Nf)/(double)(Nf+Nz);//Ty le so phan tu co it hon hoac bang 1 contact
	z()=(double) (NcZ)/(double) (Nz);//So contact trung binh cua moi hat (chi tinh voi cac hat co nhieu hon 1 contact)
	cout<<"Nombre contact:="<<NcZ<<"	Nombre particule active:="<<Nz<<"	Nombre de coordination:="<<z()<<" ratlers_= "<<ratlers_<<endl;

	//Partial connectivity
	ofstream outz("Analyse/connect/Z.txt",ios::out);
	outz<<time<<" "<<z()<<"	"<<ratlers_<<endl;
	outz.close();
	
	//Partial connectivity
	if ( exportZP)
	{
		vector <unsigned int> Npnc(NcppMax+1,0);//Number of particles with n contacts
		
		for(i=4;i<Nb;++i)
		{ 
			if ( Activate[i]==true) Npnc[Ncpp[i]]+=1;
		}
		
		ofstream out2("Analyse/connect/ZPt.txt",ios::app);
		ofstream out3("Analyse/connect/ZP.txt",ios::out);
		out2<<time<<" ";

		vector <double> Zp(NcppMax+1);
		for(i=0;i<NcppMax+1;++i)
		{ 
			Zp[i]=(double) (Npnc[i])/(double) (Nz);

			out2<<Zp[i]<<" ";
			out3<<i<<" "<<Zp[i]<<endl;
		}
		out2<<endl;
		out2.close();
		out3.close();

		Npnc.clear();
		Zp.clear();
	}
	Ncpp.clear();

	if (Ncwf==0) 
		{ cout<<" no forces "<<endl;return 0;}
	fmoy/=(double) (Ncwf);
	unsigned int Ncw=0;
	cout<<" fmoy = "<<fmoy<< " ";
		
	ofstream scalZ("Analyse/Connect/scalarZ.txt",ios::app);
	scalZ<<time<<" "<<Ncwf<<" "<<Ncf0<<endl;
	
	for ( i=0;i< Nc;++i)
	{
		ni = sys_->nwk()->clist(i);
		if( prb2.contain( sys_->nwk()->inter(ni)))
		{
			id1=sys_->nwk()->inter(ni)->first()->id();
			id2=sys_->nwk()->inter(ni)->second()->id();
			
			force=0.;
			for(unsigned r=0; r<sys_->nwk()->inter(ni)->rang(); ++r)
			{
				sys_->nwk()->inter(ni)->current()=r;
				force+=sys_->nwk()->inter(ni)->fn();
			}
			sys_->nwk()->inter(ni)->current()=0;
			
			if ( Ncpp[id1]>1)
			{
				if( force < fmoy)
				{
					Nw[id1]+=1;
					Ncw++;
				}
			}

			if ( Ncpp[id2]>1)
			{
				if( force < fmoy)
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
	pointSet weak= fracw.slidingHisto(50,.10);
	weak.xoffSet(-rmin);
	weak.xNormalise( rmax-rmin);
	weak.write("Analyse/connect/fracweak.txt");
	
	if (cvd )
	{
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
	
	//Calcul force sur un particule.


	//Connectivity as a function of granulometry
	if( calczg)
	{
		cout<<" / Zg ";
		vector < DataSet*> DZ(Nbin,NULL);
		vector < DataSet*> R(Nbin,NULL);
		vector < DataSet*> DF(Nbin,NULL);

		for(  i=0;i<Nbin;++i)
		{
			DZ[i]=new DataSet;
			R[i]=new DataSet;
			DF[i]=new DataSet;
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
				data<<sys_->spl()->body(i)->id()<<" "<<r<<" "<<Ncpp[i] <<endl;

				ZD[Ncpp[i]]->add( r );

				if(prb2.containCenter(sys_->spl()->body(i)))
				{
					jc = (unsigned int) floor ( (r-rmin)/amp);

					if ( fabs(r - rmax)< 1e-10 )  jc=Nbin-1;
					if ( fabs(r - rmin)< 1e-10 )  jc=0;

					DZ[jc]->add( Ncpp[i]);
					R [jc]->add( r );
					PDZ.add( r,Ncpp[i]);
					DF[jc]->add( Force[i]);
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
		double fmoy=0.;
		for(i=0;i<Nbin;++i)
		{
			DF[i]->extractValues();
			fmoy+=DF[i]->mean();
		}
		fmoy/=Nbin;
		for(i=0;i<Nbin;++i)
		{
			DZ[i]->extractValues();
			R[i]->extractValues();
			data2<<(R[i]->mean()-rmin)/(rmax-rmin)<<" "<<R[i]->mean()<<" "<<DZ[i]->mean()<<" "<<DF[i]->mean()<<" "<<DF[i]->mean()/fmoy<<" "<<R[i]->variance()<<" "<<DZ[i]->variance()<<endl;
		}
		data2.close();
	}

	cout<<endl;
	
	}
	return 1;
}

void Biaxial_A::granuloSpeed(unsigned int Nbin)
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

unsigned int Biaxial_A::granuloStress(unsigned int Nbinstress )
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

unsigned int Biaxial_A::granuloStress2(unsigned int Nc ,unsigned int period, unsigned int width)
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

unsigned int Biaxial_A::granuloStress3(unsigned int Nbinstress )
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

int Biaxial_A::pdfforce(bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
{
	char fichier[100];
	cout<<"	PDF  : " ;
	if (sys_->spl()->body(4)->type()==!0)
	{	
	DataSet f,fs,fd;
	double force, force2;
	for( unsigned int i=0;i< sys_->nwk()->linter().size();++i)
	{
		if( prb2.contain( sys_->nwk()->inter(i) ))
		{
			force=0.;
			force2=1.0;
			if( fn )
			{
				for (unsigned int r=0;r<sys_->nwk()->inter(i)->rang();r++)
				{
					sys_->nwk()->inter(i)->current()=r;
					force+=sys_->nwk()->inter(i)->fn();
					force2*=sys_->nwk()->inter(i)->fn();
				}
				sys_->nwk()->inter(i)->current()=0;

				if(force!= 0.) 
				{
					f.add(force);//stockage dans f des forces de contacts à traiter
				}
				
				if((sys_->nwk()->inter(i)->rang()==1) && (sys_->nwk()->inter(i)->fn()!=0.)) fs.add(sys_->nwk()->inter(i)->fn());
				if((sys_->nwk()->inter(i)->rang()==2) && (force2!=0.)) fd.add(force);
				else fs.add(force);
			}
			else
			{
				for (unsigned int r=0;r<sys_->nwk()->inter(i)->rang();r++)
				{
					sys_->nwk()->inter(i)->current()=r;
					force+=sys_->nwk()->inter(i)->ft();
					force2*=sys_->nwk()->inter(i)->ft();
				}
				sys_->nwk()->inter(i)->current()=0;

				if(force!= 0.) 
				{
					f.add(fabs(force));//stockage dans f des forces de contacts à traiter
				}
				if((sys_->nwk()->inter(i)->rang()==1) && (sys_->nwk()->inter(i)->ft()!=0.)) fs.add(sys_->nwk()->inter(i)->ft());
				if((sys_->nwk()->inter(i)->rang()==2) && (force2!=0.)) fd.add(force);
				else fs.add(force);
			}
		}
	}
	if (f.setSize() <100)
		{cout<<" no more forces (min 100 events ) : "<<f.setSize()<<endl; return 0;}

	f.extractValues();
	fs.extractValues();
	fd.extractValues();

	if( normalize )//considération des forces normalisées 
	{
		f.Normalize( f.mean() );
		fs.Normalize( fs.mean() );
		fd.Normalize( fd.mean() );
	}
	cout<<" fmoy = "<<f.mean()<<" fsmoy = "<<fs.mean()<<" fdmoy = "<<fd.mean()<<endl;
	
	f.DecreasingSort();
	fs.DecreasingSort();
	fd.DecreasingSort();
	pointSet pdf = f.kernelPdf( nbin,.01);
	pointSet pdfs = fs.kernelPdf( nbin,.01);
	pointSet pdfd = fd.kernelPdf( nbin,.01);
	
	if(fn)
		sprintf(fichier,"Analyse/pdf/pdffn.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft.txt");
	
	pdf.write(fichier);

	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffn_mob.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft_mob.txt");
	pdf.periodic()=false;
	pointSet pdfm = pdf.mobileMean2(period, width);
	pdfm.write(fichier);

	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffns_mob.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdffts_mob.txt");
	pdfs.periodic()=false;
	pointSet pdfsm = pdfs.mobileMean2(period, width);
	pdfsm.write(fichier);

	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffnd_mob.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfftd_mob.txt");
	pdfs.periodic()=false;
	pointSet pdfdm = pdfd.mobileMean2(period, width);
	pdfdm.write(fichier);

	}
	else
	{
	DataSet f;
	double force;
	for( unsigned int i=0;i< sys_->nwk()->linter().size();++i)
	{
		if( prb2.contain( sys_->nwk()->inter(i) ))
		{
			force=0.;
			if( fn )
			{	
				if(sys_->nwk()->inter(i)->fn()!=0.) f.add(sys_->nwk()->inter(i)->fn());
			}
			else
			{
				if(sys_->nwk()->inter(i)->ft()!=0.) f.add(sys_->nwk()->inter(i)->ft());
			}
		}
	}
	if (f.setSize() <100)
		{cout<<" no more forces (min 100 events ) : "<<f.setSize()<<endl; return 0;}

	f.extractValues();

	if( normalize ) f.Normalize( f.mean() );
	cout<<" fmoy = "<<f.mean()<<endl;
	
	f.DecreasingSort();
	pointSet pdf = f.kernelPdf( nbin,.01);
	
	if(fn)
		sprintf(fichier,"Analyse/pdf/pdffn.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft.txt");
	
	pdf.write(fichier);

	if( fn )
		sprintf(fichier,"Analyse/pdf/pdffn_mob.txt");
	else
		sprintf(fichier,"Analyse/pdf/pdfft_mob.txt");
	pdf.periodic()=false;
	pointSet pdfm = pdf.mobileMean2(period, width);
	pdfm.write(fichier);
	}
	return 1;
}

int Biaxial_A::granulopdf(unsigned int Nq,bool fn,int nbin,bool normalize, unsigned int period, unsigned int width) // A coupler avec la classe SET
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

void Biaxial_A::forcesMaxCorrelation( )
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

//Fonction de distribution radiale g
void Biaxial_A::rfd(unsigned int Npoint, unsigned int Nrmean)
{
	double rmean = sys_->spl()->rmoy();
	unsigned int Nbod;
	unsigned int const nvertex=64;

	double dx,dy,d;
	double xi,yi,r,x0,y0,hh,hl;
	double amp=Nrmean*rmean/(double) Npoint;
	double meandensity;
	double section;
	double meancompact=0.;
	
	vector<double> g(Npoint,0.);
	vector<double> compact(Npoint,0.);
	vector< vector<double> > Area;
	vector< vector<double> > Volume;
	vector< vector<unsigned int> > ncenter;
	
	polyg *Prb=new polyg;
		
	Nbod=sys_->spl()->lbody().size()-sys_->lctrl().size();
	meandensity= (double) (Nbod)/prb_.area();

	x0=prb_.x();
	y0=prb_.y();
	hh=prb_.hh();
	hl=prb_.hl();	
		
	gdm::vertex vertex_;
	
	Prb->x()=prb_.x();
	Prb->y()=prb_.y();
	Prb->rot()=0.;
	Prb->Vertex().clear();
	
	vertex_.x()=-prb_.hl();		vertex_.y()=-prb_.hh();		Prb->Vertex().push_back(vertex_);
	vertex_.x()=prb_.hl(); 		vertex_.y()=-prb_.hh();		Prb->Vertex().push_back(vertex_);
	vertex_.x()=prb_.hl(); 		vertex_.y()=prb_.hh();		Prb->Vertex().push_back(vertex_);
	vertex_.x()=-prb_.hl();		vertex_.y()=prb_.hh();		Prb->Vertex().push_back(vertex_);
	Prb->adjustCenter();
	Prb->Fill(1.);
	
	for( unsigned int i=sys_->lctrl().size() ; i< sys_->spl()->lbody().size() ;++i)
	meancompact+=sys_->spl()->body(i)->Area();
	meancompact/=prb_.area();
	
	cout<<"amp:="<<amp<<"	meandensity:="<<meandensity<<"	meancompact:="<<meancompact<<endl;

	vector<double> Area_;
	vector<double> Volume_;
	vector<unsigned int> ncenter_;

	for( unsigned int i=sys_->lctrl().size() ; i< sys_->spl()->lbody().size() ;++i)
	{	
		Area_.clear();
		Volume_.clear();
		ncenter_.clear();
		
		for(unsigned int k=0;k<Npoint;k++)
		{
			r=(k+1)*amp;
			
			//calcule area commun entre prb_ et un circle centre xi,yi et radial r	(Area)
			xi=sys_->spl()->body(i)->x();
			yi=sys_->spl()->body(i)->y();
			
			polyg *Cer=new polyg;
			polyg *poly=new polyg;

			Cer->x()=xi;
			Cer->y()=yi;
			Cer->rot()=0.;
			Cer->Vertex().clear();
			for (unsigned int j=0;j<nvertex;j++)
			{
				vertex_.x()=r*cos(j*2.*M_PI/nvertex);
				vertex_.y()=r*sin(j*2.*M_PI/nvertex);
				Cer->Vertex().push_back(vertex_);
			}
			Cer->adjustCenter();
			Cer->Fill(1.);
			
			if( ((xi-r)>=(x0-hl)) && ((xi+r)<=(x0+hl)) && ((yi-r)>=(y0-hh)) && ((yi+r)<=(y0+hh)) ) section=0.5*nvertex*r*r*sin(2.*M_PI/nvertex);
			else if (Inter(Prb,Cer)) 
			{
				poly=Intersection(Prb,Cer);
				if ((poly->Area()>0)) section=poly->Area();
				else
				{
					cout<<"Calcul rfd	3219:	Aire polygone negatif ! "<<poly->Area()<<endl;
					cout<<"THONG SO CUA HINH TRON"<<endl;
					cout<<"x0:="<<Cer->x()<<"	y0:="<<Cer->y()<<endl;
					for (unsigned int j=0;j<Cer->Vertex().size();j++)
					cout<<"Vertex["<<j<<"].(x):="<<Cer->Vertex(j).x()<<"	Vertex["<<j<<"].(y):="<<Cer->Vertex(j).y()<<endl;
					
					cout<<"THONG SO CUA BAO"<<endl;
					cout<<"x0:="<<Prb->x()<<"	y0:="<<Prb->y()<<endl;
					for (unsigned int j=0;j<Prb->Vertex().size();j++)
					cout<<"Vertex["<<j<<"].(x):="<<Prb->Vertex(j).x()<<"	Vertex["<<j<<"].(y):="<<Prb->Vertex(j).y()<<endl;

					cout<<"THONG SO CUA POLYGONE"<<endl;
					cout<<"x0:="<<poly->x()<<"	y0:="<<poly->y()<<endl;
					for (unsigned int j=0;j<poly->Vertex().size();j++)
					cout<<"Vertex["<<j<<"].(x):="<<poly->Vertex(j).x()<<"	Vertex["<<j<<"].(y):="<<poly->Vertex(j).y()<<endl;
					
					char fname[100];
					sprintf(fname,"Analyse/Section.ps");
					ofstream ps(fname);
					double xoffset=5.;
					double c;
					double s;

					ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
					ps<<"%%BoundingBox: 0 0 595 842"<<endl;
					ps<<"0.5 setlinewidth"<<endl;
					ps<<"0.0 setgray"<<endl;
		
					//Tracer les CERCLE
										
				c=cos(Cer->rot());
				s=sin(Cer->rot());
				
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(0).x()*c-Cer->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(0).x()*s+Cer->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<Cer->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(j).x()*c-Cer->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(j).x()*s+Cer->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}
 			 	ps<<"1 0 0 setrgbcolor"<<endl;
 			//	ps<<"gsave 1 0 0 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.3 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
		
				//Tracer BAO
				c=cos(Prb->rot());
				s=sin(Prb->rot());
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(Prb->x()+Prb->Vertex(0).x()*c-Prb->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(Prb->y()+Prb->Vertex(0).x()*s+Prb->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<Prb->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(Prb->x()+Prb->Vertex(j).x()*c-Prb->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(Prb->y()+Cer->Vertex(j).x()*s+Prb->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}	
 			 	ps<<"0 1 0 setrgbcolor"<<endl;
 				//ps<<"gsave 0 1 0 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.3 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
			
				//Tracer POLYGONE SECTION
				c=cos(poly->rot());
				s=sin(poly->rot());			
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(poly->x()+poly->Vertex(0).x()*c-poly->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(poly->y()+Prb->Vertex(0).x()*s+poly->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<poly->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(poly->x()+poly->Vertex(j).x()*c-poly->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(poly->y()+Cer->Vertex(j).x()*s+poly->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}	
 	 			ps<<"0 0 1 setrgbcolor"<<endl;
 				ps<<"gsave 0 0 1 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.3 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
				exit(0);		
				}
			}
			else {cout<<"Calcul rfd	:	il n'y a pas de section"<<endl;exit(0);}
			Area_.push_back(section);
			
			/*if((i==29)&&(k==3)) 
			{
				char fname[100];
				sprintf(fname,"Analyse/Section.ps");
				ofstream ps(fname);
				double xoffset=5.;
				double c;
				double s;

				ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
				ps<<"%%BoundingBox: 0 0 595 842"<<endl;
				ps<<"0.5 setlinewidth"<<endl;
				ps<<"0.0 setgray"<<endl;
				
				//Tracer BAO
				c=cos(Prb->rot());
				s=sin(Prb->rot());
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(Prb->x()+Prb->Vertex(0).x()*c-Prb->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(Prb->y()+Prb->Vertex(0).x()*s+Prb->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<Prb->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(Prb->x()+Prb->Vertex(j).x()*c-Prb->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(Prb->y()+Cer->Vertex(j).x()*s+Prb->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}	
 			 	ps<<"0 1 0 setrgbcolor"<<endl;
 				//ps<<"gsave 0 1 0 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.3 setlinewidth"<<endl;
				ps<<"stroke"<<endl;

				double r1=k*amp;
				double r2=(k+1)*amp;
				double section1,section2;
				Cer->x()=xi;
				Cer->y()=yi;
				Cer->rot()=0.;
				Cer->Vertex().clear();
				for (unsigned int j=0;j<nvertex;j++)
				{
					vertex_.x()=r1*cos(j*2.*M_PI/nvertex);
					vertex_.y()=r1*sin(j*2.*M_PI/nvertex);
					Cer->Vertex().push_back(vertex_);
				}
				Cer->adjustCenter();
				Cer->Fill(1.);

				if( ((xi-r1)>(x0-hl)) && ((xi+r1)<(x0+hl)) && ((yi-r1)>(y0-hh)) && ((yi+r1)<(y0+hh)) ) 
				{
					section1=0.5*nvertex*r1*r1*sin(2.*M_PI/nvertex);
					cout<<"CHECK HANG 1:="<<(xi-r1)-(x0-hl)<<"	"<< -((xi+r1)-(x0+hl))<<"	"<<((yi-r1)-(y0-hh))<<"	"<<-((yi+r1)-(y0+hh))<<endl;
				}
				else if (Inter(Prb,Cer)) 
				{
					poly=Intersection(Prb,Cer);
					if ((poly->Area()>0)) section1=poly->Area();
					//Tracer POLYGONE SECTION
					c=cos(poly->rot());
					s=sin(poly->rot());			
					ps<<"newpath"<<endl;
					ps<<xoffset+zoom_*(poly->x()+poly->Vertex(0).x()*c-poly->Vertex(0).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(poly->y()+Prb->Vertex(0).x()*s+poly->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
					for(unsigned int j=1;j<poly->Vertex().size();j++)
					{   
						ps<<xoffset+zoom_*(poly->x()+poly->Vertex(j).x()*c-poly->Vertex(j).y()*s-sys_->spl()->xmin());
						ps<<" "<<xoffset+zoom_*(poly->y()+Cer->Vertex(j).x()*s+poly->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 					}	
 	 				ps<<"1 1 1 setrgbcolor"<<endl;
 					//ps<<"gsave 0 0 1 setrgbcolor fill grestore"<<endl;
 					ps<<"closepath"<<endl;
					ps<<"1.0 setlinewidth"<<endl;
					ps<<"stroke"<<endl;				
				}
				
				//Tracer les CERCLE
										
				c=cos(Cer->rot());
				s=sin(Cer->rot());
				
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(0).x()*c-Cer->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(0).x()*s+Cer->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<Cer->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(j).x()*c-Cer->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(j).x()*s+Cer->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}
 			 	ps<<"1 0 0 setrgbcolor"<<endl;
 			//	ps<<"gsave 1 0 0 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"0.1 setlinewidth"<<endl;
				ps<<"stroke"<<endl;
				
				Cer->Vertex().clear();
				for (unsigned int j=0;j<nvertex;j++)
				{
					vertex_.x()=r2*cos(j*2.*M_PI/nvertex);
					vertex_.y()=r2*sin(j*2.*M_PI/nvertex);
					Cer->Vertex().push_back(vertex_);
				}
				Cer->adjustCenter();
				Cer->Fill(1.);
			
				c=cos(Cer->rot());
				s=sin(Cer->rot());
				
				ps<<"newpath"<<endl;
				ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(0).x()*c-Cer->Vertex(0).y()*s-sys_->spl()->xmin());
				ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(0).x()*s+Cer->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
				for(unsigned int j=1;j<Cer->Vertex().size();j++)
				{   
					ps<<xoffset+zoom_*(Cer->x()+Cer->Vertex(j).x()*c-Cer->Vertex(j).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(Cer->y()+Cer->Vertex(j).x()*s+Cer->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 				}
 			 	ps<<"1 0 0 setrgbcolor"<<endl;
 			//	ps<<"gsave 1 0 0 setrgbcolor fill grestore"<<endl;
 				ps<<"closepath"<<endl;
				ps<<"1.5 setlinewidth"<<endl;
				ps<<"stroke"<<endl;		
				
				if( ((xi-r2)>=(x0-hl)) && ((xi+r2)<=(x0+hl)) && ((yi-r2)>=(y0-hh)) && ((yi+r2)<=(y0+hh)) ) 
				{
					section2=0.5*nvertex*r2*r2*sin(2.*M_PI/nvertex);
					cout<<"CHECK HANG 2:="<<(xi-r2)-(x0-hl)<<"	"<< -((xi+r2)-(x0+hl))<<"	"<<((yi-r2)-(y0-hh))<<"	"<<-((yi+r2)-(y0+hh))<<endl;
				}
				else if (Inter(Prb,Cer)) 
				{
					poly=Intersection(Prb,Cer);
					if ((poly->Area()>0)) section2=poly->Area();
					
					//Tracer POLYGONE SECTION
					c=cos(poly->rot());
					s=sin(poly->rot());			
					ps<<"newpath"<<endl;
					ps<<xoffset+zoom_*(poly->x()+poly->Vertex(0).x()*c-poly->Vertex(0).y()*s-sys_->spl()->xmin());
					ps<<" "<<xoffset+zoom_*(poly->y()+Prb->Vertex(0).x()*s+poly->Vertex(0).y()*c-sys_->spl()->ymin())<<" moveto"<<endl;
					for(unsigned int j=1;j<poly->Vertex().size();j++)
					{   
						ps<<xoffset+zoom_*(poly->x()+poly->Vertex(j).x()*c-poly->Vertex(j).y()*s-sys_->spl()->xmin());
						ps<<" "<<xoffset+zoom_*(poly->y()+Cer->Vertex(j).x()*s+poly->Vertex(j).y()*c-sys_->spl()->ymin())<<" lineto"<<endl;
 					}	
 		 			ps<<"0 1 1 setrgbcolor"<<endl;
 					//ps<<"gsave 0 1 1 setrgbcolor fill grestore"<<endl;
 					ps<<"closepath"<<endl;
					ps<<"1.5 setlinewidth"<<endl;
					ps<<"stroke"<<endl;
				}
				cout<<"section1:="<<section1<<"		section2:="<<section2<<"	section1-section2:="<<section1-section2<<"	r1:="<<r1<<"	r2:="<<r2<<endl;
			}*/
			
			unsigned int Ncenter=0;
			for( unsigned int j = sys_->lctrl().size() ; j< sys_->spl()->lbody().size() ;++j)
			{	
				dx= xi - sys_->spl()->body(j)->x();
				dy= yi - sys_->spl()->body(j)->y();;
				d=sqrt(dx*dx+dy*dy);
				if((d>=(r-amp))&&(d<=r)&&(j!=i)) Ncenter+=1;
			}
			ncenter_.push_back(Ncenter);
			
			double volume=0.;
			for( unsigned int j=sys_->lctrl().size() ; j< sys_->spl()->lbody().size() ;++j)
			{	
				dx= xi - sys_->spl()->body(j)->x();
				dy= yi - sys_->spl()->body(j)->y();;
				d=sqrt(dx*dx+dy*dy);
				
				if((d-sys_->spl()->body(j)->sizeVerlet())>r) volume+=0.;
				else if((d+sys_->spl()->body(j)->sizeVerlet())<r) volume+=sys_->spl()->body(j)->Area();
				else
				{
					polyg  * polygon=new polyg;
					polygon=dynamic_cast<polyg*> (sys_->spl()->body(j));
				
					if (Inter(Cer,polygon)) 
					{
						poly=Intersection(Cer,polygon);
						if ((poly->Area()>0)) volume+=poly->Area();
						else 
						{
							cout<<"Aire polygone negatif !:	ligne 3489 Biaxial_A "<<poly->Area()<<endl; 
						}	
					}
				}
			}
			Volume_.push_back(volume);

			delete poly;
			delete Cer;
		}
		Area.push_back(Area_);
		Volume.push_back(Volume_);
		ncenter.push_back(ncenter_);
	}

	for(unsigned int i=0;i<Nbod;i++) g[0]+=1./double(Nbod)*ncenter[i][0]/Area[i][0]/meandensity;
	for(unsigned int i=0;i<Nbod;i++) compact[0]+=1./double(Nbod)*Volume[i][0]/Area[i][0]/meancompact;
	
	for (unsigned int k=1;k<Npoint;k++)
	{
		for(unsigned int i=0;i<Nbod;i++) 
		{
			g[k]+=1./double(Nbod)*ncenter[i][k]/(Area[i][k]-Area[i][k-1])/meandensity;
		//	if(ncenter[i][k]<0) cout<<"ncenter[i][k]:="<<ncenter[i][k]<<endl;
		//	if(Area[i][k]-Area[i][k-1]<0) cout<<"Area["<<i<<"]["<<k<<"]-Area["<<i<<"]["<<k-1<<"]:="<<Area[i][k]-Area[i][k-1]<<endl;
		}
		for(unsigned int i=0;i<Nbod;i++) compact[k]+=1./double(Nbod)*(Volume[i][k]-Volume[i][k-1])/(Area[i][k]-Area[i][k-1])/meancompact;		
	}

	delete Prb;
	
	ofstream rfd("Analyse/rfd.txt",ios::out);
	for( unsigned int i=0;i<Npoint;++i)
		rfd<<0.5*(i+0.5)*amp/rmean<<"	"<<g[i]<<"	"<<compact[i]<<endl;
	rfd.close();
}

/*void Biaxial_A::rfd( unsigned int Npoint, unsigned int Nrmean)
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

void Biaxial_A::removeBody(body2d * rm)
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

void Biaxial_A::removeRattlers()
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

void Biaxial_A::reduceRattlers()
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

	void Biaxial_A::growRattlers( double dr)
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

	void Biaxial_A::writePS( const char * fname)
	{
				
		ofstream ps(fname);
		sys_->spl()->updateBoundaries();
	sys_->spl()->radiusExtrema(1);
	
	double topy_,bottomy_,leftx_,rightx_;
	
	for (unsigned i = 0; i<sys_->spl()->lbody().size(); i++)
	{
		if (sys_->spl()->body(i)->type()==0)
		{
			topy_ = sys_->spl()->body(i)->y();
			bottomy_ = sys_->spl()->body(i)->y();
			rightx_ = sys_->spl()->body(i)->x();
			leftx_ = sys_->spl()->body(i)->x();
			break;
		}
	}
	
	for (unsigned i = 0; i<sys_->spl()->lbody().size(); i++)
	{
		if (sys_->spl()->body(i)->type()==2 && sys_->spl()->body(i)->rot() == 0)
		{
			topy_=max(topy_,sys_->spl()->body(i)->y());
			bottomy_=min(bottomy_,sys_->spl()->body(i)->y());
		}
		else if (sys_->spl()->body(i)->type()==2 && sys_->spl()->body(i)->rot() != 0)
		{
			rightx_=max(rightx_,sys_->spl()->body(i)->x());
			leftx_=min(leftx_,sys_->spl()->body(i)->x());
		}
	}
	
	
	
	double R = sys_->spl()->rmax();
	
	double Xmin = leftx_ + 1.5*R;//sys_->spl()->xmin() + R;
	double Ymin = bottomy_ + 1.5*R;//sys_->spl()->ymin() + R;
	
	double xmin_ = leftx_ - 7.*R;//sys_->spl()->xmin() - 5.*R;
	double ymin_ = bottomy_ - 7.*R;//sys_->spl()->ymin() - 5.*R;
	double xmax_ = rightx_ + 7.*R;//1.651510e-01 + 5.*R;//sys_->spl()->xmax() + 5.*R;
	double ymax_ = topy_ + 7.*R;//0.124471 + 5.*R;//sys_->spl()->ymax() + 5.*R;
	
	double zoom = zoom_;
	
	//cout << "zoom_ " << zoom_ << endl;
	
	
	double x_offset = fabs(Xmin*zoom);
	double y_offset = fabs(Ymin*zoom);
	
	//!------ Header for ps file	
	ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	//ps<<"%%BoundingBox:"<<" "<<"-30 -20 515 550"<<endl;
	ps<<"%%BoundingBox:"<<" "<<x_offset+xmin_*zoom<<" "<<y_offset+ymin_*zoom<<" "<<x_offset+xmax_*zoom<<" "<<y_offset+ymax_*zoom<<endl;
	ps<<"%%Pages: 1"<<endl;
	ps<<"0.1 setlinewidth 0. setgray "<<endl;
	
	double x_A,y_A,x_B,y_B,x_C,y_C,x_D,y_D;
	x_A = rightx_ + R;//sys_->spl()->xmax() + R;//   sys_->spl()->body(3)->x()+R;
	y_A = bottomy_ - R;//sys_->spl()->ymin() - R;//   sys_->spl()->body(0)->y()-R;
	x_B = leftx_ - R;//sys_->spl()->xmin() - R;//   sys_->spl()->body(2)->x()-R;
	y_B = bottomy_ - R;//sys_->spl()->ymin() - R;//   sys_->spl()->body(0)->y()-R;
	x_C = leftx_ - R;//sys_->spl()->xmin() - R;//   sys_->spl()->body(2)->x()-R;
	y_C = topy_ + R;//sys_->spl()->ymax() + R;//   sys_->spl()->body(1)->y()+R;
	x_D = rightx_ + R;//sys_->spl()->xmax() + R;//   sys_->spl()->body(3)->x()+R;
	y_D = topy_ + R;//sys_->spl()->ymax() + R;//   sys_->spl()->body(1)->y()+R;
	
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
	
    ps<<x_A*zoom + x_offset + 4.*R*zoom <<" "<<y_A*zoom + y_offset + 1.75*R*zoom <<" moveto "
	  <<x_B*zoom + x_offset - 4.*R*zoom <<" "<<y_B*zoom + y_offset + 1.75*R*zoom <<" lineto "
	  <<x_B*zoom + x_offset - 4.*R*zoom <<" "<<y_B*zoom + y_offset - .0*R*zoom   <<" lineto "
	  <<x_A*zoom + x_offset + 4.*R*zoom <<" "<<y_A*zoom + y_offset - .0*R*zoom   <<" lineto "
	  <<x_A*zoom + x_offset + 4.*R*zoom <<" "<<y_A*zoom + y_offset + 1.75*R*zoom <<" lineto stroke"<<endl;
	
	ps<<x_D*zoom + x_offset - 1.75*R*zoom <<" "<<y_D*zoom + y_offset + 4.*R*zoom <<" moveto "
	  <<x_A*zoom + x_offset - 1.75*R*zoom <<" "<<y_A*zoom + y_offset - 4.*R*zoom <<" lineto "
	  <<x_A*zoom + x_offset + .0*R*zoom   <<" "<<y_A*zoom + y_offset - 4.*R*zoom <<" lineto "
	  <<x_D*zoom + x_offset + .0*R*zoom   <<" "<<y_D*zoom + y_offset + 4.*R*zoom <<" lineto "
	  <<x_D*zoom + x_offset - 1.75*R*zoom <<" "<<y_D*zoom + y_offset + 4.*R*zoom <<" lineto stroke"<<endl;
	
	ps<<x_D*zoom + x_offset + 4.*R*zoom <<" "<<y_D*zoom + y_offset - 1.75*R*zoom <<" moveto "
	  <<x_C*zoom + x_offset - 4.*R*zoom <<" "<<y_C*zoom + y_offset - 1.75*R*zoom <<" lineto "
	  <<x_C*zoom + x_offset - 4.*R*zoom <<" "<<y_C*zoom + y_offset + .0*R*zoom   <<" lineto "
	  <<x_D*zoom + x_offset + 4.*R*zoom <<" "<<y_D*zoom + y_offset + .0*R*zoom   <<" lineto "
	  <<x_D*zoom + x_offset + 4.*R*zoom <<" "<<y_D*zoom + y_offset - 1.75*R*zoom <<" lineto stroke"<<endl;
	
	ps<<x_B*zoom + x_offset + 1.75*R*zoom <<" "<<y_B*zoom + y_offset - 4.*R*zoom <<" moveto "
	  <<x_C*zoom + x_offset + 1.75*R*zoom <<" "<<y_C*zoom + y_offset + 4.*R*zoom <<" lineto "
	  <<x_C*zoom + x_offset - .0*R*zoom   <<" "<<y_C*zoom + y_offset + 4.*R*zoom <<" lineto "
	  <<x_B*zoom + x_offset - .0*R*zoom   <<" "<<y_B*zoom + y_offset - 4.*R*zoom <<" lineto "
	  <<x_B*zoom + x_offset + 1.75*R*zoom <<" "<<y_B*zoom + y_offset - 4.*R*zoom <<" lineto stroke"<<endl;
	
    
	/*
	ps<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" moveto "<<x_B*zoom + x_offset<<" "<<y_B*zoom + y_offset<<" "<<" lineto "<<x_C*zoom + x_offset<<" "<<y_C*zoom + y_offset<<" lineto "<<x_D*zoom + x_offset<<" "<<y_D*zoom + y_offset<<" lineto "<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" "<<" lineto stroke"<<endl;
	ps<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" moveto "<<x_B*zoom + x_offset<<" "<<y_B*zoom + y_offset<<" "<<" lineto "<<x_C*zoom + x_offset<<" "<<y_C*zoom + y_offset<<" lineto "<<x_D*zoom + x_offset<<" "<<y_D*zoom + y_offset<<" lineto "<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" "<<" lineto stroke"<<endl;
	ps<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" moveto "<<x_B*zoom + x_offset<<" "<<y_B*zoom + y_offset<<" "<<" lineto "<<x_C*zoom + x_offset<<" "<<y_C*zoom + y_offset<<" lineto "<<x_D*zoom + x_offset<<" "<<y_D*zoom + y_offset<<" lineto "<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" "<<" lineto stroke"<<endl;
	*/
	//ps<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" moveto "<<x_B*zoom + x_offset<<" "<<y_B*zoom + y_offset<<" "<<" lineto "<<x_C*zoom + x_offset<<" "<<y_C*zoom + y_offset<<" lineto "<<x_D*zoom + x_offset<<" "<<y_D*zoom + y_offset<<" lineto "<<x_A*zoom + x_offset<<" "<<y_A*zoom + y_offset<<" "<<" lineto stroke"<<endl;
	/*
	if (sys_->spl()->body(0)->type() ==2)
		ps<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" moveto "
		<<sys_->spl()->body(0)->xmin()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmin()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymax()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymax()*zoom + y_offset<<" lineto "
		<<sys_->spl()->body(0)->xmax()*zoom + x_offset<<" "<<sys_->spl()->body(0)->ymin()*zoom + y_offset<<" lineto stroke"<<endl;
	*/
	//!------ Draw sample
	
	if(displaySample)
	{
		DataSet V, pos_W, neg_W;
		
		for(unsigned i=0 ; i<sys_->spl()->lbody().size() ; ++i)
		{
			//if (sys_->spl()->body(i)->type() ==0)
			{
				V.add(sqrt(sys_->spl()->body(i)->vx()*sys_->spl()->body(i)->vx()+sys_->spl()->body(i)->vy()*sys_->spl()->body(i)->vy()));
				if (sys_->spl()->body(i)->vrot() >= 0)
					pos_W.add(sys_->spl()->body(i)->vrot());
				else
					neg_W.add(sys_->spl()->body(i)->vrot());
			}
		}
		
		V.extractValues();
		pos_W.extractValues();
		neg_W.extractValues();
		
		//!------ Draw Particles and particles velocity
		
		for (unsigned int j=0 ; j<sys_->spl()->lbody().size() ; ++j)
		{
			if (sys_->spl()->body(j)->type() ==0)
			{
				/*
				if (sys_->spl()->body(j)->vrot()<0)
					ps<<"/cirg {" <<(sys_->spl()->body(j)->vrot()-neg_W.min())/(neg_W.max()+neg_W.min()) << " 1 0 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;//make circle 
				else
					ps<<"/cirg {1 " <<1-(sys_->spl()->body(j)->vrot()-pos_W.min())/(pos_W.max()+pos_W.min()) << " 0 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;//make circle 
				*/
				
				//if (sqrt(sys_->spl()->body(j)->vx()*sys_->spl()->body(j)->vx()+sys_->spl()->body(j)->vy()*sys_->spl()->body(j)->vy()) == 0)
					ps<<"/cirg {0.25 0.25 0.25 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;
				//else
				//	ps<<"/cirg {1 " <<(sqrt(sys_->spl()->body(j)->vx()*sys_->spl()->body(j)->vx()+sys_->spl()->body(j)->vy()*sys_->spl()->body(j)->vy())-V.min())/(V.max()+V.min()) << " 0 setrgbcolor 0.0 setlinewidth 0 360 arc gsave  fill grestore} def"<<endl;//make circle 
				
				ps<<"newpath "<<x_offset + sys_->spl()->body(j)->x()*zoom<<" "<<y_offset + sys_->spl()->body(j)->y()*zoom<<" "<<sys_->spl()->body(j)->sizeVerlet()*zoom<<" cirg"<<endl;
			}
		}
		
		neg_W.setClear();
		pos_W.setClear();
		V.setClear();
	}
	
	//!------ Draw branch vector
	/*
	if(displayBranch)
	{
		ps<<"/coul_green {1 setlinecap 0 1 0 setrgbcolor} def"<<endl;
		ps<<"% x0 y0 x1 y1 L"<<endl;
		ps<<"/L{newpath moveto lineto gsave grestore stroke}def"<<endl; 
		ps<<"0.25 setlinewidth"<<endl;
		ps<<"coul2"<<endl;
		
		for(unsigned int i=0 ; i<cla_->linterdof().size() ; ++i)
		{
			ps<<cla_->linterdof(i)->first()->mcx()*zoom<<" "<<cla_->linterdof(i)->first()->mcy()*zoom<<" "<<cla_->linterdof(i)->second()->mcx()*zoom<<" "<<cla_->linterdof(i)->second()->mcy()*zoom<<" L"<<endl;
		}
	}
	*/ 
	//!------ Draw contact forces
	if(displayForce)
	{
		DataSet NEGATIVE_FN;
		DataSet POSITIVE_FN;
		
		sys_->spl()->radiusExtrema(2);
		
		for(unsigned i=0 ; i<sys_->nwk()->clist().size() ; ++i)
		{
			if(sys_->nwk()->inter(i)->type()==0)
			{
				if(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()<0)
					NEGATIVE_FN.add(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
				else
					POSITIVE_FN.add(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn());
			}
		}
		
		
		//ps<<"/coul_red {1 setlinecap 1 0 0 setrgbcolor} def"<<endl;
		ps<<"/coul_ligth_blue {1 setlinecap 0 1 1 setrgbcolor} def"<<endl;
		//ps<<"/coul_ligth_blue {1 setlinecap 0 0 0.50 setrgbcolor} def"<<endl;
		ps<<"/coul_blue {1 setlinecap 0 0 1 setrgbcolor} def"<<endl;
		
		ps<<"/coul_blanc {1 setlinecap 1 1 1 setrgbcolor} def"<<endl;
		
		ps<<0.0*(sys_->nwk()->inter(sys_->nwk()->clist(0))->fn())<<" setlinewidth coul_blanc"<<endl;
		ps<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->x()*zoom + x_offset<<" "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->y()*zoom + y_offset<<" "<<"moveto "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->x()*zoom + x_offset<<" "
		<<sys_->nwk()->inter(sys_->nwk()->clist(0))->first()->y()*zoom + y_offset<<" lineto stroke"<<endl;
		
		
		NEGATIVE_FN.extractValues();
		POSITIVE_FN.extractValues();
		
		for(unsigned i=0 ; i<sys_->nwk()->clist().size() ; ++i)
		{
			if(sys_->nwk()->inter(sys_->nwk()->clist(i))->type()==0)
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
					ps	<<20.*(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn())/*-NEGATIVE_FN.min())/(NEGATIVE_FN.min()+NEGATIVE_FN.max())*/<<" setlinewidth coul_ligth_blue"<<endl;
					ps	<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()*zoom + x_offset<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()*zoom + y_offset<<" "<<"moveto"<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()*zoom + x_offset<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()*zoom + y_offset<<" "<<"lineto stroke"<<endl;
				}
				else 
				{
					ps	<<10*(sys_->nwk()->inter(sys_->nwk()->clist(i))->fn()-POSITIVE_FN.min())/(POSITIVE_FN.min()+POSITIVE_FN.max())<<" setlinewidth coul_blue"<<endl;
					ps	<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->x()*zoom + x_offset<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->y()*zoom + y_offset<<" "<<"moveto"<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->x()*zoom + x_offset<<" "
						<<sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->y()*zoom + y_offset<<" "<<" "<<"lineto stroke"<<endl;
				}
				
			}
		}
		
		ofstream test_("test.dat",ios::app);
		test_ << POSITIVE_FN.max() << " " << POSITIVE_FN.min() << endl;
		test_.close();
		
		NEGATIVE_FN.setClear();
		POSITIVE_FN.setClear();
	}
	

		
	ps.close();
}

//Les particules en noir pour étudier la distribution des volumes
void Biaxial_A::writePS2( const char * fname)
	{
				
		ofstream ps(fname);

		this->sys()->spl()->updateBoundaries(); 		//Ajouter le 25/11/2011 	
		double height= sys_->spl()->ymax()-sys_->spl()->ymin();
		double width = sys_->spl()->xmax()-sys_->spl()->xmin();
		double xoffset=5.;
		double zoom_=min((595-2*xoffset)/height,(842-2*xoffset)/width);
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
void Biaxial_A::display_data(char * fname)
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

bool Biaxial_A::Inter (polyg *p, polyg *q)
{
	bool Inter = false;
	double c,s,d;
	polyg *r=new polyg; //polygone qui est parti commun de p et q
	gdm::vertex q1, q2, vtex;
	vector<gdm::vertex> Psection; //Point section entre deux polygone
	
	d=sqrt((p->x()-q->x())*(p->x()-q->x()) + (p->y()-q->y())*(p->y()-q->y()));
	if ( d>(p->Rout()+q->Rout())) return false;
	else
	{
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
}


polyg * Biaxial_A::Intersection (polyg *p, polyg *q)
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
//	cout<<"size1:="<<r->Vertex().size()<<endl;
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
//	cout<<"size2:="<<r->Vertex().size()<<endl;

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
//	cout<<"size3:="<<r->Vertex().size()<<endl;

//	cout<<"ligne 3146:= "<<r->Vertex().size()<<endl;
	
	//arranger les vertex du polygone r 
	r->Arrangement();
	//Invert ordre vertex polygone if CW
	if (r->ispolysimple()<0) r->Invert();
	r->adjustCenter();
	r->Fill(1);
	return r;
}

double Biaxial_A::compactness(double &h,double &l)
{
	double compact=0;
	
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	compact+=sys_->spl()->body(i)->Area();
	compact/=(h*l);
	return compact;
}


double Biaxial_A::compactness(double &x, double &y, double &h, double &l)
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
	double d;
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	{
		d=sqrt(pow(sys_->spl()->body(i)->x()-Rec->x(),2)+pow(sys_->spl()->body(i)->y()-Rec->y(),2));
		if(d<=(sys_->spl()->body(i)->sizeVerlet()+Rec->sizeVerlet()))
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
	}
	compact/=(h*l);
	//cout<<"Compactness: "<<compact<<endl;
	return compact;
}

double Biaxial_A::compactness_disk(double &x, double &y, double &h, double &l)
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

double Biaxial_A::compactness(double &ratio)
{
	/*
	const double g=10;
	double compact=0;
	double h,w,h1;
	double vertex1_x, vertex1_y;
	double vertex2_x, vertex2_y;
	double c,s;
	double dmoy=0;
	double mmoy=0;
	double vmax=0;
	double v;
	double Ec=0;
	double R, d, alpha, area;
	unsigned int Nb;
	
	sys_->spl()->updateBoundaries();
	h = sys_->spl()->ymax()-sys_->spl()->ymin();
	w = sys_->spl()->xmax()-sys_->spl()->xmin();
	h*=ratio;
	h1=h+sys_->spl()->ymin();
	
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	{
		if (sys_->spl()->body(i)->ymax()<=h1)
		{
			compact+=sys_->spl()->body(i)->Area();
		}		
		else if (sys_->spl()->body(i)->ymin()<h1)
		if (sys_->spl()->body(i)->type()==1)
		{
			polyg * polyg_ini;
			polyg * polyg_repla = new polyg;
			polyg_ini=dynamic_cast<polyg*> (sys_->spl()->body(i));
			c=cos(polyg_ini->rot());
			s=sin(polyg_ini->rot());
			*polyg_repla=*polyg_ini;
			polyg_repla->Vertex().clear();
			gdm::vertex vertex_;
			//cout<<"polyg_ini->Vertex().size() ="<<polyg_ini->Vertex().size()<<endl;
			for (unsigned int j=0; j<polyg_ini->Vertex().size(); j++)
			{
				vertex1_x = polyg_ini->x()+polyg_ini->Vertex(j).x()*c-polyg_ini->Vertex(j).y()*s;
				vertex1_y = polyg_ini->y()+polyg_ini->Vertex(j).x()*s+polyg_ini->Vertex(j).y()*c;
				if (j<polyg_ini->Vertex().size()-1)
				{
					vertex2_x = polyg_ini->x()+polyg_ini->Vertex(j+1).x()*c-polyg_ini->Vertex(j+1).y()*s;
					vertex2_y = polyg_ini->y()+polyg_ini->Vertex(j+1).x()*s+polyg_ini->Vertex(j+1).y()*c;
				//	cout<<"j "<<j<<" vertex2_x: "<< vertex2_x <<" vertex2_y: "<<vertex2_y <<endl;
				}
				else
				{
					vertex2_x = polyg_ini->x()+polyg_ini->Vertex(0).x()*c-polyg_ini->Vertex(0).y()*s;
					vertex2_y = polyg_ini->y()+polyg_ini->Vertex(0).x()*s+polyg_ini->Vertex(0).y()*c;				
				}
				
				if ((vertex1_y<=h1) && (vertex2_y<=h1))
				{
				//	cout<<"case 1"<<endl;
					vertex_.x()=vertex1_x;
					vertex_.y()=vertex1_y;
					polyg_repla->Vertex().push_back(vertex_);
				}
				else if ((vertex1_y<=h1) && (vertex2_y>h1))
				{
				//	cout<<"case 2"<<endl;
					//Ajouter vertex j
					vertex_.x()=vertex1_x;
					vertex_.y()=vertex1_y;
					polyg_repla->Vertex().push_back(vertex_);
					
					//Ajouter section entre Vertex(j+(j+1)) et y=height
					vertex_.y()=h1;
					vertex_.x()=vertex2_x-(vertex2_x-vertex1_x)*(vertex2_y-vertex_.y())/(vertex2_y-vertex1_y);
					polyg_repla->Vertex().push_back(vertex_);
				}
				else if ((vertex1_y>h1) && (vertex2_y<=h1))
				{
				//	cout<<"case 3"<<endl;
					vertex_.y()=h1;
					vertex_.x()=vertex1_x-(vertex1_x-vertex2_x)*(vertex1_y-vertex_.y())/(vertex1_y-vertex2_y);
					polyg_repla->Vertex().push_back(vertex_);					
				}
			}	

			if (polyg_repla->Vertex().size()==0){
				compact+=0;
			}
			else
			{
				polyg_repla->adjustCenter ();
				polyg_repla->Fill(1.);
				compact+=polyg_repla->Area();
			}

			delete polyg_repla;
		}
		else if (sys_->spl()->body(i)->type()==0)
		{
			body_=sys_->spl()->body(i);
			R=body_->sizeVerlet();			
			d=h1-body_->y();
			alpha =acos(fabs(d)/R);
			area=alpha*R*R-0.5*R*R*sin(2*alpha);
			if (d>0) area=M_PI*R*R-area;			
			compact+=area;	
		}

	}
	
	compact/=(h*w);
	cout<<"Compactness: "<<compact<<endl;
	//Calcule vmax,Ec moyence
	Nb=sys_->spl()->lbody().size()-sys_->lctrl().size();
	for (unsigned int i=sys_->lctrl().size(); i<sys_->spl()->lbody().size();i++)
	{
		dmoy+=2*sys_->spl()->body(i)->sizeVerlet();
		mmoy+=sys_->spl()->body(i)->mass();
		v=sqrt(pow(sys_->spl()->body(i)->vx(),2)+pow(sys_->spl()->body(i)->vy(),2));
		vmax=(v>vmax) ? v : vmax;
		Ec+=0.5*sys_->spl()->body(i)->mass()*pow(v,2);
	}
	dmoy/=Nb;
	mmoy/=Nb;
	Ec=Ec/(Nb*mmoy*g*dmoy);
	vmax/=sqrt(g*dmoy);//vmaxcritique
	ofstream str("Analyse/compactness.txt",ios::app);
	str<<time<<"\t"<<compact<<"\t"<<vmax<<"\t"<<Ec<<endl;
	//str<<time<<"\t"<<compact<<"\t"<<sys_->spl()->xmin()<<"\t"<<sys_->spl()->xmax()<<"\t"<<sys_->spl()->ymin()<<"\t"<<sys_->spl()->ymax()<<endl;
	str.close();
	return compact;
	*/
		return 0.;
}
	
void Biaxial_A::fracture()
{
	unsigned int nfracture=0;
	unsigned int ncasse=0;
	
	for (unsigned int i=0; i<sys_->grpRel()->fragmentation()->lgrain().size();i++)
	if (sys_->grpRel()->fragmentation()->lgrain()[i]==true) nfracture++;
	
	//Calcule nombre particule casse
	
	for (unsigned int i=0; i<sys_->grpRel()->fragmentation()->Casse().size(); i++)
	if (sys_->grpRel()->fragmentation()->Casse()[i]) ncasse++;
	
	ofstream fra("Analyse/fracture.txt",ios::app);
	fra<<time<<"\t"<<epsq<<"\t"<<sys_->grpRel()->fragmentation()->numfissure()<<"\t"<<nfracture*1.0/sys_->grpRel()->fragmentation()->lgrain().size()<<"\t"<<ncasse*1.0/sys_->grpRel()->fragmentation()->lgrain().size()<<endl;
	fra.close();

	cout<<"Nombre fissure : "<<sys_->grpRel()->fragmentation()->numfissure()<<"\t Nombre particule fracture : "<<nfracture<<"\t Nombre particule casse:	"<<ncasse<<endl;
}

void Biaxial_A::granulometrique(unsigned int Nbin )
{
	unsigned int Ncluster;	// Nombre cluster
	unsigned int index;
	double volume,diametre,vtotal;
	double dmin,dmax;
	vector <double> Vcluster; //Volume cluster
	vector <double> Dcluster; // Diamètre cluster
	vector <double> Granulo;
	vector <double> Granulocumule;
	
	Ncluster=sys_->grpRel()->fragmentation()->Cluster().size();
	
	for (unsigned int i=0; i<Ncluster; i++)
	{
		volume=0;
		for (unsigned int j=0; j<sys_->grpRel()->fragmentation()->Cluster()[i].size(); j++)
		{
			index=sys_->grpRel()->fragmentation()->Cluster()[i][j]+4;
			volume+=sys_->spl()->body(index)->Area();
		}
		Vcluster.push_back(volume);
	}
	
	for (unsigned int i=0; i<Ncluster; i++)
	{
		diametre=sqrt(4*Vcluster[i]/M_PI);
		Dcluster.push_back(diametre);
	}
	
	dmin=dmax=Dcluster[0];
	for (unsigned int i=1; i<Ncluster; i++)
	{
		dmin = (dmin>Dcluster[i]) ?  Dcluster[i] : dmin;	
		dmax = (dmax<Dcluster[i]) ?  Dcluster[i] : dmax;
	}
	
	
	vtotal=0.;
	for (unsigned int i=0;i<Ncluster;i++) vtotal+=Vcluster[i];
	
	
	double amp=(dmax-dmin)/Nbin;
	for (unsigned int i=0;i<Nbin;i++)
	{
		volume=0.0;
		for (unsigned int j=0;j<Ncluster;j++)
		if ( (Dcluster[j]>=dmin+i*amp) &&(Dcluster[j]<dmin+(i+1.0)*amp)) volume+=Vcluster[j];
		Granulo.push_back(volume/vtotal);
	}
	
	for (unsigned int i=0;i<Nbin;i++)
	{
		volume=0.0;
		for (unsigned int j=0;j<Ncluster;j++)
		if ( Dcluster[j]<=(dmin+(i+1.0)*amp)) volume+=Vcluster[j];
		Granulocumule.push_back(volume/vtotal);
	}
	
		volume=0;
		for (unsigned int j=0;j<Nbin;j++) volume+=Granulo[j];
		cout<<"volume:="<<volume<<endl;

	/*ofstream gra("Analyse/granulometrique.txt",ios::app);
	for(unsigned int i=0; i<Nbin; i++)
	gra<<time<<"\t"<<i<<"\t"<<dmin+(i+0.5)*amp<<"\t"<<Granulo[i]<<"\t"<<0.1*(i+0.5)<<"\t"<<Granulocumule[i]<<endl;
	gra.close();*/

	ofstream gra("Analyse/granulometrique.txt",ios::app);
	gra<<time<<"\t"<<0<<"\t"<<dmin*1.0E+06<<"\t"<<0.0<<"\t"<<0.0<<endl;
	for(unsigned int i=0; i<Nbin; i++)
	//gra<<time<<"\t"<<i<<"\t"<<dmin+(i+0.5)*amp<<"\t"<<Granulo[i]<<"\t"<<0.1*(i+1.0)<<"\t"<<Granulocumule[i]<<endl;
	gra<<time<<"\t"<<i+1<<"\t"<<(dmin+(i+1.0)*amp)*1.0E+06<<"\t"<<0.1*(i+1.0)<<"\t"<<Granulocumule[i]<<endl;
	gra.close();
}

void Biaxial_A::heterogeneity(unsigned int Ni,unsigned int Nf)
{
	double hight,large;
	double xmin,xmax,ymin,ymax;
	DataSet* Compact;
	double compact;
	double DOH, mean, variance;// Degree of heterogeneity
	double dmoy=0.;
	unsigned int Nb;
	
	sys_->spl()->updateBoundaries();
	xmin=sys_->spl()->xmin();
	xmax=sys_->spl()->xmax();
	ymin=sys_->spl()->ymin();
	ymax=sys_->spl()->ymax();
	
	Nb=sys_->spl()->lbody().size()-sys_->lctrl().size();
	for(unsigned int i=sys_->lctrl().size();i<sys_->spl()->lbody().size();i++)
	dmoy+=2*sys_->spl()->body(i)->sizeVerlet();
	dmoy/=double(Nb);
	cout<<"dmoy:="<<dmoy<<endl;
		
	ofstream het("Analyse/heterogeneity.txt",ios::out);
	for(unsigned int l=Nf; l>Ni-1; l--)
	{
		hight=(ymax-ymin)/double(l);
		large=(xmax-xmin)/double(l);
			
		Compact = new DataSet;
		for (unsigned int j=1;j<=l;j++)
			for(unsigned int k=1;k<=l;k++)
			{
				double centrex=xmin+(j-0.5)*large;
				double centrey=ymin+(k-0.5)*hight;
				if (sys_->spl()->body(4)->type()==0) compact = compactness_disk(centrex,centrey,hight,large);
				else compact=compactness(centrex,centrey,hight,large);			
				Compact->add(compact);				
			}
		mean=Compact->Mean();
		variance=Compact->Variance();
		DOH=variance/mean*100.;
		het<<time<<"\t"<<l<<"\t"<< 0.5*(hight+large)/dmoy<<"\t"<<variance<<"\t"<<DOH<<"\t"<<Compact->Min()<<"\t"<<Compact->Max()<<"\t"<<Compact->Max()/Compact->Min()<<endl;
		delete Compact;		
	}
	het.close();
}	

void Biaxial_A::clustercontactdouble()
{
	vector <unsigned int> cluster; //cluster contient des indexs des particules qui sont liées par contact double
	vector < vector<unsigned int> > Listcluster; // List des clusters
	vector <unsigned int> Listdouble; // List des particule qui a au moins un contact double
	inter2d * oxo;
	unsigned id1,id2,id;
	
	double pc; // Proportion des particules dans les clusters;
	unsigned int Nb; // Nombre des particule dans les cluster
	
	Listdouble.clear();
	bool logic1,logic2;
	double fn1,fn2;
	for (unsigned int i=0; i<sys_->nwk()->clist().size();++i)
	{
		oxo=sys_->nwk()->inter(sys_->nwk()->clist(i));
		fn1 =0.0;
		fn2 =0.0;
		if(oxo->rang()==2)
		{
			oxo->current()=0;
			fn1=oxo->fn();
			oxo->current()=1;
			fn2=oxo->fn();
			oxo->current()=0;
		}
		
		if (oxo->rang()==2 && (fn1!=0.) && (fn2!=0.)&&(oxo->first()->id()>3)&&(oxo->second()->id()>3))
		{
			id1=oxo->first()->id();
			id2=oxo->second()->id();
			logic1=true;
			logic2=true;
			unsigned int j=0;
			while(j<Listdouble.size() and (logic1 || logic2))
			{
				if (id1==Listdouble[j]) logic1=false;
				if (id2==Listdouble[j]) logic2=false;
				j++;
			}
			if (logic1) Listdouble.push_back(id1);
			if (logic2) Listdouble.push_back(id2);
		}
	}	
	
	Nb=Listdouble.size();
	pc=double(Nb)/double(sys_->spl()->lbody().size()-sys_->lctrl().size());
	
	//cout<<"SIZE:="<<Listdouble.size()<<endl;
	Listcluster.clear();
	ClusterCD.clear();
	while(!Listdouble.empty())
	{
		//cout<<"Listdouble.size():="<<Listdouble.size()<<endl<<endl<<endl;;
		cluster.clear();
		cluster.push_back(Listdouble[0]);
		Listdouble.erase(Listdouble.begin());
		unsigned int l=0;
		while((l<cluster.size()) && (!Listdouble.empty()))
		{
			id=cluster[l];
			for(unsigned int j=0; j<sys_->nwk()->clist().size();++j)
			{
				oxo=sys_->nwk()->inter(sys_->nwk()->clist(j));
				id1=oxo->first()->id();
				id2=oxo->second()->id();
				fn1 =0.0;
				fn2 =0.0;
				if((oxo->rang()==2)&&(id1>3)&&(id2>3))
				{
					oxo->current()=0;
					fn1=oxo->fn();
					oxo->current()=1;
					fn2=oxo->fn();
					oxo->current()=0;
				}

				if ((id1==id) &&(fn1!=0.)&&(fn2!=0.)) 
				{
					unsigned int k=0;
					bool stop=true;
					while((k<Listdouble.size())&&(stop))
					{
						if (Listdouble[k]==id2) {cluster.push_back(id2);Listdouble.erase(Listdouble.begin()+k);stop=false;}
						else k++;
					}
				}
				else if ((id2==id) &&(fn1!=0.)&&(fn2!=0.)) 
				{
					unsigned int k=0;
					bool stop=true;
					while((k<Listdouble.size())&&(stop))
					{
						if (Listdouble[k]==id1) {cluster.push_back(id1);Listdouble.erase(Listdouble.begin()+k);stop=false;}
						else k++;
					}
				}
			}
			l++;
			//cout<<"cluster.size():="<<cluster.size()<<"		l:="<<l<<endl;
		}
		Listcluster.push_back(cluster);
		ClusterCD.push_back(cluster);
	}
	
	/*for(unsigned int i=0;i<Listcluster.size();i++)
	{
		for(unsigned int j=0;j<Listcluster[i].size();j++) cout<<Listcluster[i][j]<<" ";
		cout<<endl;
	}	*/
	unsigned int Nc; // Nombre des cluster
	double Nm; // Nombre moyenne des particules dans un cluster
	unsigned int Nmax=0; // Nombre maximume des particules dans un cluster 
	Nc=Listcluster.size();
	Nm=double(Nb)/double(Nc);
	
	
	
	for(unsigned int i=0;i<Listcluster.size();i++)
	{
		if (Nmax<Listcluster[i].size()) Nmax=Listcluster[i].size();
	}
	
	/*cout<<"Nombre des clusters:="<<Nc<<endl;
	cout<<"Proportion des particules dans des clusters:= "<<pc<<endl;
	cout<<"Nombre moyenne des particules dans un clusters:="<<double(Nb)/double(Nc)<<endl;
	cout<<"Nombre maximume des particules dans un clusters:="<<Nmax<<endl;*/

	ofstream clust("Analyse/cluster_number.txt",ios::app);
	clust<<time<<" "<<Nc<<" "<<pc<<" "<<Nm<<" "<<Nmax<<endl;
	clust.close();
}
/*void Biaxial_A::clustercontactdouble()
{
	vector <unsigned int> cluster; //cluster contient des indexs des particules qui sont liées par contact double
	vector < vector<unsigned int> > Listcluster; // List des clusters
	vector <unsigned int> Listdouble; // List des particule qui a au moins un contact double
	inter2d * oxo;
	unsigned id1,id2,id;
	
	double pc; // Proportion des particules dans les clusters;
	unsigned int Nb; // Nombre des particule dans les cluster
	
	Listdouble.clear();
	bool logic1,logic2;
	double fn1,fn2;
	for (unsigned int i=0; i<sys_->nwk()->clist().size();++i)
	{
		oxo=sys_->nwk()->inter(sys_->nwk()->clist(i));
		fn1 =0.0;
		fn2 =0.0;
		if(oxo->rang()==2)
		{
			oxo->current()=0;
			fn1=oxo->fn();
			oxo->current()=1;
			fn2=oxo->fn();
			oxo->current()=0;
		}
		
		if ( ((oxo->rang()==1) || (fn1*fn2==0.))&&(oxo->first()->id()>3)&&(oxo->second()->id()>3))
		{
			id1=oxo->first()->id();
			id2=oxo->second()->id();
			logic1=true;
			logic2=true;
			unsigned int j=0;
			while(j<Listdouble.size() and (logic1 || logic2))
			{
				if (id1==Listdouble[j]) logic1=false;
				if (id2==Listdouble[j]) logic2=false;
				j++;
			}
			if (logic1) Listdouble.push_back(id1);
			if (logic2) Listdouble.push_back(id2);
		}
	}	
	
	Nb=Listdouble.size();
	pc=double(Nb)/double(sys_->spl()->lbody().size()-sys_->lctrl().size());
	
	cout<<"SIZE:="<<Listdouble.size()<<endl;
	Listcluster.clear();
	ClusterCD.clear();
	while(!Listdouble.empty())
	{
		//cout<<"Listdouble.size():="<<Listdouble.size()<<endl<<endl<<endl;;
		cluster.clear();
		cluster.push_back(Listdouble[0]);
		Listdouble.erase(Listdouble.begin());
		unsigned int l=0;
		while((l<cluster.size()) && (!Listdouble.empty()))
		{
			id=cluster[l];
			for(unsigned int j=0; j<sys_->nwk()->clist().size();++j)
			{
				oxo=sys_->nwk()->inter(sys_->nwk()->clist(j));
				id1=oxo->first()->id();
				id2=oxo->second()->id();
				fn1 =0.0;
				fn2 =0.0;
				if((oxo->rang()==2)&&(id1>3)&&(id2>3))
				{
					oxo->current()=0;
					fn1=oxo->fn();
					oxo->current()=1;
					fn2=oxo->fn();
					oxo->current()=0;
				}

				if ((id1==id) &&((fn1==0.)||(fn2==0.))) 
				{
					unsigned int k=0;
					bool stop=true;
					while((k<Listdouble.size())&&(stop))
					{
						if (Listdouble[k]==id2) {cluster.push_back(id2);Listdouble.erase(Listdouble.begin()+k);stop=false;}
						else k++;
					}
				}
				else if ((id2==id) &&((fn1==0.)||(fn2==0.))) 
				{
					unsigned int k=0;
					bool stop=true;
					while((k<Listdouble.size())&&(stop))
					{
						if (Listdouble[k]==id1) {cluster.push_back(id1);Listdouble.erase(Listdouble.begin()+k);stop=false;}
						else k++;
					}
				}
			}
			l++;
			//cout<<"cluster.size():="<<cluster.size()<<"		l:="<<l<<endl;
		}
		Listcluster.push_back(cluster);
		ClusterCD.push_back(cluster);
	}
	
	for(unsigned int i=0;i<Listcluster.size();i++)
	{
		for(unsigned int j=0;j<Listcluster[i].size();j++) cout<<Listcluster[i][j]<<" ";
		cout<<endl;
	}	
	unsigned int Nc; // Nombre des cluster
	double Nm; // Nombre moyenne des particules dans un cluster
	unsigned int Nmax=0; // Nombre maximume des particules dans un cluster 
	Nc=Listcluster.size();
	Nm=double(Nb)/double(Nc);
	
	
	
	for(unsigned int i=0;i<Listcluster.size();i++)
	{
		if (Nmax<Listcluster[i].size()) Nmax=Listcluster[i].size();
	}
	
	cout<<"Nombre des clusters:="<<Nc<<endl;
	cout<<"Proportion des particules dans des clusters:= "<<pc<<endl;
	cout<<"Nombre moyenne des particules dans un clusters:="<<double(Nb)/double(Nc)<<endl;
	cout<<"Nombre maximume des particules dans un clusters:="<<Nmax<<endl;

	ofstream clust("Analyse/cluster_number.txt",ios::app);
	clust<<time<<" "<<Nc<<" "<<pc<<" "<<Nm<<" "<<Nmax<<endl;
	clust.close();
}*/
void Biaxial_A::cluster_df()
{
	unsigned int nc; // Nombre des cluster
	unsigned int nmax,nmin;
	double nmoy;
	nmax=0;
	nmin=1000;
	nmoy=0.;
	
	nc=ClusterCD.size();
	for(unsigned int i=0; i<nc;i++) 
	{
		nmax=(nmax>ClusterCD[i].size()) ? nmax :ClusterCD[i].size() ;
		nmoy+=ClusterCD[i].size();
	}
	nmoy/=nc;
	
	vector<double> ncp; // Nombre des cluster qui a p particlues.
	
	ncp.push_back(0.);
	ncp.push_back(0.);
	unsigned int np; 
	for (unsigned int i=2;i<nmax+1;i++)
	{
		np=0;
		for(unsigned int j=0;j<nc;j++) if (ClusterCD[j].size()==i) np++;
		ncp.push_back(np);
	}
	ofstream clust_df("Analyse/Cluster/cluster_df.txt",ios::out);
	for(unsigned int i=2;i<ncp.size();i++) clust_df<<i<<"	"<<ncp[i]<<"	"<<double(ncp[i])/double(nc)<<endl;
	clust_df.close();
}

void Biaxial_A::cluster_volume()
{
	vector<double> volume;
	double vol;
	double vmin, vmax, vmean, vtot;
	unsigned int nc;
	
	nc=ClusterCD.size();
	for(unsigned int i=0;i<ClusterCD.size();i++)
	{
		vol=0.;
		for (unsigned int j=0;j<ClusterCD[i].size();j++)
		vol+=sys_->spl()->body(ClusterCD[i][j])->Area();
		volume.push_back(vol);
	}
	
	vmin=1.0E+06;
	vmax=0.;
	vtot=0.;
	
	for (unsigned int i=0;i<volume.size();i++)
	{
		vmin=(vmin>volume[i]) ? volume[i] : vmin;
		vmax=(vmax<volume[i]) ? volume[i] : vmax;
		vtot+=volume[i];
	}
	vmean=vtot/nc;
	cout<<"vmin:="<<vmin<<"		vmax:="<<vmax<<"	vmean:="<<vmean<<endl;
	
	vector<double> dv;
	vector<double> dn;
	unsigned int nbin=20;
	unsigned int n;
	double amp;
	
	amp=(vmax-vmin)/nbin;
	
	dv.clear();
	dn.clear();
	for(unsigned int i=0; i<nbin;i++)
	{
		n=0;
		vol=0.;
		for(unsigned int j=0; j<nc;j++)
		if ((volume[j]>vmin+i*amp)&&(volume[j]<=vmin+(i+1)*amp))
		{
			n++;
			vol+=volume[j];
		}
		dn.push_back(double(n)/double(nc));
		dv.push_back(vol/vtot);	
	}
	
	ofstream clust_dv("Analyse/Cluster/cluster_volume.txt",ios::out);
//	clust_dv<<time<<" "<<0<<" "<<vmin<<" "<<0<<" "<<0<<endl;
	for(unsigned int i=0; i<nbin;i++)
		clust_dv<<time<<" "<<double(i+1)/double(nbin)<<" "<<vmin+(i+1)*amp<<" "<<dn[i]<<" "<<dv[i]<<"	"<<double(dn[i])/double(nc)<<"	"<<dv[i]/vtot<<"	"<<dv[i]/prb_.area()<<endl;
	
	clust_dv.close();
	
	ofstream clust_vmax_min("Analyse/Cluster/cluster_vmax_min.txt",ios::out);
		clust_vmax_min<<time<<" "<<vmin<<"	"<<vmean<<"	"<<vmax<<"	"<<vtot<<endl;
	clust_vmax_min.close();
}

//Analyse des formes et des directions des clusters
void Biaxial_A::cluster_shape()
{
	vector<double> aspect;
	vector<double> r;
	vector<double> direct;
	
	double direction;
	double xc,yc;
	double area;
	unsigned int id;
	double c,s;
	double rx,ry;
	double rmax=0.0;
	
	for( unsigned int i=sys_->lctrl().size() ; i< sys_->spl()->lbody().size() ;++i)
	rmax=max(rmax,sys_->spl()->body(i)->sizeVerlet());
	cout<<"dmax="<<rmax<<endl;
	direct.clear();
	aspect.clear();
	r.clear();
	for(unsigned int i=0; i<ClusterCD.size();i++)
	//if (ClusterCD[i].size()>2)
	{
		//Chercher barycenter du cluster
		xc=0.;yc=0;
		area=0.;
		
		for(unsigned int j=0;j<ClusterCD[i].size();j++)
		{
			id=ClusterCD[i][j];
			xc+=sys_->spl()->body(id)->x()*sys_->spl()->body(id)->Area();
			yc+=sys_->spl()->body(id)->y()*sys_->spl()->body(id)->Area();
			area+=sys_->spl()->body(ClusterCD[i][j])->Area();
		}
		xc/=area;
		yc/=area;
		
		//Déterminer la sommet qui est le plus loin de barycenter
		rx=0.;
		for(unsigned int j=0;j<ClusterCD[i].size();j++)
		{
			id=ClusterCD[i][j];
			
			polyg * polygo= new polyg;//Attention reutiliser meme mémoire
			polygo=dynamic_cast<polyg*> (sys_->spl()->body(id));	
				
			c=cos(polygo->rot());
			s=sin(polygo->rot());
			double px,py;
			for(unsigned int k=0;k<polygo->Vertex().size();k++)
			{   
				px=polygo->x()+polygo->Vertex(k).x()*c-polygo->Vertex(k).y()*s;
				py=polygo->y()+polygo->Vertex(k).x()*s+polygo->Vertex(k).y()*c;
				if(rx<sqrt(pow(px-xc,2)+pow(py-yc,2)))
				{
					rx=sqrt(pow(px-xc,2)+pow(py-yc,2));
					direction=atan((py-yc)/(px-xc));
					if (direction<0.) direction+=M_PI;
				}
				
 			}
		}
		direct.push_back(direction);
		r.push_back(rx/rmax);
		//cout<<"rx:="<<rx<<"	direction:="<<direction/M_PI*180.<<endl;
		
		//Déterminer la sommet qui est le plus loin selon la direction perpendiculaire
		ry=0.;
		for(unsigned int j=0;j<ClusterCD[i].size();j++)
		{
			id=ClusterCD[i][j];
			
			polyg * polygo= new polyg;//Attention reutiliser meme mémoire
			polygo=dynamic_cast<polyg*> (sys_->spl()->body(id));	
				
			c=cos(polygo->rot());
			s=sin(polygo->rot());
			double px,py,distance;
			for(unsigned int k=0;k<polygo->Vertex().size();k++)
			{   
				px=polygo->x()+polygo->Vertex(k).x()*c-polygo->Vertex(k).y()*s;
				py=polygo->y()+polygo->Vertex(k).x()*s+polygo->Vertex(k).y()*c;
				distance=fabs(tan(direction)*(px-xc)-(py-yc))/sqrt(1+pow(tan(direction),2));
				if(ry<distance)
				{
					ry=distance;
				}				
 			}
		}
		aspect.push_back(rx/ry);
	}
	
	double a_min,a_max,a_mean;
	a_min=1.0e+6;
	a_mean=0.;
	a_max=0.;

	ofstream clust("Analyse/Cluster/cluster.txt",ios::out);
	for(unsigned int i=0;i<aspect.size();i++) clust<<aspect[i]<<endl;
	clust.close();
	
	for(unsigned int i=0;i<aspect.size();i++)
	{
		a_min=min(a_min,aspect[i]);
		a_max=max(a_max,aspect[i]);
		a_mean+=aspect[i];
	//	cout<<"aspect:="<<aspect[i]<<"	direction:="<<direct[i]/M_PI*180.<<endl;
	}
	a_mean/=double(aspect.size());
	cout<<"a_min:="<<a_min<<"	a_max:="<<a_max<<"	a_mean:="<<a_mean<<endl;
	
	double r_min,r_max,r_mean;
	r_min=1.0e+6;
	r_mean=0.;
	r_max=0.;

	ofstream clust_r("Analyse/Cluster/cluster_r.txt",ios::out);
	for(unsigned int i=0;i<r.size();i++) clust_r<<r[i]<<endl;
	clust_r.close();
	
	for(unsigned int i=0;i<r.size();i++)
	{
		r_min=min(r_min,r[i]);
		r_max=max(r_max,r[i]);
		r_mean+=r[i];
	}
	r_mean/=double(r.size());
	cout<<"r_min:="<<r_min<<"	r_max:="<<r_max<<"	r_mean:="<<r_mean<<endl;
	ofstream clust_rmax("Analyse/Cluster/cluster_rmax.txt",ios::out);
	clust_rmax<<" "<<r_min<<" "<<r_max<<" "<<r_mean<<endl;
	clust_rmax.close();
	
//********Calcul pdf des aspects
	char  fichier[100];
	DataSet f;
	for(unsigned int i=0;i<aspect.size();i++) f.add(aspect[i]);
	f.extractValues();
//	f.Normalize( f.mean());
	f.DecreasingSort();
	pointSet pdf = f.kernelPdf( 120,.01);
	sprintf(fichier,"Analyse/Cluster/pdf_aspect.txt");
	pdf.write(fichier);
	sprintf(fichier,"Analyse/Cluster/pdf_aspect_mob.txt");
	pdf.periodic()=false;
	pointSet pdfm = pdf.mobileMean2(1, 2);
	pdfm.write(fichier);

//**************
	

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny;
	unsigned int Nf=direct.size();
	
	for(unsigned int i=0;i<Nf;i++)
	{
		nx = cos(direct[i]);
		ny = sin(direct[i]);
		F->xx() += nx*nx;
		F->xy() += nx*ny;
		F->yy() += ny*ny;
	}
	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;
	F->yx() = F->xy();
	
	double a_cluster, da_cluster;
	da_cluster=0.;
	a_cluster=0.;
	if (F != NULL ) 
	{
		F->eigenValues();
		a_cluster=2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );
		da_cluster=F->majorDirection();
		cout<<" direction " << da_cluster/M_PI*180.<<" a_cluster= "<<a_cluster<<endl;	
	}
	delete F;
	
	ofstream clust_aspect("Analyse/Cluster/cluster_aspect.txt",ios::out);
		clust_aspect<<time<<" "<<a_min<<" "<<a_mean<<" "<<a_max<<" "<<da_cluster<<" "<<a_cluster<<endl;
	clust_aspect.close();
	unsigned int Nbin=120;
	vector <unsigned int > NcN(Nbin,0);
	double amp = M_PI/(double) (Nbin);
	unsigned int rangN;

	for(unsigned int i=0;i<Nf;i++)
	{
		rangN=(unsigned int) (floor(direct[i]/amp));
		if (rangN==Nbin) rangN--;
		NcN[ rangN ]++;
	}

	pointSet p;
	for ( unsigned int i =0;i<Nbin;++i) p.add( (.5+i)*amp, (double) (NcN[i])/(double) (Nf)*(double)(Nbin)/M_PI );
	for ( unsigned int i =0;i<Nbin;++i)	p.add( (.5+i)*amp+M_PI, (double) (NcN[i])/(double) (Nf)*(double)(Nbin)/M_PI );

	pointSet pm  =p.mobileMean(2,10);
	pm.write("ANALYSE/Cluster/Ptheta_mob.txt");
}
