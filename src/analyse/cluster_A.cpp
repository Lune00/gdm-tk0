#include "cluster_A.hpp"

void cluster_A::buildListinterdof(Probe & prb)//construction de la liste des dofs

{
	linterdof_.clear();

	unsigned int N=sys_->nwk()->clist().size();
	unsigned int i, k;
	bool exist=false;
	dof* first;
	dof* second;

	for (i=0; i<N; ++i)
	{
		if(prb.contain(sys_->nwk()->inter(sys_->nwk()->clist(i)))  &&  sys_->nwk()->inter(sys_->nwk()->clist(i))->type()!= _type_dkrl)//contact rlinedisk exclu
		{
			first = sys_->nwk()->inter(sys_->nwk()->clist(i))->first()->bodyDof();//premier dof en contact
			second = sys_->nwk()->inter(sys_->nwk()->clist(i))->second()->bodyDof();//second dof en contact
			if(first==NULL || second==NULL)
			{
				continue;
		}

			exist = false;
			for (unsigned int j=0;j<linterdof_.size();++j)
			{

				if (first==linterdof_[j]->first())
				{
					if (second==linterdof_[j]->second())
					{
					exist = true;
					k=j;
					break;
					}

				}
				else if (first==linterdof_[j]->second())
				{
					if (second==linterdof_[j]->first())
					{
					exist=true;
					k=j;
					break;
					}
				}
				
				//else if (second==linterdof_[j]->first())
				//{
				//	if (first==linterdof_[j]->second())
				//	{
				//	exist = true;
				//	k=j;
				//	break;
				//	}

				//}
				//else if (second==linterdof_[j]->second())
				//{
				//	if (first==linterdof_[j]->first())
				//	{
				//	exist=true;
				//	k=j;
				//	break;
				//	}
				//}

			}		


			if (exist) 
			{
				linterdof_[k]->linter().push_back(sys_->nwk()->inter(sys_->nwk()->clist(i)));
				//cout<<" rank ins  "<<linterdof_[k]->rank()<<endl;
				
			}
			else 
			{
				linterdof_.push_back(new interdof(first,second));
				linterdof_.back()->linter().push_back(sys_->nwk()->inter(sys_->nwk()->clist(i)));
				//cout<<" rank "<<linterdof_.back()->rank()<<endl;

			}
		}
	}
	
	
	activatedDof.resize(sys_->ldof().size(),false);
	for (unsigned int j=0;j<linterdof_.size(); ++j)
	{ 
		activatedDof[linterdof_[j]->first()->id()]=true;
		activatedDof[linterdof_[j]->second()->id()]=true;
		
	}
	Ndof_=0;
	for (unsigned int j=0;j<activatedDof.size(); ++j)
	{ 
		if (activatedDof[j]) Ndof_++;
	}
	
	cout<<"Ndof = "<<Ndof_<<endl;//nombre de dof en contacts
	
	ofstream branchLengthC("Analyse/Cluster/BranchC.txt",ios::out);
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		linterdof_[i]->Frame();
		linterdof_[i]->calcFres();
		branchLengthC<<" "<<linterdof_[i]->vbranch()<<endl;
	}
	branchLengthC.close();
	
	
//debuggage par rapport aux vitesses des corps controlés par un dof
	/*for (unsigned int j=0; j<sys_->ldof().size(); ++j)
	{
		double VX,VY;
		VX=sys_->ldof(j)->mcvx();
		VY = sys_->ldof(j)->mcvy();
		cout<< " dof: " <<j<<" v_0: "<<sqrt(VX*VX+VY*VY)<<endl;
	for(unsigned int k = 0;  k<sys_->ldof(j)->lctrlBodies().size(); ++k)
		{
			if(sys_->ldof(j)->ctrlBody(k)->type() == _type_dkdk)
			{
				double VBX,VBY;
				VBX = sys_->ldof(j)->ctrlBody(k)->vx();
				VBY = sys_->ldof(j)->ctrlBody(k)->vy();
				cout<<" "<<VBX*VBX+VBY*VBY<<endl;
			}
		}
		
	}*/
}



//le 26/04/09 B.SAINT-CYR

//----------Fonction donnant la proportion de chaque type et contribution de chaque type de contacts----------------------------
vector<double>cluster_A::contactTypes()
{
	unsigned int Nmax_=10;
	vector<unsigned int>V(Nmax_,0);
	
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		V[linterdof_[i]->rank()]++;//Classe le nombre de particules qui ont k contacts

	}

	vector<double>XType_(V.size());
	
	for (unsigned int i =0 ; i<V.size() ; ++i)
	{XType_[i] = (double) V[i]/(double) linterdof_.size();}//proportion de chaque type de contacts
	
	return XType_;
}

/*vector<double> cluster_A::contactDouble()
{
	vector<interdof *>lcontactDouble_(linterdof_.size());
	//interdofk = new interdof(i_,j_);
	
	//attention il faut allouer de l'espace mémoire qu'il faudra libérer
	
	
	vector<inter2d*>Temp(2); //idem
	vector<double>XDouble_(2,0);
	
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		if(linterdof_[i]->rank() == 2)
		{
			lcontactDouble_[i] = new interdof(i_,j_);
			lcontactDouble_[i] = linterdof_[i];
			//interdofk = linterdof_[i];
			//lcontactDouble_.push_back(interdofk);
		}
		
	}
	
	for (unsigned int j =0 ; j<lcontactDouble_.size() ; ++j)
	{
		delete lcontactDouble_[j];
	}
	//delete interdofk;
	for (unsigned int i =0 ; i<lcontactDouble_.size() ; ++i)
	{
		Temp[0] = (lcontactDouble_[i]->linter())[0];
		Temp[1] = (lcontactDouble_[i]->linter())[1];
		

		if (Temp[0]->first()->id() == Temp[1]->first()->id() || Temp[0]->first()->id() == Temp[1]->second()->id() || Temp[0]->second()->id() == Temp[1]->first()->id() || Temp[0]->second()->id() == Temp[1]->second()->id())
		{
			XDouble_[0]++;//contacts Doubles
		}
		else
		{
			XDouble_[1]++;//contacts DoublesSimples formes de 2 simples
		}
	}
	XDouble_[0] /= (double) linterdof_.size();
	XDouble_[1] /= (double) linterdof_.size();
	
	return XDouble_;
}*/

void cluster_A::MeanForce()
{
	unsigned int N = 10;
	vector<DataSet>MeanRadialForce_(N);
	vector<DataSet>MeanOrthoForce_(N);
	DataSet MeanRadial_all_;
	DataSet MeanOrtho_all_;
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		MeanRadialForce_[linterdof_[i]->rank()].add(linterdof_[i]->fresNormal());//Stockage des forces de contacts par type de contacts
		MeanOrthoForce_[linterdof_[i]->rank()].add(linterdof_[i]->fresOrthoNorm());
		MeanRadial_all_.add(linterdof_[i]->fresNormal());//Stockage des forces radiales en globalité
		MeanOrtho_all_.add(linterdof_[i]->fresOrthoNorm());
	}
	for (unsigned int j=0 ; j<N ; ++j)
	{
		MeanRadialForce_[j].extractValues();
		MeanOrthoForce_[j].extractValues();
	}
		
	MeanRadial_all_.extractValues();
	MeanOrtho_all_.extractValues();
	
	for(unsigned int k=0 ; k<N ; ++k)
	{
		MeanRadial_.push_back(MeanRadialForce_[k].mean()/MeanRadial_all_.mean());
		MeanOrtho_.push_back(MeanOrthoForce_[k].mean()/MeanOrtho_all_.mean());
	}
	
}

/*void cluster_A::ContactForceRatio()
{
	double fglobal = 0;
	double fsimple = 0;
	double fdouble = 0;
	double ftriple =0;
	unsigned int N=sys_->nwk()->clist().size();

	for(unsigned int i = 0; i<N; ++i)
	{
		if(sys_->nwk()->inter(sys_->nwk()->clist(i))->type() != _type_dkrl)
		{
			fglobal += sys_->nwk()->inter(sys_->nwk()->clist(i))->fx();
		}
	}
	
	for(unsigned int i =0; i<linterdof_.size(); ++i)
	{
		if (linterdof_[i]->rank() == 2)
		{
			for( int j=0; j<linterdof_[i]->rank(); ++j)
			{
				fdouble += (linterdof_[i]->linter())[j]->fx();
			}
		}
		if (linterdof_[i]->rank() == 3)
		{
			for( int j=0; j<linterdof_[i]->rank(); ++j)
			{
				ftriple += (linterdof_[i]->linter())[j]->fx();
			}
		}
		if (linterdof_[i]->rank() == 1)
		{
			for( int j=0; j<linterdof_[i]->rank(); ++j)
			{
				fsimple += (linterdof_[i]->linter())[j]->fx();
			}
		}
	}
	
	Xforce_simple_ = fsimple;
	Xforce_double_ = fdouble;
	Xforce_triple_ = ftriple;
	
	cout<<" " <<Xforce_simple_<<" " <<Xforce_double_<<" " <<Xforce_triple_<<" " <<fglobal<<endl;
}*/

vector<double> cluster_A::connectivity()
{
	vector<unsigned int>Ncp(sys_->ldof().size(),0);
	
	for (unsigned int j=0;j<linterdof_.size();++j)
	{
		Ncp[linterdof_[j]->first()->id()]++;//compte le nombre de voisins par particules
		Ncp[linterdof_[j]->second()->id()]++;
		
	}

	

	vector<unsigned int>Npkc(30,0);
	for (unsigned int j=0;j<sys_->ldof().size();++j)
	{
		if(activatedDof[j])
		Npkc[Ncp[j]]++;//compte le nombre de particules qui ont k voisins
	}
	
	/*for(unsigned int j = 0 ; j<Npkc.size() ; ++j)
	{
		cout<<" " <<Npkc[j]<<endl;
	}*/
	vector<double> Xpkc(30);
	for (unsigned int j=0;j<Npkc.size();++j)
	{
		Xpkc[j]= (double)Npkc[j]/(double)Ndof_;//proportion des particules qui ont k voisins
	}

	return Xpkc;
}

void cluster_A::contact_connectivity()
{
	unsigned int Nlist = sys_->nwk()->clist().size();
	unsigned int Nb = sys_->ldof().size();
	vector < unsigned int> Ncpp(Nb,0);
	unsigned int id1,id2;
	for (unsigned int i =0 ; i<linterdof_.size() ; ++i)
	{
		id1=linterdof_[i]->first()->id();
		id2=linterdof_[i]->second()->id();
		
		

		
		
			Ncpp[id1]+= linterdof_[i]->rank();
			Ncpp[id2]+= linterdof_[i]->rank();
		

	}
	unsigned int Ncrl = 0;
	for (unsigned int j =0 ; j<Nlist ; ++j)
	{
		if(sys_->nwk()->inter(sys_->nwk()->clist(j))->type() == _type_dkrl)
		{
			Ncrl += Ncrl;
		}
	}
	unsigned int Nc = 0;
	unsigned int Nwc_ = 0;
	for (unsigned int k = 0 ; k<Ncpp.size() ; ++k)
	{
		if(Ncpp[k] > 1 && activatedDof[k] == true)
		{
			Nc += Ncpp[k];
			Nwc_+=1;
		}
	}
	
	Z = (double)(Nc + Ncrl)/(double)Ndof_;
	cout<<"Z:   "<<Z<<endl;//------------------------------------------------print
	
	unsigned int Nf_ = 0;
	for(unsigned int j =0; j<sys_->ldof().size(); ++j)
	{
		if(  activatedDof[j] == false && Ncpp[j]<= 1)
		Nf_ +=1;
	}
	Xratlers_ = (double)Nf_/(double)(Nf_+Nwc_);
	cout<<" ratlers: " <<Xratlers_<<endl;//----------------------------------------print
}


vector<double> cluster_A::rankStat()
{
	vector <double> Xk(10,0);
	for(unsigned int j=0;j<linterdof_.size();j++)
	{
		Xk[linterdof_[j]->rank()]+=1;
	}
	for(unsigned int j=0;j<Xk.size();j++)
	{
		Xk[j]/=(double) linterdof_.size();
	}
	
	return Xk;
}

	
/*void cluster_A::PDFFN_Cluster()
{
	char * fichier;
	
	fichier = "Analyse/Cluster/pdf/PdfFR.txt";

	DataSet Data;
	
	
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		Data.add(linterdof_[i]->fresNormal());		
	}	
		
	Data.extractValues();
	Data.Normalize( Data.mean() );
	

	pointSet pdf = Data.Rich_PDF2( Data.setSize()/500);
	pdf.write( fichier);
	
	cout <<"OK "<<endl;*/
	
//---------------------------------Partie à revoir----------------------
	
	//pointSet pdfs_all = pdf_all.mobileMean(3,1);
	//pointSet pdfs_simp = pdf_simp.mobileMean(3,1);
	//pointSet pdfs_nonsimp = pdf_nonsimp.mobileMean(3,1);
	
	//fichier_all = "Analyse/Cluster/pdf/PdffNslid_all.txt";
	//fichier_simp = "Analyse/Cluster/pdf/PdffNslid_simp.txt";
	//fichier_nonsimp = "Analyse/Cluster/pdf/PdffNslid_nonsimp.txt";

	//pdfs_all.write(fichier_all);

	//pdfs_simp.write(fichier_simp);
	//pdfs_nonsimp.write(fichier_nonsimp);

//-------------------------------------------------------------------				


/*void cluster_A::PDFFT_Cluster()
{
	char * fichier_simp, *fichier_nonsimp,*fichier_all;
	fichier_all = "Analyse/Cluster/pdf/PdffT_all.txt";
	fichier_simp = "Analyse/Cluster/pdf/PdffT_simp.txt";
	fichier_nonsimp = "Analyse/Cluster/pdf/Pdfft_nonsimp.txt";
	
	vector<DataSet>Data(2);
	DataSet Data_all;
	
	for (unsigned int i=0 ; i<linterdof_.size() ; ++i)
	{
		Data_all.add(fabs(linterdof_[i]->fresOrthoNorm()));
		if (linterdof_[i]->rank() == 1)
		{Data[0].add(fabs(linterdof_[i]->fresOrthoNorm()));}
		
		if (linterdof_[i]->rank() !=1)
		{Data[1].add(fabs(linterdof_[i]->fresOrthoNorm()));}
		
	}
	
	for (unsigned int j = 0; j<Data.size() ; ++j)
	{
		//if (Data[j].setSize()<100)
		//cout <<"No more forces : "<<endl;
		Data[j].extractValues();
		Data[j].Normalize( Data[j].mean() );
		cout<<"fmoy_tang = "<<Data[j].mean()<<endl;
	}
	
	Data_all.extractValues();
	Data_all.Normalize( Data_all.mean());
	
	pointSet pdf_all = Data_all.Rich_PDF2( Data_all.setSize()/100);
	pointSet pdf_simp = Data[0].Rich_PDF2( Data[0].setSize()/100);
	pointSet pdf_nonsimp = Data[1].Rich_PDF2( Data[1].setSize()/100);
	
	pdf_all.write( fichier_all);
	pdf_simp.write( fichier_simp);
	pdf_nonsimp.write( fichier_nonsimp);
	
//----------------------------------Partie à revoir----------------------
	pointSet pdfs_all = pdf_all.mobileMean(1,3);
	pointSet pdfs_simp = pdf_simp.mobileMean(1,3);
	pointSet pdfs_nonsimp = pdf_nonsimp.mobileMean(1,3);
	
	fichier_all = "Analyse/Cluster/pdf/PdffTslid_all.txt";
	fichier_simp = "Analyse/Cluster/pdf/PdffTslid_simp.txt";
	fichier_nonsimp = "Analyse/Cluster/pdf/PdffTslid_nonsimp.txt";

	pdfs_all.write(fichier_all);

	pdfs_simp.write(fichier_simp);


	pdfs_nonsimp.write(fichier_nonsimp);

//-------------------------------------------------------------------
}*/


/*void cluster_A::Ptheta_Cluster(unsigned int Nbin,unsigned int period,unsigned int width) 
{
	
	cout<<"	Ptheta  : ";
	//cout<<"NbinPT = "<<Nbin<<endl;

	vector <unsigned int > NcN(Nbin,0);
	vector <unsigned int > NcT(Nbin,0);

	double fresNormalmoy=0;
	double fresOrthomoy=0;
	double angN;
	double amp = M_PI/(double) (Nbin);//amplitude de chaque classe

	unsigned int NctotN=0,rangN,NctotT=0;
	
	for( unsigned int c=0;c<linterdof_.size();++c)
	{

		if( linterdof_[c]->fresNormal() !=0.)
		{
			angN= acos( linterdof_[c]->nx() ) ;//orientation des vecteurs branches 
			if(linterdof_[c]->ny()< 0.) angN = M_PI - angN;//orientation des vecteurs perpendiculaires aux vecteurs branches

			rangN=(unsigned int) ( floor( angN/amp));// calcul de la densité
			if (rangN==Nbin) rangN--;


			NcN[ rangN ]++;//remplissage de chaque classe avec son effectif
			NctotN++;
			
			fresNormalmoy+=linterdof_[c]->fresNormal();

			if( linterdof_[c]->fresOrthoNorm() != 0. )
			{
				NcT[ rangN ]++;
				fresOrthomoy += linterdof_[c]->fresOrthoNorm();
				NctotT++;
			}
		}

	
	}
	fresNormalmoy/=(double) (NctotN);
	fresOrthomoy/=(double) (NctotT);

	pointSet p;

	for ( unsigned int i =0;i<Nbin;++i)
	{
		p.add( (.5+i)*amp, (double) (NcN[i])/(double) (NctotN)*(double)(Nbin)/M_PI );

	}

	pointSet pm  =p.mobileMean(  period,width);

	
	pm.write("Analyse/Cluster/ptheta/Ptheta_mob.txt");
	
	cout<<" ok "<<endl;

}*/

void cluster_A::Stress_Cluster(Probe & prb)
{

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress");
	double dx,dy,fx,fy;
	unsigned int nc=0;
	

	
	for (unsigned int c = 0; c < linterdof_.size(); ++c)
	{
		if(linterdof_.size() != 0)
		{
		if( linterdof_[c]->fresNormal() != 0. )
		{	
			if(prb.containDofCenter(linterdof_[c]->first()) || prb.containDofCenter(linterdof_[c]->second()) )
			{
				
				dx = linterdof_[c]->vbranchx();
				dy = linterdof_[c]->vbranchy();
				
				fx = linterdof_[c]->fresNormalx();
				fy = linterdof_[c]->fresNormaly();

				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		
		}
		}
		else
		{
			cout<<"pas de listinterdof"<<endl;
		}
		
	}
	
	
	S->scalarMult( 1./prb.area() );
	S->eigenValues();
	double s1=S->l1();
	double s2=S->l2();
	sigma_=( max(s1,s2)-min(s1,s2) )/(s1+s2);
	cout<<"qoverp:  " <<sigma_<<endl;//------------------------------------------------print
	
	

	delete S;

}

void cluster_A::Stress_ClusterSimple(Probe & prb)
{

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress");
	double dx,dy,fx,fy;
	unsigned int nc=0;
	

	
	for (unsigned int c = 0; c < linterdof_.size(); ++c)
	{
		if(linterdof_.size() != 0)
		{
		if( linterdof_[c]->fresNormal() != 0. && linterdof_[c]->rank()==1.)
		{	
			if(prb.containDofCenter(linterdof_[c]->first()) || prb.containDofCenter(linterdof_[c]->second()) )
			{
				
				dx = linterdof_[c]->vbranchx();
				dy = linterdof_[c]->vbranchy();
				
				fx = linterdof_[c]->fresNormalx();
				fy = linterdof_[c]->fresNormaly();

				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		
		}
		}
		else
		{
			cout<<"pas de listinterdof"<<endl;
		}
		
	}
	
	
	S->scalarMult( 1./prb.area() );
	S->eigenValues();
	double s1=S->l1();
	double s2=S->l2();
	sigmaSimple_=( max(s1,s2)-min(s1,s2) )/(s1+s2);
	
	

	delete S;

}

void cluster_A::Stress_ClusterNsimple(Probe & prb)
{

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress");
	double dx,dy,fx,fy;
	unsigned int nc=0;
	

	
	for (unsigned int c = 0; c < linterdof_.size(); ++c)
	{
		if(linterdof_.size() != 0)
		{
		if( linterdof_[c]->fresNormal() != 0. && linterdof_[c]->rank()!=1.)
		{	
			if(prb.containDofCenter(linterdof_[c]->first()) || prb.containDofCenter(linterdof_[c]->second()) )
			{
				
				dx = linterdof_[c]->vbranchx();
				dy = linterdof_[c]->vbranchy();
				
				fx = linterdof_[c]->fresNormalx();
				fy = linterdof_[c]->fresNormaly();

				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		
		}
		}
		else
		{
			cout<<"pas de listinterdof"<<endl;
		}
		
	}
	
	
	S->scalarMult( 1./prb.area() );
	S->eigenValues();
	double s1=S->l1();
	double s2=S->l2();
	sigmaNsimple_=( max(s1,s2)-min(s1,s2) )/(s1+s2);
	
	

	delete S;

}

/*void cluster_A::Fabric_Cluster(Probe & prb)//fabrique de texture de branche
{
	
	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	

	for (unsigned int c = 0; c < linterdof_.size(); ++c)
	{		
		if(prb.containDofCenter(linterdof_[c]->first()) || prb.containDofCenter(linterdof_[c]->second()))
		{
			ang= acos( linterdof_[c]->nx() ) ;
			if(linterdof_[c]->ny()< 0.) ang=M_PI-ang;

			nx = cos(ang);
			ny = sin(ang);

			F->xx() += nx*nx;
			F->xy() += nx*ny;
			F->yy() += ny*ny;
			Nf++;
		}
		
	
	}

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();
	F->eigenValues();
	ac_ = 2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );

	delete F;

}*/

void cluster_A::frAnisoInProbe(Probe & prb)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("F");
	gdm::Tensor2x2 * G = new gdm::Tensor2x2("G");
	double nx,ny,nxt,nyt,fn,ft,ang,fmoy=0;
	unsigned int Nf=0;

	for ( unsigned int c = 0; c <  linterdof_.size(); ++c)
	{

		if (linterdof_[c]->fresNormal() != 0.0 && linterdof_[c]->fresOrthoNorm() !=0.0)
		{
			if(prb.containDofCenter(linterdof_[c]->first()) || prb.containDofCenter(linterdof_[c]->second()) )
			{
				ang= acos( linterdof_[c]->nx() ) ;
				if(linterdof_[c]->ny()< 0.) ang = M_PI - ang;


				nx = linterdof_[c]->nx();
				ny = linterdof_[c]->ny();
				
				nxt = cos(ang);
				nyt = sin(ang);
				
				fn = linterdof_[c]->fresNormal();
				ft = linterdof_[c]->fresOrthoNorm();

				fmoy+= fn;

				F->xx() += fn*nx*nx;
				F->xy() += fn*nx*ny;
				F->yx() += fn*nx*ny;
				F->yy() += fn*ny*ny;
				
				G->xx() -= ft*nyt*nxt;
				G->yx() += ft*nxt*nxt;
				G->xy() -= ft*nyt*nyt;
				G->yy() += ft*nyt*nxt;

				Nf++;
			}

		}
	
	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
		G->scalarMult(1./fmoy);
		G->scalarMult(1./(double) (Nf));
		

	}
	if (F!=NULL && G!=NULL)
	{
		F->eigenValues();
		an_ = 2.*( max(F->l2(),F->l1()) - min(F->l2(),F->l1()) )/(F->l1() + F->l2() );
		an_=an_-ac_;
		*G = *F + *G;
		if (G->xx()!=0. && G->yy()!=0.)
		{
			G->eigenValues();
			at_=2.*( max(G->l2(),G->l1()) - min(G->l2(),G->l1()) )/(G->l1() + G->l2() );
			at_-=an_+ac_;
			
		}

	}
	
	delete F;
	delete G;

}	

void  cluster_A::ContactFabric(Probe & prb) //fabrique de texture de normale au contact
{
	unsigned int Nc = sys_->nwk()->clist().size();
	if (Nc == 0)  {cout<<"liste de contacts vide"<<endl;}//return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	inter2d * interc = NULL;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(c));
		//cout<<" adresse: " <<& interc<<endl;
		if (interc->fn() != 0.0)
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;
				nx = cos(ang);
				ny = sin(ang);

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}
	
	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();
	F->eigenValues();
	ac_ = 2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );

	delete F;

}


void  cluster_A::simpleContactFabric(Probe & prb) //fabrique de texture de normale au contact
{
	unsigned int Nc = sys_->nwk()->clist().size();
	if (Nc == 0)  {cout<<"liste de contacts vide"<<endl;}//return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	inter2d * interc = NULL;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(c));
		if (interc->fn() != 0.0  && interc->rang()==1)
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;
				nx = cos(ang);
				ny = sin(ang);

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}
	
	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();
	F->eigenValues();
	acs_ = 2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );

	delete F;

}

void  cluster_A::NonsimpleContactFabric(Probe & prb) //fabrique de texture de normale au contact
{
	unsigned int Nc = sys_->nwk()->clist().size();
	if (Nc == 0)  {cout<<"liste de contacts vide"<<endl;}//return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	inter2d * interc = NULL;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = sys_->nwk()->inter(sys_->nwk()->clist(c));
		if (interc->fn() != 0.0  && interc->rang()!=1)
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;
				nx = cos(ang);
				ny = sin(ang);

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}
	
	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();
	F->eigenValues();
	acns_ = 2.*(max(F->l2(),F->l1())- min(F->l2(),F->l1()) );

	delete F;

}







