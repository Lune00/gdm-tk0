#include "shearP_CD.hpp"

void shearP_CD::read_parameters(istream & is)
{
  string token;
  pressY_=false;
    symetrical_=false;
  firstUse_=true;
useSuper_=false;
  is >> token;	
  while(is)
    {	
    if      (token == "setUnit")              is >> Unit_;		
    else if (token == "topPlateThickness")    is >> topPlateThickness_;
    else if (token == "bottomPlateThickness") is >> bottomPlateThickness_;
    else if (token == "boundariesAuto")       boundariesAuto_ = true;
	else if (token == "dverlet")			  is >> dverlet_;
        else if (token == "symetrical")			  symetrical_ =true;
	else if (token == "ChangeGroup")		  changeGroup_=true;
	else if (token == "dsuperlist")
	{
		useSuper_ = true;
		is>>dsuperList_;
	}
	else if (token == "bandwidth")
	{
		is>>bandwidth_;
	}
    else if (token == "firstUse")
		{
			is>>token;
			if (token == "YES") firstUse_=true;
			else if (token == "NO") firstUse_=false;
		}
    else if (token =="topX")
      {	
		shearRate_=false;
			is >> token;
			if(token == "FORCE")          topXmode_ = _FORCE;
			else if (token == "VELOCITY") topXmode_ = _VELOCITY;
			else if (token == "SHEARRATE") 
				{
					topXmode_ = _VELOCITY;
					shearRate_=true;
				}
			is >> topXvalue_; 
      }
		
    else if (token =="topY")
      {	
			is >> token;
			if(token == "FORCE")          topYmode_ = _FORCE;
			else if (token == "VELOCITY") topYmode_ = _VELOCITY;
			else if (token == "PRESSURE") { topYmode_ = _FORCE; pressY_=true; }
			is >> topYvalue_;
      }
	else if (token=="Nanalyze")
			is>>Nanalyze_;
    
    else if (token == "}") break;
    else cerr << "@shearP_CD::read_parameters, Unknown keyword: " << token << endl;
    
    is >> token;
    }
}

void shearP_CD::write_parameters(ostream & os)
{

}

void shearP_CD::init() 
{  

	cout<<"--------------shearP_CD::init()-----------"<<endl;
	
	double masstot=0.;
	double gx=0,gy=0.;
	cout<<" angle "<<atan( gx/gy)/M_PI*180<<endl;
	unsigned int i;
	bool isperiodic=true;
	
		spl_->radiusExtrema(0);
		this->gx()=gx;
		this->gy()=gy;
		spl_->updateBoundaries();
		
		if( useSuper_) nwk_->useSuperList()=true;
		
		//if (boundariesAuto_)  spl_->definePeriodicityCV(true);
		
		if (Unit_ == "Rmax") 
			{
				topPlateThickness_    *= spl_->rmax();
				bottomPlateThickness_ *= spl_->rmax();
				nwk_->dverlet() =  dverlet_*spl_->rmax();
				spl_->bandWidth() = bandwidth_*spl_->rmax();
				
				if ( nwk_->useSuperList() )
				{
					nwk_->dsuperList() =  dsuperList_*spl_->rmax();
					nwk_->dsuperListP()=  dsuperListP_*spl_->rmax();
					nwk_->dsuperListP()=  nwk_->dsuperList();
				}
				if (boundariesAuto_)  spl_->definePeriodicityCV(true);
				//exit(1);
			}
    
		if (Unit_ == "Rmin") 
			{
				topPlateThickness_    *= spl_->rmin();
				bottomPlateThickness_ *= spl_->rmin();
				nwk_->dverlet() =  dverlet_*spl_->rmin();
				spl_->bandWidth() = bandwidth_*spl_->rmin();
				
				if ( nwk_->useSuperList() )
				{
					nwk_->dsuperList() =  dsuperList_*spl_->rmin();
					nwk_->dsuperListP()=  dsuperListP_*spl_->rmin();
					nwk_->dsuperListP()=  nwk_->dsuperList();
				}
				
				if (boundariesAuto_)  spl_->definePeriodicityCV(true);
				//exit(1);	
			}
			
			
						
		if (pressY_)
			{
			topYvalue_ *= - spl_->boundWidth();
			}
		
		cout<<"FIRST USE : "<<firstUse_<<endl;
		if ( firstUse_)
		{
		cout<<"-------FIRST USE------"<<endl;
		spl_->SortByHeight();
		//nwk_->dverlet() = alpha * spl_->rmoy();
		for ( i=0;i<spl_->lbody().size();++i)
			{
				spl_->body(i)->id()=i;
			}
			
		ldof().push_back( new dof(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.));
		//ldof(0)->m()=0.;
		ldof(0)->id()=0;
		ldof(0)->setGravity(gx,gy);
		for (i=0;i< spl_->lbody().size();++i)
			{
				if (spl_->body(i)->y() < spl_->ymin() + bottomPlateThickness_ )
				{
					ldof(0)->plugBody(spl_->body(i));
				//	spl_->body(i)->bodyDof() = ldof(0);
				//	ldof(0)->lctrlBodies().push_back(  spl_->body(i) );
				//	ldof(0)->m() += spl_->body(i)->mass();
					//cout<<ldof(0)->m()<<endl;
				}
			}
		cout<<endl;
		defineTopPlate2();
		ldof().push_back( new dof(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.));
		//ldof().push_back( new dof(topXmode_, _FORCE, _VELOCITY, topXvalue_, topYvalue_, 0.));//_______________MODIF
		//ldof(1)->m()=0.;
		ldof(1)->id()=1;
		ldof(1)->setGravity(gx,gy);
		for ( i=0;i<topPlate_.size();++i)
			{
				ldof(1)->plugBody(spl_->body(topPlate(i)));
				
				//spl_->body(topPlate(i))->bodyDof() = ldof(1);
				//ldof(1)->lctrlBodies().push_back(  spl_->body(topPlate(i)) );
				//ldof(1)->m()+= spl_->body(topPlate(i))->mass();
				//cout<<ldof(1)->m()<<endl;
			//cout<<ldof(1)->ctrlBody(i)->id()<<"  "<<ldof(1)->ctrlBody(i)<<"   "<<spl_->body(topPlate(i))->bodyDof()<<endl;
			}
		for( i=0;i< ldof().size();++i)
		{
			//ldof(i)->exportBodyId( spl_->includeFrom());
		//	ldof(i)->computeMassCenter();
		//	ldof(i)->computeMoment();
			ldof(i)->isPeriodic()=isperiodic;
		
		//	ldof(i)->print();
		}
		
		}
		else  //firstUse==false
		{
			cout<<"--------NOT FIRST USE-------"<<endl;
		/*	for ( i=0; i<spl_->lbody().size() ;++i)
			{
				spl_->body(i)->Fill( 2650 );
				masstot+=spl_->body(i)->mass();

				if (spl_->body(i)->id() != i ) cout<<"Pb id "<<i<<endl;
			}
		*/
		//	vector <unsigned int> idBodies;
			cout<<" taille dof "<<ldof().size()<<endl;
		
			ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.);
			ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
			
		/*	
				ldof(i)->m()=0.;
				ldof(i)->id()=i;
				ldof(i)->setGravity(gx,gy);
				
			//	ldof(i)->plugBody(i,spl_->includeFrom(),idBodies);
			//	ldof(i)->lctrlBodies().clear();
				for (unsigned int j=0 ;j < ldof(i)->lctrlBodies().size();++j)
				{
					//spl_->body(idBodies[j])->bodyDof() = ldof(i);
					//ldof(i)->lctrlBodies().push_back(  spl_->body(idBodies[j]) );
					ldof(i)->m() += ldof(i)->ctrlBody(j)->mass();
					//cout<<ldof(i)->ctrlBody(j)->id()<<endl;
				}
			*/
				for (i=0;i<ldof().size();++i)
				{
				//idBodies.clear();
			//	ldof(i)->computeMassCenter();
			//	ldof(i)->computeMoment();
				ldof(i)->isPeriodic()=isperiodic;
			//	ldof(i)->print();
			}
		}
		spl_->updateBands();//Id based
						
		//nwk_->dverlet() =  dverlet_*spl_->rmax();
		
		body2d * Rmin=NULL,*Rmax=NULL;
		
		for(unsigned int i=0;i<spl_->lbody().size();++i)
		{
			if ( spl_->body(i)->sizeVerlet() == spl_->rmin()) Rmin=spl_->body(i);
			if ( spl_->body(i)->sizeVerlet() == spl_->rmax()) Rmax=spl_->body(i);
			//cout<<spl_->body(i)->mass()<<endl;
		}
		if( Rmin==NULL) cout<<" ---------- rmin non retrouve"<<endl;
		if( Rmax==NULL) cout<<" ---------- rmax non retrouve"<<endl;
		
		cout<<endl<<"Periode ="<<spl_->boundWidth()<<endl;
		cout<<"Limite gauche = "<<spl_->leftBoundary()<<" lim droite ="<<spl_->rightBoundary()<<endl;
		cout<<scientific<<"Rmin = "<<spl_->rmin()<<" Rmax ="<<spl_->rmax()<<endl;
		cout<<"Rmax/Rmin = "<<spl_->rmax()/spl_->rmin()<<endl;
		cout<<"Taille de la liste de controle = "<<lctrl().size()<<endl;
    
        if(symetrical_ && shearRate_)
       {
        cout<<"Imposed symetrical shear rate: "<<topXvalue_<<endl;
           topXvalue_ *= ldof(1)->mcy()-ldof(0)->mcy();
           topXvalue_ *= 0.5;
           
        ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
        ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY,-topXvalue_, 0.,  0.);
        
       }
    
    
 //Impose une masse a la paroi libre egale à celle du système:
    
    double Mparoi=0.;
    
    for(unsigned int i=0;i<spl_->lbody().size();++i){
        
        if (spl_->body(i)->bodyDof() == NULL)  Mparoi += spl_->body(i)->mass();
           
        }
        
    ldof(1)->ComputeImposedMass(Mparoi);
    
		cout<<scientific<<"		Speed of topPlate : "<<topXvalue_<<endl;
        
        cout<<scientific<<"		Mass of topPlate : "<<Mparoi<<endl;
		
		cout<<"Parametre d'inertie:  I = "<< topXvalue_/( spl_->ymax()-spl_->ymin())*sqrt(-1.*masstot/topYvalue_/ spl_->boundWidth())<<endl;
		cout<<" I + petite particule : "<<topXvalue_/( spl_->ymax()-spl_->ymin())*sqrt(Rmin->mass()/(-1.*topYvalue_))<<endl;
		
		cout<<" I + grosse particule : "<<topXvalue_/( spl_->ymax()-spl_->ymin())*sqrt(Rmax->mass()/(-1.*topYvalue_))<<endl;

		if (topXvalue_ != 0. )
		cout<<"Durée de l'essais pour 50% de Xdef = "<<spl_->boundWidth()*.5/topXvalue_<<endl;
	//cout<<" taille echant "<<spl_->lbody().size()<<endl;
		
}

void shearP_CD::drive() 
{
	
			
}


void shearP_CD::defineTopPlate() 
{
	unsigned long i,j,Ndep=0,N,Nsuiv;
	vector <unsigned int> temp;
	vector <unsigned int> voisin;
	
	bool ChercheSuivant;
	double x,y,r,rot,rotmin;
	double d,dalert=spl_->rmin()*2.;//////////////////////////////WARNING=    pour le tres polydisp, a travailler
	
		for (i=0;i< spl_->lbody().size();++i)
		{
			if (spl_->body(i)->y() > spl_->ymax() - topPlateThickness_ )
			{
				//topPlate_.push_back(i);
				temp.push_back(i);
				if (spl_->body(i)->ymax()==spl_->ymax())   Ndep=i;
				cout<<"i ="<<i<<endl;
			}
		}
		cout<<"Particule de depart "<<Ndep<<endl;
		
		N=Ndep;
		topPlate_.push_back(Ndep);
		Nsuiv=1;
		ChercheSuivant=true;
		
		while( (Nsuiv!=Ndep) && ChercheSuivant)
		{
		x=spl_->body(N)->x();
		y=spl_->body(N)->y();
		r=spl_->body(N)->xmax()-spl_->body(N)->x();
		//cout<<"rapport taille "<< r/spl_->rmax()<<"  "<<r/spl_->rmin()<<endl;
		//cout<<"Voisin de "<<N<<" :"<<endl;
		for (i=0;i<temp.size();++i)
		{
			j=temp[i];
			d=sqrt( pow(spl_->body(j)->y()-y,2)+pow(spl_->body(j)->x()-x,2) )-r-(spl_->body(j)->xmax()-spl_->body(j)->x());
			//cout<<"j= "<<j<<" d= "<<d<<" dalert "<<dalert<<endl;
			if ( ( d< dalert) && (j!=N) && (spl_->body(j)->x()>=x))
			{
			voisin.push_back(j);
			//cout<<"		particule "<<j<<"  d= "<<d<<endl;
			}
		}
		
		if ( x> spl_->xmax()-5.*spl_->rmax() )
		{
			
			for (i=0;i<temp.size();++i)
			{
				j=temp[i];
				d=sqrt( pow(spl_->body(j)->y()-y,2)+pow(spl_->body(j)->x()+spl_->boundWidth()-x,2) )-r-(spl_->body(j)->xmax()-spl_->body(j)->x());
				//cout<<"j= "<<j<<" d= "<<d<<endl;
					if ( ( d< dalert) && (j!=N) && (spl_->body(j)->x()+spl_->boundWidth() >=x))
					{
					voisin.push_back(j);
					//cout<<"		particule periodique "<<j<<"  d= "<<d<<endl;
					}
			}
			
		}
		
		//Nsuiv=0;
		if (voisin.size()==1)
		{
			
			Nsuiv=voisin[0];
			//cout<<"		"<<Nsuiv<<" est le seul voisin de "<<N<<endl;
		}
		else 
		{
			rotmin=200;
			for (i=0;i<voisin.size();++i)
			{
				j=voisin[i];
				rot= atan( (spl_->body(j)->y()-y)/(spl_->body(j)->x()-x)) - M_PI/2.;
				//cout<<"rot avec "<<j<<"  "<< rot<<endl;
				if ( rot<rotmin)
				{
						Nsuiv=j;
				}
				//cout<<" Particule suivante "<<Nsuiv<<endl;
				
			}
		}
		voisin.clear();
		if (Nsuiv!=Ndep)
			{	
			topPlate_.push_back(Nsuiv);
			
			}
		else
			ChercheSuivant=false;
		cout<<"		Particule "<<Nsuiv<<" dans la surface  ///  Ndep = "<< Ndep<<endl;
		getchar();
		N=Nsuiv;
		}//while
	cout<<"		Fin defineTopPlate "<<endl;
}

void shearP_CD::defineTopPlate2() 
{
	unsigned long i;
	
		for (i=0;i< spl_->lbody().size();++i)
		{
			if (spl_->body(i)->y() > spl_->ymax() - topPlateThickness_ )
			{
				//topPlate_.push_back(i);
				topPlate_.push_back(i);
				
			}
		}
	//cout<<"		Fin defineTopPlate :"<<topPlate_.size()<<endl;
}


void shearP_CD::trans() 
{

}

void shearP_CD::share()
{
	/*
	double vref=30.*(pow(spl_->body(topPlate_.back())->vy(),2) + pow(spl_->body(topPlate_.back())->vx(),2));
	unsigned int Nstop=0;
	double k,v2;
	for( unsigned int i =0; i< spl_->lbody().size();++i)
	{
		if( spl_->body(i)->bodyDof()==NULL)
		{ 
			v2=spl_->body(i)->vy()*spl_->body(i)->vy()+spl_->body(i)->vx()*spl_->body(i)->vx();
		if ( v2 > vref)
		{
			k=vref/v2;
			spl_->body(i)->vx()*=k;
			spl_->body(i)->vy()*=k;
			//spl_->body(i)->vrot()=0.;
			Nstop++;
			//cout<<"arret de "<<i<<endl;
		}
		}
	}
	cout<<" @shearP_CD::share : Nstop = "<<Nstop<<endl; 
	*/
}

int shearP_CD::check()
{
  return 1;
}

/*void shearP_CD::little_analyse(double time)
{
	
}*/

void shearP_CD::stress_strain()
{}
