#include "dev_cluster.hpp"

void dev_cluster::read_parameters(istream & is)
{
	string token;
	pressY_=false;
	firstUse_=true;

	is >> token;	
	while(is)
	{	
		if      (token == "setUnit")              is >> Unit_;		
		else if (token == "topPlateThickness")    is >> topPlateThickness_;
		else if (token == "bottomPlateThickness") is >> bottomPlateThickness_;
		else if (token == "boundariesAuto")       boundariesAuto_ = true;
		else if (token == "kverlet")              is >> kverlet_;
		else if (token == "dverlet")			  is >> dverlet_;
		else if (token == "GRAVITY")			  {is>> grav;is >> angle_;}
		else if (token == "dsuperlist")
		{
			nwk_->useSuperList()=true;
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

			is >> token;
			if(token == "FORCE")          topXmode_ = _FORCE;
			else if (token == "VELOCITY") topXmode_ = _VELOCITY;
			is >> topXvalue_; 
		}
		else if (token =="topR")
		{	

			is >> token;
			if(token == "FORCE")          topRmode_ = _FORCE;
			else if (token == "VELOCITY") topRmode_ = _VELOCITY;
			is >> topRvalue_; 
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
		else cerr << "@dev_cluster::read_parameters, Unknown keyword: " << token << endl;

		is >> token;
	}
}

void dev_cluster::write_parameters(ostream & os)
{

}

void dev_cluster::init() 
{  

	cout<<"--------------dev_cluster::init()-----------"<<endl;

	double masstot=0.;

	double gx,gy;
	cout<<angle_<<endl;
	if( grav > 0)
	{
		gx= grav*sin(angle_/180.*M_PI);
		gy=-grav* cos(angle_/180.*M_PI);
	}
	else
		{gx=gy=0.;}
	unsigned int i;
	spl_->radiusExtrema(0);
	this->gx()=gx;
	this->gy()=gy;
		//this->setGravity(gx,gy);
	spl_->updateBoundaries();

	if (boundariesAuto_)  spl_->definePeriodicityCV(true);

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
			nwk_->dsuperListP()=nwk_->dsuperList();
		}
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
			nwk_->dsuperListP()=nwk_->dsuperList();
		}	
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
			spl_->body(i)->Fill( 2650 );
			masstot+=spl_->body(i)->mass();
		}

		ldof(0)->x()=_VELOCITY;
		ldof(0)->y()=_VELOCITY;
		ldof(0)->rot()=_VELOCITY;
		ldof(0)->xval()=0.;
		ldof(0)->yval()=0.;
		ldof(0)->rotval()=0.;		

		ldof(1)->x()   =topXmode_;
		ldof(1)->y()   =topYmode_;
		ldof(1)->rot() =topRmode_;
		ldof(1)->xval()=topXvalue_;
		ldof(1)->yval()=topYvalue_;
		ldof(1)->rotval()=topRvalue_;
		
		for( i=0;i< ldof().size();++i)
		{
			ldof(i)->m()=0;
			for(unsigned int j=0;j< ldof(i)->lctrlBodies().size();++j)
			{
				ldof(i)->m()+= ldof(i)->ctrlBody(j)->mass();
			}
			ldof(i)->id()=i;
			ldof(i)->setGravity(gx,gy);
			//ldof(i)->exportBodyId( spl_->includeFrom());
			ldof(i)->computeMassCenter();
			ldof(i)->computeMoment();
			ldof(i)->isPeriodic()=false;
			//ldof(i)->print();
			//ldof(i)->print();
		}

	}
	else  //firstUse==false
	{
		cout<<"--------NOT FIRST USE-------"<<endl;
		for ( i=0; i<spl_->lbody().size() ;++i)
		{
			spl_->body(i)->Fill( 2650 );
			masstot+=spl_->body(i)->mass();

			if (spl_->body(i)->id() != i ) cout<<"Pb id "<<i<<endl;
		}
		vector <unsigned int> idBodies;

			//ldof().push_back( new dof(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.));
		ldof(0)->x()=_VELOCITY;
		ldof(0)->y()=_VELOCITY;
		ldof(0)->rot()=_VELOCITY;
		ldof(0)->xval()=0.;
		ldof(0)->yval()=0.;
		ldof(0)->rotval()=0.;
			//ldof().push_back( new dof(topXmode_, topYmode_, topRmode_ , topXvalue_, topYvalue_, topRvalue_));
		ldof(1)->x()   =topXmode_;
		ldof(1)->y()   =topYmode_;
		ldof(1)->rot() =topRmode_;
		ldof(1)->xval()=topXvalue_;
		ldof(1)->yval()=topYvalue_;
		ldof(1)->rotval()=topRvalue_;

		for (i=0;i<ldof().size();++i)
		{
			ldof(i)->m()=0.;
			ldof(i)->id()=i;
			ldof(i)->setGravity(gx,gy);

				//ldof(i)->plugBody(i,spl_->includeFrom(),idBodies);

			for (unsigned int j=0 ;j < idBodies.size();++j)
			{
					//spl_->body(idBodies[j])->bodyDof() = ldof(i);
					//ldof(i)->lctrlBodies().push_back(  spl_->body(idBodies[j]) );
				ldof(i)->m() += spl_->body(idBodies[j])->mass();
					//cout<<spl_->body(idBodies[j])->id()<<endl;
			}
				//idBodies.clear();

		}
	}
	
	spl_->updateBands();//Id based

		//nwk_->dverlet() =  dverlet_*spl_->rmax();
	cout<<" taille dof "<<ldof().size()<<endl;

	cout<<endl<<"Periode ="<<spl_->boundWidth()<<endl;
	cout<<"Limite gauche = "<<spl_->leftBoundary()<<" lim droite ="<<spl_->rightBoundary()<<endl;
	cout<<scientific<<"Rmin = "<<spl_->rmin()<<" Rmax ="<<spl_->rmax()<<endl;
	cout<<"Rmax/Rmin = "<<spl_->rmax()/spl_->rmin()<<endl;
	cout<<"Taille de la liste de controle = "<<lctrl().size()<<endl;

	cout<<"Parametre d'inertie:  I = "<< topXvalue_/( spl_->ymax()-spl_->ymin())*sqrt(-1.*masstot/topYvalue_/ spl_->boundWidth())<<endl;

}

void dev_cluster::drive() 
{


}


void dev_cluster::defineTopPlate() 
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

void dev_cluster::defineTopPlate2() 
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


void dev_cluster::trans() 
{

}

void dev_cluster::share()
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
	cout<<" @dev_cluster::share : Nstop = "<<Nstop<<endl; 
	*/
}

int dev_cluster::check()
{
	return 1;
}

void dev_cluster::stress_strain()
{}

