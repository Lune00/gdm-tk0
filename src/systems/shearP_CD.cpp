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
		else if (token == "gravite")
		{
			gravite_ =true;
			is >> multig_;
		}
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

	cout<<"------oooo------shearP_CD::init()------oooo----"<<endl;

	//Rajouter gravite dans les options
	//gravite value
	double gx, gy ;
	if(gravite_)
	{
		gx = 0. ;
		gy = -9.81 * multig_ ;
	}
	else
	{
		gx = 0. ;
		gy = 0. ;
	}

	//Calcul rmin_, rmax_ et rmoy_ de l'échantillon
	spl_->radiusExtrema(0);
	this->gx()=gx;
	this->gy()=gy;

	//Rescale limit indicated by default (4rmax), rmin or rmax
	SetUnity();

	//Définition des limites et des conditions périodiques (taille et bande)
	spl_->updateBoundaries(); // min,max

	//Il faudra rajouter éventuellement la possibilité de voir xmin_ et xmax_ défini par l'utilisateur
	//Pour l'instant c'est automatiquement défini a partir des limites

	if (boundariesAuto_)
	{
		spl_->definePeriodicity(bandwidth_);
	}
	else
	{
		spl_->definePeriodicity(bandwidth_);
	}
	//spl_->updateBands();//Id based

	if (useSuper_) nwk_->useSuperList()=true;
	if (pressY_)   topYvalue_ *= - spl_->boundWidth();


	// firstUse_ TRUE à TESTER (pas encore tester)
	if ( firstUse_)
	{
		cout<<".Génération des dofs : FIRST USE -> génération automatique des dofs"<<endl;
		spl_->SortByHeight();
		for ( unsigned int i=0;i<spl_->lbody().size();++i)
		{
			spl_->body(i)->id()=i;
		}

		ldof().push_back( new dof(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.));
		ldof(0)->id()=0;
		ldof(0)->setGravity(gx,gy);

		for (unsigned int i=0;i< spl_->lbody().size();++i)
		{
			if (spl_->body(i)->y() < spl_->ymin() + bottomPlateThickness_ )
			{
				ldof(0)->plugBody(spl_->body(i));
			}
		}
		defineTopPlate2();
		ldof().push_back( new dof(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.));

		ldof(1)->id()=1;
		ldof(1)->setGravity(gx,gy);
		for (unsigned int i=0;i<topPlate_.size();++i)
		{
			ldof(1)->plugBody(spl_->body(topPlate(i)));
		}
		for(unsigned int i=0;i< ldof().size();++i)
		{
			ldof(i)->isPeriodic()=true;
		}

	}
	else  //firstUse==false
	{
		cout<<".Génération des dofs : NOT FIRST USE -> dof pré-générés"<<endl;
		cout<<".Nombre de dofs : "<<ldof().size()<<endl;

	}

	//Imposer shear rate : H particule plus haute - particule plus basse
	if (shearRate_)
	{
		//Store initial shearRate
		shearRate_init=topXvalue_;
		cout<<".Imposed shear rate : "<<topXvalue_<<endl;

		double h = ldof(1)->lowerBody()->y() - ldof(0)->lowerBody()->y();
		// h*=-1.;
		cout<<".Sample thickness="<<h<<endl;
		topXvalue_ *= h;//ldof(1)->mcy()-ldof(0)->mcy();
		//On stocke en mémoire le taux de déformation imposé (que l'on pourra remettre a jour si l'épaisseur de l'échantillon varie) 
	}
	//Les deux plaques sont en mouvement, le taux de déformation est donc égal a 2v/h, on divise par deux pour avoir le taux de def utilisateur donné 
	if(symetrical_)
	{
		cout<<".Imposed symetrical shear rate: "<<topXvalue_<<endl;
		//topXvalue_ *= ldof(1)->mcy()-ldof(0)->mcy();
		topXvalue_ *= 0.5;

		ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
		ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY,-topXvalue_, 0.,  0.);
	}
	else
	{

		cout<<".Imposed non-symetrical shear rate: "<<topXvalue_<<endl;
		ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
		ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.);
	}


	//Periodise les dofs:
	for (unsigned int i=0;i<ldof().size();++i)
	{
		ldof(i)->isPeriodic()=true;
	}
	spl_->updateBands();//Id based

	cout<<endl<<".Taille zone périodique (d)="<<spl_->boundWidth()/(spl_->rmin()+spl_->rmax())<<endl;
	cout<<".Taille des bandes periodiques (d) : "<<spl_->bandWidth()/(spl_->rmin()+spl_->rmax())<<endl;
	cout<<scientific<<"Rmin = "<<spl_->rmin()<<" Rmax ="<<spl_->rmax()<<endl;
	cout<<"Rmax/Rmin = "<<spl_->rmax()/spl_->rmin()<<endl;

	//Impose une masse a la paroi libre egale à celle du système:

	double Mparoi=0.;

	for(unsigned int i=0 ;i<spl_->lbody().size();i++){

		if (spl_->body(i)->bodyDof() == NULL)  Mparoi += spl_->body(i)->mass();

	}

	ldof(1)->ComputeImposedMass(Mparoi);

	cout<<scientific<<".Speed of topPlate : "<<topXvalue_<<endl;

	cout<<scientific<<".Mass of topPlate : "<<Mparoi<<endl;

	cout<<".Nombre de pas (pas de temps 0.0001) pour l'essais pour 100% de Xdef = "<<spl_->boundWidth()/(topXvalue_*0.0001)<<endl;

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

void shearP_CD::SetUnity()
{
	if(Unit_ == "default")
	{
		if (topPlateThickness_ == 0.)		{ topPlateThickness_ = 2.*spl_->rmax();		}
		if (bottomPlateThickness_ == 0.)	{ bottomPlateThickness_ = 2.*spl_->rmax();	}
		if (bandwidth_ == 0.)				{ bandwidth_ = 2.*spl_->rmax();				}
		if (dverlet_ == 0.)
		{
			dverlet_ = spl_->rmax();
			if (useSuper_)
			{
				if (dsuperList_ == 0. || dsuperListP_ == 0.)
				{
					dsuperList_		= 2.*spl_->rmax();
					dsuperListP_	= 2.*spl_->rmax();
				}
			}
		}
	}
	else if (Unit_ == "Rmax" || Unit_ == "Rmin")
	{
		int a = topPlateThickness_ ;
		cerr << "DBG TTTTTTTTTTTTTTTTTTTTTTTTTTTT a:" << a << endl ;
		if (topPlateThickness_ == 0.)		{ topPlateThickness_ = 2.;		}
		if (bottomPlateThickness_ == 0.)	{ bottomPlateThickness_ = 2.;	}
		if (bandwidth_ == 0.)				{ bandwidth_ = 2.;				}
		if (dverlet_ == 0.)
		{
			dverlet_ = 1.;
			if (useSuper_)
			{
				if (dsuperList_ == 0. || dsuperListP_ == 0.)
				{
					dsuperList_		= 2.;
					dsuperListP_	= 2.;
				}
			}
		}
	}

	if (Unit_ == "Rmax")
	{
		topPlateThickness_    *= spl_->rmax();
		bottomPlateThickness_ *= spl_->rmax();
		nwk_->dverlet() =  dverlet_*spl_->rmax();
		spl_->bandWidth() = bandwidth_*spl_->rmax();

		if ( nwk_->useSuperList() )
		{
			nwk_->dsuperList()	=  dsuperList_*spl_->rmax();
			nwk_->dsuperListP()	=  dsuperListP_*spl_->rmax();
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
		}
	}
}

void shearP_CD::perturbation(){
	cout<<"Je perturbe dans un cisaillement periodique"<<endl;
}
void shearP_CD::updateShear(){


	if(!shearRate_) return ;
	cout<<". Taux de deformation initial imposé : "<<shearRate_init<<endl;

	topXvalue_ = shearRate_init ;
	//On recalcule l epaisseur
	double h = ldof(1)->lowerBody()->y() - ldof(0)->lowerBody()->y();
	// On re applique le taux de cisaillement initial:
	topXvalue_ *= h ;


	if(symetrical_)
	{
		cout<<".Imposed symetrical shear rate: "<<topXvalue_<<endl;
		//topXvalue_ *= ldof(1)->mcy()-ldof(0)->mcy();
		topXvalue_ *= 0.5;

		ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
		ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY,-topXvalue_, 0.,  0.);
	}
	else
	{

		cout<<".Imposed non-symetrical shear rate: "<<topXvalue_<<endl;
		ldof(1)->affect(topXmode_, topYmode_, _VELOCITY, topXvalue_, topYvalue_, 0.);
		ldof(0)->affect(_VELOCITY, _VELOCITY, _VELOCITY, 0., 0., 0.);
	}
}
void shearP_CD::printMetrics()
{
	double xmin,xmax,ymin,ymax;
	ofstream metrics("Analyse/metrics.txt");

	xmin = spl()->leftBoundary();
	xmax = spl()->rightBoundary();
	ymin = ldof(0)->lowerBody()->y();
	ymax = ldof(1)->lowerBody()->y();
	metrics << "xmin "<<xmin<<endl;
	metrics << "xmax "<<xmax<<endl;
	metrics << "ymin "<<ymin<<endl;
	metrics << "ymax "<<ymax<<endl;
	metrics.close();
}


