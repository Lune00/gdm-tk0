#include "pgpg.hpp"

void pgpg::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
	is >> *Id1 >> *Id2 >> rank_
	  >> nx_ >> ny_
	  >> fn_ >> ft_ >> frot_;
	if ( rank_ > 1 ) 
		{ is >> fn2_ >> ft2_ ;}
}

void pgpg::write(ostream & os)
{
	os.precision(15);
	os << "pgpg " << i_->id() << ' ' << j_->id() 
	  << ' ' << rank_ << ' '
	  << nx_ << ' ' << ny_ << ' ' 
	<< fn_ << ' ' << ft_ << ' ' << frot_ ;
	
	if ( rank_ > 1 )
		os << ' ' << fn2_ << ' ' << ft2_ << endl << flush;
	else 
	    os << endl << flush;
}

bool pgpg::Activate()
{
	/*gdm::vertex c1, c2, d1,d2;
	
	c1.x()=i_->x();
	c1.y()=i_->y();
	c2.x()=j_->x();
	c2.y()=j_->y();

	vector<gdm::vertex> v1,v2;
	v1=i_->SectionSegment(c1,c2);
	v2=j_->SectionSegment(c1,c2);
	d1.x()=v1[0].x();
	d1.y()=v1[0].y();
	d2.x()=v2[0].x();
	d2.y()=v2[0].y();
	dist=sqrt((d1.x()-d2.x())*(d1.x()-d2.x())+(d1.y()-d2.y())*(d1.y()-d2.y()));
	if ((dist<this->dAct())||(rank_>0)) return  true;
	else return false;*/
	
	if (rank_>0) return true;
	else return false;
}

double pgpg::Dist()
{
//cerr << "rang = " << rank_ s< endl;
	if (rank_ == 0)
		return 1.;
	else
		return -1.e-12;
}


// adaptation de la procedure anapaire(...) de J.J. Moreau
void pgpg::Frame()
{
	unsigned int v = 0;
	double c,s;
	unsigned int iCriticalVertex = 0, jCriticalVertex = 0;
	unsigned int iCriticalVertex0, jCriticalVertex0;
	int sens = 0;
	int obj = 0;
	unsigned int iprevVertex, jprevVertex;
	unsigned int inextVertex, jnextVertex;
	double iScalarProduct = -1.0E+10, jScalarProduct = 1.0E+10, vScalarProduct;
	double sepx, sepy, norm;
	double sepx0, sepy0;
	double delx, dely;
	double over0;

	double centrx,centry;
	double tanx=0,tany=0;
	double dist1,dist2;
	double abs1=0,abs2=0;
	double base1=0,base2=0;
	
	double diametre=this->first()->Rout() > this->second()->Rout() ? this->second()->Rout() : this->first()->Rout();
	double epsilon=1e-2*diametre;
	
	unsigned int nVertexi = i_->Vertex().size();
	unsigned int nVertexj = j_->Vertex().size();

	gdm::vertex * iVertex[nVertexi];
	gdm::vertex * jVertex[nVertexj];
	gdm::vertex * iNorm[nVertexi];
	gdm::vertex * jNorm[nVertexj];

// Global vertex coordinates and face normals
// Polygon i
	c = cos(i_->rot());
	s = sin(i_->rot());
	for(v = 0; v < nVertexi; ++v)
	{
	// Vertexes
		iVertex[v]->x() = i_->x() + c * i_->Vertex(v).x() - s * i_->Vertex(v).y();
		iVertex[v]->y() = i_->y() + s * i_->Vertex(v).x() + c * i_->Vertex(v).y();

	// Normales
		iNorm[v]->x() = c * i_->Normal(v).x() - s * i_->Normal(v).y();
		iNorm[v]->y() = s * i_->Normal(v).x() + c * i_->Normal(v).y();
	}

// Polygon j
	c = cos(j_->rot());
	s = sin(j_->rot());	
	for(v = 0; v < nVertexj ; ++v)
	{
	// Vertexes
		jVertex[v]->x() = j_->x() + c * j_->Vertex(v).x() - s * j_->Vertex(v).y();
		jVertex[v]->y() = j_->y() + s * j_->Vertex(v).x() + c * j_->Vertex(v).y();

	// Normales
		jNorm[v]->x() = c * j_->Normal(v).x() - s * j_->Normal(v).y();
		jNorm[v]->y() = s * j_->Normal(v).x() + c * j_->Normal(v).y();
	}  

// First estimate of the separation vector (i.e. unit vector of the intercenter line)
	sepx = j_->x() - i_->x();
	sepy = j_->y() - i_->y();
	sepx *= (norm = 1.0 / sqrt(sepx * sepx + sepy * sepy));
	sepy *= norm;

// Research of the critical vertexes and of the corresponding overlap 
	for(v = 0 ; v < nVertexi ; ++v)
	{
	// Max for body i
		if((vScalarProduct = sepx * iVertex[v]->x() + sepy * iVertex[v]->y()) > iScalarProduct)	
		{
			iCriticalVertex = v ;
			iScalarProduct = vScalarProduct;
		}
	}
	for(v = 0 ; v < nVertexj ; ++v)
	{
	// Min for body j
		if((vScalarProduct = sepx * jVertex[v]->x() + sepy * jVertex[v]->y()) < jScalarProduct)
		{
			jCriticalVertex = v ;
			jScalarProduct = vScalarProduct ;
		}
	}
	over0 = iScalarProduct - jScalarProduct;
	
	//cout<<"over0:= "<<over0<<endl;
	if(over0 <= -epsilon)
	{
		rank_ = 0;
		overlap_=-1;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
		nx_ = -sepx;
		ny_ = -sepy;
		tx_ = -ny_;
		ty_ =  nx_;
		return;
	}

	iCriticalVertex0 = iCriticalVertex;
	jCriticalVertex0 = jCriticalVertex;
	sepx0 = sepx;
	sepy0 = sepy;

// vecteur joignant les sommets critiques, dit segment actif
	delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x();
	dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y();

	sens = ((sepx * dely - sepy * delx) > 0.0) ? -1 : 1;

	switch(sens)
	{
		case -1:
		while(sepx * dely - sepy * delx > 0.0)
		{
		// retenir la moins ample des deux rotations négatives possibles
			iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);	
			jprevVertex = (jCriticalVertex == 0) ? (nVertexj - 1) : (jCriticalVertex - 1);
			if(   sepx * iNorm[iprevVertex]->x() + sepy * iNorm[iprevVertex]->y() 
				> -sepx * jNorm[jprevVertex]->x() - sepy * jNorm[jprevVertex]->y() )
			{
				sepx = iNorm[iprevVertex]->x();
				sepy = iNorm[iprevVertex]->y();
				iCriticalVertex = iprevVertex;	// jCriticalVertex est inchange
				obj = 1;
			}
			else	// i.e. B actif
			{
				sepx = -jNorm[jprevVertex]->x();
				sepy = -jNorm[jprevVertex]->y();
				jCriticalVertex = jprevVertex;	// iCriticalVertex est inchangé
				obj = -1 ;		
			}
		// actualiser 
			delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x() ;
			dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y() ;
			over0 = sepx * delx + sepy * dely ;
		}
		break ;

		case 1:
		while(sepx * dely - sepy * delx < 0.0)
		{
		// retenir la plus petite des deux rotations positives possibles
			if(   sepx * iNorm[iCriticalVertex]->x() + sepy * iNorm[iCriticalVertex]->y() 
				> -sepx * jNorm[jCriticalVertex]->x() - sepy * jNorm[jCriticalVertex]->y() )
			{
				sepx = iNorm[iCriticalVertex]->x() ;
				sepy = iNorm[iCriticalVertex]->y() ;
				iCriticalVertex = iCriticalVertex == nVertexi - 1 ? 0 : iCriticalVertex + 1 ;  // jCriticalVertex est inchangé
				obj = 1 ;
			}
			else	// i.e. B actif
			{
				sepx = -jNorm[jCriticalVertex]->x() ;
				sepy = -jNorm[jCriticalVertex]->y() ;
				jCriticalVertex = jCriticalVertex == nVertexj - 1 ? 0 : jCriticalVertex + 1 ;	// iCriticalVertex est inchangé
				obj = -1 ;		
			}
		// actualiser le segment actif
			delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x() ;
			dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y() ;
			over0 = sepx * delx + sepy * dely ;
		}
		break ;
	}
//	cout<<"over0:= "<<over0<<endl;	
	if(over0 <= -epsilon)	
	{
		rank_ = 0;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
		nx_ = -sepx;
		ny_ = -sepy;
		tx_ = -ny_;
		ty_ =  nx_;
		overlap_ = -1;
		return;
	}

// l'overlap est calculé; identification des points de contact
	centrx = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
	centry = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
	x_ = y_ = x2_ = y2_ = 0.0;
	nx_ = -sepx;	// sep est une normale sortante de l'objet i_ ou entrante de l'objet j_
	ny_ = -sepy;
	tx_ = -ny_;
	ty_ =  nx_;

	rank_ = 0;
	//rankco_=0;
	//overlap_ = 1;

	switch(obj)
	{
		case -1:	// B est actif
			// point principal:
		x_ = iVertex[iCriticalVertex]->x() - 0.5 * over0 * sepx;
		y_ = iVertex[iCriticalVertex]->y() - 0.5 * over0 * sepy;
		overlap_ = over0;

	// voir si d'autres sommets de A sont en contact
		iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
		dist1 = (iVertex[iprevVertex]->x() - centrx) * sepx + (iVertex[iprevVertex]->y() - centry) * sepy;
	//	cout<<"dist11:= "<<dist1<<endl;
		if(dist1 > -epsilon)
		{
			x2_ = iVertex[iprevVertex]->x() - dist1 * sepx;
			y2_ = iVertex[iprevVertex]->y() - dist1 * sepy;
			rank_ = 1;
		//	rankco_=1;
		}
		inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);
		dist2 = (iVertex[inextVertex]->x() - centrx) * sepx + (iVertex[inextVertex]->y() - centry) * sepy;
	//	cout<<"dist12:= "<<dist2<<endl;
		if(dist2 > -epsilon && dist2 > dist1)
		{
			x2_ = iVertex[inextVertex]->x() - dist2 * sepx;
			y2_ = iVertex[inextVertex]->y() - dist2 * sepy;
			rank_ = 1;
		//	rankco_=1;
		}

		if(rank_ == 1)	// cas de 2 points de contact
		{
			
			tanx = -sepy ;
			tany = sepx ;
			abs1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
			abs2 = (x2_ - centrx) * tanx + (y2_ - centry) * tany;
			switch(sens)
			{
				case -1:
				jnextVertex = (jCriticalVertex == nVertexj - 1) ? 0 : (jCriticalVertex + 1);
				base1 = (jVertex[jnextVertex]->x() - centrx) * tanx + (jVertex[jnextVertex]->y() - centry) * tany;
				base2 = (jVertex[jCriticalVertex]->x() - centrx) * tanx + (jVertex[jCriticalVertex]->y() - centry) * tany;
				break ;
				case 1:
				base1 = (jVertex[jCriticalVertex]->x() - centrx) * tanx + (jVertex[jCriticalVertex]->y() - centry) * tany;
				jprevVertex = jCriticalVertex == 0 ? nVertexj - 1 : jCriticalVertex - 1 ;
				base2 = (jVertex[jprevVertex]->x() - centrx) * tanx + (jVertex[jprevVertex]->y() - centry) * tany;
				break ;
			}

		}        
		break;


		case 1:	// A est actif
	// point principal:
		x_ = jVertex[jCriticalVertex]->x() + 0.5 * over0 * sepx;
		y_ = jVertex[jCriticalVertex]->y() + 0.5 * over0 * sepy;
		overlap_ = over0;		
	// voir si d'autres sommets de B sont en contact
		jprevVertex = (jCriticalVertex == 0) ? (nVertexj - 1) : (jCriticalVertex - 1);
		dist1 = (jVertex[jprevVertex]->x() - centrx) * sepx + (jVertex[jprevVertex]->y() - centry) * sepy;
	//	cout<<"dist21:= "<<dist1<<endl;
		if(dist1 < epsilon)
		{
			x2_ = jVertex[jprevVertex]->x() - dist1 * sepx;
			y2_ = jVertex[jprevVertex]->y() - dist1 * sepy;
			rank_ = 1;
		//	rankco_=1;
		}
		jnextVertex = (jCriticalVertex == nVertexj - 1) ? 0 : (jCriticalVertex + 1);
		dist2 = (jVertex[jnextVertex]->x() - centrx) * sepx + (jVertex[jnextVertex]->y() - centry) * sepy;
	//	cout<<"dist22:= "<<dist2<<endl;
		
		if(dist2 < epsilon && dist2 < dist1)
		{
			x2_ = jVertex[jnextVertex]->x() - dist2 * sepx;
			y2_ = jVertex[jnextVertex]->y() - dist2 * sepy;
			rank_ = 1;
		//	rankco_=1;
		}

		if(rank_ == 1)	// cas de 2 points de contact
		{
			
			tanx = -sepy;
			tany = sepx;
			abs1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
			abs2 = (x2_ - centrx) * tanx + (y2_ - centry) * tany;
			switch(sens)
			{
				case -1:
				base1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
				inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);
				base2 = (iVertex[inextVertex]->x() - centrx) * tanx + (iVertex[inextVertex]->y() - centry) * tany;	
				break ;
				case 1:
				iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
				base1 = (iVertex[iprevVertex]->x() - centrx) * tanx + (iVertex[iprevVertex]->y() - centry) * tany;
				base2 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
				break ;
			}
		}	
		break ;
	}

	if(abs2 > base2)  abs2 = base2;
	if(abs2 < base1)  abs2 = base1;
	x2_ = centrx + abs2 * tanx;
	y2_ = centry + abs2 * tany;
	rank_ += 1;	
	//rankco_+=1;
	
	if (rank_==2) longeur_=sqrt(pow(x2_-x_,2)+pow(y2_-y_,2));
	else longeur_=0;
	
	if ( (rank_==2) && (longeur_<0.05*diametre) ) 
	{
		rank_ = 0;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
		nx_ = -sepx;
		ny_ = -sepy;
		tx_ = -ny_;
		ty_ =  nx_;
		overlap_ = -1;
	}
	
	/*if (rank_==2) 
	{
		rank_ = 1;
		x_ = x2_ = 0.5 * (x_+x2_);
		y_ = y2_ = 0.5 * (y_+y2_);
	}*/
	
	if(rank_==1) fn2_=ft2_=0;
	

} // END pgpg::Frame

void pgpg::Frame2()
{
	unsigned int v = 0;
	double c,s;
	unsigned int iCriticalVertex = 0, jCriticalVertex = 0;
	unsigned int iCriticalVertex0, jCriticalVertex0;
	int sens = 0;
	int obj = 0;
	unsigned int iprevVertex, jprevVertex;
	unsigned int inextVertex, jnextVertex;
	double iScalarProduct = -1.0E+10, jScalarProduct = 1.0E+10, vScalarProduct;
	double sepx, sepy, norm;
	double sepx0, sepy0;
	double delx, dely;
	double over0;

	double centrx,centry;
	double tanx=0,tany=0;
	double dist1,dist2;
	double abs1=0,abs2=0;
	double base1=0,base2=0;
	
	double diametre=this->first()->Rout() > this->second()->Rout() ? this->second()->Rout() : this->first()->Rout();
	double epsilon=1e-2*diametre;
	
	unsigned int nVertexi = i_->Vertex().size();
	unsigned int nVertexj = j_->Vertex().size();

	gdm::vertex * iVertex[nVertexi];
	gdm::vertex * jVertex[nVertexj];
	gdm::vertex * iNorm[nVertexi];
	gdm::vertex * jNorm[nVertexj];

// Global vertex coordinates and face normals
// Polygon i
	c = cos(i_->rot());
	s = sin(i_->rot());
	for(v = 0; v < nVertexi; ++v)
	{
	// Vertexes
		iVertex[v]->x() = i_->x() + c * i_->Vertex(v).x() - s * i_->Vertex(v).y();
		iVertex[v]->y() = i_->y() + s * i_->Vertex(v).x() + c * i_->Vertex(v).y();

	// Normales
		iNorm[v]->x() = c * i_->Normal(v).x() - s * i_->Normal(v).y();
		iNorm[v]->y() = s * i_->Normal(v).x() + c * i_->Normal(v).y();
	}

// Polygon j
	c = cos(j_->rot());
	s = sin(j_->rot());	
	for(v = 0; v < nVertexj ; ++v)
	{
	// Vertexes
		jVertex[v]->x() = j_->x() + c * j_->Vertex(v).x() - s * j_->Vertex(v).y();
		jVertex[v]->y() = j_->y() + s * j_->Vertex(v).x() + c * j_->Vertex(v).y();

	// Normales
		jNorm[v]->x() = c * j_->Normal(v).x() - s * j_->Normal(v).y();
		jNorm[v]->y() = s * j_->Normal(v).x() + c * j_->Normal(v).y();
	}  

// First estimate of the separation vector (i.e. unit vector of the intercenter line)
	sepx = j_->x() - i_->x();
	sepy = j_->y() - i_->y();
	sepx *= (norm = 1.0 / sqrt(sepx * sepx + sepy * sepy));
	sepy *= norm;

// Research of the critical vertexes and of the corresponding overlap 
	for(v = 0 ; v < nVertexi ; ++v)
	{
	// Max for body i
		if((vScalarProduct = sepx * iVertex[v]->x() + sepy * iVertex[v]->y()) > iScalarProduct)	
		{
			iCriticalVertex = v ;
			iScalarProduct = vScalarProduct;
		}
	}
	for(v = 0 ; v < nVertexj ; ++v)
	{
	// Min for body j
		if((vScalarProduct = sepx * jVertex[v]->x() + sepy * jVertex[v]->y()) < jScalarProduct)
		{
			jCriticalVertex = v ;
			jScalarProduct = vScalarProduct ;
		}
	}
	over0 = iScalarProduct - jScalarProduct;
	
	//cout<<"over0:= "<<over0<<endl;
	if(over0 <= -epsilon)
	{
		cout<<"	over:="<<over0<<"	rang:="<<rank_<<endl;
		//rank_ = 0;
		overlap_=-epsilon;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
		/*nx_ = -sepx;
		ny_ = -sepy;*/
		tx_ = -ny_;
		ty_ =  nx_;
		return;
	}

	iCriticalVertex0 = iCriticalVertex;
	jCriticalVertex0 = jCriticalVertex;
	sepx0 = sepx;
	sepy0 = sepy;

// vecteur joignant les sommets critiques, dit segment actif
	delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x();
	dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y();

	sens = ((sepx * dely - sepy * delx) > 0.0) ? -1 : 1;

	switch(sens)
	{
		case -1:
		while(sepx * dely - sepy * delx > 0.0)
		{
		// retenir la moins ample des deux rotations négatives possibles
			iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);	
			jprevVertex = (jCriticalVertex == 0) ? (nVertexj - 1) : (jCriticalVertex - 1);
			if(   sepx * iNorm[iprevVertex]->x() + sepy * iNorm[iprevVertex]->y() 
				> -sepx * jNorm[jprevVertex]->x() - sepy * jNorm[jprevVertex]->y() )
			{
				sepx = iNorm[iprevVertex]->x();
				sepy = iNorm[iprevVertex]->y();
				iCriticalVertex = iprevVertex;	// jCriticalVertex est inchange
				obj = 1;
			}
			else	// i.e. B actif
			{
				sepx = -jNorm[jprevVertex]->x();
				sepy = -jNorm[jprevVertex]->y();
				jCriticalVertex = jprevVertex;	// iCriticalVertex est inchangé
				obj = -1 ;		
			}
		// actualiser 
			delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x() ;
			dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y() ;
			over0 = sepx * delx + sepy * dely ;
		}
		break ;

		case 1:
		while(sepx * dely - sepy * delx < 0.0)
		{
		// retenir la plus petite des deux rotations positives possibles
			if(   sepx * iNorm[iCriticalVertex]->x() + sepy * iNorm[iCriticalVertex]->y() 
				> -sepx * jNorm[jCriticalVertex]->x() - sepy * jNorm[jCriticalVertex]->y() )
			{
				sepx = iNorm[iCriticalVertex]->x() ;
				sepy = iNorm[iCriticalVertex]->y() ;
				iCriticalVertex = iCriticalVertex == nVertexi - 1 ? 0 : iCriticalVertex + 1 ;  // jCriticalVertex est inchangé
				obj = 1 ;
			}
			else	// i.e. B actif
			{
				sepx = -jNorm[jCriticalVertex]->x() ;
				sepy = -jNorm[jCriticalVertex]->y() ;
				jCriticalVertex = jCriticalVertex == nVertexj - 1 ? 0 : jCriticalVertex + 1 ;	// iCriticalVertex est inchangé
				obj = -1 ;		
			}
		// actualiser le segment actif
			delx = iVertex[iCriticalVertex]->x() - jVertex[jCriticalVertex]->x() ;
			dely = iVertex[iCriticalVertex]->y() - jVertex[jCriticalVertex]->y() ;
			over0 = sepx * delx + sepy * dely ;
		}
		break ;
	}
//	cout<<"over0:= "<<over0<<endl;	
	if(over0 <= -epsilon)	
	{
		cout<<"	over:="<<over0<<"	rang:="<<rank_<<endl;
		//rank_ = 0;
		overlap_=-epsilon;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
		/*nx_ = -sepx;
		ny_ = -sepy;*/
		tx_ = -ny_;
		ty_ =  nx_;
		return;
	}

// l'overlap est calculé; identification des points de contact
	centrx = 0.5 * (iVertex[iCriticalVertex]->x() + jVertex[jCriticalVertex]->x());
	centry = 0.5 * (iVertex[iCriticalVertex]->y() + jVertex[jCriticalVertex]->y());
	x_ = y_ = x2_ = y2_ = 0.0;
	/*nx_ = -sepx;	// sep est une normale sortante de l'objet i_ ou entrante de l'objet j_
	ny_ = -sepy;*/
	tx_ = -ny_;
	ty_ =  nx_;

	//rank_ = 0;
	
	switch(obj)
	{
		case -1:	// B est actif
			// point principal:
		x_ = iVertex[iCriticalVertex]->x() - 0.5 * over0 * sepx;
		y_ = iVertex[iCriticalVertex]->y() - 0.5 * over0 * sepy;
		overlap_ = over0;

	// voir si d'autres sommets de A sont en contact
		iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
		dist1 = (iVertex[iprevVertex]->x() - centrx) * sepx + (iVertex[iprevVertex]->y() - centry) * sepy;
	//	cout<<"dist11:= "<<dist1<<endl;
		if(dist1 > -epsilon)
		{
			x2_ = iVertex[iprevVertex]->x() - dist1 * sepx;
			y2_ = iVertex[iprevVertex]->y() - dist1 * sepy;
		//	rank_ = 1;
		}
		inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);
		dist2 = (iVertex[inextVertex]->x() - centrx) * sepx + (iVertex[inextVertex]->y() - centry) * sepy;
	//	cout<<"dist12:= "<<dist2<<endl;
		if(dist2 > -epsilon && dist2 > dist1)
		{
			x2_ = iVertex[inextVertex]->x() - dist2 * sepx;
			y2_ = iVertex[inextVertex]->y() - dist2 * sepy;
		//	rank_ = 1;
		}
		
	//	if(rank_ == 1)	// cas de 2 points de contact
		if(rank_ == 2)	// cas de 2 points de contact
		{
			
			tanx = -sepy ;
			tany = sepx ;
			abs1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
			abs2 = (x2_ - centrx) * tanx + (y2_ - centry) * tany;
			switch(sens)
			{
				case -1:
				jnextVertex = (jCriticalVertex == nVertexj - 1) ? 0 : (jCriticalVertex + 1);
				base1 = (jVertex[jnextVertex]->x() - centrx) * tanx + (jVertex[jnextVertex]->y() - centry) * tany;
				base2 = (jVertex[jCriticalVertex]->x() - centrx) * tanx + (jVertex[jCriticalVertex]->y() - centry) * tany;
				break ;
				case 1:
				base1 = (jVertex[jCriticalVertex]->x() - centrx) * tanx + (jVertex[jCriticalVertex]->y() - centry) * tany;
				jprevVertex = jCriticalVertex == 0 ? nVertexj - 1 : jCriticalVertex - 1 ;
				base2 = (jVertex[jprevVertex]->x() - centrx) * tanx + (jVertex[jprevVertex]->y() - centry) * tany;
				break ;
			}

		}        
		break;


		case 1:	// A est actif
	// point principal:
		x_ = jVertex[jCriticalVertex]->x() + 0.5 * over0 * sepx;
		y_ = jVertex[jCriticalVertex]->y() + 0.5 * over0 * sepy;
		overlap_ = over0;		
	// voir si d'autres sommets de B sont en contact
		jprevVertex = (jCriticalVertex == 0) ? (nVertexj - 1) : (jCriticalVertex - 1);
		dist1 = (jVertex[jprevVertex]->x() - centrx) * sepx + (jVertex[jprevVertex]->y() - centry) * sepy;
	//	cout<<"dist21:= "<<dist1<<endl;
		if(dist1 < epsilon)
		{
			x2_ = jVertex[jprevVertex]->x() - dist1 * sepx;
			y2_ = jVertex[jprevVertex]->y() - dist1 * sepy;
		//	rank_ = 1;
		//	rankco_=1;
		}
		jnextVertex = (jCriticalVertex == nVertexj - 1) ? 0 : (jCriticalVertex + 1);
		dist2 = (jVertex[jnextVertex]->x() - centrx) * sepx + (jVertex[jnextVertex]->y() - centry) * sepy;
	//	cout<<"dist22:= "<<dist2<<endl;
		
		if(dist2 < epsilon && dist2 < dist1)
		{
			x2_ = jVertex[jnextVertex]->x() - dist2 * sepx;
			y2_ = jVertex[jnextVertex]->y() - dist2 * sepy;
		//	rank_ = 1;
		//	rankco_=1;
		}

		if(rank_ == 2)	// cas de 2 points de contact
		{
			
			tanx = -sepy;
			tany = sepx;
			abs1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
			abs2 = (x2_ - centrx) * tanx + (y2_ - centry) * tany;
			switch(sens)
			{
				case -1:
				base1 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
				inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);
				base2 = (iVertex[inextVertex]->x() - centrx) * tanx + (iVertex[inextVertex]->y() - centry) * tany;	
				break ;
				case 1:
				iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
				base1 = (iVertex[iprevVertex]->x() - centrx) * tanx + (iVertex[iprevVertex]->y() - centry) * tany;
				base2 = (iVertex[iCriticalVertex]->x() - centrx) * tanx + (iVertex[iCriticalVertex]->y() - centry) * tany;
				break ;
			}
		}	
		break ;
	}

	if(abs2 > base2)  abs2 = base2;
	if(abs2 < base1)  abs2 = base1;
	x2_ = centrx + abs2 * tanx;
	y2_ = centry + abs2 * tany;
	//rank_ += 1;	
		
	if (rank_==2) longeur_=sqrt(pow(x2_-x_,2)+pow(y2_-y_,2));
	else longeur_=0;
	
	/*if ( (rank_==2) && (longeur_<0.05*diametre) ) 
	{
		rank_ = 0;
		x_ = x2_ = 0.5 * (iVertex[iCriticalVertex].x() + jVertex[jCriticalVertex].x());
		y_ = y2_ = 0.5 * (iVertex[iCriticalVertex].y() + jVertex[jCriticalVertex].y());
		nx_ = -sepx;
		ny_ = -sepy;
		tx_ = -ny_;
		ty_ =  nx_;
		overlap_ = -1;
	}*/
	
	if(rank_==1) fn2_=ft2_=0;
	

} // END pgpg::Frame


void pgpg::Kin()
{  
	//original calcul ok
	/*double cix = x_ - i_->x();
	double ciy = y_ - i_->y(); 

	double cjx = x_ - j_->x();
	double cjy = y_ - j_->y();

	double vx = i_->vx() - j_->vx()
		-ciy * i_->vrot() + cjy * j_->vrot();
	double vy = i_->vy() - j_->vy()
		+cix * i_->vrot() - cjx * j_->vrot();

	vn_ = vx * nx_ + vy * ny_;
	vt_ = vx * tx_ + vy * ty_;

	if (rank_ == 2  )
	{
		cix = x2_ - i_->x();
		ciy = y2_ - i_->y();

		cjx = x2_ - j_->x();
		cjy = y2_ - j_->y();

		vx = i_->vx() - j_->vx()
			-ciy * i_->vrot() + cjy * j_->vrot();
		vy = i_->vy() - j_->vy()
			+cix * i_->vrot() - cjx * j_->vrot();

		vn2_ = vx * nx_ + vy * ny_;
		vt2_ = vx * tx_ + vy * ty_;    
		}*/

	if(i_->bodyDof()== NULL && j_->bodyDof()== NULL)
		{
			double cix = x_ - i_->x();
			double ciy = y_ - i_->y(); 

			double cjx = x_ - j_->x();
			double cjy = y_ - j_->y();

			double vx = i_->vx() - j_->vx()
				-ciy * i_->vrot() + cjy * j_->vrot();
			double vy = i_->vy() - j_->vy()
				+cix * i_->vrot() - cjx * j_->vrot();

			vn_ = vx * nx_ + vy * ny_;
			vt_ = vx * tx_ + vy * ty_;

			if (rank_ == 2  )
			{
				cix = x2_ - i_->x();
				ciy = y2_ - i_->y();

				cjx = x2_ - j_->x();
				cjy = y2_ - j_->y();

				vx = i_->vx() - j_->vx()
					-ciy * i_->vrot() + cjy * j_->vrot();
				vy = i_->vy() - j_->vy()
					+cix * i_->vrot() - cjx * j_->vrot();

				vn2_ = vx * nx_ + vy * ny_;
				vt2_ = vx * tx_ + vy * ty_;    
			}
			return;
		}
		else if ( i_->bodyDof() != NULL && j_->bodyDof()== NULL )
		{
			double cix = x_ - i_->bodyDof()->mcx();
			double ciy = y_ - i_->bodyDof()->mcy(); 

			double cjx = x_ - j_->x();
			double cjy = y_ - j_->y();

			double vx = i_->bodyDof()->mcvx() - j_->vx()
				-ciy * i_->bodyDof()->mcvrot() + cjy * j_->vrot();
			double vy = i_->bodyDof()->mcvy() - j_->vy()
				+cix * i_->bodyDof()->mcvrot() - cjx * j_->vrot();

			vn_ = vx * nx_ + vy * ny_;
			vt_ = vx * tx_ + vy * ty_;

			if (rank_ == 2  )
			{
				cix = x2_ - i_->bodyDof()->mcx();
				ciy = y2_ - i_->bodyDof()->mcy();

				cjx = x2_ - j_->x();
				cjy = y2_ - j_->y();

				vx = i_->bodyDof()->mcvx() - j_->vx()
					-ciy * i_->bodyDof()->mcvrot() + cjy * j_->vrot();
				vy = i_->bodyDof()->mcvy() - j_->vy()
					+cix * i_->bodyDof()->mcvrot() - cjx * j_->vrot();

				vn2_ = vx * nx_ + vy * ny_;
				vt2_ = vx * tx_ + vy * ty_;    
			}
			return;
		}
		else if( i_->bodyDof()== NULL && j_->bodyDof() != NULL ) 
		{
			double cix = x_ - i_->x();
			double ciy = y_ - i_->y(); 

			double cjx = x_ - j_->bodyDof()->mcx();
			double cjy = y_ - j_->bodyDof()->mcy();

			double vx = i_->vx() - j_->bodyDof()->mcvx()
				-ciy * i_->vrot() + cjy * j_->bodyDof()->mcvrot();
			double vy = i_->vy() - j_->bodyDof()->mcvy()
				+cix * i_->vrot() - cjx * j_->bodyDof()->mcvrot();

			vn_ = vx * nx_ + vy * ny_;
			vt_ = vx * tx_ + vy * ty_;

			if (rank_ == 2  )
			{
				cix = x2_ - i_->x();
				ciy = y2_ - i_->y();

				cjx = x2_ - j_->bodyDof()->mcx();
				cjy = y2_ - j_->bodyDof()->mcy();

				vx = i_->vx() - j_->bodyDof()->mcvx()
					-ciy * i_->vrot() + cjy * j_->bodyDof()->mcvrot();
				vy = i_->vy() - j_->bodyDof()->mcvy()
					+cix * i_->vrot() - cjx * j_->bodyDof()->mcvrot();

				vn2_ = vx * nx_ + vy * ny_;
				vt2_ = vx * tx_ + vy * ty_;    
			}
			return;
		}

		else if(i_->bodyDof() != NULL && j_->bodyDof() != NULL) 
		{
			double cix = x_ -i_->bodyDof()->mcx();
			double ciy = y_ -i_->bodyDof()->mcy();

			double cjx = x_ -j_->bodyDof()->mcx();
			double cjy = y_ -j_->bodyDof()->mcy();

			double vx = i_->bodyDof()->mcvx() - j_->bodyDof()->mcvx() - ciy * i_->bodyDof()->mcvrot() + cjy * j_->bodyDof()->mcvrot();
			double vy = i_->bodyDof()->mcvy() - j_->bodyDof()->mcvy() + cix * i_->bodyDof()->mcvrot() - cjx * j_->bodyDof()->mcvrot();

			vn_ = vx*nx_ + vy*ny_;
			vt_ = vx*tx_ + vy*ty_ ;
	//	vrot_ =  i_->bodyDof()->mcvrot() - j_->bodyDof()->mcvrot();
			if (rank_ == 2  )
			{
				cix = x2_ - i_->bodyDof()->mcx();
				ciy = y2_ - i_->bodyDof()->mcy();

				cjx = x2_ - j_->bodyDof()->mcx();
				cjy = y2_ - j_->bodyDof()->mcy();

				vx = i_->bodyDof()->mcvx() - j_->bodyDof()->mcvx() - ciy * i_->bodyDof()->mcvrot() + cjy * j_->bodyDof()->mcvrot();
				vy = i_->bodyDof()->mcvy() - j_->bodyDof()->mcvy() + cix * i_->bodyDof()->mcvrot() - cjx * j_->bodyDof()->mcvrot();

				vn2_ = vx * nx_ + vy * ny_;
				vt2_ = vx * tx_ + vy * ty_;    
			}
			return;
		}

	}

	void pgpg::Res()
	{
		double fx  = fn_ * nx_ + ft_ * tx_;
		double fy  = fn_ * ny_ + ft_ * ty_;
		double cix = x_ - i_->x();
		double ciy = y_ - i_->y();

		double cjx = x_ - j_->x();
		double cjy = y_ - j_->y();

		if (i_->bodyDof()==NULL)
		{
			i_->fx()   += fx;
			i_->fy()   += fy;
			i_->frot() += (cix * fy - ciy * fx); 
		}
		else
		{
			i_->bodyDof()->resDof(fx,fy,0,x_,y_);
		} 

		if (j_->bodyDof()==NULL)
		{
			j_->fx()   -= fx;
			j_->fy()   -= fy;
			j_->frot() += (-cjx * fy + cjy * fx);
		}
		else
		{
			j_->bodyDof()->resDof(-fx,-fy,0,x_,y_);
		}

		if (rank_ == 2 )
		{
			fx  = fn2_ * nx_ + ft2_ * tx_;
			fy  = fn2_ * ny_ + ft2_ * ty_;
			cix = x2_ - i_->x();
			ciy = y2_ - i_->y(); 
			cjx = x2_ - j_->x();
			cjy = y2_ - j_->y();

			if (i_->bodyDof()==NULL)
			{
				i_->fx()   += fx;
				i_->fy()   += fy;
				i_->frot() += (cix * fy - ciy * fx); 
			}
			else
			{
				i_->bodyDof()->resDof(fx,fy,0,x2_,y2_);
			} 

			if (j_->bodyDof()==NULL)
			{
				j_->fx()   -= fx;
				j_->fy()   -= fy;
				j_->frot() += (-cjx * fy + cjy * fx);
			}
			else
			{
				j_->bodyDof()->resDof(-fx,-fy,0,x2_,y2_);
			}
		}
	}

	void pgpg::Res(const double dfn, const double dft, const double dfs)
	{	
		double dfx = dfn * nx_ + dft * tx_;
		double dfy = dfn * ny_ + dft * ty_;
		double cix, ciy, cjx, cjy;

		if (rank_ == 2 && current_ == 1)
		{
			cix = x2_ - i_->x();
			ciy = y2_ - i_->y(); 

			cjx = x2_ - j_->x();
			cjy = y2_ - j_->y();

			if (i_->bodyDof()==NULL)
			{
				i_->fx()   += dfx;
				i_->fy()   += dfy;
				i_->frot() += (cix * dfy - ciy * dfx) + dfs;
			}
			else
			{
				i_->bodyDof()->resDof(dfx,dfy,0,x2_,y2_);
			}  

			if (j_->bodyDof()==NULL)
			{
				j_->fx()   -= dfx;
				j_->fy()   -= dfy;
				j_->frot() += (-cjx * dfy + cjy * dfx) - dfs;
			}
			else
			{
				j_->bodyDof()->resDof(-dfx,-dfy,0,x2_,y2_);	
			}

			return;
		}

		cix = x_ - i_->x();
		ciy = y_ - i_->y(); 

		cjx = x_ - j_->x();
		cjy = y_ - j_->y();

		if (i_->bodyDof()==NULL)
		{
			i_->fx()   += dfx;
			i_->fy()   += dfy;
			i_->frot() += (cix * dfy - ciy * dfx);  
		}
		else
		{
			i_->bodyDof()->resDof(dfx,dfy,0,x_,y_);
		}

		if (j_->bodyDof()==NULL)
		{
			j_->fx()   -= dfx;
			j_->fy()   -= dfy;
			j_->frot() += (-cjx * dfy + cjy * dfx);
		}
		else
		{
			j_->bodyDof()->resDof(-dfx,-dfy,0,x_,y_);

		}
	}

	void pgpg::CDcoeff()
	{
		double mi,mj,momi,momj;
		double cin,cit,cjn,cjt;

		double mn,mt;

// temporaire
		double en = 0.0;
		double et = 0.0;

		if(i_->bodyDof() != NULL) 
		{
			mi   = i_->bodyDof()->m();
			momi = i_->bodyDof()->mom();

			double cix = x_ - i_->bodyDof()->mcx();
			double ciy = y_ - i_->bodyDof()->mcy();
			cin = cix * nx_ + ciy * ny_;
			cit = cix * tx_ + ciy * ty_;
		}
		else 
		{
			mi   = i_->mass();
			momi = i_->mom();
			double cix = x_ - i_->x();
			double ciy = y_ - i_->y();
			cin = cix * nx_ + ciy * ny_;
			cit = cix * tx_ + ciy * ty_;
		}

		if(j_->bodyDof() != NULL) 
		{
			mj   = j_->bodyDof()->m();
			momj = j_->bodyDof()->mom();

			double cjx = x_ -j_->bodyDof()->mcx();
			double cjy = y_ -j_->bodyDof()->mcy();
			cjn = cjx * nx_ + cjy * ny_;
			cjt = cjx * tx_ + cjy * ty_;
		}
		else 
		{
			mj   = j_->mass();
			momj = j_->mom();
			double cjx = x_ - j_->x();
			double cjy = y_ - j_->y();
			cjn = cjx * nx_ + cjy * ny_;
			cjt = cjx * tx_ + cjy * ty_;
		}

		mn = 1.0/(1.0/mi+1.0/mj
			+ (cit*cit)/momi + (cjt*cjt)/momj);
		facn0_ = (1.0 + en) * mn;
		facn1_ = mn / mi;
		facn2_ = mn / mj;
		facn3_ = mn*cit/momi;
		facn4_ = mn*cjt/momj;

		mt     = 1.0 / (1.0 / mi + 1.0 / mj 
			+ (cin*cin)/momi + (cjn*cjn)/momj);
		fact0_ = (1.0 + et) * mt;
		fact1_ = mt / mi;
		fact2_ = mt / mj;
		fact3_ = mt * cin / momi;
		fact4_ = mt * cjn / momj;

/*	ms     = 1.0 / (1.0 / momi + 1.0 / momj);
		facs0_ = (1.0 + et) * ms;
		facs1_ = ms / momi;
		facs2_ = ms / momj;*/


/*
			double mn = 1.0/(1.0/i_->mass()+1.0/j_->mass()
			+ (cit*cit)/i_->mom() + (cjt*cjt)/j_->mom());
		facn0_ = (1.0 + en)*mn;
		facn1_ = mn/i_->mass();
		facn2_ = mn/j_->mass();
		facn3_ = mn*cit/i_->mom();
		facn4_ = mn*cjt/j_->mom();

		double mt = 1.0 / (1.0 / i_->mass() + 1.0 / j_->mass() 
			+ (cin*cin)/i_->mom() + (cjn*cjn)/j_->mom());
		fact0_ = (1.0 + et) * mt;
		fact1_ = mt / i_->mass();
		fact2_ = mt / j_->mass();
		fact3_ = mt * cin / i_->mom();
		fact4_ = mt * cjn / j_->mom();
	*/

		if (rank_ < 2) return;

		if(i_->bodyDof() != NULL) 
		{
			mi   = i_->bodyDof()->m();
			momi = i_->bodyDof()->mom();

			double cix = x2_ - i_->bodyDof()->mcx();
			double ciy = y2_ - i_->bodyDof()->mcy();
			cin = cix * nx_ + ciy * ny_;
			cit = cix * tx_ + ciy * ty_;
		}
		else 
		{
			mi   = i_->mass();
			momi = i_->mom();
			double cix = x2_ - i_->x();
			double ciy = y2_ - i_->y();
			cin = cix * nx_ + ciy * ny_;
			cit = cix * tx_ + ciy * ty_;
		}

		if(j_->bodyDof() != NULL) 
		{
			mj   = j_->bodyDof()->m();
			momj = j_->bodyDof()->mom();

			double cjx = x2_ -j_->bodyDof()->mcx();
			double cjy = y2_ -j_->bodyDof()->mcy();
			cjn = cjx * nx_ + cjy * ny_;
			cjt = cjx * tx_ + cjy * ty_;
		}
		else 
		{
			mj   = j_->mass();
			momj = j_->mom();
			double cjx = x2_ - j_->x();
			double cjy = y2_ - j_->y();
			cjn = cjx * nx_ + cjy * ny_;
			cjt = cjx * tx_ + cjy * ty_;
		}

		mn = 1.0/(1.0/mi+1.0/mj
			+ (cit*cit)/momi + (cjt*cjt)/momj);
		facn0_2_ = (1.0 + en) * mn;
		facn1_2_ = mn / mi;
		facn2_2_ = mn / mj;
		facn3_2_ = mn*cit/momi;
		facn4_2_ = mn*cjt/momj;

		mt     = 1.0 / (1.0 / mi + 1.0 / mj 
			+ (cin*cin)/momi + (cjn*cjn)/momj);
		fact0_2_ = (1.0 + et) * mt;
		fact1_2_ = mt / mi;
		fact2_2_ = mt / mj;
		fact3_2_ = mt * cin / momi;
		fact4_2_ = mt * cjn / momj;

// calculer les _facs

	}

	double pgpg::An(const double dt)
	{
		if( i_->bodyDof()==NULL && j_->bodyDof()==NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- facn0_2_ * vn2_/dt + fn2_
					- facn1_2_ * (i_->fx()*nx_ + i_->fy()*ny_)
					+ facn2_2_ * (j_->fx()*nx_ + j_->fy()*ny_)
					+ facn3_2_ * i_->frot()
					- facn4_2_ * j_->frot() 
					);     
			}

			return (
				- facn0_ * vn_/dt + fn_
				- facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
				+ facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
				+ facn3_ * i_->frot()
				- facn4_ * j_->frot() 
				);  
		}

		else if( i_->bodyDof() !=NULL && j_->bodyDof()==NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- facn0_2_ * vn2_/dt + fn2_
					- facn1_2_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
					+ facn2_2_ * (j_->fx()*nx_ + j_->fy()*ny_)
					+ facn3_2_ * i_->bodyDof()->mcfrot()
					- facn4_2_ * j_->frot() 
					);     
			}

			return (
				- facn0_ * vn_/dt + fn_
				- facn1_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
				+ facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
				+ facn3_ * i_->bodyDof()->mcfrot()
				- facn4_ * j_->frot() 
				);  
		}
		else if( i_->bodyDof()==NULL && j_->bodyDof()!=NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- facn0_2_ * vn2_/dt + fn2_
					- facn1_2_ * (i_->fx()*nx_ + i_->fy()*ny_)
					+ facn2_2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
					+ facn3_2_ * i_->frot()
					- facn4_2_ * j_->bodyDof()->mcfrot() 
					);     
			}

			return (
				- facn0_ * vn_/dt + fn_
				- facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
				+ facn2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
				+ facn3_ * i_->frot()
				- facn4_ * j_->bodyDof()->mcfrot() 
				);  
		}
		else if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- facn0_2_ * vn2_/dt + fn2_
					- facn1_2_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
					+ facn2_2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
					+ facn3_2_ * i_->bodyDof()->mcfrot()
					- facn4_2_ * j_->bodyDof()->mcfrot() 
					);     
			}

			return (
				- facn0_ * vn_/dt + fn_
				- facn1_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
				+ facn2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
				+ facn3_ * i_->bodyDof()->mcfrot()
				- facn4_ * j_->bodyDof()->mcfrot() 
				);  
		}
		else 
		{
			cout<<" probleme pgpg"<<endl;
			return 0;
		}

	}

	double pgpg::At(const double dt)
	{
		if( i_->bodyDof()==NULL && j_->bodyDof()==NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- fact0_2_ * vt2_/dt + ft2_
					- fact1_2_ * (i_->fx()*tx_ + i_->fy()*ty_)
					+ fact2_2_ * (j_->fx()*tx_ + j_->fy()*ty_)
					- fact3_2_ * i_->frot()
					+ fact4_2_ * j_->frot() 
					);     
			}

			return (
				- fact0_ * vt_/dt + ft_
				- fact1_ * (i_->fx()*tx_ + i_->fy()*ty_)
				+ fact2_ * (j_->fx()*tx_ + j_->fy()*ty_)
				- fact3_ * i_->frot()
				+ fact4_ * j_->frot() 
				);
		}
		else if( i_->bodyDof() !=NULL && j_->bodyDof()==NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- fact0_2_ * vt2_/dt + ft2_
					- fact1_2_ * (i_->bodyDof()->mcfx()*tx_ + i_->bodyDof()->mcfy()*ty_)
					+ fact2_2_ * (j_->fx()*tx_ + j_->fy()*ty_)
					- fact3_2_ * i_->bodyDof()->mcfrot()
					+ fact4_2_ * j_->frot() 
					);     
			}

			return (
				- fact0_ * vt_/dt + ft_
				- fact1_ * (i_->bodyDof()->mcfx()*tx_ + i_->bodyDof()->mcfy()*ty_)
				+ fact2_ * (j_->fx()*tx_ + j_->fy()*ty_)
				- fact3_ * i_->bodyDof()->mcfrot()
				+ fact4_ * j_->frot() 
				);
		}
		else if ( i_->bodyDof() ==NULL && j_->bodyDof() !=NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- fact0_2_ * vt2_/dt + ft2_
					- fact1_2_ * (i_->fx()*tx_ + i_->fy()*ty_)
					+ fact2_2_ * (j_->bodyDof()->mcfx()*tx_ + j_->bodyDof()->mcfy()*ty_)
					- fact3_2_ * i_->frot()
					+ fact4_2_ * j_->bodyDof()->mcfrot() 
					);     
			}

			return (
				- fact0_ * vt_/dt + ft_
				- fact1_ * (i_->fx()*tx_ + i_->fy()*ty_)
				+ fact2_ * (j_->bodyDof()->mcfx()*tx_ + j_->bodyDof()->mcfy()*ty_)
				- fact3_ * i_->frot()
				+ fact4_ * j_->bodyDof()->mcfrot() 
				);
		}
		else if ( i_->bodyDof() !=NULL && j_->bodyDof() !=NULL )
		{
			if (rank_ == 2 && current_ == 1)
			{
				return (
					- fact0_2_ * vt2_/dt + ft2_
					- fact1_2_ * (i_->bodyDof()->mcfx()*tx_ + i_->bodyDof()->mcfy()*ty_)
					+ fact2_2_ * (j_->bodyDof()->mcfx()*tx_ + j_->bodyDof()->mcfy()*ty_)
					- fact3_2_ * i_->bodyDof()->mcfrot()
					+ fact4_2_ * j_->bodyDof()->mcfrot() 
					);     
			}

			return (
				- fact0_ * vt_/dt + ft_
				- fact1_ * (i_->bodyDof()->mcfx()*tx_ + i_->bodyDof()->mcfy()*ty_)
				+ fact2_ * (j_->bodyDof()->mcfx()*tx_ + j_->bodyDof()->mcfy()*ty_)
				- fact3_ * i_->bodyDof()->mcfrot()
				+ fact4_ * j_->bodyDof()->mcfrot() 
				);
		}
		else 
		{
			cout<<" pb pgpg at"<<endl;
			return 0;
		}
	}

	double pgpg::As(const double dt)
	{
		return (
			- facs0_ * vrot_/dt + frot_
			- facs1_ * i_->frot()
			+ facs2_ * j_->frot()
			);
	}

	void pgpg::writeMGP(ostream & os)
	{
// NOT YET IMPLEMENTED IN MGPOST
	}
