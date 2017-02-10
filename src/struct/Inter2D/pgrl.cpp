#include "pgrl.hpp"

void pgrl::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
	is >> *Id1 >> *Id2 >> rank_
	  >> nx_ >> ny_
	  >> fn_ >> ft_ >> frot_;
	if ( rank_ > 1 ) 
		{ is >> fn2_ >> ft2_ ;}
}

void pgrl::write(ostream & os)
{
	os << "pgrl " << i_->id() << ' ' << j_->id() 
	  << ' ' << rank_ << ' '
	  << nx_ << ' ' << ny_ << ' ' 
	<< fn_ << ' ' << ft_ << ' ' << frot_ ;
	
	if ( rank_ > 1 )
		os << ' ' << fn2_ << ' ' << ft2_ << endl << flush;
	else 
	    os << endl << flush;
}

bool pgrl::Activate()
{
	if (rank_ > 0) return true;
	return false;
}

double pgrl::Dist()
{
// temporairement pour la dynamique des contacts
// ça devrait suffire ??

//cerr << "rang = " << rank_ << endl;
	if (rank_ == 0)
		return 1.;
	else
		return -1.e-12;
}


// adaptation de la procedure anapaire(...) de J.J. Moreau
void pgrl::Frame()
{
	
	unsigned int v = 0;
	double c,s;
	unsigned int iCriticalVertex = 0;
	unsigned int iprevVertex;
	unsigned int inextVertex;  
	
	double lx ;
	double ly ;
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());
	
	
	double h  = 0.5 * j_->L();
	double lt , la ;
	double prod,over=1.0E+6;

	unsigned int nVertexi = i_->Vertex().size();
//	cout<<" nver "<<nVertexi<<endl;

	gdm::vertex * iVertex[nVertexi];
	gdm::vertex * iNorm[nVertexi];

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


// Research of the critical vertexes and of the corresponding overlap 
	for(v = 0 ; v < nVertexi ; ++v)
	{
		lx= iVertex[v]->x() - j_->x();
		ly= iVertex[v]->y() - j_->y();
		lt= lx * ay - ly * ax;

		if( fabs( lt )<= over)//Contact le plus proche de l'axe du rline
		{
			over = fabs(lt);
			iCriticalVertex = v;
		}

	}
	//cout<<" icrit "<<iCriticalVertex<<endl;
	
	iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
	inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);

//Contact selon la longueur? double contact?
	lx= iVertex[iCriticalVertex]->x() - j_->x();
	ly= iVertex[iCriticalVertex]->y() - j_->y();
	la= lx * ax + ly * ay;

	if (fabs(la)<= h )//Contact sur la longueur du rline
	{
		lt = lx * ay - ly * ax;
		if(lt > 0.0) 
		{
			nx_ = ay;
			ny_ = -ax;
			tx_ = -ny_;
			ty_ = nx_;
		}
		else
		{
			nx_ = -ay;
			ny_ = ax;
			tx_ = -ny_;
			ty_ = nx_;
		}
		
		x_ = j_->x() + la*ax + j_->R()*nx_;
		y_ = j_->y() + la*ay + j_->R()*ny_;
		
		if (fabs(lt) <  j_->R() ) rank_ = 1;
		else 
		{
			rank_ = 0;
			overlap_=-1.;
			return;
		}
		// deuxieme point de contact avec iprev ou inext

		lx= iVertex[iprevVertex]->x() - x_;
		ly= iVertex[iprevVertex]->y() - y_;
		prod= lx*nx_ + ly*ny_;
		
		overlap_= j_->R()- fabs(lt);
		
		if ( prod < 0 )//Deuxieme contact
		{
			x2_ = x_ + lx;
			y2_ = y_ + ly;
			rank_+=1;
			//return;
		}

		lx= iVertex[inextVertex]->x() - x_;
		ly= iVertex[inextVertex]->y() - y_;
		prod= lx*nx_ + ly*ny_;
	
		if ( prod < 0 )//Deuxieme contact
		{
			x2_ = x_ + lx;
			y2_ = y_ + ly;
			rank_+=1;
			//return;
		}
	}
	else //Contact avec le bout
	{

	}
	
	if(rank_==1) fn2_=ft2_=0;
	
} // END pgrl::Frame

void pgrl::Frame2()
{
	
	unsigned int v = 0;
	double c,s;
	unsigned int iCriticalVertex = 0;
	unsigned int iprevVertex;
	unsigned int inextVertex;  
	
	double lx ;
	double ly ;
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());
	
	
	double h  = 0.5 * j_->L();
	double lt , la ;
	double prod,over=1.0E+6;

	unsigned int nVertexi = i_->Vertex().size();
//	cout<<" nver "<<nVertexi<<endl;

	gdm::vertex * iVertex[nVertexi];
	gdm::vertex * iNorm[nVertexi];

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


// Research of the critical vertexes and of the corresponding overlap 
	for(v = 0 ; v < nVertexi ; ++v)
	{
		lx= iVertex[v]->x() - j_->x();
		ly= iVertex[v]->y() - j_->y();
		lt= lx * ay - ly * ax;

		if( fabs( lt )<= over)//Contact le plus proche de l'axe du rline
		{
			over = fabs(lt);
			iCriticalVertex = v;
		}

	}
	//cout<<" icrit "<<iCriticalVertex<<endl;
	
	iprevVertex = (iCriticalVertex == 0) ? (nVertexi - 1) : (iCriticalVertex - 1);
	inextVertex = (iCriticalVertex == nVertexi - 1) ? 0 : (iCriticalVertex + 1);

//Contact selon la longueur? double contact?
	lx= iVertex[iCriticalVertex]->x() - j_->x();
	ly= iVertex[iCriticalVertex]->y() - j_->y();
	la= lx * ax + ly * ay;

	if (fabs(la)<= h )//Contact sur la longueur du rline
	{
		lt = lx * ay - ly * ax;
		if(lt > 0.0) 
		{
			nx_ = ay;
			ny_ = -ax;
			tx_ = -ny_;
			ty_ = nx_;
		}
		else
		{
			nx_ = -ay;
			ny_ = ax;
			tx_ = -ny_;
			ty_ = nx_;
		}
		
		x_ = j_->x() + la*ax + j_->R()*nx_;
		y_ = j_->y() + la*ay + j_->R()*ny_;
		
	/*	if (fabs(lt) <  j_->R() ) rank_ = 1;
		else 
		{
			rank_ = 0;
			overlap_=-1.;
			return;
		}*/
		// deuxieme point de contact avec iprev ou inext

		lx= iVertex[iprevVertex]->x() - x_;
		ly= iVertex[iprevVertex]->y() - y_;
		prod= lx*nx_ + ly*ny_;
		
		overlap_= j_->R()- fabs(lt);
		
		if ( prod < 0 )//Deuxieme contact
		{
			x2_ = x_ + lx;
			y2_ = y_ + ly;
		//	rank_+=1;
			//return;
		}

		lx= iVertex[inextVertex]->x() - x_;
		ly= iVertex[inextVertex]->y() - y_;
		prod= lx*nx_ + ly*ny_;
	
		if ( prod < 0 )//Deuxieme contact
		{
			x2_ = x_ + lx;
			y2_ = y_ + ly;
		//	rank_+=1;
			//return;
		}
	}
	else //Contact avec le bout
	{

	}
	
	if(rank_==1) fn2_=ft2_=0;
	
} // END pgrl::Frame

void pgrl::Kin()
{  
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

	if (rank_ == 2)
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

void pgrl::Res()
{
	/*double fx  = fn_ * nx_ + ft_ * tx_;
	double fy  = fn_ * ny_ + ft_ * ty_;
	double cix = x_ - i_->x();
	double ciy = y_ - i_->y();
	double cjx = x_ - j_->x();
	double cjy = y_ - j_->y();

	i_->fx()   += fx;
	i_->fy()   += fy;
	i_->frot() += (cix * fy - ciy * fx);  

	j_->fx()   -= fx;
	j_->fy()   -= fy;
	j_->frot() += (-cjx * fy + cjy * fx);

	if (rank_ == 2)
	{
		fx  = fn2_ * nx_ + ft2_ * nx_;
		fy  = fn2_ * ny_ + ft2_ * ny_;
		cix = x2_ - i_->x();
		ciy = y2_ - i_->y(); 
		cjx = x2_ - j_->x();
		cjy = y2_ - j_->y();

		i_->fx()   += fx;
		i_->fy()   += fy;
		i_->frot() += (cix * fy - ciy * fx);  

		j_->fx()   -= fx;
		j_->fy()   -= fy;
		j_->frot() += (-cjx * fy + cjy * fx); 
	}*/
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

void pgrl::Res(const double dfn, const double dft, const double dfs)
{	
	/*double dfx = dfn * nx_ + dft * tx_;
	double dfy = dfn * ny_ + dft * ty_;
	double cix, ciy, cjx, cjy;

	if (rank_ == 2 && current_ == 1)
	{
		cix = x2_ - i_->x();
		ciy = y2_ - i_->y(); 
		cjx = x2_ - j_->x();
		cjy = y2_ - j_->y();

		i_->fx()   += dfx;
		i_->fy()   += dfy;
		i_->frot() += (cix * dfy - ciy * dfx) + dfs;  

		j_->fx()   -= dfx;
		j_->fy()   -= dfy;
		j_->frot() += (-cjx * dfy + cjy * dfx) - dfs;

		return;
	}

	cix = x_ - i_->x();
	ciy = y_ - i_->y(); 
	cjx = x_ - j_->x();
	cjy = y_ - j_->y();

	i_->fx()   += dfx;
	i_->fy()   += dfy;
	i_->frot() += (cix * dfy - ciy * dfx);  

	j_->fx()   -= dfx;
	j_->fy()   -= dfy;
	j_->frot() += (-cjx * dfy + cjy * dfx);*/
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

void pgrl::CDcoeff()
{
	/*double cix = x_ - i_->x();
	double ciy = y_ - i_->y();
	double cit = cix * tx_ + ciy * ty_;
	double cin = cix * nx_ + ciy * ny_;

	double cjx = x_ - j_->x();
	double cjy = y_ - j_->y();
	double cjt = cjx * tx_ + cjy * ty_;
	double cjn = cjx * nx_ + cjy * ny_; 

// temporaire
	double en = 0.0;
    double et = 0.0;

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
	fact4_ = mt * cjn / i_->mom();

	if (rank_ < 2) return;

	cix = x2_ - i_->x();
	ciy = y2_ - i_->y();
	cit = cix*tx_ + ciy*ty_;
	cin = cix * nx_ + ciy * ny_;

	cjx = x2_ - j_->x();
	cjy = y2_ - j_->y();
	cjt = cjx*tx_ + cjy*ty_; 
	cjn = cjx * nx_ + cjy * ny_; 
	
	mn = 1.0/(1.0/i_->mass()+1.0/j_->mass()
		+ (cit*cit)/i_->mom() + (cjt*cjt)/j_->mom());
	facn0_2_ = (1.0 + en)*mn;
	facn1_2_ = mn/i_->mass();
	facn2_2_ = mn/j_->mass();
	facn3_2_ = mn*cit/i_->mom();
	facn4_2_ = mn*cjt/j_->mom();
	
	mt = 1.0 / (1.0 / i_->mass() + 1.0 / j_->mass() 
		+ (cin*cin)/i_->mom() + (cjn*cjn)/j_->mom());
	fact0_2_ = (1.0 + et) * mt;
	fact1_2_ = mt / i_->mass();
	fact2_2_ = mt / j_->mass();
	fact3_2_ = mt * cin / i_->mom();
	fact4_2_ = mt * cjn / i_->mom();
	*/
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

double pgrl::An(const double dt)
{

	/*if (rank_ == 2 && current_ == 1)
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
		*/
		
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
		cout<<" probleme pgrl an"<<endl;
		return 0;
	}
	
}

double pgrl::At(const double dt)
{
	/*if (rank_ == 2 && current_ == 1)
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
	*/
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
			cout<<" pb pgrl at"<<endl;
			return 0;
		}
}

double pgrl::As(const double dt)
{
	return (
		- facs0_ * vrot_/dt + frot_
		- facs1_ * i_->frot()
		+ facs2_ * j_->frot()
		);
}

void pgrl::writeMGP(ostream & os)
{
// NOT YET IMPLEMENTED IN MGPOST
}
