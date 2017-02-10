#include "dkrl.hpp"

void dkrl::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
	is >> *Id1 >> *Id2
	  >> nx_ >> ny_
	  >> fn_ >> ft_ >> frot_;
}

void dkrl::write(ostream & os)
{
	os << "dkrl " << i_->id() << ' ' << j_->id() << ' '
	  << nx_ << ' ' << ny_ << ' ' 
	  << fn_ << ' ' << ft_ << ' ' << frot_
	  << endl << flush;
	
}

double dkrl::Dist()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double h  = 0.5 * j_->L();
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());

	double la = lx*ax + ly*ay;

	if(fabs(la)<=h)
		return fabs(lx*ay-ly*ax) - i_->R() - j_->R();

	if (la>h)
	{
		lx -= h*ax;
		ly -= h*ay;
		return sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
	}
	else
	{
		lx += h*ax;
		ly += h*ay;
		return sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
	}
}

double dkrl::Overlap()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double h  = 0.5 * j_->L();
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());

	double la = lx*ax + ly*ay;
	double d;

	if(fabs(la)<=h)
		d= fabs(lx*ay-ly*ax) - i_->R() - j_->R();

	else
	{if (la>h)
	{
		lx -= h*ax;
		ly -= h*ay;
		d= sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
	}
	else
	{
		lx += h*ax;
		ly += h*ay;
		d= sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
	}
	}
	if (d>0.) return -1.;
	else return -d;
}

bool dkrl::Activate()
{
	//if (rank_ > 0) return true;
	//return false;
	return activate();
}

void dkrl::Frame()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());
	double h  = 0.5 * j_->L();

	double la = lx*ax + ly*ay;

	if (fabs(la)<=h)
	{
		double prod = lx * ay - ly * ax;
		if(prod > 0.0) 
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
		
		if (fabs(prod) - this->dAct() < (i_->R() + j_->R())) 
		{
		this->activate()=true;
			//rankco_=1;
			//overlap_ = fabs(prod) - (i_->R() + j_->R());
		}
		else
		{
			this->activate()=false;
			//overlap_=1.;
		}
		
		if (fabs(prod) < (i_->R() + j_->R())) 
		{
			rank_ = 1;
			//rankco_=1;
			//overlap_ = fabs(prod) - (i_->R() + j_->R());
		}
		else
		{
			rank_=0;
			//overlap_=1.;
		}
	}
	else 
	{
		//cout<<"detection bout"<<endl;
		if (la>h)
		{
			la  = h;
			lx -= h*ax;
			ly -= h*ay;
			double l = sqrt(lx*lx+ly*ly);
			double invl = 1.0 / l;

			nx_ = lx*invl;
			ny_ = ly*invl;
			tx_ = -ny_;
			ty_ = nx_;

			if (l < (i_->R() + j_->R()))
			{
				rank_ = 1;
			//	rankco_=1;
				//overlap_ = l - (i_->R() + j_->R());
			}
			else
			{
				rank_=0;
				//overlap_=1.;
			}
		}
		else
		{
			la  = -h;
			lx += h*ax;
			ly += h*ay;
			double l = sqrt(lx*lx+ly*ly);
			double invl = 1.0 / l;

			nx_ = lx*invl;
			ny_ = ly*invl;
			tx_ = -ny_;
			ty_ = nx_;	

			if (l < (i_->R() + j_->R()))
			{
				rank_ = 1;
				//rankco_=1;
				//overlap_ = l - (i_->R() + j_->R());
			}
			else
			{
				rank_=0;
				//overlap_=1.;
			}
		}
	}

	x_ = j_->x() + la*ax + j_->R()*nx_;
	y_ = j_->y() + la*ay + j_->R()*ny_;	
}

void dkrl::Frame2()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double ax = cos(j_->rot());
	double ay = sin(j_->rot());
	double h  = 0.5 * j_->L();

	double la = lx*ax + ly*ay;

	if (fabs(la)<=h)
	{
		double prod = lx * ay - ly * ax;
		if(prod > 0.0) 
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
		
		if (fabs(prod) - this->dAct() < (i_->R() + j_->R())) 
		{
		this->activate()=true;
			//rankco_=1;
			//overlap_ = fabs(prod) - (i_->R() + j_->R());
		}
		else
		{
			this->activate()=false;
			//overlap_=1.;
		}
		
		if (fabs(prod) < (i_->R() + j_->R())) 
		{
			rank_ = 1;
			//rankco_=1;
			//overlap_ = fabs(prod) - (i_->R() + j_->R());
		}
		else
		{
			rank_=0;
			//overlap_=1.;
		}
	}
	else 
	{
		//cout<<"detection bout"<<endl;
		if (la>h)
		{
			la  = h;
			lx -= h*ax;
			ly -= h*ay;
			double l = sqrt(lx*lx+ly*ly);
			double invl = 1.0 / l;

			nx_ = lx*invl;
			ny_ = ly*invl;
			tx_ = -ny_;
			ty_ = nx_;

			if (l < (i_->R() + j_->R()))
			{
				rank_ = 1;
			//	rankco_=1;
				//overlap_ = l - (i_->R() + j_->R());
			}
			else
			{
				rank_=0;
				//overlap_=1.;
			}
		}
		else
		{
			la  = -h;
			lx += h*ax;
			ly += h*ay;
			double l = sqrt(lx*lx+ly*ly);
			double invl = 1.0 / l;

			nx_ = lx*invl;
			ny_ = ly*invl;
			tx_ = -ny_;
			ty_ = nx_;	

			if (l < (i_->R() + j_->R()))
			{
				rank_ = 1;
				//rankco_=1;
				//overlap_ = l - (i_->R() + j_->R());
			}
			else
			{
				rank_=0;
				//overlap_=1.;
			}
		}
	}

	x_ = j_->x() + la*ax + j_->R()*nx_;
	y_ = j_->y() + la*ay + j_->R()*ny_;	
}

void dkrl::Kin()
{
	//WARNING :  vrot n'est pas calculée----------------------
	if(i_->bodyDof()== NULL && j_->bodyDof()== NULL)
	{	
		double lx = i_->x() - j_->x();
		double ly = i_->y() - j_->y();
		double ax = cos(j_->rot());
		double ay = sin(j_->rot());

		double la = lx*ax + ly*ay;

		double Rvrot = i_->R() * i_->vrot() + j_->R() * j_->vrot();
		//On peut calculer les vitesse a partir de la position du point de contact
		double vx    = i_->vx() - j_->vx()
			+ Rvrot*ny_ + la*ay * j_->vrot();//---------faux cv 29/01
			//-i->vrot()*i->R()*nx_ - j_->R() * j_->vrot()*ny_
		double vy    = i_->vy() - j_->vy()
			- Rvrot*nx_ - la*ax * j_->vrot();//---------faux
			//-i->vrot()*i->R()*ny_ - j_->R() * j_->vrot()*nx_


		vn_ = vx*nx_ + vy*ny_;
		vt_ = vx*tx_ + vy*ty_;
		return;
	}
	else if ( i_->bodyDof() != NULL && j_->bodyDof()== NULL )
	{
		double cix = x_ - i_->bodyDof()->mcx();
		double ciy = y_ - i_->bodyDof()->mcy();

		double cjx = x_ - j_->x();
		double cjy = y_ - j_->y();

		double vx = i_->bodyDof()->mcvx() - j_->vx() - ciy * i_->bodyDof()->mcvrot() + cjy * j_->vrot();
		double vy = i_->bodyDof()->mcvy() - j_->vy() + cix * i_->bodyDof()->mcvrot() - cjx * j_->vrot();

		vn_ = vx*nx_ + vy*ny_;
		vt_ = vx*tx_ + vy*ty_ ;
		vrot_ =  i_->bodyDof()->mcvrot() - j_->vrot();
		return;
	}
	else if( i_->bodyDof()== NULL && j_->bodyDof() != NULL ) 
	{
		double cix = x_ -i_->x();
		double ciy = y_ -i_->y();

		double cjx = x_ -j_->bodyDof()->mcx();
		double cjy = y_ -j_->bodyDof()->mcy();

		double vx = i_->vx() - j_->bodyDof()->mcvx() - ciy* i_->vrot() + cjy * j_->bodyDof()->mcvrot();
		double vy = i_->vy() - j_->bodyDof()->mcvy() + cix* i_->vrot() - cjx * j_->bodyDof()->mcvrot();

		vn_ = vx*nx_ + vy*ny_;
		vt_ = vx*tx_ + vy*ty_ ;
		vrot_ =  i_->vrot() - j_->bodyDof()->mcvrot();
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
		vrot_ =  i_->bodyDof()->mcvrot() - j_->bodyDof()->mcvrot();
		return;
	}

}

void dkrl::Res()
{

	double fx  = fn_*nx_ + ft_*tx_;
	double fy  = fn_*ny_ + ft_*ty_;
	double cjx = x_ - j_->x();
	double cjy = y_ - j_->y();

	if (i_->bodyDof()==NULL)
	{
		i_->fx()   += fx;
		i_->fy()   += fy;
		i_->frot() -=  i_->R() * ft_;
	}
	else
	{
		i_->bodyDof()->resDof(fx,fy,-i_->R() * ft_ + frot_,x_,y_);
	}
	if (j_->bodyDof()==NULL)
	{
		j_->fx()   -= fx;
		j_->fy()   -= fy;
		j_->frot() += (-cjx*fy + cjy*fx);

	}
	else
	{
		j_->bodyDof()->resDof(-fx,-fy ,-j_->R() * ft_ - frot_,x_,y_);
	}

}

void dkrl::Res(const double dfn, const double dft, const double dfs)
{
	double dfx = dfn*nx_ + dft*tx_;
	double dfy = dfn*ny_ + dft*ty_;
	double cjx = x_ - j_->x();
	double cjy = y_ - j_->y();
	if (i_->bodyDof()==NULL)
	{
		i_->fx()   += dfx;
		i_->fy()   += dfy;
		i_->frot() -= i_->R() * dft - dfs;
	}
	else
	{
		i_->bodyDof()->resDof(dfx,dfy,-i_->R() * dft + dfs,x_,y_);
	}
	if (j_->bodyDof()==NULL)
	{
		j_->fx()   -= dfx;
		j_->fy()   -= dfy;
		j_->frot() += (-cjx*dfy + cjy*dfx) - dfs;
	}
	else
	{
		j_->bodyDof()->resDof(-dfx,-dfy,-j_->R() * dft - dfs,x_,y_);
	}
}

void dkrl::CDcoeff()
{
	double mn, mt, ms;
	double en = 0.0;
	double et = 0.0;
	double mi,mj,momi,momj;

	double cin,cit,cjn,cjt;

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
		cin  = - i_->R();
		cit=0;
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
		double cjx = x_ - j_->x();
		double cjy = y_ - j_->y();
		mj   = j_->mass();
		momj = j_->mom();
		cjn  = cjx*nx_ + cjy*ny_;
		cjt  = cjx*tx_ + cjy*ty_;

	}

	mn = 1.0/(1.0/mi+1.0/mj
		+ (cit*cit)/momi + (cjt*cjt)/momj);
	facn0_ = (1.0 + en) * mn;
	facn1_ = mn / mi;
	facn2_ = mn / mj;
	facn3_ = mn*cit/momi;
	facn4_ = mn*cjt/momj;

	//cout<<mn<<" "<<facn0_<<" "<<facn1_<<" "<<facn2_<<endl;
	mt     = 1.0 / (1.0 / mi + 1.0 / mj 
		+ (cin*cin)/momi + (cjn*cjn)/momj);
	fact0_ = (1.0 + et) * mt;
	fact1_ = mt / mi;
	fact2_ = mt / mj;
	fact3_ = mt * cin / momi;
	fact4_ = mt * cjn / momj;

	ms     = 1.0 / (1.0 / momi + 1.0 / momj);
	facs0_ = (1.0 + et) * ms;
	facs1_ = ms / momi;
	facs2_ = ms / momj;

}

double dkrl::An(const double dt) 
{
	if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
		return (
		- facn0_ * vn_/dt + fn_
		- facn1_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
		+ facn2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
		+ facn3_ * i_->bodyDof()->mcfrot()
		- facn4_ * j_->bodyDof()->mcfrot()
		);
	else if(i_->bodyDof()!=NULL && j_->bodyDof()==NULL)
		return (
		- facn0_ * vn_/dt + fn_ 
		- facn1_ * (i_->bodyDof()->mcfx()*nx_ + i_->bodyDof()->mcfy()*ny_)
		+ facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
		+ facn3_ * i_->bodyDof()->mcfrot()
		);
	else if(i_->bodyDof()==NULL && j_->bodyDof()!=NULL)
		return (
		- facn0_ * vn_/dt + fn_
		- facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
		+ facn2_ * (j_->bodyDof()->mcfx()*nx_ + j_->bodyDof()->mcfy()*ny_)
		- facn4_ * j_->bodyDof()->mcfrot()
		);

	else
		return (
		- facn0_ * vn_/dt + fn_
		- facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
		+ facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
		- facn4_ * j_->frot() 
		);

}

double dkrl::At(const double dt) 
{
	if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
		return (- fact0_ * vt_/dt + ft_
		- fact1_ * (i_->bodyDof()->mcfx()*tx_+i_->bodyDof()->mcfy()*ty_)
		+ fact2_ * (j_->bodyDof()->mcfx()*tx_+j_->bodyDof()->mcfy()*ty_)
		- fact3_ * i_->bodyDof()->mcfrot()
		+ fact4_ * j_->bodyDof()->mcfrot()//+ -> -
		);

	else if(i_->bodyDof()!=NULL && j_->bodyDof()==NULL )
		return (
		- fact0_ * vt_/dt + ft_
		- fact1_ * (i_->bodyDof()->mcfx()*tx_+i_->bodyDof()->mcfy()*ty_)
		+ fact2_ * (j_->fx()*tx_+ j_->fy()*ty_)
		- fact3_ * i_->bodyDof()->mcfrot()
		+ fact4_ * j_->frot()//+ -> -
		);

	else if(i_->bodyDof()==NULL && j_->bodyDof()!=NULL)
		return (
		- fact0_ * vt_/dt + ft_
		- fact1_ * (i_->fx()*tx_+ i_->fy()*ty_)
		+ fact2_ * (j_->bodyDof()->mcfx()*tx_+j_->bodyDof()->mcfy()*ty_)
		- fact3_ * i_->frot()
		+ fact4_ * j_->bodyDof()->mcfrot()//+ -> -
		);
	else
		return (
		- fact0_ * vt_/dt + ft_
		- fact1_ * (i_->fx()*tx_+i_->fy()*ty_)
		+ fact2_ * (j_->fx()*tx_+j_->fy()*ty_)
		- fact3_ * i_->frot()
		+ fact4_ * j_->frot() 
		);

}

double dkrl::As(const double dt)
{
	if( i_->bodyDof()==NULL && j_->bodyDof()==NULL  )
		return (
          - facs0_ * vrot_/dt + frot_
          - facs1_ * i_->frot()
          + facs2_ * j_->frot()
          );
	else if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
		return (
		- facs0_ * vrot_/dt + frot_
          - facs1_ * i_->bodyDof()->mcfrot()
          + facs2_ * j_->bodyDof()->mcfrot()
		);

	else if(i_->bodyDof()!=NULL && j_->bodyDof()==NULL )
		return (
		- facs0_ * vrot_/dt + frot_
          - facs1_ * i_->bodyDof()->mcfrot()
          + facs2_ * j_->frot()
		);

	else if(i_->bodyDof()==NULL && j_->bodyDof()!=NULL)
		return (
		- facs0_ * vrot_/dt + frot_
          - facs1_ * i_->frot()
          + facs2_ * j_->bodyDof()->mcfrot()
		);
	else
		{cout<<" pb dkrl"<<endl;return 0;}
}

void dkrl::writeMGP (ostream & os)
{
// Not yet used in mgpost
}
