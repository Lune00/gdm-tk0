#include "dkdk.hpp"

void dkdk::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
	is >> *Id1 >> *Id2
		>> nx_ >> ny_
		>> fn_ >> ft_ >> frot_;
}

void dkdk::write(ostream & os)
{
	os << "dkdk " << i_->id() << ' ' << j_->id() << ' '
		<< nx_ << ' ' << ny_ << ' ' 
		<< fn_ << ' ' << ft_ << ' ' << frot_
		<< endl << flush;
}

//! Compute the distance
double dkdk::Dist()
{
	double lx = j_->x() - i_->x();
	double ly = j_->y() - i_->y();

	return sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
}

double dkdk::Overlap()
{
	double lx = j_->x() - i_->x();
	double ly = j_->y() - i_->y();
	double d= sqrt(lx*lx+ly*ly) - i_->R() - j_->R();
	if( d > 0 ) return -1.;
	else return -d;
}

bool dkdk::Activate()
{
	//if (rank_ > 0) return true;
	//return false;
	return this->activate();
}

//! Build the local frame
void dkdk::Frame()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double l  = sqrt(lx * lx + ly * ly);
	double invl = 1.0 / l;  

	nx_ = lx * invl;
	ny_ = ly * invl;
	tx_ = -ny_;
	ty_ =  nx_;

	x_ = i_->x() - i_->R() * nx_;
	y_ = i_->y() - i_->R() * ny_;

	if( l - this->dAct() < (i_->R() + j_->R()))
	{
		this->activate()=true;
	}
	else
	{
		this->activate()=false;
	}

	if (l < (i_->R() + j_->R())) //ajouter dAct() le 2 Janvier 2012
	{
		rank_ = 1; // pour les actions a distance ????
	//rankco_=1;
	}
	else rank_ = 0;

}

void dkdk::Frame2()
{
	double lx = i_->x() - j_->x();
	double ly = i_->y() - j_->y();
	double l  = sqrt(lx * lx + ly * ly);
	double invl = 1.0 / l;  

	nx_ = lx * invl;
	ny_ = ly * invl;
	tx_ = -ny_;
	ty_ =  nx_;

	x_ = i_->x() - i_->R() * nx_;
	y_ = i_->y() - i_->R() * ny_;

	if( l - this->dAct() < (i_->R() + j_->R()))
	{
		this->activate()=true;
	}
	else
	{
		this->activate()=false;
	}

	if (l < (i_->R() + j_->R())) //ajouter dAct() le 2 Janvier 2012
	{
		rank_ = 1; // pour les actions a distance ????
	//rankco_=1;
	}
	else rank_ = 0;

}

//! Compute the relative velocity
void dkdk::Kin()
{
	if(i_->bodyDof()== NULL && j_->bodyDof()== NULL)
	{
		double vx = i_->vx() - j_->vx();
		double vy = i_->vy() - j_->vy();

		vn_ = vx*nx_ + vy*ny_;
		vt_ = vx*tx_ + vy*ty_ 
			- i_->R() * i_->vrot() - j_->R() * j_->vrot();
		vrot_ = i_->vrot() - j_->vrot();
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

void dkdk::Res()
{		
	double fx = fn_*nx_ + ft_*tx_;
	double fy = fn_*ny_ + ft_*ty_;
	if (i_->bodyDof()==NULL)
	{
		i_->fx()   += fx;
		i_->fy()   += fy;
		i_->frot() -=  i_->R() * ft_ - frot_;
	}
	else
		i_->bodyDof()->resDof(fx,fy,-i_->R() * ft_ + frot_,x_,y_);//On peut supprimer le passage du moment, utile seulement pour les monocorps


	if (j_->bodyDof()==NULL)
	{
		j_->fx()   -= fx;
		j_->fy()   -= fy;
		j_->frot() -= j_->R() * ft_ + frot_;
	}
	else
		j_->bodyDof()->resDof(-fx,-fy ,-j_->R() * ft_ - frot_,x_,y_);

}

//! Increment the force resultants in the mass center of each body
void dkdk::Res(const double dfn, const double dft, const double dfs)
{		
	double dfx = dfn*nx_ + dft*tx_;
	double dfy = dfn*ny_ + dft*ty_;

	if (i_->bodyDof()==NULL)
	{
		i_->fx()   += dfx;
		i_->fy()   += dfy;
		i_->frot() -=  i_->R() * dft - dfs;
	}
	else
	{
		i_->bodyDof()->resDof(dfx,dfy,-i_->R() * dft + dfs,x_,y_);
	//cout<<"ifx = "<<i_->bodyDof()->fx()<<" ify = "<<i_->bodyDof()->fy()<<endl;
	}

	if (j_->bodyDof()==NULL)
	{
		j_->fx()   -= dfx;
		j_->fy()   -= dfy;
		j_->frot() -= j_->R() * dft + dfs;
	}
	else
	{
		j_->bodyDof()->resDof(-dfx,-dfy,-j_->R() * dft - dfs,x_,y_);
//cout<<"jfx = "<<j_->bodyDof()->fx()<<" jfy = "<<j_->bodyDof()->fy()<<endl;
	}
}

//WIP

void dkdk::CDcoeff(GroupRelationData * gR) // CDcoeff(double en,double et) ou CDcoeff(grpRel&)
{
	double mn, mt, ms;

//temporaire
	//double en = 1.0;
	//double et = 0.0;
	int g1 = i_->grp();
	int g2 = j_->grp();
	//Recuperer les parametres
	double mi,mj,momi,momj;

	double en = gR->getParameter("en",g1,g2);
	double et = gR->getParameter("et",g1,g2);
	//double mu = gR->getParameter("mu",g1,g2);

//	cout<<"en = "<<en<<" ----- et = "<<et<<endl;
//	cout<<"mu = "<<mu<<endl;
/*	if(g1 == 2 || g2 ==2){
		cout<<"Particule temoin"<<endl;
		cout<<"mu "<<g1<<"/"<<g2<<" = "<<gR->getParameter("mu",g1,g2)<<endl;
	}
*/

	double cin,cit,cjn,cjt;

	//cin=cit=cjn=cjt=0;

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
		mj   = j_->mass();
		momj = j_->mom();
		cjn  = j_->R();
		cjt=0;

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
//Original
void dkdk::CDcoeff() // CDcoeff(double en,double et) ou CDcoeff(grpRel&)
{
	double mn, mt, ms;

//temporaire
	double en = 0.0;
	double et = 0.0;
	//unsigned int g1 = i_->grp();
	//unsigned int g2 = j_->grp();
	//Recuperer les parametres
	double mi,mj,momi,momj;

	//double en = getParameter("en",g1,g2);
	//double et = getParameter("et",g1,g2);
	//cout<<"en = "<<en<<" ----- et = "<<et<<endl;

	double cin,cit,cjn,cjt;

	//cin=cit=cjn=cjt=0;

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
		mj   = j_->mass();
		momj = j_->mom();
		cjn  = j_->R();
		cjt=0;

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
// Force normal
double dkdk::An(const double dt) //équation 4.18 page 88 Charles
{
	if( i_->bodyDof()==NULL && j_->bodyDof()==NULL  )
		return (
		- facn0_ * vn_/dt + fn_
		- facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
		+ facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
		);
	else if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
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
		{cout<<" pb dkdk"<<endl;return 0;}

}
//Force friction
double dkdk::At(const double dt)
{
	if( i_->bodyDof()==NULL && j_->bodyDof()==NULL  )
		return (
		- fact0_ * vt_/dt + ft_
		- fact1_ * (i_->fx()*tx_+i_->fy()*ty_)
		+ fact2_ * (j_->fx()*tx_+j_->fy()*ty_)
		- fact3_ * i_->frot()
		+ fact4_ * j_->frot() 
		);
	else if( i_->bodyDof()!=NULL && j_->bodyDof()!=NULL )
		return (
		- fact0_ * vt_/dt + ft_
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
		{cout<<" pb dkdk"<<endl;return 0;}
}

double dkdk::As(const double dt)
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
		{cout<<" pb dkdk"<<endl;return 0;}

}

void dkdk::writeMGP(ostream & os)
{
	os << "    <DKDKx antac=\"" << j_->id()+1
		<< "\" rn=\"" << fn_
		<< "\" rt=\"" << ft_ 
		<< "\" M=\""  << frot_ << "\"/>" << endl;
}
