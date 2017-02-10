#include "oedometric.hpp"

void Oedometric::read_parameters(istream & is)
{
	string token;
	string mode;

	is >> token;	
	while(is)
	{	
		if (token == "yval")    { is >> yval_; }
		else if (token == "ymod")    { is >> mode; ymod_ = loadingMode(mode); }
		else if (token == "adjust")  adjust_ = true;
		else if (token == "gravity") gravity_= true;
		else if (token == "increment") {is >> incr_;}
		
		else if (token == "}")       break;
		else cerr << "@Biaxial::read_parameters, Unknown parameter: " << token << endl;

		is >> token;
	}
}

void Oedometric::write_parameters(ostream & os)
{
	//os << "biaxial" << endl;
	//os << "xval   " << xval_ << endl;
	//os << "yval   " << yval_ << endl;
	// TODO: use gdm::modeToString(...) instead
	//os << "xmod   " << loadingMode(xmod_) << endl;
	//os << "ymod   " << loadingMode(ymod_) << endl;
}

// TODO: init est appelé 2 fois !!!!! a voir...
void Oedometric::init()
	{  
	k_=0;

// Check the type of first four bodies
	for (unsigned int i=0;i<2;++i) 
		if (typeid(*(spl_->body(i))) != typeid(rline))
		gdm::fatal("@Biaxial::init, the first two bodies must be rlines!");
	spl_->radiusExtrema(2);
	
	cout<<"appel init"<<endl;
	
	
	double x,xmin,xmax;
	double y,ymin,ymax;
	xmin = xmax = spl_->body(2)->x();
	ymin = ymax = spl_->body(2)->y();
	for (unsigned int i=3;i<spl_->lbody().size();++i)
	{
		x = spl_->body(i)->xmin();
		xmin = (xmin > x) ? x : xmin;
		
		y = spl_->body(i)->ymin();
		ymin = (ymin > y) ? y : ymin;
			
		x = spl_->body(i)->xmax();
		xmax = (xmax < x) ? x : xmax;
		
		y = spl_->body(i)->ymax();
		ymax = (ymax < y) ? y : ymax;

	} 
	
	xmax_ = xmax;
	xmin_ = xmin;
	
	ymax_ = ymax;
	ymin_ = ymin;
	
	lx_ = xmax - xmin;
	ly_ = ymax - ymin;
	
	
	//ofstream stress_strain_out("stress_strain.dat",ios::out);
	//stress_strain_out<<"# sigma_yy  epsilon_yy"<<endl;
	//stress_strain_out.close();
	
	ofstream stress_strain_out("stress_strain.dat",ios::out);
	stress_strain_out<<"# p   F   epsilon_yy"<<endl;
	stress_strain_out.close();

// Defining the position of the boundary rlines
	rline * rl[2];
	rl[0] = dynamic_cast<rline*>(spl_->body(0));
	rl[1] = dynamic_cast<rline*>(spl_->body(1));
	
	for(unsigned int i=0;i<spl_->lbody().size();++i)
	spl_->body(i)->id()=i;
				
	bottom_ = rl[0];
	for (unsigned int i=1;i<2;++i) 	
	{	if( rl[i]->rot() == 0)
		bottom_ = (rl[i]->y() < bottom_->y()) ? rl[i] : bottom_;
	}	
	top_ = rl[0];
	for (unsigned int i=1;i<2;++i) 
	{	if( rl[i]->rot() == 0)
		top_    = (rl[i]->y() > top_->y())    ? rl[i] : top_;
	}
	
	idBottom_ = bottom_->id();
	idTop_    = top_->id();
	
	y0_ = top_->y();
	defy_ = 0.;

	
	
// We create 2 controls
	for (unsigned int i=0;i<2;++i) lctrl_.push_back(control());

// bottom is fixed
	lctrl_[idBottom_].x()      = _VELOCITY;
	lctrl_[idBottom_].y()      = _VELOCITY;
	lctrl_[idBottom_].rot()    = _VELOCITY;	
	lctrl_[idBottom_].xval()   = 0.0;
	lctrl_[idBottom_].yval()   = 0.0;
	lctrl_[idBottom_].rotval() = 0.0;

// top 
	lctrl_[idTop_].x()         = _VELOCITY;
	lctrl_[idTop_].y()         = ymod_;
	lctrl_[idTop_].rot()       = _VELOCITY;	
	lctrl_[idTop_].xval()      = 0.0;
	lctrl_[idTop_].yval()      = -yval_;
	lctrl_[idTop_].rotval()    = 0.0;
	
	


// When a pressure is applied, this is in fact a force from the algorithm point of view !!
	if (lctrl_[idTop_].y()   == _PRESSURE) lctrl_[idTop_].y()   = _FORCE;


// If mode is PRESSURE, the forces have to be computed 
	drive();
	

// Adjust length and positions of rlines
	if (adjust_)
	{
		const double rratio = 0.01;
		const double lratio = 2.;
		
		double x,y,xmin,xmax,ymin,ymax;
		xmin = xmax = spl_->body(2)->x();
		ymin = ymax = spl_->body(2)->y();
		for (unsigned int i=3;i<spl_->lbody().size();++i)
		{
			x = spl_->body(i)->xmin();
			y = spl_->body(i)->ymin();
			xmin = (xmin > x) ? x : xmin;
			ymin = (ymin > y) ? y : ymin;
			
			x = spl_->body(i)->xmax();
			y = spl_->body(i)->ymax();
			xmax = (xmax < x) ? x : xmax;
			ymax = (ymax < y) ? y : ymax;
		} 
		
		double lx = xmax - xmin;
		double ly = ymax - ymin;
		double lmin = (lx > ly) ? ly : lx;
		double r = rratio * lmin;
		
		bottom_->x()   = xmin + 0.5 * lx;
		bottom_->y()   = ymin - r;
		bottom_->rot() = 0.0;
		bottom_->R()   = r;
		bottom_->L()   = lratio * lx;
		
		
		top_->x()   = xmin + 0.5 * lx;
		top_->y()   = ymax + r;
		top_->rot() = 0.0;
		top_->R()   = r;
		top_->L()   = lratio * lx;		
	}
//dof

	//if( ymod_==_VELOCITY)
	//{
	//cout<<"Shear rate init = "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())<<endl;
	//double masstot=0.;
	//double m =0.;
	//spl_->radiusExtrema(4);
	
	//body2d * Rmin=NULL,*Rmax=NULL;

	//for (unsigned int i=4; i< spl_->lbody().size();++i)
	/*{
		
		if(spl_->body(i)->bodyDof() != NULL)
		{
			if(spl_->body(i)->bodyDof()->multiBody())
			masstot = spl_->body(i)->bodyDof()->m();
		}
		else 
		{
			masstot = spl_->body(i)->mass();
		}
		if ( spl_->body(i)->sizeVerlet() == spl_->rmin()) Rmin=spl_->body(i);
		if ( spl_->body(i)->sizeVerlet() > .99* spl_->rmax()) Rmax=spl_->body(i);
	}*/
	
	//if( Rmin==NULL) cout<<" ---------- rmin non retrouve"<<endl;
	//if( Rmax==NULL) cout<<" ---------- rmax non retrouve"<<endl;
	
	/*cout<<" Inertial number : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(masstot/xval_)<<endl;
	cout<<" I smallest particle : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(Rmin->mass()/(xval_))<<endl;
	cout<<" I bigger particle : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(Rmax->mass()/(xval_))<<endl;*/

	//}

	//cout<<"Biax : init ok"<<endl;
	//cout<<" taille dof = "<<ldof().size()<<endl;
	
//	lctrl_[idRight_].print();
//	lctrl_[idTop_].print();
					  
	
}

void Oedometric::drive() 
{
	k_++;
	if (ymod_ == _PRESSURE)  
	{
		lctrl_[idTop_].yval() = -yval_ * (xmax_-xmin_) +k_*(-yval_ * (xmax_-xmin_))*incr_;
		
	}
		
	stress_strain();
		

}

void Oedometric::trans() 
{

}

void Oedometric::share() 
{
	
	
	
}

int Oedometric::check()
{
	return 1;
}

void Oedometric::stress_strain() // stress_strain data writing
{
	
	rline * rl[2];
	rl[0] = dynamic_cast<rline*>(spl_->body(0));
	rl[1] = dynamic_cast<rline*>(spl_->body(1));
	
	bottom_ = rl[0];
	for (unsigned int i=1;i<2;++i) 	
	{	if( rl[i]->rot() == 0)
		bottom_ = (rl[i]->y() < bottom_->y()) ? rl[i] : bottom_;
	}	
	top_ = rl[0];
	for (unsigned int i=1;i<2;++i) 
	{	if( rl[i]->rot() == 0)
		top_    = (rl[i]->y() > top_->y())    ? rl[i] : top_;
	}
	
	

	defy_=(y0_-top_->y())/y0_;


	ofstream stress_strain_app("stress_strain.dat",ios::app);
	stress_strain_app<<-1.0*lctrl_[idTop_].yval()/(xmax_-xmin_)<<"  "<<-1.0*lctrl_[idTop_].yval()<<"   "<<defy_<<endl;
	stress_strain_app.close();

}
