#include "biaxial.hpp"

void Biaxial::read_parameters(istream & is)
{
	string token;
	string mode;

	is >> token;	
	while(is)
	{	
		if      (token == "xval")    { is >> xval_; }
		else if (token == "yval")    { is >> yval_; }
		else if (token == "xmod")    { is >> mode; xmod_ = loadingMode(mode); } 
		// TODO use gdm::stringToMode(...) instead
		else if (token == "ymod")    { is >> mode; ymod_ = loadingMode(mode); }
		else if (token == "adjust")  adjust_ = true;
		else if (token == "gravity")
		{
			gravity_= true;
		}
		else if (token == "increment") {is >> incr_; increment_ = true;}
		
// Modif du 28/08/09 
		//else if (token == "massFraction") {is >> massFraction_;}
		else if (token == "}")       break;
		else cerr << "@Biaxial::read_parameters, Unknown parameter: " << token << endl;

		is >> token;
	}
}

void Biaxial::write_parameters(ostream & os)
{
	os << "biaxial" << endl;
	os << "xval   " << xval_ << endl;
	os << "yval   " << yval_ << endl;
	// TODO: use gdm::modeToString(...) instead
	os << "xmod   " << loadingMode(xmod_) << endl;
	os << "ymod   " << loadingMode(ymod_) << endl;
}

// TODO: init est appelé 2 fois !!!!! a voir...
void Biaxial::init()
	{  
	k_=0;
// Check the type of first four bodies
	for (unsigned int i=0;i<4;++i) 
		if (typeid(*(spl_->body(i))) != typeid(rline))
		gdm::fatal("@Biaxial::init, the first four bodies must be rlines!");
	spl_->radiusExtrema(4); // Calcul Rmin, Rmax, Rmoy des diks


	ofstream stress_strain_out("stress_strain.dat",ios::out);
	stress_strain_out<<"# p   F    defy"<<endl;
	stress_strain_out.close();

		//gy()=-10.;
// Defining the position of the boundary rlines
	rline * rl[4];
	rl[0] = dynamic_cast<rline*>(spl_->body(0));
	rl[1] = dynamic_cast<rline*>(spl_->body(1));
	rl[2] = dynamic_cast<rline*>(spl_->body(2));
	rl[3] = dynamic_cast<rline*>(spl_->body(3));
	
	for(unsigned int i=0;i<spl_->lbody().size();++i)
	spl_->body(i)->id()=i;
	
	if ( gravity_) this->gy()=-9.81;
				
	bottom_ = rl[0];
	for (unsigned int i=1;i<4;++i) 	
	{	if( rl[i]->rot() == 0)
		bottom_ = (rl[i]->y() < bottom_->y()) ? rl[i] : bottom_;
	}	
	top_ = rl[0];
	for (unsigned int i=1;i<4;++i) 
	{	if( rl[i]->rot() == 0)
		top_    = (rl[i]->y() > top_->y())    ? rl[i] : top_;
	}
	left_ = rl[0];
	for (unsigned int i=1;i<4;++i) 
	{	
		if( rl[i]->rot() != 0)
		left_   = (rl[i]->x() < left_->x())   ? rl[i] : left_;
	}
	right_ = rl[0];
	for (unsigned int i=1;i<4;++i) 
	{	
		if( rl[i]->rot() != 0)
		right_  = (rl[i]->x() > right_->x())  ? rl[i] : right_;
	}
	
	idBottom_ = bottom_->id();
	idTop_    = top_->id();
	idLeft_   = left_->id();
	idRight_  = right_->id();
	
	y0_ = top_->y();
	defy_ = 0.;

	
	
// We create 4 controls
	for (unsigned int i=0;i<4;++i) lctrl_.push_back(control());
//	ldof().push_back(new dof(xmod_, _VELOCITY, _VELOCITY, -xval_, 0., 0.));

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

// left is fixed
	lctrl_[idLeft_].x()        = _VELOCITY;
	lctrl_[idLeft_].y()        = _VELOCITY;
	lctrl_[idLeft_].rot()      = _VELOCITY;
	lctrl_[idLeft_].xval()     = 0.0;
	lctrl_[idLeft_].yval()     = 0.0;
	lctrl_[idLeft_].rotval()   = 0.0;

// right
	//ldof(0)->plugBody(spl_->body(idRight_));
	
	lctrl_[idRight_].x()       = xmod_;
	lctrl_[idRight_].y()       = _VELOCITY;
	lctrl_[idRight_].rot()     = _VELOCITY;	
	lctrl_[idRight_].xval()    = -xval_;
	lctrl_[idRight_].yval()    = 0.0;
	lctrl_[idRight_].rotval()  = 0.0;
	

// When a pressure is applied, this is in fact a force from the algorithm point of view !!
	if (lctrl_[idRight_].x() == _PRESSURE) lctrl_[idRight_].x() = _FORCE;
	if (lctrl_[idTop_].y()   == _PRESSURE) lctrl_[idTop_].y()   = _FORCE;

// If mode is PRESSURE, the forces have to be computed 
	drive();


		const double rratio = 0.1;
		const double lratio = 3.;
		
		double x,y,xmin,xmax,ymin,ymax;
		xmin = spl_->body(4)->xmin();
		xmax = spl_->body(4)->xmax();
		ymin = spl_->body(4)->ymin();
		ymax = spl_->body(4)->ymax();

		for (unsigned int i=5;i<spl_->lbody().size();++i)
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
		cout<<"xmin:="<<xmin<<"  xmax:= "<<xmax<<" ymin:= "<<ymin<<"  ymax:= "<<ymax<<endl;
		cout<<"lx:= "<<lx<<"  "<<"ly:= "<<ly<<endl;

// Adjust length and positions of rlines
		
	if (adjust_)
	{	
		double lmin = (lx > ly) ? ly : lx;
		double r = rratio * lmin;
		bottom_->x()   = xmin + 0.5 * lx;
		bottom_->y()   = ymin - r;
		bottom_->rot() = 0.0;
		bottom_->R()   = r;
		bottom_->L()   = lratio * lx;
		
		left_->x()   = xmin - r;
		left_->y()   = ymin + 0.5 * ly;
		left_->rot() = 0.5 * M_PI;
		left_->R()   = r;
		left_->L()   = lratio * ly;
		
		right_->x()   = xmax + r;
		right_->y()   = ymin + 0.5 * ly;
		right_->rot() = 0.5 * M_PI;
		right_->R()   = r;
		right_->L()   = lratio * ly;
		
		top_->x()   = xmin + 0.5 * lx;
		top_->y()   = ymax + r;
		top_->rot() = 0.0;
		top_->R()   = r;
		top_->L()   = lratio * lx;	
	}
//dof

	if( ymod_==_VELOCITY)
	{
	cout<<"Shear rate init = "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())<<endl;
	double masstot=0.;
	//double m =0.;
	//spl_->radiusExtrema(4);
	
	body2d * Rmin=NULL,*Rmax=NULL;

	for (unsigned int i=4; i< spl_->lbody().size();++i)
	{
		
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
	}
	
	if( Rmin==NULL) cout<<" ---------- rmin non retrouve"<<endl;
	if( Rmax==NULL) cout<<" ---------- rmax non retrouve"<<endl;
	
	cout<<" Inertial number : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(masstot/xval_)<<endl;
	cout<<" I smallest particle : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(Rmin->mass()/(xval_))<<endl;
	cout<<" I bigger particle : "<<fabs(lctrl_[idTop_].yval())/(top_->y()-bottom_->y())*sqrt(Rmax->mass()/(xval_))<<endl;

	}
	
	cout<<"Biax : init ok"<<endl;
	cout<<" taille dof = "<<ldof().size()<<endl;
}

void Biaxial::drive() 
{
	k_++;
	
	if (xmod_ == _PRESSURE)
	{	
		lctrl_[idRight_].xval() = -xval_ * (top_->ymin()-bottom_->ymax());
	}
	if (ymod_ == _PRESSURE)  
	{
		lctrl_[idTop_].yval() = -yval_ * (right_->xmin()-left_->xmax());
		if(increment_)
		{
			lctrl_[idTop_].yval() = -yval_ * (right_->xmin()-left_->xmax()) +k_*(-yval_ * (right_->xmin()-left_->xmax()))*incr_;
		}
	}
		
	if ( bottom_->x() +.5* bottom_->L() < right_->x()+3.*right_->R() )
	{
		bottom_->L()*=1.5;
		top_->L()*=1.5;
	}
		
if(increment_)
{stress_strain();}
}

void Biaxial::trans() 
{

}

void Biaxial::share() 
{
	
}

int Biaxial::check()
{
// check could be done during initilization!
// Maybe this virtual function could be removed...
	return 1;
}


void Biaxial::stress_strain()
{
	ofstream stress_strain_app("stress_strain.dat",ios::app);
	stress_strain_app<<-1.0*lctrl_[idTop_].yval()<<"   "<<-1.0*lctrl_[idTop_].yval()/(right_->xmin()-left_->xmax())<<"  "<<defy_<<endl;
	stress_strain_app.close();
}
