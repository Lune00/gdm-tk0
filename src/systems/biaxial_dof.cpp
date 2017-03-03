#include "biaxial_dof.hpp"

void Biaxial_dof::read_parameters(istream & is)
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
		else if (token == "gravity") gravity_= true;
		else if (token == "}")       break;
		else cerr << "@Biaxial_dof::read_parameters, Unknown parameter: " << token << endl;
		
		is >> token;
	}
}

void Biaxial_dof::write_parameters(ostream & os)
{
	os << "Biaxial_dof" << endl;
	os << "xval   " << xval_ << endl;
	os << "yval   " << yval_ << endl;
	// TODO: use gdm::modeToString(...) instead
	os << "xmod   " << loadingMode(xmod_) << endl;
	os << "ymod   " << loadingMode(ymod_) << endl;
}

// TODO: init est appelé 2 fois !!!!! a voir...
void Biaxial_dof::init()
	{  
// Check the type of first four bodies
	for (unsigned int i=0;i<4;++i) 
		if ((*(spl_->body(i))).type() != _type_rline)
		gdm::fatal("@Biaxial_dof::init, the first four bodies must be rlines!");
	spl_->radiusExtrema(4);
	cout<<" biax dof init "<<endl;

		//gy()=-10.;
// Defining the position of the boundary rlines
	rline * rl[4];
	rl[0] = dynamic_cast<rline*>(spl_->body(0));
	rl[1] = dynamic_cast<rline*>(spl_->body(1));
	rl[2] = dynamic_cast<rline*>(spl_->body(2));
	rl[3] = dynamic_cast<rline*>(spl_->body(3));
	
	for(unsigned int i=0;i<spl_->lbody().size();++i)
	spl_->body(i)->id()=i;
	
	if ( gravity_) this->gy()=-10.;
	
	//double R= rl[0]->R();
			
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
	
	cout<<" bas "<<idBottom_<<" haut "<<idTop_<<endl;
	cout<<" droite "<<idRight_<<" gauche "<<idLeft_<<endl;
	
// We create 4 dof
	//for (unsigned int i=0;i<4;++i) lctrl_.push_back(control());
//	ldof().push_back(new dof(xmod_, _VELOCITY, _VELOCITY, -xval_, 0., 0.));
	cout<<" biax dof init "<<endl;

// bottom is fixed
	ldof().push_back( new dof(_VELOCITY,_VELOCITY,_VELOCITY,0.,0.,0.));
	ldof().back()->plugBody(spl_->body(idBottom_));
	ldof().back()->id()=0;
	

// top 
	ldof().push_back( new dof(_VELOCITY,ymod_,_VELOCITY,0.,-yval_,0.));
	ldof().back()->plugBody(spl_->body(idTop_));
	ldof().back()->id()=1;
	
	
	/*lctrl_[idTop_].x()         = _VELOCITY;
	lctrl_[idTop_].y()         = ymod_;
	lctrl_[idTop_].rot()       = _VELOCITY;	
	lctrl_[idTop_].xval()      = 0.0;
	lctrl_[idTop_].yval()      = -yval_;
	lctrl_[idTop_].rotval()    = 0.0;*/

// left is fixed
	ldof().push_back( new dof(_VELOCITY,_VELOCITY,_VELOCITY,0.,0.,0.));
	ldof().back()->plugBody(spl_->body(idLeft_));
	ldof().back()->id()=2;
	ldof().back()->mcrot()=spl_->body(idLeft_)->rot();
	cout<<spl_->body(idLeft_)->rot()<<endl;
	
	/*lctrl_[idLeft_].x()        = _VELOCITY;
	lctrl_[idLeft_].y()        = _VELOCITY;
	lctrl_[idLeft_].rot()      = _VELOCITY;
	lctrl_[idLeft_].xval()     = 0.0;
	lctrl_[idLeft_].yval()     = 0.0;
	lctrl_[idLeft_].rotval()   = 0.0;*/

// right
	//ldof(0)->plugBody(spl_->body(idRight_));
	ldof().push_back( new dof(xmod_,_VELOCITY,_VELOCITY,-xval_,0.,0.));
	ldof().back()->plugBody(spl_->body(idRight_));
	ldof().back()->id()=3;
	ldof().back()->mcrot()=spl_->body(idRight_)->rot();
	
	for(unsigned int i=0;i<ldof().size();++i)
	{
		ldof(i)->multiBody()=false;
	}
/*	lctrl_[idRight_].x()       = xmod_;
	lctrl_[idRight_].y()       = _VELOCITY;
	lctrl_[idRight_].rot()     = _VELOCITY;	
	lctrl_[idRight_].xval()    = -xval_;
	lctrl_[idRight_].yval()    = 0.0;
	lctrl_[idRight_].rotval()  = 0.0;*/
	

// When a pressure is applied, this is in fact a force from the algorithm point of view !!
	if (spl_->body(idRight_)->bodyDof()->x() == _PRESSURE) spl_->body(idRight_)->bodyDof()->x() = _FORCE;
    if (spl_->body(idTop_)->bodyDof()->y() == _PRESSURE) spl_->body(idTop_)->bodyDof()->y() = _FORCE;
	

// If mode is PRESSURE, the forces have to be computed 
	drive();

// Adjust length and positions of rlines
	if (adjust_)
	{
		const double rratio = 0.01;
		const double lratio = 2.;
		
		double x,y,xmin,xmax,ymin,ymax;
		xmin = xmax = spl_->body(4)->x();
		ymin = ymax = spl_->body(4)->y();
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
	cout<<"Shear rate init = "<<fabs(spl_->body(idTop_)->bodyDof()->yval())/(top_->y()-bottom_->y())<<endl;
	double masstot=0.;
	//spl_->radiusExtrema(4);
	
	body2d * Rmin=NULL,*Rmax=NULL;
	for (unsigned int i=4; i< spl_->lbody().size();++i)
	{
		masstot+=spl_->body(i)->mass();
		if ( spl_->body(i)->sizeVerlet() == spl_->rmin()) Rmin=spl_->body(i);
		if ( spl_->body(i)->sizeVerlet() > .99* spl_->rmax()) Rmax=spl_->body(i);
	}
	if( Rmin==NULL) cout<<" ---------- rmin non retrouve"<<endl;
	if( Rmax==NULL) cout<<" ---------- rmax non retrouve"<<endl;
	
	cout<<" Inertial number : "<<fabs(spl_->body(idTop_)->bodyDof()->yval())/(top_->y()-bottom_->y())*sqrt(masstot/xval_)<<endl;
	cout<<" I smallest particle : "<<fabs(spl_->body(idTop_)->bodyDof()->yval())/(top_->y()-bottom_->y())*sqrt(Rmin->mass()/(xval_))<<endl;
	cout<<" I bigger particle : "<<fabs(spl_->body(idTop_)->bodyDof()->yval())/(top_->y()-bottom_->y())*sqrt(Rmax->mass()/(xval_))<<endl;
	
	}
	
	cout<<"Biax : init ok"<<endl;
	cout<<" taille dof = "<<ldof().size()<<endl;
	
//	lctrl_[idRight_].print();
//	lctrl_[idTop_].print();
/*	for(unsigned int i=0;i<ldof().size();++i)
	{
		ldof(i)->print();
	}
	*/				  
	
}

void Biaxial_dof::drive() 
{
	//cout<<"drive "<<endl;
	
	if (xmod_ == _PRESSURE)
		spl_->body(idRight_)->bodyDof()->xval() = -xval_ * (top_->ymin()-bottom_->ymax());

	if (ymod_ == _PRESSURE)  
		spl_->body(idTop_)->bodyDof()->yval() = -yval_ * (right_->xmin()-left_->xmax());
		
	if ( bottom_->x() +.5* bottom_->L() < right_->x()+3.*right_->R() )
	{
		bottom_->L()*=1.5;
		top_->L()*=1.5;
	}
		
	//cout<<"drive ok"<<endl;

}

void Biaxial_dof::trans() 
{

}

void Biaxial_dof::share() 
{
	
	
	
}

int Biaxial_dof::check()
{
// check could be done during initilization!
// Maybe this virtual function could be removed...
	return 1;
}


void Biaxial_dof::stress_strain()
{}
