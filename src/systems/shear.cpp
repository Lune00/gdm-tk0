#include "shear.hpp"

void Shear::read_parameters(istream & is)
{
	string token;
	string mode;

	is >> token;	
	while(is)
	{	
		if      (token == "pressure") { is >> pressure_; }
		else if (token == "rate")     { is >> rate_; }
		else if (token == "adjust")  adjust_ = true;
		else if (token == "}")       break;
		else cerr << "@Shear::read_parameters, Unknown parameter: " << token << endl;
        cout<<"void shear"<<endl;

		is >> token;
	}
}

void Shear::write_parameters(ostream & os)
{
	os << "shear" << endl;
	os << "pressure " << pressure_ << endl;
	os << "rate     " << rate_ << endl;
}

// TODO: init est appelé 2 fois !!!!! a voir...
void Shear::init()
{  
// Some checks
	for (unsigned int i=0;i<4;++i) 
		if ((*(spl_->body(i))).type() != _type_rline)
		gdm::fatal("@Shear::init, the first four bodies must be rlines!");

// y de left et right doivent etre les memes

// Defining the position of the boundary rlines
// WARNING: these tests are simplistic
	rline * rl[4];
	rl[0] = dynamic_cast<rline*>(spl_->body(0));
	rl[1] = dynamic_cast<rline*>(spl_->body(1));
	rl[2] = dynamic_cast<rline*>(spl_->body(2));
	rl[3] = dynamic_cast<rline*>(spl_->body(3));
			
	bottom_ = rl[0];
	for (unsigned int i=1;i<4;++i) bottom_ = (rl[i]->y() < bottom_->y()) ? rl[i] : bottom_;
	top_ = rl[0];
	for (unsigned int i=1;i<4;++i) top_    = (rl[i]->y() > top_->y()) ? rl[i] : top_;
	left_ = rl[0];
	for (unsigned int i=1;i<4;++i) left_   = (rl[i]->x() < left_->x()) ? rl[i] : left_;
	right_ = rl[0];
	for (unsigned int i=1;i<4;++i) right_  = (rl[i]->x() > right_->x()) ? rl[i] : right_;
	
	idBottom_ = bottom_->id();
	idTop_    = top_->id();
	idLeft_   = left_->id();
	idRight_  = right_->id();
	
// We create 4 controls
	for (unsigned int i=0 ; i<4 ; ++i) lctrl_.push_back(control());

// Bottom is fixed
	lctrl_[idBottom_].x()      = _VELOCITY;
	lctrl_[idBottom_].y()      = _VELOCITY;
	lctrl_[idBottom_].rot()    = _VELOCITY;	
	lctrl_[idBottom_].xval()   = 0.0;
	lctrl_[idBottom_].yval()   = 0.0;
	lctrl_[idBottom_].rotval() = 0.0;

// Pressure on the top 
	lctrl_[idTop_].x()         = _VELOCITY;
	lctrl_[idTop_].y()         = _FORCE;
	lctrl_[idTop_].rot()       = _VELOCITY;	
	lctrl_[idTop_].xval()      = 0.0;
	lctrl_[idTop_].yval()      = -pressure_ * (right_->x()-left_->x()-right_->R()-left_->R());
	lctrl_[idTop_].rotval()    = 0.0;

// left 
	lctrl_[idLeft_].x()        = _VELOCITY;
	lctrl_[idLeft_].y()        = _VELOCITY;
	lctrl_[idLeft_].rot()      = _VELOCITY;
	lctrl_[idLeft_].xval()     = 0.0;
	lctrl_[idLeft_].yval()     = 0.0;
	lctrl_[idLeft_].rotval()   = -rate_;

// right
	lctrl_[idRight_].x()       = _VELOCITY;
	lctrl_[idRight_].y()       = _VELOCITY;
	lctrl_[idRight_].rot()     = _VELOCITY;	
	lctrl_[idRight_].xval()    = 0.0;
	lctrl_[idRight_].yval()    = 0.0;
	lctrl_[idRight_].rotval()  = -rate_;

// Rotation centers
	xc0_ = left_->x() - (left_->y()-bottom_->y()) * tan(0.5*M_PI-left_->rot()); 
	yc0_ = bottom_->y();
	

}

void Shear::drive() // drive(double dt)
{
	// en fait il faudrait recuperer dt !!

	double dt = 0.0001;
	
	double c = cos(-rate_ * dt);
	double s = sin(-rate_ * dt);
	double X0 = left_->x() - xc0_;
	double Y0 = left_->y() - yc0_;
	double dx = (c * X0 - s * Y0) - X0;
	double dy = (s * X0 + c * Y0) - Y0;
	double dxdt = dx / dt;
	double dydt = dy / dt;
	
	//left_->x()  += dx;
	//left_->y()  += dy;
	left_->vx() += dxdt;
	left_->vy() += dydt;
	
	//right_->x() += dx;
	//right_->y() += dy;
	right_->vx() += dxdt;
	right_->vy() += dydt;
}

void Shear::trans() 
{

}

void Shear::share() 
{

}

int Shear::check()
{
// check could be done during initilization!
// Maybe this virtual function could be removed...
	return 1;
}

void Shear::stress_strain()
{}
