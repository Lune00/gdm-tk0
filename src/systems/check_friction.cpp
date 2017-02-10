#include "check_friction.hpp"

void Check_friction::read_parameters(istream & is)
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
		else if (token == "gravity") gravity_= true;
		else if (token == "}")       break;
		else cerr << "@Check_friction::read_parameters, Unknown parameter: " << token << endl;

		is >> token;
	}
}

void Check_friction::write_parameters(ostream & os)
{
	os << "Check_friction" << endl;
	os << "xval   " << xval_ << endl;
	os << "yval   " << yval_ << endl;
	// TODO: use gdm::modeToString(...) instead
	os << "xmod   " << loadingMode(xmod_) << endl;
	os << "ymod   " << loadingMode(ymod_) << endl;
}

// TODO: init est appelé 2 fois !!!!! a voir...
void Check_friction::init()
	{  
// Check the type of first four bodies
	for (unsigned int i=0;i<2;++i) 
		if (typeid(*(spl_->body(i))) != typeid(polyg))
		gdm::fatal("@Check_friction::init, the first four bodies must be polyg!");
	//spl_->radiusExtrema(4);
	
	one_ = dynamic_cast<polyg*>(spl_->body(0));
	two_ = dynamic_cast<polyg*>(spl_->body(1));


	
	for(unsigned int i=0;i<spl_->lbody().size();++i)
	spl_->body(i)->id()=i;
	
	if ( gravity_) 
	{
		this->gy()=-10.;
		cout<<" gravity activated "<<endl;
	}
	
	//double R= rl[0]->R();
			
	
	
	//cout<<" bas "<<idBottom_<<" haut "<<idTop_<<endl;
	//cout<<" droite "<<idRight_<<" gauche "<<idLeft_<<endl;
	
// We create 4 controls
	for (unsigned int i=0;i<2;++i) lctrl_.push_back(control());

// bottom is fixed
	lctrl_[0].x()      = _VELOCITY;
	lctrl_[0].y()      = _VELOCITY;
	lctrl_[0].rot()    = _VELOCITY;	
	lctrl_[0].xval()   = 0.0;
	lctrl_[0].yval()   = 0.0;
	lctrl_[0].rotval() = 0.0;

// top 
	lctrl_[1].x()         = xmod_;;
	lctrl_[1].y()         = ymod_;
	lctrl_[1].rot()       = _VELOCITY;	
	lctrl_[1].xval()      = -xval_;
	lctrl_[1].yval()      = -yval_;
	lctrl_[1].rotval()    = 0.0;



// When a pressure is applied, this is in fact a force from the algorithm point of view !!
//	if (lctrl_[idRight_].x() == _PRESSURE) lctrl_[idRight_].x() = _FORCE;
//	if (lctrl_[idTop_].y()   == _PRESSURE) lctrl_[idTop_].y()   = _FORCE;
	

// If mode is PRESSURE, the forces have to be computed 
	drive();


//dof


	cout<<"Check Friction: init ok"<<endl;
	
//	lctrl_[idRight_].print();
//	lctrl_[idTop_].print();
					  
	
}

void Check_friction::drive() 
{
	/*if (xmod_ == _PRESSURE)
		lctrl_[idRight_].xval() = -xval_ * (top_->ymin()-bottom_->ymax());

	if (ymod_ == _PRESSURE)  
		lctrl_[idTop_].yval() = -yval_ * (right_->xmin()-left_->xmax());
		
	if ( bottom_->x() +.5* bottom_->L() < right_->x()+3.*right_->R() )
	{
		bottom_->L()*=1.5;
		top_->L()*=1.5;
	}
	*/	
	//pgpg* inter;
	inter2d * inter;
	for( unsigned int i = 0 ; i< nwk_->clist().size();++i)
	{
		inter = nwk_->inter(nwk_->clist(i));
		cout<<inter->fx()<<" "<<inter->fy()<<" "<<inter->second()->vx()<<" "<<inter->second()->vy()<<" "<<
		inter->second()->fx()<<" "<<inter->second()->fy()<<endl;
	}

}

void Check_friction::trans() 
{

}

void Check_friction::share() 
{
	
	
	
}

int Check_friction::check()
{
// check could be done during initilization!
// Maybe this virtual function could be removed...
	return 1;
}


void Check_friction::stress_strain()
{}
