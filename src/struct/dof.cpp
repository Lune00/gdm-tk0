#include "dof.hpp"

void dof::check()
{
	if (x_ == _FORCE && y_ == _FORCE && rot_ == _FORCE
		&& xval_ == 0.0 && yval_ == 0.0 && rotval_ == 0.0)
	{
		gdm::warning("@control::check, controlled body is free !");
	}
}

void dof::print()
{
	cout<< "dof id : "<<id_<<endl;
	cout <<"x  : ";
	if (x_ == _FORCE)    cout << "force    ";
	if (x_ == _VELOCITY) cout << "velocity ";
	cout << xval_ << endl;

	cout << "y  : ";
	if (y_ == _FORCE)    cout << "force    ";
	if (y_ == _VELOCITY) cout << "velocity ";
	cout << yval_ << endl; 

	cout << "rot: ";
	if (rot_ == _FORCE)    cout << "force    ";
	if (rot_ == _VELOCITY) cout << "velocity ";
	cout << rotval_ << endl; 
	cout<<"Current values = "<<endl;
	cout<<" mcx: "<<mcx_<<endl;
	cout<<" mcy: "<<mcy_<<endl;
	cout<<" mcrot: "<<mcrot_<<endl;
	cout<<"Gravity : xG= "<<xGravity_<<"  yG= "<<yGravity_<<endl;
	cout<<"Mass    : "<<m_<<endl;
	cout<<"Controled bodies : ";
	for (unsigned int i=0;i<ctrlBodies_.size();++i)
	{
		cout<<ctrlBodies_[i]->id()<<" ";
	}
	cout<<endl<<endl;
}

void dof::exportBodyId(string includeFile)
{
	ostringstream name;
	string fname;
	name<<"./dof_"<<this->id()<<"_"<<includeFile<<".dof";
	fname=name.str();
	cout<<"Ecriture de "<<fname<<endl;
	ofstream output ( fname.c_str() , ios::out | ios::trunc);
	//output.open( fname.c_str() , ios::out | ios::in);
	//Conversion de ostringstream en string puis en const char*
	if ( ! output ) cout<<"erreur creation de "<<fname.c_str()<<endl;
	else
		for (unsigned int i=0;i<lctrlBodies().size();++i)
	{
		output<< ctrlBody(i)->id()<<" ";
	}

	output.close();
}

unsigned int dof::writeSpl( ostream & os)
{	
	unsigned int i;
	
	for(i=0;i<lctrlBodies().size();++i)
	{
		ctrlBody(i)->write(os);
	}
	if (area_ != -1)
		os<<"area "<<area_<<endl;
	return (i-1);
}

void dof::computeMassCenter()//fonction utilisee dans la fonction membre insert_cluster_in_disk; classe alteration.hpp
{
//	cout<<" computeMass ";
	//cout<<mcx_<<" "<<mcy_<<endl;

	if( multiBody_)//plusieurs bodies controles par un dof
	{
		double xt=0.,yt=0.;
		for( unsigned int i=0;i<lctrlBodies().size();++i)//lctrlBodies = nombre de bodies controlés par un dof
		{
			xt += ctrlBody(i)->x() * ctrlBody(i)->mass();// calcul des coordonnées barycentriques du centre de masse
			yt += ctrlBody(i)->y() * ctrlBody(i)->mass();
		}
		mcx_ = xt/m_;
		mcy_ = yt/m_;
		computeVertex();
	}
	else
	{
		mcx_=ctrlBody(0)->x();
		mcy_=ctrlBody(0)->y();
		mcrot_=ctrlBody(0)->rot();
	}
	//cout<<"Mass center "<<mcx_<<" "<<mcy_<<endl;
}

void dof::fill( double density)//il faut faire la fonction fill() fonctionner pour le cas multibody
{

	m_=0.;
	if( area_ == -1 )
		for( unsigned int i=0;i<lctrlBodies().size();++i)
	{
		m_ += ctrlBody(i)->mass();
	}
	else
	{
		m_ = area_ * density;
		for( unsigned int i=0;i<lctrlBodies().size();++i)
		{
			ctrlBody(i)->mass() = m_ / lctrlBodies().size();
		}
	}

}
void dof::computeVertex()
{
	vertex_.clear();
	if( multiBody_)
	{
		for( unsigned int i=0;i<lctrlBodies().size();++i)
		{
			vertex_.push_back( gdm::vertex(ctrlBody(i)->x() - mcx_, ctrlBody(i)->y() - mcy_) );
	//	if( id_==1) cout<<"Vertex "<<vertex_.back().x()<<" "<<vertex_.back().y()<<endl;
		}
	}
	//cout<<" mcx_ = "<<mcx_<<" mcy_ = "<<mcy_<<" mom_ = "<<mom_<<endl;

}

void dof::ComputeImposedMass(double spl_mass_tot_)
{
    
	double surf_ = 0.;
	for( unsigned int i=0;i<lctrlBodies().size();++i)
    surf_ += ctrlBody(i)->Area();
    
	for( unsigned int i=0;i<lctrlBodies().size();++i)
	{
		ctrlBody(i)->mass() = spl_mass_tot_*ctrlBody(i)->Area()/surf_;
		ctrlBody(i)->Fill(ctrlBody(i)->mass()/ctrlBody(i)->Area());
	}
	
	m_ = spl_mass_tot_;
	computeMassCenter();
	computeVertex();
	computeMoment();
    
}

void dof::computeMoment()
{
	mom_=0.;
	for( unsigned int i=0;i<lctrlBodies().size();++i)
	{
		mom_ += ctrlBody(i)->mom() + 
			ctrlBody(i)->mass() * ( pow(ctrlBody(i)->x() -mcx_,2)+ pow(ctrlBody(i)->y() -mcy_,2) );
	}
}

body2d * dof::rightBody()
{
	double xmax=ctrlBody(0)->xmax();
	body2d * temp=ctrlBody(0);
	for( unsigned int i=1;i<lctrlBodies().size();++i)
	{
		if( ctrlBody(i)->xmax() > xmax )
		{
			xmax=ctrlBody(i)->xmax();
			temp=ctrlBody(i);
		} 
	}
	return temp;
}

body2d * dof::leftBody()
{
	double xmin=ctrlBody(0)->xmin();
	body2d * temp=ctrlBody(0);
	for( unsigned int i=1;i<lctrlBodies().size();++i)
	{
		if( ctrlBody(i)->xmin() < xmin )
		{
			xmin=ctrlBody(i)->xmin();
			temp=ctrlBody(i);
		} 
	}
	return temp;
}

void dof::translateBodies( double rb, double lb, double P )
{
	//cout<<"translate";
	for( unsigned int i=0;i<lctrlBodies().size();++i)
	{
		if( ctrlBody(i)->x() > rb)
			ctrlBody(i)->x() -=P ;

		if( ctrlBody(i)->x() < lb)
			ctrlBody(i)->x() +=P ;
	}
	computeMassCenter();
	computeVertex();
	computeMoment();
}
void dof::translate( double P )
{
	//cout<<"translate dof"<<endl;
	for( unsigned int i=0;i<lctrlBodies().size();++i)
	{
		ctrlBody(i)->x() +=P ;

	}
	mcx_ +=P;
	//computeMassCenter();
}

body2d* dof::lowerBody()
{
	body2d* temp = ctrlBodies_[0];
	for ( unsigned int i=1;i<ctrlBodies_.size();++i)
	{
		if(  ctrlBodies_[i]->ymin() < temp->ymin() ) temp = ctrlBodies_[i];
	}
	return (temp);	

}


void dof::plugBody(unsigned int N,string includeFile,vector <unsigned int> & idBodies)
{
	ostringstream name;
	string fname;
	unsigned int Ntemp;

	name<<"./dof_"<<N<<"_"<<includeFile<<".dof";
	fname=name.str();
	cout<<"Ouverture de "<<fname<<endl;
	ifstream input ( fname.c_str() , ios::in);
	if ( ! input ) cout<<" *****************************ouverture impossible "<<fname<<endl;
	else
	{
		while( true )
		{
			input>>Ntemp;
			idBodies.push_back(Ntemp);
			//cout<<Ntemp<<endl;
			if ( input.eof() ) break;
		}
		//cout<<endl;
		idBodies.pop_back();
	}
	input.close();
	cout<<"lecture dof OK "<<endl;
}

void dof::plugBody( body2d * bod)
{
	this->lctrlBodies().push_back( bod);
	bod->bodyDof()=this;
	m_+=bod->mass();
}
void dof::affect(unsigned int x,unsigned int y,unsigned int rot,double xval, double yval, double rotval)
{
	x_   = x;
	y_   = y;
	rot_ = rot;

	xval_   = xval;
	yval_   = yval;
	rotval_ = rotval;
}

void dof::resDof(double dfx, double dfy, double dfrot,double xi,double yi)
{
	mcfx_ +=  dfx;
	mcfy_ +=  dfy;

	double dx = xi - mcx_;
	double dy = yi - mcy_;
	//cout<<dx<<endl;

	mcfrot_ += dfy*dx - dfx*dy;
}


void dof::imposeForce()
{

	if( x() == _FORCE)
		mcfx_   = xval_ + xGravity_ * m_;

	if( y() == _FORCE)
		mcfy_   = yval_ + yGravity_ * m_;//GRAVITY

	if( rot() == _FORCE)
		mcfrot_   = rotval_ ;

}
void dof::imposeVelocity()
{
	if( x()   == _VELOCITY)     mcvx_   = xval_;
	if( y()   == _VELOCITY)     mcvy_   = yval_ ;   
	if( rot() == _VELOCITY)     mcvrot_ = rotval_ ;

	if( multiBody_)
	{
		double dx,dy;
		for( unsigned int i=0;i < lctrlBodies().size();++i)
		{
			dx= ctrlBody(i)->x() - mcx_ ;
			dy= ctrlBody(i)->y() - mcy_ ;

			ctrlBody(i)->vx() = mcvx_ - dy*mcvrot_;//doute..............
			ctrlBody(i)->vy() = mcvy_ + dx*mcvrot_;
			ctrlBody(i)->vrot() = 0.;
		}
	}

}
void dof::move( double dt)
{
	mcx_   += mcvx_ * dt;
	mcy_   += mcvy_ * dt;
	mcrot_ += mcvrot_ * dt;

	if( multiBody_)
	{
		double c = cos(mcrot_);
		double s = sin(mcrot_);
		for (unsigned int i = 0 ; i < lctrlBodies().size() ; ++i)	
		{		
			ctrlBody(i)->x()   = mcx_ + c * vertex_[i].x() - s * vertex_[i].y();
			ctrlBody(i)->y()   = mcy_ + s * vertex_[i].x() + c * vertex_[i].y();
			ctrlBody(i)->rot() = mcrot_;
		}
	}
	else
	{
		ctrlBody(0)->x()   = mcx_ ;
		ctrlBody(0)->y()   = mcy_ ;
		ctrlBody(0)->rot() = mcrot_;
	}

	//if( isPeriodic_ ) computeMassCenter();
}
//Increment de vitesse du au chargement du dof
void dof::imposeVelocityOfForce(double dt)
{
	if( x()   == _FORCE) mcvx_   += dt * mcfx_ / m_;
	if( y()   == _FORCE) mcvy_   += dt * mcfy_ / m_;
	if( rot() == _FORCE) mcvrot_ += dt * mcfrot_ /mom_;

	if( multiBody_)
	{
		double dx,dy;
		for( unsigned int i=0;i<lctrlBodies().size();++i)
		{
			dx= ctrlBody(i)->x() - mcx_ ;
			dy= ctrlBody(i)->y() - mcy_ ;

			ctrlBody(i)->vx() = mcvx_ - dy*mcvrot_;
			ctrlBody(i)->vy() = mcvy_ + dx*mcvrot_;
			ctrlBody(i)->vrot() = 0.;
		}
	}


}


void dof::imposeForceOfVelocity()
{

	if( x()   == _VELOCITY) mcfx_   = 0.0;
	if( y()   == _VELOCITY) mcfy_   = 0.0;
	if( rot() == _VELOCITY) mcfrot_ = 0.0;	

}


