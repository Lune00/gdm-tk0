#include "tensor.hpp"

void gdm::Tensor2x2::print()
      {
      cout << name_ << " =" << endl;
      cout << '[' << endl;
      cout << ' ' << xx_ << '\t' << xy_ << endl;
      cout << ' ' << yx_ << '\t' << yy_ << endl;
      cout << ']' << endl;
	cout <<" v1 = "<<l1_<< " v2 = "<<l2_<<endl;
      }

void gdm::Tensor2x2::symmetrize()
      {
      xy_ = 0.5 * (xy_ + yx_);
      yx_ = xy_;
      }


void gdm::Tensor2x2::scalarMult( double k)
      {
      xx_ *= k;
      xy_ *= k;
	  yy_ *= k;
	  yx_ *= k;
      }

void gdm::Tensor2x2::eigenValues()
{
	double b= -yy_ - xx_;
	double c=-xy_*yx_+xx_*yy_;
	double d=b*b-4*c;
	//cout<<"delta "<<d<<endl;
	l1_=(-b-sqrt(d))/2.;
	l2_=(-b+sqrt(d))/2.;
	//cout<<"l1 = "<<l1_<<"  l2_= "<<l2_<<endl;
	
}
void gdm::Tensor2x2::eigenVectors()
{
	this->eigenValues();
	
	double norme;
	if (l1_==0 && l2_==0) 
	{//cout<<" val propre nulle"<<endl;
	}
	else
	{
	vx1_=xy_ / ( - xx_+ l1_);
	vy1_=1.;
	norme=sqrt(vx1_*vx1_ + vy1_*vy1_);
	vx1_/=norme;
	vy1_/=norme;
	
	vx2_=xy_/( - xx_ + l2_);
	vy2_=1.;
	norme=sqrt(vx2_*vx2_ + vy2_*vy2_);
	vx2_/=norme;
	vy2_/=norme;
	}
}

double gdm::Tensor2x2::majorDirection( )

{
	this->eigenVectors();
	//cout<<endl<<" l1 "<<l1_<<" l2 "<<l2_<<endl;
	//cout<<" v1 "<<(vy1_ < 0. ? 2.*M_PI-acos(vx1_): acos(vx1_))/M_PI*180.<<endl;
	//cout<<" v2 "<<(  vy2_ < 0. ? 2.*M_PI-acos(vx2_): acos(vx2_))/M_PI*180.<<endl;
	if ( l1_ >= l2_)
		return (  vy1_ < 0. ? 2.*M_PI-acos(vx1_): acos(vx1_)); 
	else 
		return (  vy2_ < 0. ? 2.*M_PI-acos(vx2_): acos(vx2_)); 
}

