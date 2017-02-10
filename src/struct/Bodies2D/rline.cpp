#include "rline.hpp"

void rline::read (istream & is)
{
  is >> grp_ >> R_ >> L_ 
     >> x_  >> y_  >> rot_ 
     >> vx_ >> vy_ >> vrot_ ;
  x0_=x_;
  y0_=y_;
  rot0_=rot_;

}

void rline::write (ostream & os)
{
  os << "rline " << grp_ << ' ' << R_ << ' ' << L_ << ' '
     << x_  << ' ' << y_  << ' ' << rot_  << ' ' 
     << vx_ << ' ' << vy_ << ' ' << vrot_ 
     << endl << flush;
}

void rline::writeMGP (ostream & os)
{
  os << "    <JONCx  ax1=\"" <<  0.5 * L_ << "\" ax2=\"" << R_ << "\">" << endl
     << "     <position x=\"" << x_ 
	 << "\" y=\"" << y_ << "\" rot=\"" << rot_ << "\"/>" << endl 
     << "    </JONCx>" << endl << flush;
}

void rline::writePS (ostream & os)
{
  //os << "";
}

void rline::writeM(ostream &os)
{
    
}

rline* rline::duplicate()
{
  rline* copy = new rline;
  *copy = *this;
  return copy;
}

void rline::Fill(double density)
{
  double R2   = R_ * R_;
  double Area = 2.0 * L_ * R_ + M_PI * R2;
  mass_ = Area * density;
  
  

  //_mom = density*(M_PI*R2*(0.25*L2 + 0.5*R2)
	//	+ _L*_R*(0.25*L2+R2))/Area; 
  
  // la valeur de _mom est COMPLETEMENT FAUSSE
  // A CALCULER ...
  mom_ = mass_/3.*(R2+L_*L_) //inertia of rectangle
       + 2.*( (M_PI/4.-8./(9.*M_PI))*R2*R2 // inertia of half circle in mass center of half circle
	   + M_PI*R2/2. * pow( L_+ 4*R_/(3*M_PI),2));// moment transport formula
}
