#include "massPoint.hpp"

void massPoint::read (istream & is)
{
  is >> grp_ >> R_ 
     >> x_  >> y_
     >> vx_ >> vy_;
}

void massPoint::write (ostream & os)
{
  os << "massPoint " << grp_ << ' ' << R_ << ' '
     << x_  << ' ' <<  y_ << ' ' 
     << vx_ << ' ' << vy_
     << endl << flush;
}

void massPoint::writeMGP (ostream & os)
{
  os << "    <DISKx id=\"" << id_ << "\" r=\"" << R_ << "\">" << endl
     << "     <position x=\"" << x_ << "\" y=\"" << y_ << "\"/>" << endl
     << "     <velocity x=\"" << vx_ << "\" y=\"" << vy_ << "\"/>" << endl	 
     << "    </DISKx>" << endl << flush;
}

void massPoint::writePS (ostream & os)
{
  os << "";
}

void massPoint::writeM(ostream &os)
{
    
}


massPoint* massPoint::duplicate()
{
  massPoint* copy = new massPoint;
  *copy = *this; 
  return copy;
}

void massPoint::Fill(double density)
{
  double Area = M_PI * R_ * R_;

  mass_ = Area * density;
  mom_ = 0.5 * R_ * R_ * mass_;
}
