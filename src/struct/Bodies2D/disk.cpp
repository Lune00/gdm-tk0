#include "disk.hpp"

void disk::read (istream & is)
{
  is >> grp_ >> R_ 
     >> x_  >> y_  >> rot_
     >> vx_ >> vy_ >> vrot_;
  x0_=x_;
  y0_=y_;
  rot0_=rot_;
}

void disk::write (ostream & os)
{
  os << "disk " << grp_ << ' ' << R_ << ' '
     << x_  << ' ' <<  y_  << ' ' << rot_  << ' ' 
     << vx_ << ' ' << vy_  << ' ' << vrot_
     << endl << flush;
}

void disk::writeM (ostream & os)
{
    os << R_ << ' '<< x_  << ' ' <<  y_  << ' ' << rot_  << ' '
    << vx_ << ' ' << vy_ << vrot_
    << endl << flush;
}

void disk::writeMGP (ostream & os)
{
  os << "    <DISKx id=\"" << id_ << "\" r=\"" << R_ << "\">" << endl
     << "     <position x=\"" << x_ 
     << "\" y=\"" << y_ << "\" rot=\"" << rot_ << "\"/>" << endl
     << "     <velocity x=\"" << vx_ 
     << "\" y=\"" << vy_ << "\" rot=\"" << vrot_ << "\"/>" << endl	 
     << "    </DISKx>" << endl << flush;
}

void disk::writePS (ostream & os)
{
  os << "";
}

disk* disk::duplicate()
{
  disk* copy = new disk;
  *copy = *this; 
  return copy;
}

void disk::Fill(double density)
{
  double Area = M_PI * R_ * R_;

  mass_ = Area * density;
  mom_ = 0.5 * R_ * R_ * mass_;
}
