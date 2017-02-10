#include "mpmp.hpp"

void mpmp::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
  is >> *Id1 >> *Id2
  >> nx_ >> ny_
  >> fn_ >> ft_;
}

void mpmp::write(ostream & os)
{
  os << "mpmp " << i_->id() << ' ' << j_->id() << ' '
  << nx_ << ' ' << ny_ << ' ' 
  << fn_ << ' ' << ft_
  << endl << flush;
}

double mpmp::Dist()
{
  double lx = j_->x() - i_->x();
  double ly = j_->y() - i_->y();

  return sqrt(lx*lx+ly*ly);
}

bool mpmp::Activate()
{
  return true;
}

void mpmp::Frame()
{
  double lx = i_->x() - j_->x();
  double ly = i_->y() - j_->y();
  double l  = sqrt(lx * lx + ly * ly);
  double invl = 1.0 / l;  
	
  nx_ = lx * invl;
  ny_ = ly * invl;
  tx_ = -ny_;
  ty_ =  nx_;
  
  x_ = 0.5 * (i_->x() + j_->x());
  y_ = 0.5 * (i_->y() + j_->y());
}

void mpmp::Kin()
{
  double vx = i_->vx() - j_->vx();
  double vy = i_->vy() - j_->vy();
  
  vn_ = vx*nx_ + vy*ny_;
  vt_ = vy*tx_ + vy*ty_;
}

void mpmp::Res()
{		
  double fx = fn_*nx_ + ft_*nx_;
  double fy = fn_*ny_ + ft_*ny_;
		
  i_->fx() += fx;
  i_->fy() += fy;

  j_->fx() -= fx;
  j_->fy() -= fy;
}

void mpmp::Res(const double dfn, const double dft, const double dfs)
{		
  double dfx = dfn*nx_ + dft*nx_;
  double dfy = dfn*ny_ + dft*ny_;
		
  i_->fx() += dfx;
  i_->fy() += dfy;

  j_->fx() -= dfx;
  j_->fy() -= dfy;
}

void mpmp::CDcoeff()
{
  
}

double mpmp::An(const double dt)
{
  return 0.0;
}

double mpmp::At(const double dt)
{
  return 0.0;
}

double mpmp::As(const double dt)
{
  return 0.0;
}

void mpmp::writeMGP(ostream & os)
{
  os << "    <DKDKx antac=\"" << j_->id()+1
     << "\" rn=\"" << fn_
     << "\" rt=\"" << ft_ 
     << "\"/>" << endl;
}
