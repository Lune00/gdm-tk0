#include "control.hpp"

void control::check()
{
  if (x_ == _FORCE && y_ == _FORCE && rot_ == _FORCE
      && xval_ == 0.0 && yval_ == 0.0 && rotval_ == 0.0)
    {
    gdm::warning("@control::check, controlled body is free !");
    }
}

void control::print()
{
  cout << "x: ";
  if (x_ == _FORCE)    cout << "force    ";
  if (x_ == _VELOCITY) cout << "velocity ";
  cout << xval_ << endl;
  
  cout << "y: ";
  if (y_ == _FORCE)    cout << "force    ";
  if (y_ == _VELOCITY) cout << "velocity ";
  cout << yval_ << endl; 
  
  cout << "rot: ";
  if (rot_ == _FORCE)    cout << "force    ";
  if (rot_ == _VELOCITY) cout << "velocity ";
  cout << rotval_ << endl;  
}



