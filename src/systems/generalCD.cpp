#include "generalCD.hpp"

// todo rename as 'Basic'

void generalCD::read_parameters(istream & is)
{
  string token;
  unsigned int xmod = _VELOCITY, ymod = _VELOCITY, rotmod = _VELOCITY;
  double yvalue = 0.0, xvalue = 0.0, rotvalue = 0.0;
  unsigned int currentControl = 0;
  unsigned int nextControl    = 0;
  
  is >> token;	
  while(is)
    {	
    
    if (token == "Gravity")
      {
      double g_mod;
      is >> g_mod;
      gx() = 0.0;
      gy() = -g_mod;
      }
    else if (token == "..")
      {
      is >> nextControl;
      
      // Read 3 commands
      for (unsigned int i=0 ; i<3 ; ++i)
        {
        is >> token;
        
        if (token == "blocx")
          { xmod = _VELOCITY ; xvalue = 0.0; }
        else if (token == "blocy")
          { ymod = _VELOCITY ; yvalue = 0.0; }
        else if (token == "blocrot")
          { rotmod = _VELOCITY ; rotvalue = 0.0; }
        else if (token == "freex")
          { xmod = _FORCE ; xvalue = 0.0; }
        else if (token == "freey")
          { ymod = _FORCE ; yvalue = 0.0; }
        else if (token == "freerot")
          { rotmod = _FORCE ; rotvalue = 0.0; }
        else if (token == "vx")
          { xmod = _VELOCITY ; is >> xvalue; }
        else if (token == "vy")
          { ymod = _VELOCITY ; is >> yvalue; }
        else if (token == "vrot")
          { rotmod = _VELOCITY ; is >> rotvalue; }
        else if (token == "fx")
          { xmod = _FORCE ; is >> xvalue; }
        else if (token == "fy")
          { ymod = _FORCE ; is >> yvalue; }
        else if (token == "frot")
          { rotmod = _FORCE ; is >> rotvalue; }
        else
          { cerr << "@generalCD::read_parameters, Unknown keyword: " << token << endl; }
        }
      
      for(unsigned int i=currentControl ; i<nextControl ; ++i)
        lctrl_.push_back(control(xmod,ymod,rotmod,xvalue,yvalue,rotvalue));
      
      currentControl = nextControl;
      }
    else if (token == "}") break;
    else cerr << "@generalCD::read_parameters, Unknown keyword: " << token << endl;
    
    is >> token;
    }
}

void generalCD::write_parameters(ostream & os)
{
  os << "generalCD" << endl;
  os << "Gravity " << gy() << endl;
  for (unsigned int i=0 ; i<lctrl_.size() ; ++i)
    {
    os << ".. " << i+1 << endl; 
    if (lctrl_[i].x() == _FORCE) os << "fx ";
    else os << "vx ";
    os << lctrl_[i].xval() << " ";
    if (lctrl_[i].y() == _FORCE) os << "fy ";
    else os << "vy ";
    os << lctrl_[i].yval() << " ";
    if (lctrl_[i].rot() == _FORCE) os << "frot ";
    else os << "vrot ";
    os << lctrl_[i].rotval() << endl;
    }
  
  /*
  if(lctrl_[1].x() == _FORCE && _xpres == 0.0)
    os << "forceX    " << -_xval  << endl;
  if(lctrl_[3].y() == _FORCE && _ypres == 0.0)
    os << "forceY    " << -_yval  << endl;	
  
  if(lctrl_[1].x() == _FORCE && _xpres > 0.0)
    os << "pressureX " << _xpres << endl;
  if(lctrl_[3].y() == _FORCE && _ypres > 0.0)
    os << "pressureY " << _ypres << endl;
  
  if(lctrl_[1].x() == _VELOCITY)
    os << "velocityX " << -_xval  << endl;
  if(lctrl_[3].y() == _VELOCITY)
    os << "velocityY " << -_yval  << endl;
   */
}

void generalCD::init() 
{  

}

void generalCD::drive() 
{
  
}

void generalCD::trans() 
{

}

void generalCD::share() 
{
    
}

int generalCD::check()
{
  return 1;
}


void generalCD::stress_strain()
{}
