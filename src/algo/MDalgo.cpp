/// EN MEGA TRAVAUX !!!!!!! 


// stand
//    step --> update step num., update positions, drive & share functions of the system,
//             update system velocities
//    look --> contact,fres,mass matrix & array of coefficients,
//             iter, velocities update
//    hand --> increment time, update verlet list 
//
//
#include "MDalgo.hpp"

using namespace std;

void MDalgo::read_parameters(istream & is)
{
  string token;
  
  is >> token;	
  while(is)
  {	
    if      (token == "dt")      is >> dt_;
    else if (token == "ns")      is >> ns_;
    else if (token == "nsi")     is >> nsi_;
    else if (token == "nsf")     is >> nsf_;
    else if (token == "nwpost")  is >> nwpost_;
    else if (token == "nver")    is >> nver_;
    else if (token == "dverlet") is >> nwk_->dverlet();
    else if (token == "}")       break;
    else cerr << "@MDalgo::read_parameters, Unknown parameter: " << token << endl;
    
    is >> token;
  }
}

void MDalgo::write_parameters(ostream & os)
{
  os << "dt      " << dt_             << endl;
  os << "ns      " << ns_             << endl;
  os << "nsi     " << nsi_            << endl;
  os << "nsf     " << nsf_            << endl;
  os << "nwpost  " << nwpost_         << endl;
  os << "nver    " << nver_           << endl;
  os << "dverlet " << nwk_->dverlet() << endl;
}

// Take care of contacts and forces
void MDalgo::velocityVerletStep(unsigned int ns)
{
  ns_=ns;
  
  unsigned int N     = spl_->lbody().size();
  unsigned int Ncont = sys_->lctrl().size();
  
  // a rajouter a MDalgo
  double dt2_2 = 0.5 * (dt_ * dt_);
  double dt_2  = 0.5 * dt_;
  
  // .......... Update positions (t + dt) and velocities (t + dt/2)
  
  // Imposed DOF
  for (unsigned int i=0 ; i<Ncont ; ++i)
    {
    switch (sys_->ctrl(i).x()) 
      {
      case _VELOCITY : 
        spl_->body(i)->fx() = 0.0;
        spl_->body(i)->x() += (spl_->body(i)->vx() = sys_->ctrl(i).xval()) * dt_2;
        break;
        
      case _FORCE :
        spl_->body(i)->fx() = sys_->ctrl(i).xval();
        spl_->body(i)->x()  += spl_->body(i)->vx() * dt_ + spl_->body(i)->ax() * dt2_2; 
        spl_->body(i)->vx() += spl_->body(i)->ax() * dt_2;
        break;
      }

    switch (sys_->ctrl(i).y()) 
      {
      case _VELOCITY : 
        spl_->body(i)->fy() = 0.0;        
        spl_->body(i)->y() += (spl_->body(i)->vy() = sys_->ctrl(i).yval()) * dt_2;
        break;
        
      case _FORCE :
        spl_->body(i)->fy() = sys_->ctrl(i).yval();
        spl_->body(i)->y()  += spl_->body(i)->vy() * dt_ + spl_->body(i)->ay() * dt2_2; 
        spl_->body(i)->vy() += spl_->body(i)->ay() * dt_2;
        break;
      }
    
    // TODO rotations
    }	
  
  // Free bodies
  for (unsigned int i=Ncont ; i<N ; ++i)
    {
    spl_->body(i)->x()    += spl_->body(i)->vx()   * dt_ + spl_->body(i)->ax()   * dt2_2;
    spl_->body(i)->y()    += spl_->body(i)->vy()   * dt_ + spl_->body(i)->ay()   * dt2_2;
    spl_->body(i)->rot()  += spl_->body(i)->vrot() * dt_ + spl_->body(i)->arot() * dt2_2;
    
    spl_->body(i)->vx()   += spl_->body(i)->ax()   * dt_2;
    spl_->body(i)->vy()   += spl_->body(i)->ay()   * dt_2;
    spl_->body(i)->vrot() += spl_->body(i)->arot() * dt_2;
    }
   
  
  // ........... Compute forces
  forces();
  
  // ........... Update velocities (t + dt)
  
  // Imposed DOF
  double invmassi;
  
  for (unsigned int i=0 ; i<Ncont ; ++i)	
  {
    invmassi = 1.0 / spl_->body(i)->mass();
    
    if(sys_->ctrl(i).x() == _FORCE)
      {
      spl_->body(i)->ax() = spl_->body(i)->fx() * invmassi;
      spl_->body(i)->vx() += dt_2 * spl_->body(i)->ax();
      }

    if(sys_->ctrl(i).y() == _FORCE)
      {
      spl_->body(i)->ay() = spl_->body(i)->fy() * invmassi;
      spl_->body(i)->vy() += dt_2 * spl_->body(i)->ay();
      }
    
    // if(sys_->ctrl(i).rot() == _FORCE) TODO

  }
  
  // Free bodies
  for (unsigned int i=Ncont ; i<N ; ++i)
  {
    invmassi = 1.0 / spl_->body(i)->mass();
    
    spl_->body(i)->ax()   = spl_->body(i)->fx()   * invmassi;
    spl_->body(i)->ay()   = spl_->body(i)->fy()   * invmassi;
    spl_->body(i)->arot() = spl_->body(i)->frot() / spl_->body(i)->mom();    
        
    spl_->body(i)->vx()   += dt_2 * spl_->body(i)->ax();
    spl_->body(i)->vy()   += dt_2 * spl_->body(i)->ay();
    spl_->body(i)->vrot() += dt_2 * spl_->body(i)->arot();

  }	
}

void MDalgo::forces()
{
  if (nwk_->linter().empty()) return;
  
  for (unsigned int k=0 ; k < nwk_->linter().size() ; ++k)
    {		
    nwk_->inter(k)->Kin();
    law_->computeForces(k);
    //law_->fn(k);
    //law_->ft(k);
    //law_->frot(k);
    }  
}


// Increment time, Update Verlet list (if necessary)
void MDalgo::hand()
{
  tm_ += dt_;
  
  if (ns_%nver_ == 0)
  {
    nwk_->verlet(spl_, grpRel_);
  }
}


// Initialize the computation
void MDalgo::stand()
{
  // check existence of the required parameters
  law_->checkRequiredParameters();
  
  // Compute mass of bodies
  if(grpDat_->exist("density"))
    {
    double density = 2500.0;
    for (unsigned int i=0 ; i<spl_->lbody().size() ; ++i)
      {
      density = grpDat_->getParameter("density",spl_->body(i)->grp());
      spl_->body(i)->Fill(density);
      }
    }
  else
    {
    spl_->fill(2500.0);
    dtk::warning("@MDalgo::stand, GroupData density not defined (so it is set to 2500)");
    }
  
  // Set up the initial contact list and the force list
  nwk_->verlet(spl_, grpRel_);
  forces();
}
