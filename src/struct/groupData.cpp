#include "groupData.hpp"

GroupData::GroupData() { npar_ = 0 ; ngrp_ = 0; }

GroupData::GroupData(unsigned int ngrp) : ngrp_(ngrp)
{ 
  npar_ = 0;
}

GroupData::~GroupData() 
{
  // Free memory for each parameter
  for(unsigned int p=0 ; p<npar_ ; ++p)
  {
    if(lpar_[p] != 0) 
    {
      delete [] lpar_[p]; 
      lpar_[p] = 0;
    }
  }
}

bool GroupData::exist(string name)
{
  mapIdParam::const_iterator ip = idParam_.find(name);
  if (ip == idParam_.end()) return false;
  else return true;
}

unsigned int GroupData::getId(string name)
{
  mapIdParam::const_iterator ip = idParam_.find(name);
  if (ip != idParam_.end())
  {
    return (unsigned int) ip->second;
  }
  else
  {
    return ngrp_;
  }    
}

double GroupData::getParameter(string name, unsigned int g) const
{
  // Retrieve the group identifier
  unsigned int idPar;
  mapIdParam::const_iterator ip = idParam_.find(name);
  if (ip != idParam_.end())
  {
    idPar = ip->second;
  }
  else
  {
    cerr << "@GroupData::getParameter, parameter " << name << " not found" << endl;
    return 0.0;
  }

  // Retrieve the value
  if (g < ngrp_)
  {
    return lpar_[idPar][g];
  }
  else
  {
    cerr << "@GroupData::getParameter, bad group number" << endl;
    cerr << "see 'ngrp' in section GroupData" << endl;
  }

  return 0.0;    
}


double GroupData::getParameterQuickly(unsigned int idPar, unsigned int g) const
{
  // Here we have some confidence in the user!
  return lpar_[idPar][g];    
}

void GroupData::setParameter(string name, unsigned int g, double value)
{
  // Retrieve the parameter identifier
  unsigned int idPar;
  mapIdParam::const_iterator ip = idParam_.find(name);
  if (ip != idParam_.end())
  {
    idPar = ip->second;
  }
  else
  {
    cerr << "@GroupData::setParameter, unknown parameter" << name << endl;
    return;
  }

  // Affect the value
  if (g < ngrp_)
  {
    lpar_[idPar][g] = value;
  }
  else
  {
    cerr << "@GroupData::setParameter, bad group number" << endl;
  }

  return;    
}

void GroupData::addParameter(string name)
{
  double * p = 0;

  idParam_[name] = npar_++;

  p = new double [ngrp_];

  for (unsigned int i=0 ; i<ngrp_ ; ++i)
    p[i] = 0.0;

  lpar_.push_back(p);
}

void GroupData::read(istream & is)
{
  string token;

  is >> token;	
  while(is)
  {	
    if      (token == "ngrp")      
    {
      is >> ngrp_;
      if (ngrp_ == 0) cerr << "GroupData::read, ngrp can not be 0" << endl;
    }
    else if (token == "parameter")
    { 
      if (ngrp_ == 0) cerr << "GroupData::read, ngrp can not be 0" << endl;
      string name;
      is >> name;
      addParameter(name);
    }
    else if (token == "setall")
    { 
      string parName;
      double value;

      is >> parName >> value;

      for (unsigned int g=0;g<ngrp_;++g)
      {
	setParameter(parName,g,value);
      }
    }
    else if (token == "set")
    { 
      string parName;
      unsigned int g;
      double value;

      is >> parName >> g >> value;

      setParameter(parName,g,value);
    }
    else if (token == "}") break;
    else cerr << "@GroupData::read, Unknown token: " << token << endl;

    is >> token;
  }
}
//Lecture a partir de group.ini s'il existe
void GroupData::read(const char * fname)
{
  cerr<<"On est dans GroupData::read(const char*)"<<endl;
  ifstream datafile(fname);
  if(!datafile)
  {
    cerr<<"@GroupData::read, cannot open file "<<fname<<endl;
  }
  string token;
  datafile >> token;
  while(datafile)
  {
    if(token=="GroupData{")
    {
      while(token!="}")
      {
	datafile >> token ;

	if(token == "ngrp")      
	{
	  datafile >> ngrp_;
	  if (ngrp_ == 0) cerr << "GroupData::read, ngrp can not be 0" << endl;
	}
	else if (token == "parameter")
	{ 
	  if (ngrp_ == 0) cerr << "GroupData::read, ngrp can not be 0" << endl;
	  string name;
	  datafile >> name;
	  addParameter(name);
	  cerr<<"set "<<name<<endl;
	}
	else if (token == "setall")
	{ 
	  string parName;
	  double value;
	  datafile >> parName >> value;
	  cerr<<"setall : "<<parName<<" "<<ngrp_<<endl;
	  for (unsigned int g=0;g<ngrp_;++g)
	  {
	    setParameter(parName,g,value);
	  }
	}
	else if (token == "set")
	{ 
	  string parName;
	  unsigned int g;
	  double value;
	  datafile >> parName >> g >> value;
	  setParameter(parName,g,value);
	}
      }
      cerr<<"GroupData defined."<<endl;
    }
      datafile >> token;
  }
}

void GroupData::write(ostream & os)
{
  os << "ngrp " << ngrp_ << endl;

  mapIdParam::iterator imap = idParam_.begin();
  string parName;
  while (imap != idParam_.end())
  {
    parName = imap->first;
    os << "parameter " << parName << endl;

    for (unsigned int g=0;g<ngrp_;++g)
    {
      if (lpar_[imap->second][g]) 
      {
	os << "set " << parName << " " 
	  << g << " " 
	  << lpar_[imap->second][g] << endl;
      }
    }

    ++imap;
  }
}
