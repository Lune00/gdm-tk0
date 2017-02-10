#include "groupRelationData.hpp"

GroupRelationData::GroupRelationData() { npar_ = 0 ; ngrp_ = 0; fragmentation_=false; frag_ = 0 ;}

GroupRelationData::GroupRelationData(unsigned int ngrp) : ngrp_(ngrp)
{ 
  npar_ = 0;
  initActivator();
}

GroupRelationData::~GroupRelationData() 
{
  // Free memory for the action table
  if (act_ != 0)
    {
    for(unsigned int i=0 ; i<ngrp_ ; ++i)
      {
      delete [] act_[i];
      }
    delete [] act_;
    act_ = 0;
    }
  
  // Free memory for each parameter table
  for(unsigned int p=0 ; p<npar_ ; ++p)
    {
    if (lpar_[p] != 0)
      {
      for(unsigned int i=0 ; i<ngrp_ ; ++i)
        {
        delete [] lpar_[p][i];
        }
      delete [] lpar_[p];
      lpar_[p] = 0;
      }
    }
}

bool GroupRelationData::act(unsigned int g1, unsigned int g2) const
{
  if (g1 < ngrp_ && g2 < ngrp_)
    return act_[g1][g2];
  else
    cerr << "@GroupRelationData::act, bad group number g1= "<<g1<<" g2= "<<g2 << endl;
  
  return false;
}

void GroupRelationData::activate(unsigned int g1, unsigned int g2)
{
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    act_[g1][g2] = true;
    act_[g2][g1] = true;
    }
  else
    cerr << "@GroupRelationData::activate, bad group number" << endl;
}

void GroupRelationData::deactivate(unsigned int g1, unsigned int g2)
{
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    act_[g1][g2] = false;
    act_[g2][g1] = false;
    }
  else
    cerr << "@GroupRelationData::deactivate, bad group number" << endl;
}

bool GroupRelationData::exist(string name)
{
  map<string, unsigned int >::const_iterator ip = idParam_.find(name);
  if (ip == idParam_.end()) return false;
  else return true;
}
/*
bool GroupRelationData::existLaw(string name)
{
  map<string, unsigned int >::const_iterator ip = idLaw_.find(name);
  if (ip == idLaw_.end()) return false;
  else return true;
}
*/
unsigned int GroupRelationData::getId(string name)
{
  map<string, unsigned int >::const_iterator ip = idParam_.find(name);
  if (ip != idParam_.end())
    {
    return (unsigned int) ip->second;
    }
  else
    {
    return npar_;
    }    
}
/*
unsigned int GroupRelationData::getIdLaw(string name)
{
  map<string, unsigned int >::const_iterator ip = idLaw_.find(name);
  if (ip != idLaw_.end())
    {
    return (unsigned int) ip->second;
    }
  else
    {
    return npar_;
    }    
}
*/
double GroupRelationData::getParameter(string name, unsigned int g1, unsigned int g2) const
{
  // Retrieve the parameter identifier
  unsigned int idPar;
  map<string, unsigned int >::const_iterator ip = idParam_.find(name);
  if (ip != idParam_.end())
    {
    idPar = ip->second;
    }
  else
    {
    cerr << "@GroupRelationData::getParameter, parameter " << name << " not found" << endl;
    return 0.0;
    }
  
  // Retrieve the value
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    return lpar_[idPar][g1][g2];
    }
  else
    {
    cerr << "@GroupRelationData::getParameter, bad group number" << endl;
    cerr << "see 'ngrp' in section GroupRelationData in your data" << endl;
    }
  
  return 0.0;    
}

cohesionLaw * GroupRelationData::getLaw( unsigned int g1, unsigned int g2) const
{
  // Retrieve the parameter identifier
  /*unsigned int idPar;
  map<string, unsigned int >::const_iterator ip = idLaw_.find(name);
  if (ip != idLaw_.end())
    {
    idPar = ip->second;
    }
  else
    {
    cerr << "@GroupRelationData::getLaw, law " << name << " not found" << endl;
    return NULL;
    }
  
  // Retrieve the value
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    return  llaw_[idPar][g1][g2];
    }
  else
    {
    cerr << "@GroupRelationData::getLaw, bad group number" << endl;
    cerr << "see 'ngrp' in section GroupRelationData in your data" << endl;
    }
  */
  return lco_[ g1 * (ngrp_-1) + g2 ];    
}

double GroupRelationData::getParameterQuickly(unsigned int idPar, unsigned int g1, unsigned int g2) const
{
  // Here we have some confidence in the user !
  return lpar_[idPar][g1][g2];    
}
/*
cohesionLaw * GroupRelationData::getLawQuickly(unsigned int idPar, unsigned int g1, unsigned int g2) const
{
  // Here we have some confidence in the user !
  return  llaw_[idPar][g1][g2];    
}*/

void GroupRelationData::setParameter(string name, unsigned int g1, unsigned int g2, double value)
{
  // Retrieve the parameter identifier
  unsigned int idPar;
  map<string, unsigned int >::const_iterator ip = idParam_.find(name);
/*for(ip=idParam_.begin();ip!=idParam_.end();++ip)
{
	cout<<ip->first<<" p  "<<ip->second<<endl;
}*/
  if (ip != idParam_.end())
    {
    idPar = ip->second;
    }
  else
    {
    cerr << "@GroupRelationData::setParameter, unknown parameter" << name << endl;
    return;
    }
  
  // Affect the value
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    lpar_[idPar][g1][g2] = value;
    lpar_[idPar][g2][g1] = value;
    }
  else
    {
    cerr << "@GroupRelationData::setParameter, bad group number" << endl;
    }
  
  return;    
}

void GroupRelationData::setLaw(string name, unsigned int g1, unsigned int g2,istream & is)
{
  // Retrieve the parameter identifier
	//cout<<" setLaw:::"<<endl;
	lco_[ g1 * (ngrp_-1) + g2 ] = cohesionLaw::factory(name);
	
	lco_[ g2 * (ngrp_-1) + g1 ] -> read(is);
	lco_[ g2 * (ngrp_-1) + g1 ] = lco_[ g1 * (ngrp_-1) + g2 ];
	//cohesionLaw * colaw = 
  //unsigned int idLaw;
 /* map<string, unsigned int >::const_iterator il = idLaw_.find(name);
for(il=idLaw_.begin();il!=idLaw_.end();++il)
{
	cout<<il->first<<" "<<il->second<<endl;
}
  if (il != idLaw_.end())
    {
    idLaw = il->second;
    }
  else
    {
    cerr << "@GroupRelationData::setLaw, unknown law " << name << endl;
    return;
    }
	//cout<<" idLAw = "<<idLaw<<endl;
  // Affect the law
  if (g1 < ngrp_ && g2 < ngrp_)
    {
    (llaw_[idLaw][g1][g2])->read(is);
    llaw_[idLaw][g2][g1] = llaw_[idLaw][g1][g2];
    }
  else
    {
    cerr << "@GroupRelationData::setLaw, bad group number " << endl;
    }
  */
  return;    
}

void GroupRelationData::addParameter(string name)
{
  double ** p = 0;
  
  //idParam_[name] = npar_++;
  idParam_.insert(mapIdParam::value_type(name,npar_++));
  
  p = new double * [ngrp_];
  if (p != 0)
    {
    for(unsigned int i=0;i<ngrp_;i++)
      {
      p[i] = new double [ngrp_];
      }
    }
  else
    cerr << "@GroupRelationData::addParameter, allocation problem" << endl;
  
  for (unsigned int i=0 ; i<ngrp_ ; ++i)
    for (unsigned int j=0 ; j<ngrp_ ; ++j)
      p[i][j] = 0.0;
  
  lpar_.push_back(p);
}
/*
void GroupRelationData::addLaw(string name)
{
  cohesionLaw *** p = 0;
  
  //idParam_[name] = npar_++;
  idLaw_.insert(map<string ,unsigned int >::value_type(name,nlaw_++));
  
  p = new cohesionLaw ** [ngrp_];
  if (p != 0)
    {
    for(unsigned int i=0;i<ngrp_;i++)
      {
      p[i] = new cohesionLaw * [ngrp_];
      }
    }
  else
    cerr << "@GroupRelationData::addLaw, allocation problem" << endl;
  
  cohesionLaw * colaw = cohesionLaw::factory(name);
  for (unsigned int i=0 ; i<ngrp_ ; ++i)
    for (unsigned int j=0 ; j<ngrp_ ; ++j)
      p[i][j] = colaw;
  
  llaw_.push_back(p);
}
*/
void GroupRelationData::initActivator()
{
  act_ = new bool * [ngrp_];
  for(unsigned int i=0;i<ngrp_;i++)
    {
    act_[i] = new bool [ngrp_];
	for(unsigned int j=0;j<ngrp_;++j)
		lco_.push_back(cohesionLaw::factory("nocohesion"));
    }
//	cout<<" taille lco "<<lco_.size()<<endl;
  // By default, all bodies can act on all other bodies 
  for (unsigned int i=0;i<ngrp_;i++)
    for (unsigned int j=0;j<ngrp_;j++)
      act_[i][j] = true;
}

void GroupRelationData::read(istream & is)
{
  string token;
  
  is >> token;	
  while(is)
    {	
    if      (token == "ngrp")      
      {
      is >> ngrp_;
      if (ngrp_ == 0) cerr << "GroupRelationData::read, ngrp can not be 0" << endl;
      else 
		{
			initActivator();
		}
		
      }
    else if (token == "noact")
      { 
      unsigned int i,j; 
      is >> i ; is >> j;
      deactivate(i,j);
      }
    else if (token == "parameter")
      { 
      if (ngrp_ == 0) cerr << "GroupRelationData::read, ngrp can not be 0" << endl;
      string name;
      is >> name;
      addParameter(name);
      }
    else if (token == "setall")
      { 
      string parName;
      double value;
      
      is >> parName >> value;
      
      for (unsigned int g1=0;g1<ngrp_;++g1)
        for (unsigned int g2=0;g2<ngrp_;++g2)
          setParameter(parName,g1,g2,value);
      }
    else if (token == "set")
      { 
      string parName;
      unsigned int g1,g2;
      double value;
      
      is >> parName >> g1 >> g2 >> value;
      
      setParameter(parName,g1,g2,value);
      }
	// cohesion laws
	else if (token == "cohesionLaw")
      { 
	//	cout<<"set law "<<endl;
      string lawName;
      unsigned int g1,g2;
      //double value;
      is >> lawName >> g1 >> g2;
	  //cout<<" nom loi = "<<lawName<<endl;

      setLaw(lawName , g1 , g2, is);
      }
    else if (token=="Fragmentation")
      {
      	fragmentation_=true;
      	if (frag_ == 0) frag_ = new Fragmentation();
      	frag_->read(is);
      }
    else if (token == "}") break;
    else cerr << "@GroupRelationData::read, Unknown token: " << token << endl;
    
    is >> token;
    }
}

void GroupRelationData::write(ostream & os)
{
  os << "ngrp " << ngrp_ << endl;
  
  for (unsigned int g1=0;g1<ngrp_;++g1)
    {
    for (unsigned int g2=0;g2<ngrp_;++g2)  
      {
      if (!act_[g1][g2]) os << "noact " << g1 << " " << g2 << endl;
      }
    }
  
  mapIdParam::iterator imap = idParam_.begin();
  string parName;
  while (imap != idParam_.end())
    {
    parName = imap->first;
    os << "parameter " << parName << endl;
    
    for (unsigned int g1=0;g1<ngrp_;++g1)
      {
      for (unsigned int g2=0;g2<ngrp_;++g2)  
        {
        if (lpar_[imap->second][g1][g2]) 
          {
          os << "set " << parName << " " 
          << g1 << " " << g2 << " " 
          << lpar_[imap->second][g1][g2] << endl;
          }
        }
      }
    
    ++imap;
    }
  
}


