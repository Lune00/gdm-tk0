#ifndef _groupdata_hpp
#define _groupdata_hpp

#include<fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

//! \brief Define the groups and their parameters
//! \author V. Richefeu
class GroupData
{
  typedef map<string, unsigned int > mapIdParam;
  
public:
  
  GroupData();
  GroupData(unsigned int ngrp);
  ~GroupData();
  
  //! Return true if the parameter exist
  //! \param name  Parameter name 
  bool exist(string name);
  
  //! Get a parameter identifier from his name
  //! \param name  Parameter name 
  unsigned int getId(string name);
  
  //! Get a parameter value from his name
  //! \param name  Parameter name 
  //! \param g     Group number
  double getParameter(string name, unsigned int g) const;
  
  //! Get quickly a parameter value from an identifier
  //! \param idPar Parameter identifier (see function getId) 
  //! \param g     Group number
  double getParameterQuickly(unsigned int idPar, unsigned int g) const;
  
  //! Set a parameter value
  //! \param name  Name of the parameter 
  //! \param g     Group number
  //! \param value Parameter value
  void setParameter(string name, unsigned int g, double value);
  
  void addParameter(string name);
  void read(istream & is);
  void read(const char *);
  void write(ostream & os);
  
private:
  
  unsigned int ngrp_;                  //< Number of groups
  unsigned int npar_;                  //< Number of parameters
  map<string ,unsigned int > idParam_; //< A map to relate the parameter identifiers to their name
  vector<double *> lpar_;              //< The list of parameters
  
};

#endif // _groupdata_hpp
