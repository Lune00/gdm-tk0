#ifndef _groupRelationData_hpp
#define _groupRelationData_hpp

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "cohesionLaw.hpp"
#include "Fragmentation.hpp"

using namespace std;

//! \brief Define the visibility and the parameters between groups
//! \author V. Richefeu
class GroupRelationData
{
  typedef map<string ,unsigned int > mapIdParam;
  //typedef map<string ,unsigned int > mapIdLaw;
  
public:
  
  GroupRelationData();
  GroupRelationData(unsigned int ngrp); 
  ~GroupRelationData();
  
  bool act(unsigned int g1, unsigned int g2) const;
  
  void activate(unsigned int g1, unsigned int g2);
  
  void deactivate(unsigned int g1, unsigned int g2);
  
  //! Return true if the parameter exist
  //! \param name  Parameter name 
  bool exist(string name);
 
  //! Get a parameter identifier from his name
  //! \param name  Parameter name
  unsigned int getId(string name);
  
  //! Get a parameter value from his name
  //! \param name  Parameter name 
  //! \param g1    First group number
  //! \param g2    Second group number
  double getParameter(string name, unsigned int g1, unsigned int g2) const;

  //! Get quickly a parameter value from an identifier
  //! \param idPar Parameter identifier (see function getId) 
  //! \param g1    First group number
  //! \param g2    Second group number
  double getParameterQuickly(unsigned int idPar, unsigned int g1, unsigned int g2) const;
  
  //! Set a parameter value
  //! \param name  Name of the parameter 
  //! \param g1    First group number
  //! \param g2    Second group number
  //! \param value Parameter value
  void setParameter(string name, unsigned int g1, unsigned int g2, double value);

 void addParameter(string name);


 /* //! Return true if the law exist
  //! \param name  Law name 
  bool existLaw(string name);

  //! Get a law identifier from his name
  //! \param name  Law name
  unsigned int getIdLaw(string name);

  //! Get a parameter value from his name
  //! \param name  Parameter name 
  //! \param g1    First group number
  //! \param g2    Second group number
  cohesionLaw * getLaw(string name, unsigned int g1, unsigned int g2) const;

  //! Get quickly a parameter value from an identifier
  //! \param idPar Parameter identifier (see function getId) 
  //! \param g1    First group number
  //! \param g2    Second group number
  cohesionLaw * getLawQuickly(unsigned int idPar, unsigned int g1, unsigned int g2) const;
*/
  //! Set a parameter value
  //! \param name  Name of the parameter 
  //! \param g1    First group number
  //! \param g2    Second group number
  //! \param value Parameter value
  void setLaw(string name, unsigned int g1, unsigned int g2, istream & is);

//! Get a parameter value from his name
  //! \param name  Parameter name 
  //! \param g1    First group number
  //! \param g2    Second group number
  cohesionLaw * getLaw( unsigned int g1, unsigned int g2) const;
  
  bool existfragmentation() const {return fragmentation_;}
  bool & existfragmentation() {return fragmentation_;}

  Fragmentation* fragmentation() const {return frag_;}
  Fragmentation* & fragmentation() {return frag_;}
//  void addLaw(string name);
  
  

  
  void initActivator();
  
  void read(istream & is);
  
  void write(ostream & os);
  
private:

  unsigned int ngrp_;                  //
  unsigned int npar_;                  //
  bool ** act_;                        // 
  bool fragmentation_;
  map<string ,unsigned int > idParam_; //
  vector<double **> lpar_;             //

  //map<string ,unsigned int > idLaw_; //
  //vector<cohesionLaw *** >llaw_;
  vector <cohesionLaw* > lco_;
  //unsigned int nlaw_;
  Fragmentation * frag_;
};

#endif // _GroupRelationData_hpp
