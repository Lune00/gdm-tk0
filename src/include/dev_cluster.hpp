#ifndef _dev_cluster_h
#define _dev_cluster_h

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>
#include "io.hpp"
#include "system.hpp"



class dev_cluster : public System
{

 
	
 vector<unsigned int> topPlate_;
 double bottomPlateThickness_,topPlateThickness_;
 string Unit_;
 bool boundariesAuto_;
 bool pressY_;
 bool firstUse_;
 double dverlet_;
 double dsuperList_;
 double dsuperListP_;
 double bandwidth_;
 double kverlet_;
 unsigned int topXmode_ , topYmode_ , topRmode_;
 double		  topXvalue_, topYvalue_, topRvalue_;
double grav,angle_,rotation_;
  unsigned int Nanalyze_;//????????????


 protected:
 	
 public:
  
  void read_parameters (istream&);
  void write_parameters(ostream&);
  void init ();
  void drive();
 
  void trans();
  void share();
  int  check();
  
  void stress_strain();

  
  //void analyze(double);
  //void write_analyze(FILE *,double);------dans l'analyseur

  ~dev_cluster() { }
   dev_cluster(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel) { }

   dev_cluster() : System() { }
   
   vector<unsigned int> & ltopPlate()       { return topPlate_; }
   vector<unsigned int>   ltopPlate() const { return topPlate_; }
   
   unsigned int & topPlate(unsigned int i)       { return topPlate_[i]; }
   unsigned int   topPlate(unsigned int i) const { return topPlate_[i]; }
	
   void defineTopPlate();
   void defineTopPlate2();
   
   //unsigned int & Nanalyze()          { return Nanalyze_; }
   //unsigned int   Nanalyze()    const { return Nanalyze_; }
  
   string & Unit()       { return Unit_;}
   string   Unit() const { return Unit_;}
   bool & firstUse()       { return firstUse_;}
   bool   firstUse() const { return firstUse_;}
   bool & pressY()       { return pressY_;}
   bool   pressY() const { return pressY_;}
   bool & boundariesAuto()       { return boundariesAuto_;}
   bool   boundariesAuto() const { return boundariesAuto_;}
   unsigned int & topXmode()          { return topXmode_; }
   unsigned int   topXmode()    const { return topXmode_; }
   unsigned int & topYmode()          { return topYmode_; }
   unsigned int   topYmode()    const { return topYmode_; }
   double & topXvalue()          { return topXvalue_; }
   double   topXvalue()    const { return topXvalue_; }
   double & topYvalue()          { return topYvalue_; }
   double   topYvalue()    const { return topYvalue_; }

   
};

#endif // _dev_cluster_h
