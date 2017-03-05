#ifndef _shearP_CD_h
#define _shearP_CD_h

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>
#include "io.hpp"
#include "system.hpp"


class shearP_CD : public System
{

 vector<unsigned int> topPlate_;
 double bottomPlateThickness_,topPlateThickness_;
 string Unit_;
 bool boundariesAuto_;
 bool pressY_;
 bool firstUse_;
 bool changeGroup_;
 bool shearRate_;
 double shearRate_init;
 double dverlet_;
 bool symetrical_;
 bool reverse_;
 bool useSuper_; 
 bool gravite_;
 unsigned int multig_; // multiple de la gravite g0=-9.81
 double dsuperList_;
 double dsuperListP_;
 double bandwidth_; // normalement associée a sample, mais permet l'entrée utilisateur
 unsigned int topXmode_ , topYmode_ ;
 double		  topXvalue_, topYvalue_;
 unsigned int Nanalyze_;


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
  void SetUnity();
  void perturbation();
  ~shearP_CD() {cerr<<"ShearP_CD destructeur."<<endl; }
    shearP_CD(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel) {  AllFalse(); }
    shearP_CD() : System() { AllFalse(); }
   
   vector<unsigned int> & ltopPlate()       { return topPlate_; }
   vector<unsigned int>   ltopPlate() const { return topPlate_; }
   
   unsigned int & topPlate(unsigned int i)       { return topPlate_[i]; }
   unsigned int   topPlate(unsigned int i) const { return topPlate_[i]; }
	
   void defineTopPlate();
   void defineTopPlate2();
  
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

    
    void AllFalse()
    {
        topPlate_.clear();
        Unit_ = "default";
	topXmode_ = 1;
	topYmode_ = 0;
        symetrical_=false;
        reverse_ =false;
        bandwidth_=0.;
    }
    
   
};

#endif // _shearP_CD_h
