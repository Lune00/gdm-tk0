#ifndef _system_h
#define _system_h

#include <iostream>
#include <string>
#include <vector>
#include "control.hpp"
#include "sample.hpp"
#include "network.hpp"
#include "dof.hpp"
#include "groupRelationData.hpp"
//#include "probe.hpp"

using namespace std;

//! \brief System DOF (Virtual class)
//! \author V. Richefeu
class System
{
	protected:

		Sample  * spl_;
		Network * nwk_;
		GroupRelationData * grpRel_;
		vector<control>  lctrl_;
		vector <dof*> ldof_;
		double gx_,gy_;
	public:

		virtual void read_parameters (istream&) = 0;
		virtual void write_parameters(ostream&) = 0;
		virtual void init () = 0;
		virtual void drive() = 0;
		virtual void trans() = 0;
		virtual void share() = 0;
		virtual int  check() = 0;
		virtual void stress_strain() = 0; 

		virtual ~System() {cerr<<"System Destructor called"<<endl; }
		System(Sample* spl, Network* nwk, GroupRelationData * grpRel) : spl_(spl), nwk_(nwk), grpRel_(grpRel) { }
		System() { 
			spl_ = NULL;
			nwk_ = NULL;
			grpRel_ = NULL;
			gx_ = 0.0 ;
			gy_ = 0.0; 
		}

		static System* factory(string type);

		void plugSample   (Sample  * spl) { spl_ = spl; }
		void plugNetwork  (Network * nwk) { nwk_ = nwk; }
		void plugGroupRelationData (GroupRelationData * grpRel) { grpRel_ = grpRel; }

		vector<control> & lctrl()       { return lctrl_; }
		vector<control>   lctrl() const { return lctrl_; }

		control & ctrl(unsigned int i)       { return lctrl_[i]; }
		control   ctrl(unsigned int i) const { return lctrl_[i]; }

		vector <dof*> & ldof()       { return ldof_; }
		vector <dof*>   ldof() const { return ldof_; }
		dof* ldof(unsigned int i)  {return ldof_[i];}

		Sample  * spl() const {return spl_;}
		Network * nwk() const {return nwk_;}
		GroupRelationData * grpRel() const {return grpRel_;}

		void computeAccField() { }

		double & gx()       { return gx_; }
		double   gx() const { return gx_; }
		double & gy()       { return gy_; }
		double   gy() const { return gy_; }

};

#endif // _system_h
