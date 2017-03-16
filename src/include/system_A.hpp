#ifndef _system_analyze_h
#define _system_analyze_h

#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>
#include "system.hpp"

using namespace std;

//! \brief System analysis  (Virtual class)
//! \author C. Voivret
class System_A
{

	protected:

		System * sys_;

	public:

		//System_analyze(){}
		//System_analyze(System *sys_a): sys_(sys_a) {}
		virtual ~System_A() {}

		virtual void initAnalyse()=0;
		virtual void analyse(double,unsigned int, unsigned int)=0;
		virtual void read_parameters(istream &)=0;
		virtual void allFalse()=0;
		virtual void plugRef() =0;

		System * & sys()       { return sys_; }
		System *   sys() const { return sys_; } 

		static System_A * factory( string type);

		void plugSystem   (System  * sys) { sys_ = sys; }



};

#endif // _system_analyze_h



