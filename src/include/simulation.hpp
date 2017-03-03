#ifndef _simulation_hpp
#define _simulation_hpp

#include <iostream>
#include <fstream>
#include <assert.h>
#include <iomanip> 
#include "sys/types.h"
#include "signal.h"
#include "algo.hpp"
#include "io.hpp"
#include "purge.h"
#include "dof.hpp"
#include "system_A.hpp"


// FONCTION EXTERNE A LA CLASS Simulation = PROBLEME 
// POUR LA SAUVEGARDE
// .. a voir
// si Simulation est un singleton peut-etre que ca s'arrange ...? 

void gdm_signal_handler(int);

//! \brief 
//! \author V. Richefeu
class Simulation
{	

	protected:

		Algo              * algo_;
		Sample            * spl_;
		Network           * nwk_;
		GroupData         * grpDat_;
		GroupRelationData * grpRel_;
		System            * sys_;
		System_A          * sysA_;

		struct sigaction sigact_;

		bool doAnalyse_;

		unsigned int speakLevel_;
		unsigned int nSpeak_;
		bool twoFilesHist_;
		bool historyNetwork_;
		bool compactHist_;
		unsigned int nHist_;
		unsigned int numFileHist_;
		unsigned int nAnalyse_;

		unsigned int  nsi_;     //!< Initial step number
		unsigned int  nsf_;     //!< Final step number
		unsigned int  ns_;      //!< Current step number
		double        time_;    //!< Cumulative time

	public:

		Simulation() 
		{
			algo_   = NULL;
			sys_    = NULL;
			spl_    = NULL;
			nwk_    = NULL;
			grpDat_ = NULL;
			grpRel_ = NULL;
			sysA_   = NULL;

			speakLevel_     = 1;
			nSpeak_         = 1;
			nHist_          = 1;
			nAnalyse_       = 1;
			historyNetwork_ = false;
			numFileHist_    = 0;
			twoFilesHist_   = false;
			doAnalyse_      = false;
			compactHist_    = false;

			ns_  = 0;
			nsi_ = 0;
			nsf_ = 0;
			time_ = 0.0;

			ofstream time("time.txt",ios::out);
			time.close();
		}

		// virtual ~Simulation() { }
		~Simulation() 
		{
			delete algo_; 
			cerr<<"algo deleted."<<endl;
			delete sys_ ;
			cerr<<"sys deleted."<<endl;
			delete nwk_;
			cerr<<"nwk deleted."<<endl;
			delete grpDat_;
			delete grpRel_;
			delete sysA_;
			cerr<<"sysA deleted."<<endl;
		}

		void read_data(const char* fname);
		void save_data(const char* fname);
		void load_history ( const char * fname);

		void init();
		void run();
		void speak();
		//	static void quit() { cout<<"je quite "<<endl;}

		Algo*   algo()     const { return algo_; }
		Algo* & algo()           { return algo_; }

		System*   sys()     const { return sys_; }
		System*  &sys()           { return sys_; }

		Network*  nwk()     const { return nwk_; }
		Network* & nwk()      { return nwk_; }

		Sample*   spl()     const { return spl_; }
		Sample*   &spl()      { return spl_; }

		System_A* sysA()    const { return sysA_; }
		System_A* &sysA()     { return sysA_; }

		unsigned int    nsi()     const { return nsi_; }
		unsigned int  & nsi()           { return nsi_; }
		unsigned int    ns()      const { return ns_;  }
		unsigned int  & ns()            { return ns_;  }
		unsigned int    nsf()     const { return nsf_; }
		unsigned int  & nsf()           { return nsf_; }
		unsigned int    nHist()     const { return nHist_; }
		unsigned int  & nHist()           { return nHist_; }
		unsigned int    nAnalyse()     const { return nAnalyse_; }
		unsigned int  & nAnalyse()           { return nAnalyse_; }
		double          time()    const { return time_; }
		double        & time()          { return time_; }
		bool          doAnalyse()    const { return doAnalyse_; }
		bool        & doAnalyse()          { return doAnalyse_; }
		// ...
};



#endif // _simulation_hpp








