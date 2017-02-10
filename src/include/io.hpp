#ifndef _io_hpp
#define _io_hpp

//! \file io.hpp
//! Some routines for data Input/Output.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "algo.hpp"
#include "sample.hpp"
#include "network.hpp"
#include "system.hpp"
#include "simulation.hpp"


 
using namespace std;

//! \brief Check whether a file exists
//! \author V. Richefeu
bool fileExists(const std::string& fileName);

//! \brief Write a mgpost file (see http://www.lmgc.univ-montp2.fr/~richefeu) 
//! \author V. Richefeu
void write_mgpost(const char * fname, Sample&, Network&, int = 0, double = 0.0);

//! \brief Write postscript file 
//! \author V. Richefeu
void write_ps(const char * fname, Sample&, Network&, int = 0, double = 0.0);

//! \brief Read data from file 
//! \author V. Richefeu
void read_data(const char* fname, Algo&, Sample&, Network&, System&, 
               GroupData&, GroupRelationData&);

//! \brief Save in file the current data concerning algorithm, 
//!        sample, network, system, and referrer in the native format
//! \author V. Richefeu 
void save_data(const char* fname, Algo&, Sample&, Network&, System&, 
               GroupData&, GroupRelationData&, bool saveNetwork = false);

//! \brief Write a history file 
//! \author V. Richefeu
void history_write(unsigned int number, Sample&, Network&, GroupRelationData&, bool TwoFiles = false, bool saveNetwork = false, bool compact = false);

//! \brief Read a history file 
//! \author V. Richefeu
void history_read(const char * name, Sample& spl, Network& nwk);



//! \brief Read a sample 
//! \author V. Richefeu
void sample_read(const char * name, Sample& spl);

//! \brief Write a sample 
//! \author V. Richefeu
void sample_write(const char * name, Sample& spl);

//! \brief Write a network 
//! \author V. Richefeu
void network_write(const char * name, Network& nwk);

//! \brief Convert sample for tapio-K 
//! \author V. Richefeu
void sample_write_tapioK(const char * name, Sample& spl);




#endif //_io_hpp
