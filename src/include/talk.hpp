#ifndef _gdm_talk_hpp
#define _gdm_talk_hpp

/** \mainpage Source References of Granular Dynamics Methods ToolKit
*  
*  gdm-tk is a toolkit designed to discrete element simulations...
*  
*  To report bugs, suggest enhancements, etc. to the Authors, contact
*  V. Richefeu.
*
*  richefeu@lmgc.univ-montp2.fr
*/

#include <iostream>
#include <fstream>
#include <cstdlib>

namespace gdm
{
  //! \brief Print a message
  void message(const char * txt);
  
  //! \brief Print a informative message and stop the computation
  void fatal(const char * txt);
  
  //! \brief Print a warning message
  void warning(const char * txt);
  
  //! \brief Print a log message in an outgoing stream
  void log(const char *txt, std::ostream & os);
  
  //! \brief A fun function that write 'GDM'
  void fun_GDM(std::ostream & os = std::cout);
  
  // 
  void debugMarker(const char*);
}

#endif //_gdm_talk_hpp
