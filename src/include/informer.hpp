#ifndef _informer_h
#define _informer_h

#include <iostream>
#include <string>
#include <vector>
#include "control.hpp"
#include "sample.hpp"
#include "network.hpp"

using namespace std;

//! \brief Display information on screen (Virtual class)
//! \author V. Richefeu
// TODO: Transform this class for processing inquiries...
class Informer
{

 protected:

  unsigned int countNumber_;
  vector<string> inquiries_;
	
 public:
  
  virtual ~Informer() { }
  Informer() { countNumber_ = 1;}
  Informer(unsigned int n) { countNumber_ = n; }
    
  virtual void read(istream&) = 0;
  virtual void write(ostream&) = 0;
  virtual void talk() = 0;

  bool mustTalk(unsigned int step) { return (step % countNumber_ == 0); }
  void setFrequency(unsigned int n) { countNumber_ = n; }
  void addInquiry(string inquiry) { inquiries_.push_back(inquiry); }
  
};

#endif // _informer_h



