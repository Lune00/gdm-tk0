#ifndef _CDinformer_h
#define _CDinformer_h

#include "informer.hpp"
#include "CDalgo.hpp"

using namespace std;

//! \brief Display information relative the Contact Dynamics on screen
//! \author V. Richefeu
class CDinformer : public Informer
{

 protected:
  
  CDalgo * algo_;
  
	
 public:

  virtual ~CDinformer() { }
  CDinformer() : Informer() { algo_ = 0; }
  CDinformer(CDalgo* algo) : Informer(1) { algo_ = algo; }
  CDinformer(CDalgo* algo, unsigned int n) : Informer(n) { algo_ = algo; }
  
  void read(istream & is)
  {
    string token;
    while(is)
      {
      is >> token;
      if (token == "nbIterations") inquiries_.push_back(token);
      else if (token == "Step") inquiries_.push_back(token);
      else if (token == "}") break;
      else cerr << "@CDinformer: unknown inquiry (" << token << ")" << endl;
      }
  }
  
  void write(ostream&) { }
  
  void talk()
  {
    unsigned int N = inquiries_.size();
    if (N == 0) return;
    
    cout << "---[ CD algo ]" << endl;
    for (unsigned int i=0;i<inquiries_.size();++i)
      {
      if (inquiries_[i] == "nbIterations") 
        cout << "nbIterations = " << algo_->niter() << endl;
      
      if (inquiries_[i] == "Step")
        cout << "Step = " << algo_->ns() << " / " << algo_->nsf() << endl;
      }
    cout << endl << flush;
  }

};

#endif // _CDinformer_h



