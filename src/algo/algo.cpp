#include "algo.hpp"
#include "CDalgo.hpp"
#include "talk.hpp"
//#include "MDalgo.hpp"
//#include "LBalgo.hpp"

Algo* Algo::factory(string type)
{
  if(type == "CD")   return new CDalgo;
  //if(type == "MD")   return new MDalgo;
  //if(type == "LB")   return new LBalgo;
	 
  
  cerr << "@Algo::factory, Unknown type of algorithm: " << type << endl << flush;
  gdm::fatal("No Algorothm is defined!");
  return 0;
}

void Algo::algoFill()
{
	// Compute mass of bodies
	double density = 2500.0;
cout<<"algo.cpp fill density"<<endl;	
cout<<"grpDat_"<<grpDat_->exist("density")<<endl;
		if(grpDat_->exist("density"))
		{
			cout<<"spl_->lbody().size()="<<spl_->lbody().size()<<endl;
			for (unsigned int i=0;i<spl_->lbody().size();++i)
			{
				density = grpDat_->getParameter("density",spl_->body(i)->grp());
				spl_->body(i)->Fill(density);
			}
			
		}
		else
		{
			cout<<"Else fill"<<endl;
			spl_->fill(2500.0);
			gdm::warning("@algo::algoFill, GroupData density not defined (so it is set to 2500)");
		}
cout<<"algo.cpp compute mass,gravity..."<<endl;
		for(unsigned int i=0;i<sys_->ldof().size();++i)
		{
			sys_->ldof(i)->fill(density);
			sys_->ldof(i)->setGravity(sys_->gx(),sys_->gy() );
			//sys_->ldof(i)->computeVertex();
			sys_->ldof(i)->computeMassCenter();
			sys_->ldof(i)->computeMoment();
			//sys_->ldof(i)->print();
		}
}


