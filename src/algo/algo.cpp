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
  gdm::fatal("No Algorithm is defined!");
  return 0;
}

void Algo::algoFill()
{
	// Compute mass of bodies
	double density = 2500.0;
		if(grpDat_->exist("density"))
		{
			cout<<"algoFill, taille systeme charge : "<<spl_->lbody().size()<<endl;
			for (unsigned int i=0;i<spl_->lbody().size();++i)
			{

				density = grpDat_->getParameter("density",spl_->body(i)->grp());
				spl_->body(i)->Fill(density);
			}
		cout<<"on sort de la boucle sur lbody"<<endl;	
		}
		else
		{
			spl_->fill(2500.0);
			gdm::warning("@algo::algoFill, GroupData density not defined (so it is set to 2500)");
		}
		for(unsigned int i=0;i<sys_->ldof().size();++i)
		{
			sys_->ldof(i)->fill(density);
			sys_->ldof(i)->setGravity(sys_->gx(),sys_->gy() );
			sys_->ldof(i)->computeMassCenter();
			sys_->ldof(i)->computeMoment();
		}
		cout<<"Fill ok"<<endl;
}


