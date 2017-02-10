#include "stress.hpp"

gdm::Tensor2x2 * StressInProbe(Probe & prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	cout<<"nc = "<<Nc<<endl;
	if (Nc == 0) 
	{
		cout<<"clist  is empty";
		return NULL;
	}

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress");
	double dx,dy,fx,fy;
	unsigned int nc=0;
	
	inter2d * interc;
	
	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->type() != 2 && interc->type() != 5)//Pkoi les contacts avec les rlines sont exclus
		{	
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				dx = interc->Vbranchx();
				dy = interc->Vbranchy();
				fx=0.;
				fy=0.;
				for(unsigned int r=0;r<interc->rang();r++)
				{
					interc->current()=r;
					fx += interc->fx();
					fy += interc->fy();
				}
				interc->current()=0;
				//if((fabs(fx)>1.0E+10)||(fabs(fy)>1.0E+10)) cout<<"fx:="<<fx<<"	fy:="<<fy<<"	rang:="<<interc->rang()<<"	id1:="<<interc->first()->id()<<"	id2:="<<interc->second()->id()<<endl;
				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		}

//cout<<" S1 et S1 "<<S->l1()<<" "<<S->l2()<<endl;

	}
	cout<<" nc = "<<nc<<endl;
	S->scalarMult( 1./prb.area() );
	return S;
}

gdm::Tensor2x2 * StressInProbe_sv(Probe & prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) 
	{
		cout<<"clist  is empty";
		return NULL;
	}

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress in simple contact");
	double dx,dy,fx,fy,fn1,fn2;
	unsigned int nc=0;
	
	inter2d * interc;
	
	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->type() != 2 && interc->type() != 5)//Pkoi les contacts avec les rlines sont exclus
		{	
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				dx = interc->Vbranchx();
				dy = interc->Vbranchy();
				fx=0.;
				fy=0.;
				if(interc->rang()==1)
				{
					fx += interc->fx();
					fy += interc->fy();
				}
				else if (interc->rang()==2) 
				{
					interc->current()=0;
					fn1=interc->fn();
					interc->current()=1;
					fn2=interc->fn();
					interc->current()=0;
					if(fn1*fn2==0.)
					{
						for(unsigned int r=0;r<interc->rang();r++)
						{
							interc->current()=r;
							fx += interc->fx();
							fy += interc->fy();
						}
						interc->current()=0;
					}
				}
							
				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		}
	}
	S->scalarMult( 1./prb.area() );
	return S;
}

gdm::Tensor2x2 * StressInProbe_ss(Probe & prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) 
	{
		cout<<"clist  is empty";
		return NULL;
	}

	gdm::Tensor2x2 * S = new gdm::Tensor2x2("stress in simple double");
	double dx,dy,fx,fy,fn1,fn2;
	unsigned int nc=0;
	
	inter2d * interc;
	
	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if (interc->type() != 2 && interc->type() != 5 && (interc->rang()==2))//Pkoi les contacts avec les rlines sont exclus
		{	
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				dx = interc->Vbranchx();
				dy = interc->Vbranchy();
				fx=0.;
				fy=0.;
				
				interc->current()=0;
				fn1=interc->fn();
				interc->current()=1;
				fn2=interc->fn();
				interc->current()=0;
				if(fn1*fn2!=0.)
				{
					for(unsigned int r=0;r<interc->rang();r++)
					{
						interc->current()=r;
						fx += interc->fx();
						fy += interc->fy();
					}
					interc->current()=0;
				}					
				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
			}
		}
	}
	S->scalarMult( 1./prb.area() );
	return S;
}

gdm::Tensor2x2 * PartialLengthStressInProbe(Probe & prb, Sample& spl, Network& nwk,unsigned int Npoint,double lmoy,double rmax,double pressure)
{
	unsigned int Nc = nwk.clist().size();
	inter2d * interc;
	
	//cout<<"nc = "<<Nc<<endl;
	if (Nc == 0) 
	{
		cout<<"clist  is empty";
		return NULL;
	}
	//Sorting network by length of branch vectors...
	vector<inter2d*> contact;
	
	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if ( interc->fn() != 0. && interc->type() != 2 && interc->type() != 5)
		{	
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
					contact.push_back(nwk.inter(nwk.clist(c)));
			}
		}
	}

	std::stable_sort(contact.begin(), contact.end(), inter2d::compareLength() );


	gdm::Tensor2x2 * S = new gdm::Tensor2x2("PartialStress");

	
	double dx,dy,fx,fy,s1,s2,ds;
	double qop,q,p;
	unsigned int nc=0;
	unsigned int freq= (unsigned int) floor(contact.size()/Npoint);
	unsigned int clim=contact.size()-1;
	ofstream GS_out("Analyse/PartialLengthstress.txt",ios::out);
	
	for (unsigned int c = 0; c < contact.size(); ++c)
	{
		interc = contact[c];

				dx = interc->Vbranchx();
				dy = interc->Vbranchy();

				fx = interc->fx();
				fy = interc->fy();

				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
				
				if (c % freq  == 0 ||c==clim)
				{
				S->scalarMult( 1./prb.area() );
				S->eigenValues();
				s1=S->l1();
				s2=S->l2();
				p=.5*(s1+s2);
				q= .5*(max(s1,s2)-min(s1,s2));
				qop=q/p;
				ds=S->majorDirection();
				GS_out<<c<<" "<<sqrt(dx*dx+dy*dy)/lmoy<<" "<<sqrt(dx*dx+dy*dy)/rmax<<" "<<p<<" "<<q<<" "<<q/p<<" "<<q/pressure<<" "<<ds<<endl;
				S->scalarMult( prb.area() );
				
				}
	
//cout<<" S1 et S1 "<<S->l1()<<" "<<S->l2()<<endl;

	}
	GS_out.close();
	
	cout<<" nc = "<<nc<<endl;
	return S;
	
}

gdm::Tensor2x2 * PartialNormalForceStressInProbe(Probe & prb, Sample& spl, Network& nwk,unsigned int Npoint,double fnmoy,double pressure)
{
	unsigned int Nc = nwk.clist().size();
	inter2d * interc;
	
	//cout<<"nc = "<<Nc<<endl;
	if (Nc == 0) 
	{
		cout<<"clist  is empty";
		return NULL;
	}
	//Sorting network by length of branch vectors...
	vector<inter2d*> contact;
	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if ( interc->fn() != 0. && interc->type() != 2 && interc->type() != 5)
		{	
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
					contact.push_back(nwk.inter(nwk.clist(c)));
			}
		}
	}

	std::stable_sort(contact.begin(), contact.end(), inter2d::compareNormalForce() );


	gdm::Tensor2x2 * S = new gdm::Tensor2x2("PartialStress");

	
	double dx,dy,fx,fy,s1,s2,ds;
	double qop,q,p;
	unsigned int nc=0;
	unsigned int freq= (unsigned int) floor(contact.size()/Npoint);
	unsigned int clim=contact.size()-1;
	ofstream GS_out("Analyse/PartialNormalForcestress.txt",ios::out);
	
	for (unsigned int c = 0; c < contact.size(); ++c)
	{
		interc = contact[c];

				dx = interc->Vbranchx();
				dy = interc->Vbranchy();

				fx = interc->fx();
				fy = interc->fy();

				S->xx() += fx*dx;
				S->xy() += fx*dy;
				S->yy() += fy*dy;
				S->yx() += fy*dx;
				nc++;
				
				if (c % freq  == 0 ||c==clim )
				{
				S->scalarMult( 1./prb.area() );
				S->eigenValues();
				s1=S->l1();
				s2=S->l2();
				p=.5*(s1+s2);
				q= .5*(max(s1,s2)-min(s1,s2));
				qop=q/p;
				ds=S->majorDirection();
				
				GS_out<<c<<" "<<interc->fn()<<" "<<interc->fn()/fnmoy<<" "<<p<<" "<<q<<" "<<qop<<" "<<q/pressure<<" "<<ds<<endl;
				S->scalarMult( prb.area() );
				}
	
	}
	GS_out.close();
	
	cout<<" nc = "<<nc<<endl;
	return S;
	
}

void IMT_Body(Sample &spl, Network &nwk, vector <gdm::Tensor2x2  > & bodyStress,double rho)
//Internal moment tensor of each body ( divide by the area )
{

	unsigned int tab1,tab2;
	double dx,dy,fx,fy;
	double periode=0.;
	//double angle;
	bool NOdkdkP=true;
	dkdkP * interP;
	inter2d * interc;
	for( unsigned int i=0; i< nwk.clist().size(); ++i)
	{
		interc=nwk.inter(nwk.clist(i));

	/*	for( unsigned int k=0; k< spl.lbody().size() ; ++k)
		{
			if ( spl.body(k)->id() == nwk.inter(i)->first()->id() ) tab1=k;
			if ( spl.body(k)->id() == nwk.inter(i)->second()->id() ) tab2=k;
			} */

				if ( interc->fn() == 0.) continue;

			tab1 = interc->first() -> id();
			tab2 = interc->second()-> id();
	//	cout<<"tab1 = "<<tab1<<" tab2 = "<<tab2<<endl;

			if( NOdkdkP && interc->type() == _type_dkdkP )
			{
				interP=dynamic_cast <dkdkP *>(interc);
				NOdkdkP=false;
				periode= interP->P();
			//cout<<"Periode = "<<periode<<endl;
			}

			if ( interc->type() == _type_dkdk || interc->type() == _type_dkdkP)
			{	

				fx = interc->fn() * interc->nx() + interc->ft() * interc->tx();
				fy = interc->fn() * interc->ny() + interc->ft() * interc->ty();

				dx = interc->x() - interc->first()->x();
				dy = interc->y() - interc->first()->y();

				bodyStress[tab1].xx() += fx*dx;
				bodyStress[tab1].xy() += fx*dy;
				bodyStress[tab1].yy() += fy*dy;
				bodyStress[tab1].yx() += fy*dx;

				if( interc->type() == _type_dkdkP )
					dx = interc->x() - interc->second()->x() - periode ;
				else
					dx = interc->x() - interc->second()->x();

				dy = interc->y() - interc->second()->y();

				bodyStress[tab2].xx() -= fx*dx;
				bodyStress[tab2].xy() -= fx*dy;
				bodyStress[tab2].yy() -= fy*dy;
				bodyStress[tab2].yx() -= fy*dx;


			}
			else cout<<" @StressBody : run only with dkdk and dkdkP "<<endl;

		}

		for( unsigned int k=0; k< spl.lbody().size() ; ++k)
		{

			bodyStress[k].scalarMult(rho/spl.body(k)->Area());
			//cout<<spl.body(k)->Area()<<" "<<bodyStress[k].yy()<<endl;
		}
		//cout<<" TEST    S1 et S1 "<<test.l1()<<" "<<test.l2()<<endl;

		//cout<<"fin bodystress"<<endl;


	}


