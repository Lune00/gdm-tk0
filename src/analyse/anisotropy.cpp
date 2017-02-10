#include "anisotropy.hpp"

gdm::Tensor2x2 * Fabric(Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.linter().size();
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny;

	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(c);
		if (interc->fn() != 0.0)
		{
			nx = interc->nx();
			ny = interc->ny();

			F->xx() += nx*nx;
			F->xy() += nx*ny;
			F->yy() += ny*ny;
		}
	}

	F->xx() /= (double)Nc;
	F->xy() /= (double)Nc;
	F->yy() /= (double)Nc;

	F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * FabricInProbe(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn,nx,ny,ang;
	unsigned int Nf=0;
	
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;

				nx = cos(ang);
				ny = sin(ang);
			//if( nx*nx + ny*ny > 1.) cout<<"alert fab "<<nx*nx + ny*ny<<endl;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}
	//cout<<" number of contact "<<Nf<<endl;

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * FabricInProbe_ss(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn,nx,ny,ang;
	unsigned int Nf=0;
	bool simple;
	double fn1,fn2;
	
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			fn1=interc->fn();
			interc->current()=1;
			fn2=interc->fn();
			interc->current()=0;
			if (fn1*fn2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if(fn!=0.0)
		{
			Nf++;
			if (!simple)
			{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;

				nx = cos(ang);
				ny = sin(ang);
			//if( nx*nx + ny*ny > 1.) cout<<"alert fab "<<nx*nx + ny*ny<<endl;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
			}
		}
		}
	}
	//cout<<" number of contact "<<Nf<<endl;

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * FabricInProbe_sv(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn,nx,ny,ang;
	unsigned int Nf=0;
	bool simple;
	double fn1,fn2;
		
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			fn1=interc->fn();
			interc->current()=1;
			fn2=interc->fn();
			interc->current()=0;
			if (fn1*fn2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			Nf++;
			if(simple)
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;

				nx = cos(ang);
				ny = sin(ang);
				//if( nx*nx + ny*ny > 1.) cout<<"alert fab "<<nx*nx + ny*ny<<endl;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
			}
		}
	}
	//cout<<" number of contact "<<Nf<<endl;

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

//Calcule selon direction de la ligne qui lie des centres des particules
gdm::Tensor2x2 * Fabric2InProbe(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn2,fn,ft,nx,ny,l,nx2,ny2;
	unsigned int Nf=0;
	
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}
	//cout<<" number of contact "<<Nf<<endl;

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * Fabric2InProbe_ss(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn2,fn,ft,nx,ny,l,nx2,ny2;
	unsigned int Nf=0, Nss=00;
	double f1,f2;
	bool simple;
	
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		
		if (fn2 != 0.0 && (prb.containCenter(interc->first()) || prb.containCenter(interc->second())))
		{
			Nf++;
			if(!simple)
			{
				Nss++;
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
			}
		}
		
	}

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * Fabric2InProbe_sv(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double fn2,fn,ft,nx,ny,l,nx2,ny2;
	unsigned int Nf=0, Nsv=0;
	double f1,f2;
	bool simple;
	
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		
		if (fn2 != 0.0 && (prb.containCenter(interc->first()) || prb.containCenter(interc->second())))
		{
			Nf++;
			if(simple)
			{
				Nsv++;
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
			}
		}
		
	}

	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;

	F->yx() = F->xy();

	return F;
}

//-------------------------------------------------------------------------
gdm::Tensor2x2 * SimpleContactFabricInProbe(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if (interc->fn() != 0.0  && interc->rang()==1)
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;
				nx = cos(ang);
				ny = sin(ang);

				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}

//	cout<<" number of simple "<<Nf<<endl;

/*	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;*/

		F->yx() = F->xy();

	return F;
}

gdm::Tensor2x2 * DoubleContactFabricInProbe(Probe& prb, Sample& spl, Network& nwk)
{
	unsigned int Nc = nwk.clist().size();
	if (Nc == 0) return 0;

	gdm::Tensor2x2 * F = new gdm::Tensor2x2;
	double nx,ny,ang;
	unsigned int Nf=0;
	inter2d * interc;

	for (unsigned int c = 0; c < Nc; ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if (interc->fn() != 0.0 && interc->rang()==2)
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				ang= acos( interc->nx() ) ;
				if(interc->ny()< 0.) ang=M_PI-ang;
				nx = cos(ang);
				ny = sin(ang);


				F->xx() += nx*nx;
				F->xy() += nx*ny;
				F->yy() += ny*ny;
				Nf++;
			}
		}
	}

//	cout<<" number of double "<<Nf<<endl;

/*	F->xx() /= (double) Nf;
	F->xy() /= (double) Nf;
	F->yy() /= (double) Nf;*/

		F->yx() = F->xy();

	return F;
}



gdm::Tensor2x2 * fnAnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn");
	double nx,ny,fn,angN,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				fmoy+= fn;

				F->xx() += fn*nx*nx;
				F->xy() += fn*nx*ny;
				F->yx() += fn*nx*ny;
				F->yy() += fn*ny*ny;
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * fnAnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn");
	double nx,ny,fn,angN,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	double f1,f2;
	bool simple;

	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				fmoy+= fn;
				if(!simple)
				{
					F->xx() += fn*nx*nx;
					F->xy() += fn*nx*ny;
					F->yx() += fn*nx*ny;
					F->yy() += fn*ny*ny;
				}
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * fnAnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn");
	double nx,ny,fn,angN,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	double f1,f2;
	bool simple;

	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				fmoy+= fn;
				if(simple)
				{
					F->xx() += fn*nx*nx;
					F->xy() += fn*nx*ny;
					F->yx() += fn*nx*ny;
					F->yy() += fn*ny*ny;
				}
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}


gdm::Tensor2x2 * fn2AnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn2");
	double nx,ny,fn,ft,fn2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;

				F->xx() += fn2*nx2*nx2;
				F->xy() += fn2*nx2*ny2;
				F->yx() += fn2*ny2*nx2;
				F->yy() += fn2*ny2*ny2;
				Nf++;
			}
		}
	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	return F;
}

gdm::Tensor2x2 * fn2AnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn2");
	double nx,ny,fn,ft,fn2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;
				Nf++;
				if(!simple)
				{
					F->xx() += fn2*nx2*nx2;
					F->xy() += fn2*nx2*ny2;
					F->yx() += fn2*ny2*nx2;
					F->yy() += fn2*ny2*ny2;
				}
			}
		}

	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	return F;
}

gdm::Tensor2x2 * fn2AnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn2");
	double nx,ny,fn,ft,fn2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		ft=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
			ft+=interc->ft();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;
		nx = interc->nx();
		ny = interc->ny();
		fn2=fn*(nx*nx2+ny*ny2)+ft*(-ny*nx2+nx*ny2);
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;
				Nf++;
				if(simple)
				{
					F->xx() += fn2*nx2*nx2;
					F->xy() += fn2*nx2*ny2;
					F->yx() += fn2*ny2*nx2;
					F->yy() += fn2*ny2*ny2;
				}
			}
		}

	}

	if(Nf != 0 ) 
	{
		fmoy/=(double) (Nf);
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	return F;
}


gdm::Tensor2x2 * ftAnisoInProbe(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft");

	inter2d * interc;
	double nx,ny,ft,fn,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;

	double angT;
	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		if (ft != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angT= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angT=M_PI-angT;
				
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				fmoy+= fn;

				F->xx() -= ft*ny*nx;
				F->yx() -= ft*ny*ny;
				F->xy() += ft*nx*nx;
				F->yy() += ft*ny*nx;

				Nf++;
			}		
		}
	}
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * ftAnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft");

	inter2d * interc;
	double nx,ny,ft,fn,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;

	double angT;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (ft != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn;
				Nf++;
				if(!simple)
				{
					angT= acos( interc->nx() ) ;
					if(interc->ny()< 0.) angT=M_PI-angT;				
					nx = interc->nx();
					ny = interc->ny();//<0. ? -interc->ny():interc->ny();

					F->xx() -= ft*ny*nx;
					F->yx() -= ft*ny*ny;
					F->xy() += ft*nx*nx;
					F->yy() += ft*ny*nx;
				}
			}		
		}
	}
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * ftAnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft");

	inter2d * interc;
	double nx,ny,ft,fn,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;

	double angT;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		
		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (ft != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn;
				Nf++;
				if(simple)
				{
					angT= acos( interc->nx() ) ;
					if(interc->ny()< 0.) angT=M_PI-angT;				
					nx = interc->nx();
					ny = interc->ny();//<0. ? -interc->ny():interc->ny();

					F->xx() -= ft*ny*nx;
					F->yx() -= ft*ny*ny;
					F->xy() += ft*nx*nx;
					F->yy() += ft*ny*nx;
				}
			}		
		}
	}
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * ft2AnisoInProbe(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft2");

	inter2d * interc;
	double nx,ny,fn,ft,fn2,ft2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;

		nx = interc->nx();
		ny = interc->ny();
		fn2=nx2*(fn*nx-ft*ny)+ny2*(fn*ny+ft*nx);
		ft2=-ny2*(fn*nx-ft*ny)+nx2*(fn*ny+ft*nx);
		
		if (ft2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;
				F->xx() -= ft2*ny2*nx2;
				F->yx() += ft2*nx2*nx2;
				F->xy() -= ft2*ny2*ny2;
				F->yy() += ft2*ny2*nx2;
				Nf++;
			}		
		}
	}
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * ft2AnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft2");

	inter2d * interc;
	double nx,ny,fn,ft,fn2,ft2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;

		nx = interc->nx();
		ny = interc->ny();
		fn2=nx2*(fn*nx-ft*ny)+ny2*(fn*ny+ft*nx);
		ft2=-ny2*(fn*nx-ft*ny)+nx2*(fn*ny+ft*nx);
		
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;
				Nf++;
				if(!simple)
				{
					F->xx() -= ft2*ny2*nx2;
					F->yx() += ft2*nx2*nx2;
					F->xy() -= ft2*ny2*ny2;
					F->yy() += ft2*ny2*nx2;
				}
			}		
		}
	}
	
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * ft2AnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft2");

	inter2d * interc;
	double nx,ny,fn,ft,fn2,ft2,nx2,ny2,l,fmoy=0;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		ft=0.;
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			ft+=interc->ft();
			fn+=interc->fn();
		}
		interc->current()=0;
		
		l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
		nx2=interc->Vbranchx()/l;
		ny2=interc->Vbranchy()/l;

		nx = interc->nx();
		ny = interc->ny();
		fn2=nx2*(fn*nx-ft*ny)+ny2*(fn*ny+ft*nx);
		ft2=-ny2*(fn*nx-ft*ny)+nx2*(fn*ny+ft*nx);
		
		if (fn2 != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				fmoy+= fn2;
				Nf++;
				if(simple)
				{
					F->xx() -= ft2*ny2*nx2;
					F->yx() += ft2*nx2*nx2;
					F->xy() -= ft2*ny2*ny2;
					F->yy() += ft2*ny2*nx2;
				}
			}		
		}
	}
	
	fmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./fmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}

gdm::Tensor2x2 * lengthAnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("length");
	double fn,nx,ny,l,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				lmoy+= l;

				F->xx() += l*nx*nx;
				F->xy() += l*nx*ny;
				F->yx() += l*nx*ny;
				F->yy() += l*ny*ny;
				Nf++;
			}

		}

	}
	lmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * length2AnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("length 2");
	double fn,nx,ny,l,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;
				lmoy+= l;

				F->xx() += l*nx*nx;
				F->xy() += l*nx*ny;
				F->yx() += l*nx*ny;
				F->yy() += l*ny*ny;
				Nf++;
			}
		}
	}
	lmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * length2AnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("length 2");
	double fn,nx,ny,l,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;
				lmoy+= l;
				Nf++;
				if(!simple)
				{
					F->xx() += l*nx*nx;
					F->xy() += l*nx*ny;
					F->yx() += l*nx*ny;
					F->yy() += l*ny*ny;
				}
			}
		}

	}
	lmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * length2AnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("length 2");
	double fn,nx,ny,l,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				l=sqrt(pow(interc->Vbranchx(),2)+pow(interc->Vbranchy(),2));
				nx=interc->Vbranchx()/l;
				ny=interc->Vbranchy()/l;
				lmoy+= l;
				Nf++;
				if(simple)
				{
					F->xx() += l*nx*nx;
					F->xy() += l*nx*ny;
					F->yx() += l*nx*ny;
					F->yy() += l*ny*ny;
				}
			}
		}

	}
	lmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * lnAnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("ln");
	double fn,nx,ny,ln,lnmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				ln  = interc->Vbranchx()*nx+interc->Vbranchy()*ny;

				lnmoy+= ln;

				F->xx() += ln*nx*nx;
				F->xy() += ln*nx*ny;
				F->yx() += ln*nx*ny;
				F->yy() += ln*ny*ny;
				Nf++;
			}

		}

	}
	lnmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lnmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * lnAnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("ln");
	double fn,nx,ny,ln,lnmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				ln  = interc->Vbranchx()*nx+interc->Vbranchy()*ny;

				lnmoy+= ln;
				if(!simple)
				{
					F->xx() += ln*nx*nx;
					F->xy() += ln*nx*ny;
					F->yx() += ln*nx*ny;
					F->yy() += ln*ny*ny;
				}
				Nf++;
			}

		}
		
	}
	lnmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lnmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * lnAnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("ln");
	double fn,nx,ny,ln,lnmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				ln  = interc->Vbranchx()*nx+interc->Vbranchy()*ny;

				lnmoy+= ln;
				if(simple)
				{
					F->xx() += ln*nx*nx;
					F->xy() += ln*nx*ny;
					F->yx() += ln*nx*ny;
					F->yy() += ln*ny*ny;
				}
				Nf++;
			}

		}
		
	}
	lnmoy/=(double) (Nf);

	if(Nf != 0 ) 
	{
		F->scalarMult(1./lnmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * ltAnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("lt");
	double fn,tx,ty,nx,ny,lt,ln,lmoy=0.;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				tx = interc->tx();
				ty = interc->ty();//<0. ? -interc->ny():interc->ny();
				ln=interc->Vbranchx()*nx+interc->Vbranchy()*ny;
				lt=interc->Vbranchx()*tx+interc->Vbranchy()*ty;

				lmoy+= ln;

				F->xx() += lt*nx*tx;
				F->xy() += lt*nx*ty;
				F->yx() += lt*ny*tx;
				F->yy() += lt*ny*ty;
				Nf++;
			}
		}

	}
	lmoy/=(double) (Nf);
	
	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	
	return F;

}

gdm::Tensor2x2 * ltAnisoInProbe_ss(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("lt");
	double fn,tx,ty,nx,ny,lt,ln,lmoy=0.;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				tx = interc->tx();
				ty = interc->ty();//<0. ? -interc->ny():interc->ny();
				ln=interc->Vbranchx()*nx+interc->Vbranchy()*ny;
				lt=interc->Vbranchx()*tx+interc->Vbranchy()*ty;

				lmoy+= ln;
				if(!simple)
				{
					F->xx() += lt*nx*tx;
					F->xy() += lt*nx*ty;
					F->yx() += lt*ny*tx;
					F->yy() += lt*ny*ty;
				}
				Nf++;
			}
		}
	}
	lmoy/=(double) (Nf);
	
	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	
	return F;

}

gdm::Tensor2x2 * ltAnisoInProbe_sv(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("lt");
	double fn,tx,ty,nx,ny,lt,ln,lmoy=0.;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;
	double f1,f2;
	bool simple;

	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if(interc->rang()==2) 
		{
			interc->current()=0;
			f1=interc->fn();
			interc->current()=1;
			f2=interc->fn();
			interc->current()=0;
			if (f1*f2!=0.) simple=false;
			else simple=true;
		}
		else simple=true;
		

		fn=0.;
		for(unsigned int r=0;r<interc->rang();r++)
		{
			interc->current()=r;
			fn+=interc->fn();
		}
		interc->current()=0;
		
		if (fn != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();

				tx = interc->tx();
				ty = interc->ty();//<0. ? -interc->ny():interc->ny();
				ln=interc->Vbranchx()*nx+interc->Vbranchy()*ny;
				lt=interc->Vbranchx()*tx+interc->Vbranchy()*ty;

				lmoy+= ln;
				if(simple)
				{
					F->xx() += lt*nx*tx;
					F->xy() += lt*nx*ty;
					F->yx() += lt*ny*tx;
					F->yy() += lt*ny*ty;
				}
				Nf++;
			}
		}
	}
	lmoy/=(double) (Nf);
	
	if(Nf != 0 ) 
	{
		F->scalarMult(1./lmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;
	
	return F;

}

gdm::Tensor2x2 * fnlAnisoInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fnl");
	double nx,ny,fn,angN,flmoy=0,l;
	unsigned int Nf=0;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;

				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				fn = interc->fn();

				if ( interc->rang()==2)
				{
					interc->current()=1;
					fn+=interc->fn();
					interc->current()=0;
				}

				flmoy+= l*fn;

				F->xx() += l*fn*nx*nx;
				F->xy() += l*fn*nx*ny;
				F->yx() += l*fn*nx*ny;
				F->yy() += l*fn*ny*ny;
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		flmoy/=(double) (Nf);
		F->scalarMult(1./flmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;

}

gdm::Tensor2x2 * ftlAnisoInProbe(Probe& prb, Sample& spl, Network& nwk, unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ftl");

	inter2d * interc;
	double nx,ny,ft,fn,flmoy=0,l;
	unsigned int Nf=0;
	if (nwk.clist().empty()) return NULL;

	double angT;
	for (unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->ft() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angT= acos( interc->nx() ) ;
			//if (angT<0.) cout<<angT<<endl;
				if(interc->ny()< 0.) angT=M_PI-angT;


				nx = cos(angT);//interc->nx();
				ny = sin(angT);//interc->ny()<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				
				ft = interc->ft();
				fn = interc->fn();

				if ( interc->rang()==2)
				{
					interc->current()=1;
					ft+=interc->ft();
					fn+=interc->fn();
					interc->current()=0;
				}
				flmoy+= l*fn;

				F->xx() -= l*ft*ny*nx;
				F->yx() += l*ft*nx*nx;
				F->xy() -= l*ft*ny*ny;
				F->yy() += l*ft*ny*nx;

				Nf++;
			}		
		}
	}
	flmoy/=(double) (Nf);
	if(Nf != 0 ) 
	{
		F->scalarMult(1./flmoy);
		F->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	return F;
}


pointSet  PthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	inter2d * interc;
	DataSet p;
	double angN;
	for( unsigned int i=0;i<nwk.clist().size();++i)
	{
		interc = nwk.inter(nwk.clist(i));

		if( interc->fn() !=0. && prb.contain(interc))
		{
			angN= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angN = M_PI - angN;

			p.add( angN);

		}
	}
	return (p.kernelPdf(nc,1.1/(double)(nc)));
}

pointSet  SimpleContactPthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	inter2d * interc;
	DataSet p;
	double angN;
	for( unsigned int i=0;i<nwk.clist().size();++i)
	{
		interc = nwk.inter(nwk.clist(i));

		if( interc->fn() !=0. && prb.contain(interc) && interc->rang()==1)
		{
			angN= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angN = M_PI - angN;

			p.add( angN);

		}


	}
	return (p.kernelPdf(nc,1.1/(double)(nc)));
}

pointSet  DoubleContactPthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	inter2d * interc;
	DataSet p;
	double angN;
	for( unsigned int i=0;i<nwk.clist().size();++i)
	{
		interc = nwk.inter(nwk.clist(i));

		if( interc->fn() !=0. && prb.contain(interc) && interc->rang()==2)
		{
			angN= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angN = M_PI - angN;

			p.add( angN);

		}


	}
	return (p.kernelPdf(nc,1.1/(double)(nc)));
}


gdm::Tensor2x2 * fnAnisoInProbe2(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Fn");
	double nx,ny;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	vector <unsigned int > NcN(Nbin,0);
	vector <double> FN(Nbin,0.);
	double angN;
	double amp = M_PI/(double) (Nbin);
	unsigned int NctotN=0,rangN;

	double fmoy=0.;
	unsigned int c;

	for ( c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if ( interc->fn() != 0.0 && prb.contain(interc))
		{
			angN= acos( interc->nx() ) ;
			//if(interc->ny()< 0.) angN = M_PI - angN;

			rangN = (unsigned int) ( floor( angN/amp) );
			if(rangN == Nbin) rangN--;

			//force moyenne fonction de theta
			FN[ rangN ] += interc->fn();

			if ( interc->rang()==2)
			{
				interc->current()=1;
				FN[ rangN ]+=interc->fn();
				fmoy += interc->fn();
				interc->current()=0;
			}

			fmoy += interc->fn();			
			NcN[ rangN ]+=1;
			NctotN++;	
		}

	}
	fmoy/=(double) (NctotN);

	ofstream F_out("fnT.txt",ios::out);
	cout<<" fn moy aniso "<<fmoy<<endl;

	for( unsigned int i=0;i<Nbin;++i)
	{
		FN[i] /= ( (double) (NcN[i]) );
		F_out<<setw(10)<<(.5+i)*amp<<" "<<NcN[i]/(double) NctotN*Nbin/M_PI<<" "<<FN[i]<<" "<<fmoy<<endl;
	}

	F_out.close();

	/*for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 && prb.contain(interc))
		{
			angN= acos( interc->nx() ) ;
				//if(interc->ny()< 0.) angN = M_PI - angN;

			rangN = (unsigned int) ( floor( angN/amp) );
			if(rangN == Nbin) rangN--;

			nx = interc->nx();
			ny = interc->ny();

				//cout<<FN[rangN]<<" "<<nx<<" "<<ny<<endl;

			F->xx() += FN[ rangN ]*nx*nx;
			F->xy() += FN[ rangN ]*nx*ny;
			F->yx() += FN[ rangN ]*nx*ny;
			F->yy() += FN[ rangN ]*ny*ny;
		}

	}
	*/
	for( unsigned int i=0;i<Nbin;++i)
	{
		nx = cos ( (i+.5) *amp);
		ny = sin ( (i+.5) *amp);
		F->xx() += FN[ i ]*nx*nx;
		F->xy() += FN[ i ]*nx*ny;
		F->yx() += FN[ i ]*nx*ny;
		F->yy() += FN[ i ]*ny*ny;

	}
	return F;

}

gdm::Tensor2x2 * ftAnisoInProbe2(Probe& prb, Sample& spl, Network& nwk,unsigned int Nbin)
{
	gdm::Tensor2x2 * F = new gdm::Tensor2x2("Ft");
	double nx,ny;
	inter2d * interc;
	if ( nwk.clist().empty()) return NULL;

	vector <unsigned int > NcT(Nbin,0);
	vector <double> FT(Nbin,0.);
	double angT;
	double amp = M_PI/(double) (Nbin);
	unsigned int NctotT=0,rangN;

	double fmoy=0.;
	unsigned int c;

	for ( c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));
		if ( interc->ft() != 0.0 && prb.contain(interc))
		{
			angT= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angT = M_PI - angT;

			rangN = (unsigned int) ( floor( angT/amp) );
			if(rangN == Nbin) rangN--;

			//force moyenne fonction de theta
			FT[ rangN ] += interc->ft();

			if ( interc->rang()==2)
			{
				interc->current()=1;
				FT[ rangN ]+=interc->ft();
				fmoy += interc->ft();
				interc->current()=0;
			}

			fmoy += interc->ft();			
			NcT[ rangN ]+=1;
			NctotT++;	
		}

	}
	fmoy/=(double) (NctotT);

	ofstream F_out("ftT.txt",ios::out);


	for( unsigned int i=0;i<Nbin;++i)
	{

		FT[i] /= ( (double) (NcT[i]) );
		//FN[i]=1.;
		F_out<<setw(10)<<(.5+i)*amp<<" "<<NcT[i]/(double) NctotT*Nbin/M_PI<<" "<<FT[i]<<" "<<fmoy<<endl;

	}

	F_out.close();

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->ft() != 0.0 && prb.contain(interc))
		{
			angT= acos( interc->nx() ) ;
			if(interc->ny()< 0.) angT = M_PI - angT;

			rangN = (unsigned int) ( floor( angT/amp) );
			if(rangN == Nbin) rangN--;

			nx = interc->nx();
			ny = interc->ny();

			F->xx() -= FT[ rangN ]*ny*nx;
			F->xy() += FT[ rangN ]*nx*nx;
			F->yx() -= FT[ rangN ]*ny*ny;
			F->yy() += FT[ rangN ]*ny*ny;


		}

	}


	return F;

}

gdm::Tensor2x2 * fnlthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nt,unsigned int nl)
{
	//if ( nwk.clist().empty()) return NULL;

	vector <pointSet> fnlt;
	vector <pointSet> Nfnlt;

	double ampt,ampl;
	ampt=M_PI/(double)(nt);
	pointSet temp;
	
	for(unsigned int i=0;i<nl;++i)
	{
		fnlt.push_back(temp);
		Nfnlt.push_back(temp);
	}

	for(unsigned int j=0;j<nl;++j)
	{
	for(unsigned int i=0;i<nt;++i)
	{
		fnlt[j].add((.5+i)*ampt,0);
		Nfnlt[j].add((.5+i)*ampt,0);
	//	cout<<(.5+i)*ampt/M_PI*180<<endl;

	}
	}
	//cout<<"alloc ok"<<endl;
	double lmax=0,lmin=1000000;
	double nx,ny,l,angN,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				if (l<lmin) lmin=l;
				if (l>lmax) lmax=l;

				lmoy+= l;
				Nf++;
			}

		}

	}
	lmoy/=(double) (Nf);

	ampl=(lmax-lmin)/(double)(nl);
	unsigned int rangt,rangl;
//	cout<<"ampl "<<ampl<<endl;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );
				rangl=(unsigned int) floor( (l-lmin)/ampl );

				if(rangt==nt) rangt--;
				if(rangl==nl) rangl--;
				if(rangl<0) rangl=0;

				fnlt[rangl].py(rangt)+=interc->fn();
				Nfnlt[rangl].py(rangt)+=1;
			}

		}

	}
	for(unsigned int j=0;j<nl;++j)
	{
		for(unsigned int i=0;i<nt;++i)
		{
			if(Nfnlt[j].py(i)>0 )
			fnlt[j].py(i)/=Nfnlt[j].py(i);
			else
				fnlt[j].py(i)=0;
		}
	}
	
	ofstream out1("Analyse/angDistrib/fnlt.txt",ios::out);
	ofstream out2("Analyse/angDistrib/test.txt",ios::out);
	
	for(unsigned int i=0;i<nl;++i)
	{
		for(unsigned int j=0;j<nt;++j)
		{
			//cout<<fnlt[i].px(j)/M_PI*180<<endl;
			
			out1<<(.5+i)*ampl*cos(fnlt[i].px(j))<<" "<<(.5+i)*ampl*sin(fnlt[i].px(j))<<" "<<fnlt[i].py(j)<<endl;
		}
		out1<<endl<<endl;
	}
	
	pointSet fnltmoy;
	double moy=0,fnmoy=0;
		for(unsigned int j=0;j<nt;++j)
		{
			//cout<<fnlt[i].px(j)/M_PI*180<<endl;
			moy=0.;
			for(unsigned int i=0;i<nl;++i)
			{
				if(fnlt[i].py(j)>0 )
				moy+= fnlt[i].py(j) * (.5+i)*ampl;
			}
			moy/=(double)(nl);
			fnltmoy.add(fnlt[0].px(j),moy);
			fnmoy+=moy;
		//	out2<<fnlt[0].px(j)<<" "<<moy<<endl;
		}
		fnmoy/=(double)(nt);
		
		for(unsigned int j=0;j<nt;++j)
		{
	
			out2<<fnlt[0].px(j)<<" "<<fnltmoy.py(j)/fnmoy<<endl;
		}
	
	out1.close();
	out2.close();
	
	double fn;
	cout<<"fnmoy = "<<fnmoy<<endl;
	gdm::Tensor2x2 * Fn = new gdm::Tensor2x2("Fn");
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;
				rangt=(unsigned int) floor( angN/ampt );
				if(rangt==nt) rangt--;
				

				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				fn = interc->fn();

				if ( interc->rang()==2)
				{
					interc->current()=1;
					fn+=interc->fn();
					interc->current()=0;
				}


				Fn->xx() += fnltmoy.py(rangt)*nx*nx;
				Fn->xy() += fnltmoy.py(rangt)*nx*ny;
				Fn->yx() += fnltmoy.py(rangt)*nx*ny;
				Fn->yy() += fnltmoy.py(rangt)*ny*ny;
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		//fmoy/=(double) (Nf);
		Fn->scalarMult(1./fnmoy);
		Fn->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	
	return Fn;

}

pointSet  CorrelationthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	// C=<fl>/(<f><l>)-1
	inter2d * interc;
	DataSet temp;
	double angN;
	vector <DataSet>  fnlt,fnt,lt;
	double nx,ny,l,ampt=M_PI/(double)(nc);
	unsigned int rangt;
	
	for(unsigned int i=0;i<nc;++i)
	{
		//Ct.push_back(temp);
		fnlt.push_back(temp);	
		fnt.push_back(temp);
		lt.push_back(temp);
	}

	
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );

				if(rangt==nc) rangt--;

				fnlt[rangt].add(interc->fn()*l);
				fnt[rangt].add(interc->fn());
				lt[rangt].add(l);
			}

		}

	}
	//cout<<"boucle ok "<<endl;
	pointSet Ct;
	for(unsigned int i=0;i<nc;++i)
	{
		fnlt[i].extractValues();
		fnt[i].extractValues();
		lt[i].extractValues();
	//	cout<<i<<" "<<fnlt[i].setSize()<<endl;
		Ct.add((.5+i)*ampt, fnlt[i].mean()/(fnt[i].mean()*lt[i].mean())-1.);
	}
	//Ct.extractValues();
	//Ct.yNormalise(Ct.ymean());
	return (Ct);
	
}

gdm::Tensor2x2 * CorrelationAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	// C=<fl>/(<f><l>)-1
	inter2d * interc;
	DataSet temp;
	double angN;
	vector <DataSet>  fnlt,fnt,lt;
	double nx,ny,l,ampt=M_PI/(double)(nc);
	unsigned int rangt;
	
	for(unsigned int i=0;i<nc;++i)
	{
		//Ct.push_back(temp);
		fnlt.push_back(temp);	
		fnt.push_back(temp);
		lt.push_back(temp);
	}

	//distribution angulaire du coefficient de correlation
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );

				if(rangt==nc) rangt--;

				fnlt[rangt].add(interc->fn()*l);
				fnt[rangt].add(interc->fn());
				lt[rangt].add(l);
			}

		}

	}
	//cout<<"boucle ok "<<endl;
	pointSet Ct;
	for(unsigned int i=0;i<nc;++i)
	{
		fnlt[i].extractValues();
		fnt[i].extractValues();
		lt[i].extractValues();
	//	cout<<i<<" "<<fnlt[i].setSize()<<endl;
		Ct.add((.5+i)*ampt, fnlt[i].mean()/(fnt[i].mean()*lt[i].mean())-1.);
	}
	
	gdm::Tensor2x2 * Cn = new gdm::Tensor2x2("Cn");
	if ( nwk.clist().empty()) return NULL;
	unsigned int Nf=0;
	double Cmoy=0;
	//Calcul anisotropie
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;
				rangt=(unsigned int) floor( angN/ampt );
				if(rangt==nc) rangt--;
				

				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				

				Cn->xx() += Ct.py(rangt)*nx*nx;
				Cn->xy() += Ct.py(rangt)*nx*ny;
				Cn->yx() += Ct.py(rangt)*nx*ny;
				Cn->yy() += Ct.py(rangt)*ny*ny;
				Nf++;
				Cmoy+=Ct.py(rangt);
			}

		}

	}

	if(Nf != 0 ) 
	{
		Cmoy/=(double) (Nf);
		Cn->scalarMult(1./Cmoy);
		Cn->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	
	
	
	
	return (Cn);
}
gdm::Tensor2x2 * fnlcAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	//if ( nwk.clist().empty()) return NULL;

	inter2d * interc;
	DataSet temp;
	double angN;
	vector <DataSet>  fnlt,fnt,lt;
	double nx,ny,l,ampt=M_PI/(double)(nc);
	double fn;
	unsigned int rangt,Nf=0;
	double ctfnlmoy=0,prod=0;
	
	for(unsigned int i=0;i<nc;++i)
	{
		//Ct.push_back(temp);
		fnlt.push_back(temp);	
		fnt.push_back(temp);
		lt.push_back(temp);
	}

	
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );

				if(rangt==nc) rangt--;

				fnlt[rangt].add(interc->fn()*l);
				fnt[rangt].add(interc->fn());
				lt[rangt].add(l);
			}

		}

	}
	//cout<<"boucle ok "<<endl;
	pointSet Ct;
	for(unsigned int i=0;i<nc;++i)
	{
		fnlt[i].extractValues();
		fnt[i].extractValues();
		lt[i].extractValues();
	//	cout<<i<<" "<<fnlt[i].setSize()<<endl;
		Ct.add((.5+i)*ampt, fnlt[i].mean()/(fnt[i].mean()*lt[i].mean())-1.);
	}
	

	gdm::Tensor2x2 * Fn = new gdm::Tensor2x2("Fn");
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;
				rangt=(unsigned int) floor( angN/ampt );
				if(rangt==nc) rangt--;
				

				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				fn = interc->fn();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				
				if ( interc->rang()==2)
				{
					interc->current()=1;
					fn+=interc->fn();
					interc->current()=0;
				}

				prod = Ct.py(rangt)*fn*l;
				
				Fn->xx() += prod*nx*nx;
				Fn->xy() += prod*nx*ny;
				Fn->yx() += prod*nx*ny;
				Fn->yy() += prod*ny*ny;
				ctfnlmoy+=prod;
				Nf++;
			}

		}

	}
	ctfnlmoy/=(double)(Nf);
	//Fn->print();
	if(Nf != 0 ) 
	{
		//fmoy/=(double) (Nf);
		//Fn->scalarMult(1./ctfnlmoy);
		//Fn->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	
	return Fn;

}

gdm::Tensor2x2 * ftlcAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc)
{
	//if ( nwk.clist().empty()) return NULL;

	inter2d * interc;
	DataSet temp;
	double angN;
	vector <DataSet>  fnlt,fnt,lt;
	double nx,ny,l,ampt=M_PI/(double)(nc);
	double ft;
	unsigned int rangt,Nf=0;
	double ctftlmoy=0,prod=0;
	
	for(unsigned int i=0;i<nc;++i)
	{
		//Ct.push_back(temp);
		fnlt.push_back(temp);	
		fnt.push_back(temp);
		lt.push_back(temp);
	}

	
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );

				if(rangt==nc) rangt--;

				fnlt[rangt].add(interc->ft()*l);
				fnt[rangt].add(interc->ft());
				lt[rangt].add(l);
			}

		}

	}
	//cout<<"boucle ok "<<endl;
	pointSet Ct;
	for(unsigned int i=0;i<nc;++i)
	{
		fnlt[i].extractValues();
		fnt[i].extractValues();
		lt[i].extractValues();
	//	cout<<i<<" "<<fnlt[i].setSize()<<endl;
		Ct.add((.5+i)*ampt, fnlt[i].mean()/(fnt[i].mean()*lt[i].mean())-1.);
	}
	

	gdm::Tensor2x2 * Ft = new gdm::Tensor2x2("Ft");
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->fn() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;
				rangt=(unsigned int) floor( angN/ampt );
				if(rangt==nc) rangt--;
				

				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				ft = interc->ft();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				
				if ( interc->rang()==2)
				{
					interc->current()=1;
					ft+=interc->ft();
					interc->current()=0;
				}

				prod = Ct.py(rangt)*ft*l;
				
				Ft->xx() += prod*nx*nx;
				Ft->xy() += prod*nx*ny;
				Ft->yx() += prod*nx*ny;
				Ft->yy() += prod*ny*ny;
				ctftlmoy+=prod;
				Nf++;
			}

		}

	}
	ctftlmoy/=(double)(Nf);
	//Fn->print();
	if(Nf != 0 ) 
	{
		//fmoy/=(double) (Nf);
		//Fn->scalarMult(1./ctfnlmoy);
		//Fn->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	
	return Ft;

}

gdm::Tensor2x2 * ftlthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nt,unsigned int nl)
{
	//if ( nwk.clist().empty()) return NULL;

	vector <pointSet> ftlt;
	vector <pointSet> Nftlt;

	double ampt,ampl;
	ampt=M_PI/(double)(nt);
	pointSet temp;
	
	for(unsigned int i=0;i<nl;++i)
	{
		ftlt.push_back(temp);
		Nftlt.push_back(temp);
	}

	for(unsigned int j=0;j<nl;++j)
	{
	for(unsigned int i=0;i<nt;++i)
	{
		ftlt[j].add((.5+i)*ampt,0);
		Nftlt[j].add((.5+i)*ampt,0);
	//	cout<<(.5+i)*ampt/M_PI*180<<endl;

	}
	}
	//cout<<"alloc ok"<<endl;
	double lmax=0,lmin=1000000;
	double nx,ny,l,angN,lmoy=0;
	unsigned int Nf=0;
	inter2d * interc;
	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->ft() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());
				if (l<lmin) lmin=l;
				if (l>lmax) lmax=l;

				lmoy+= l;
				Nf++;
			}

		}

	}
	lmoy/=(double) (Nf);

	ampl=(lmax-lmin)/(double)(nl);
	unsigned int rangl,rangt;
//	cout<<"ampl "<<ampl<<endl;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->ft() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;


				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				l  = sqrt(interc->Vbranchx()*interc->Vbranchx()+interc->Vbranchy()*interc->Vbranchy());

				rangt=(unsigned int) floor( angN/ampt );
				rangl=(unsigned int) floor( (l-lmin)/ampl );

				if(rangt==nt) rangt--;
				if(rangl==nl) rangl--;
				if(rangl<0) rangl=0;

				ftlt[rangl].py(rangt)+=fabs(interc->ft());
				Nftlt[rangl].py(rangt)+=1;
			}

		}

	}
	for(unsigned int j=0;j<nl;++j)
	{
		for(unsigned int i=0;i<nt;++i)
		{
			if(Nftlt[j].py(i)>0 )
			ftlt[j].py(i)/=Nftlt[j].py(i);
			else
				ftlt[j].py(i)=0;
		}
	}
	
	ofstream out1("Analyse/angDistrib/ftlt.txt",ios::out);
	ofstream out2("Analyse/angDistrib/testT.txt",ios::out);
	
	for(unsigned int i=0;i<nl;++i)
	{
		for(unsigned int j=0;j<nt;++j)
		{
			//cout<<fnlt[i].px(j)/M_PI*180<<endl;
			
			out1<<(.5+i)*ampl*cos(ftlt[i].px(j))<<" "<<(.5+i)*ampl*sin(ftlt[i].px(j))<<" "<<ftlt[i].py(j)<<endl;
		}
		out1<<endl<<endl;
	}
	
	pointSet ftltmoy;
	double moy=0,fnmoy=0;
		for(unsigned int j=0;j<nt;++j)
		{
			//cout<<fnlt[i].px(j)/M_PI*180<<endl;
			moy=0.;
			for(unsigned int i=0;i<nl;++i)
			{
				if(ftlt[i].py(j)>0 )
				moy+= ftlt[i].py(j) * (.5+i)*ampl;
			}
			moy/=(double)(nl);
			ftltmoy.add(ftlt[0].px(j),moy);
			fnmoy+=moy;
		//	out2<<fnlt[0].px(j)<<" "<<moy<<endl;
		}
		fnmoy/=(double)(nt);
		
		for(unsigned int j=0;j<nt;++j)
		{
	
			out2<<ftlt[0].px(j)<<" "<<ftltmoy.py(j)/fnmoy<<endl;
		}
	
	out1.close();
	out2.close();
	
	double fn;
	cout<<"ftmoy = "<<fnmoy<<endl;
	gdm::Tensor2x2 * Ft = new gdm::Tensor2x2("Fn");
	if ( nwk.clist().empty()) return NULL;

	for ( unsigned int c = 0; c <  nwk.clist().size(); ++c)
	{
		interc = nwk.inter(nwk.clist(c));

		if (interc->ft() != 0.0 )
		{
			if(prb.containCenter(interc->first()) || prb.containCenter(interc->second()) )
			{
			//	cout<<"contact"<<endl;
				angN= acos( interc->nx() ) ;
				if(interc->ny()< 0.) angN = M_PI - angN;
				rangt=(unsigned int) floor( angN/ampt );
				if(rangt==nt) rangt--;
				

				nx = interc->nx();
				ny = interc->ny();//<0. ? -interc->ny():interc->ny();
				fn = interc->fn();

				if ( interc->rang()==2)
				{
					interc->current()=1;
					fn+=interc->fn();
					interc->current()=0;
				}


				Ft->xx() -= ftltmoy.py(rangt)*nx*nx;
				Ft->xy() += ftltmoy.py(rangt)*nx*ny;
				Ft->yx() -= ftltmoy.py(rangt)*nx*ny;
				Ft->yy() += ftltmoy.py(rangt)*ny*ny;
								
				Nf++;
			}

		}

	}

	if(Nf != 0 ) 
	{
		//fmoy/=(double) (Nf);
		Ft->scalarMult(1./fnmoy);
		Ft->scalarMult(1./(double) (Nf));
	}
	else return NULL;

	
	return Ft;

}


