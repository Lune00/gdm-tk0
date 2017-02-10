#include "solidfraction.hpp"


// For disk only
double solidFraction( heightProbe & prb, Sample& spl,Network & nwk)
{
	//cout<<" Solod fraction in heightProbe "<<endl;
	
	unsigned int Nb = spl.lbody().size();
	if(Nb==0) return 0;
	
	heightProbe prb2( prb.h2()-spl.rmax(), prb.h2()+spl.rmax());
	
	double middle= ( prb.h2() + prb.h1() )*0.5;
	body2d * ghost;
	double Sarea=0,A,h,d,r,sf=0.;
	//cout<<" prb.h1 "<<prb.h1()<<" prb.h2 "<<prb.h2()<<endl;
	for (unsigned int c = 0; c < Nb; ++c)
    {
		
		ghost = spl.body(c)->duplicate();
		ghost->y()= middle + fabs( ghost->y() - middle);
		
		if ( prb.containCenter ( ghost) )
		{
			h= ghost->ymax()-prb.h2();
			if ( prb2.containCenter(ghost) && h>0)
			{
				d=prb.h2()- ghost->y();
				r=ghost->ymax()-ghost->y();
				A=r*r*acos( (r-h)/r)- (r-h)*sqrt(2*r*h-h*h);
				Sarea+= ghost->Area()-A;
				//cout<<" dessous depasse "<< (ghost->Area()-A)/ghost->Area()<<endl;
			}
			else
			{
			Sarea +=ghost->Area();
			}
			
		}
		else if ( prb2.containCenter( ghost) )
		{
			h= prb.h2()-ghost->ymin();
			if ( h>0)
			{
				d=prb.h2()-ghost->y();
				r=ghost->ymax()-ghost->y();
				A=r*r*acos( (r-h)/r)- (r-h)*sqrt(2*r*h-h*h);
				Sarea += A;
				//cout<<" dessus depasse "<< (A)/ghost->Area()<<endl;
			}
			
		}
		
		delete ghost;
	}
	
	double lx,ly;
	unsigned int nc=0 ;
	double r1,r2;
	//check for interpenetration
	for (unsigned int i=0;i< nwk.clist().size();++i)
	{
		nc = nwk.clist(i);
		
		if ( nwk.inter(nc)->rang()>0 && prb.contain( nwk.inter(nc))) 
		{
			lx= nwk.inter(nc)->first()->x() - nwk.inter(nc)->second()->x();
			ly= nwk.inter(nc)->first()->y() - nwk.inter(nc)->second()->y();
			d=sqrt(lx*lx + ly*ly);
			//Run only with disks.....
			
			r1= nwk.inter(nc)->first()->sizeVerlet();
			r2= nwk.inter(nc)->second()->sizeVerlet();
			
			if( d < r1 +r2 )
			{
				//cout<< (d*d + r1*r1 - r2*r2)/(2.*d*r1)<<" "<<( r1*r1 *acos( (d*d + r1*r1 - r2*r2)/(2.*d*r1))
				//	+ r2*r2 * acos ( d*d - r1*r1 + r2*r2)/(2.*d*r2)
				//	-.5*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2)) )<<endl;
				
					
			Sarea -= ( 
				  r1*r1 * acos( (d*d + r1*r1 - r2*r2)/(2.*d*r1))
				+ r2*r2 * acos( (d*d - r1*r1 + r2*r2)/(2.*d*r2))
				-.5*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
						);
			}
			
		}
	
	}
	if (spl.isMonoPeriodic() ) sf=Sarea/( ( prb.h2() -prb.h1())*spl.boundWidth() );
	else cout<<"@SolidFractionInHeightProbe, run only with X periodiq sample "<<endl;
	
	return sf;
}


// For disk only
double solidFraction( circularProbe & prb, Sample& spl,Network & nwk)
{
	//cout<<" Solid fraction in circularProbe "<<endl;
	unsigned int Nb = spl.lbody().size();
	if(Nb==0) return 0;
	
	double lx,ly;
	unsigned int nc=0 ;
	double r1,r2;
	
	double d,sf=0.;	
	
	//Run only with disk
	for( unsigned int i=0; i< Nb ; ++i)
	{
		if( spl.body(i)->sizeVerlet()==0.) continue;
		
		/*if ( prb.probeInsideBody( spl.body(i)))
		{
			//sf+=prb.area();
			cout<<" insideBody "<<endl;
			return 1.;
		}*/
		
		if ( prb.intersection( spl.body(i) ) || prb.containCenter( spl.body(i) ) )
		{
		//	cout<<" intersection  ";
			 if( prb.containEntireBody(  spl.body(i) ) )
			{

			//	cout<<" EntireBody "<<spl.body(i)->Area()<<endl;
			//	cout<<spl.body(i)->x()<<" "<<spl.body(i)->y()<<" "<<spl.body(i)->sizeVerlet()<<endl;
			//	cout<<prb.x()<< " "<<prb.y()<<" "<<prb.R()<<endl;
				sf += spl.body(i)->Area();
				//cout<<" entire "<<endl;
			}
			else
			{
			lx= prb.x() - spl.body(i)->x();
			ly= prb.y() - spl.body(i)->y();
			d=sqrt(lx*lx + ly*ly);
			//Run only with disks.....
			
			r1= prb.R();
			r2= spl.body(i)->sizeVerlet();
		//	cout<<d<<" "<<r1<<" "<<r2<<endl;
			sf+= ( 
				  r1*r1 * acos( (d*d + r1*r1 - r2*r2)/(2.*d*r1))
				+ r2*r2 * acos( (d*d - r1*r1 + r2*r2)/(2.*d*r2))
				-.5*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
						);
					//cout<<endl;
			}
		}
	//	cout<<sf<<endl;
		
		
	}
	//cout<<sf<<" sf  "<<flush;
	//Check for overlap	
	for (unsigned int i=0;i< nwk.clist().size();++i)
	{
		nc = nwk.clist(i);
		
		if ( nwk.inter(nc)->rang()>0 && prb.contain( nwk.inter(nc)) && nwk.inter(nc)->type()==_type_dkdk ) 
		{
			
			lx= nwk.inter(nc)->first()->x() - nwk.inter(nc)->second()->x();
			ly= nwk.inter(nc)->first()->y() - nwk.inter(nc)->second()->y();
			d=sqrt(lx*lx + ly*ly);
			//Run only with disks.....
			
			r1= nwk.inter(nc)->first()->sizeVerlet();
			r2= nwk.inter(nc)->second()->sizeVerlet();
			
			if( r1==0 || r2==0 ) continue;
			//cout<<d-r1-r2<<" "<<(d*d + r1*r1 - r2*r2)/(2.*d*r1)<<endl;
			if( d-r1-r2 < -spl.rmin()*.01)
			sf-= ( 
				  r1*r1 * acos( (d*d + r1*r1 - r2*r2)/(2.*d*r1))
				+ r2*r2 * acos( (d*d - r1*r1 + r2*r2)/(2.*d*r2))
				-.5*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
						);
					//cout<<i<<" "<<sf<<endl;

		}
	
	}
	//cout<<sf<<" ** "<<sf/prb.area()<<endl;
	sf /= prb.area(); 
	
	return( sf );
}

//for disk only
double solidFraction( rectangularProbe & prb, Sample& spl,Network & nwk)
{
	//cout<<" Solid fraction in circularProbe "<<endl;
	unsigned int Nb = spl.lbody().size();
	if(Nb==0) return 0;
	
	vector <dof *> usedDof;
 	vector <dof *> ::iterator it;
	bool used = false;
	
	double lx,ly;
	unsigned int nc=0 ;
	double r1,r2;
	
	double d,sf=0.;	
	
	//Run only with disk
	for( unsigned int i=0; i< Nb ; ++i)
	{
		if( spl.body(i)->sizeVerlet()==0.) continue;
		
		if( spl.body(i)->type()==2){ continue;}
		/*if ( prb.probeInsideBody( spl.body(i)))
		{
			//sf+=prb.area();
			cout<<" insideBody "<<endl;
			return 1.;
		}*/
		//cout<<"corps i= "<<i<<" type = "<<spl.body(i)->type()<<endl;
	
		//	cout<<" intersection  ";
			if( prb.containCenter( spl.body(i)) || prb.containEntireBody(  spl.body(i) ) )
			{

			//	cout<<" EntireBody "<<spl.body(i)->Area()<<endl;
			//	cout<<spl.body(i)->x()<<" "<<spl.body(i)->y()<<" "<<spl.body(i)->sizeVerlet()<<endl;
			//	cout<<prb.x()<< " "<<prb.y()<<" "<<prb.R()<<endl;
			//	cout<<"	contenu "<<endl;
				if( spl.body(i)->bodyDof()!=NULL)
				{
					used=false;
					for(it=usedDof.begin();it!=usedDof.end();it++)
					{
						if( *it==spl.body(i)->bodyDof())
						{
							used=true;
							break;
						}
					}
					
					if( ! used )
					{
						/*cout<<"	corps pris en compte ";
						for( unsigned int j=0;j<spl.body(i)->bodyDof()->lctrlBodies().size();++j)
							cout<<spl.body(i)->bodyDof()->ctrlBody(j)->id()<<" ";
						cout<<endl;	
						*/
						usedDof.push_back( spl.body(i)->bodyDof());
						sf+= spl.body(i)->bodyDof()->Area();
					}
				}
				else
				sf += spl.body(i)->Area();
				//cout<<" entire "<<endl;
			}
			else
			{
				/*cout<<"intersection ?? "<<prb.intersection( spl.body(i))<<endl;
				cout<< spl.body(i)->xmax() <<" < ? "<<prb.x()+prb.hl()<<endl;
				cout<< spl.body(i)->xmin() <<" > ? "<<prb.x()-prb.hl()<<endl;
				cout<< spl.body(i)->ymax() <<" < ? "<<prb.y()+prb.hh()<<endl;
				cout<< spl.body(i)->ymin() <<" > ? "<<prb.y()-prb.hh()<<endl;
				*/
					
				
			}
			
		
		//cout<<sf<<endl;
		
		
	}
	cout<<sf<<" sf  "<<flush;
	//Check for overlap	
	for (unsigned int i=0;i< nwk.clist().size();++i)
	{
		nc = nwk.clist(i);
		
		if ( nwk.inter(nc)->rang()>0 && prb.contain( nwk.inter(nc)) && nwk.inter(nc)->type()==_type_dkdk ) 
		{
			
			lx= nwk.inter(nc)->first()->x() - nwk.inter(nc)->second()->x();
			ly= nwk.inter(nc)->first()->y() - nwk.inter(nc)->second()->y();
			d=sqrt(lx*lx + ly*ly);
			//Run only with disks.....
			
			r1= nwk.inter(nc)->first()->sizeVerlet();
			r2= nwk.inter(nc)->second()->sizeVerlet();
			
			if( (d-r1-r2)/ r1 < -1e-9) continue;
			if( (d-r1-r2)/ r2 < -1e-9) continue;
			
			if( r1==0 || r2==0 ) continue;
			//cout<<d-r1-r2<<" "<<(d*d + r1*r1 - r2*r2)/(2.*d*r1)<<endl;
			//if( d-r1-r2 < -spl.rmin()*.01)
			sf-= ( 
				  r1*r1 * acos( (d*d + r1*r1 - r2*r2)/(2.*d*r1))
				+ r2*r2 * acos( (d*d - r1*r1 + r2*r2)/(2.*d*r2))
				-.5*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
						);
					//cout<<i<<" "<<sf<<endl;

		}
	
	}
	cout<<sf<<" ** "<<sf/prb.area()<<endl;
	sf /= prb.area(); 
	
	return( sf );
}


