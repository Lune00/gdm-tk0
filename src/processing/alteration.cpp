#include "alteration.hpp"

// Not yet tested
void limit_speed(Sample& spl, double vlim)
{
	double vx,vy,vel2;
	double limitor;

	for (unsigned int i=0 ; i<spl.lbody().size() ; ++i)
	{
		vx = spl.body(i)->vx();
		vy = spl.body(i)->vy();
		vel2 = vx*vx + vy*vy;
		if (vel2 > (vlim*vlim)) 
		{
			limitor = vlim / sqrt(vel2);
			spl.body(i)->vx() *= limitor;
			spl.body(i)->vy() *= limitor;
		}
	}
}

// Not yet tested
void homothetie(Sample& spl, double h, double xc, double yc)
{
	if (h == 1.0) return;

	for (unsigned int i=0 ; i<spl.lbody().size() ; ++i)
	{
		spl.body(i)->x() = xc + (spl.body(i)->x() - xc)*h;
		spl.body(i)->y() = yc + (spl.body(i)->y() - yc)*h;
	}
}

void insert_polyg_in_disk(Sample& spl, unsigned int nVertexMin, unsigned int nVertexMax,double dR, double aspectRatio)
{
	body2d * B = 0;
	double Rmax;

	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		if(typeid(*(spl.body(i))) == typeid(disk))
		{
			string s("polyg");
			B = body2d::factory(s);
			B->x() = spl.body(i)->x();
			B->y() = spl.body(i)->y();
			B->rot() = 0.0;
			Rmax = (dynamic_cast<disk*>(spl.body(i)))->R();
			mutate_polyg(B, nVertexMin, nVertexMax, Rmax, dR, aspectRatio,false);
			//B->adjustCenter();
			spl.substituteBody(i,*B);
		}
	}
}

void insert_regular_polyg_in_disk(Sample& spl, unsigned int nVertex)
{
	body2d * B = 0;
	double Rmax=1;

	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		if(typeid(*(spl.body(i))) == typeid(disk))
		{
			string s("polyg");
			B = body2d::factory(s);
			B->x() = spl.body(i)->x();
			B->y() = spl.body(i)->y();
			B->rot() = 0.0;
			Rmax = (dynamic_cast<disk*>(spl.body(i)))->R();
			mutate_regular_polyg(B, nVertex, Rmax);
			//B->adjustCenter();
			spl.substituteBody(i,*B);
			//cout<<"aire calc = "<<.5*Rmax*Rmax<<endl;
			//cout<<"aire "<<B->Area()<<endl;
			//B->write( cout );
			//cout<<endl;
		}
	}
}


void insert_polyg_in_disk_CEGEO(Sample& spl, unsigned int nVertexMin, unsigned int nVertexMax, double eta)
{
	body2d * B = 0;
	double Rmax;
	cout<<" poly in circle CEGEO"<<endl;
	unsigned int Vmin=(unsigned int) ceil( 3.28/acos(1-eta));
	vector <double> rsave;
	if( nVertexMin< Vmin) nVertexMin=Vmin;

	spl.updateBoundaries();

	double height= spl.ymax()-spl.ymin();
	double width = spl.xmax()-spl.xmin();

	cout<<"Vmax "<<nVertexMax<<" Vmin "<<nVertexMin<<endl;
	cout<<" eta "<<eta<<endl;

	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		if(typeid(*(spl.body(i))) == typeid(disk))
		{
			string s("polyg");
			B = body2d::factory(s);
			B->x() = spl.body(i)->x();
			B->y() = spl.body(i)->y();
			B->rot() = 0.0;
			Rmax = (dynamic_cast<disk*>(spl.body(i)))->R();
			rsave.push_back(Rmax);
			mutate_polyg(B, nVertexMin, nVertexMax, Rmax, eta*Rmax, 1.,true);
			spl.substituteBody(i,*B);
		}
	}
	

	cout<<" Ecriture PS"<<endl;
	ofstream ps("verif.ps",ios::out);

	double zoom=min( 842/height, 595/width);
	cout<<"zoom "<<zoom<<endl;

//if (removeVisu ) removeRattlers();

	ps<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	ps<<"%%BoundingBox: 0 0 595 842"<<endl;
	ps<<"%%Pages: 1"<<endl;
	ps<<"1 setlinecap"<<endl;
	ps<<"% Procedure C3: trace un cercle vide de centre xy"<<endl;
	ps<<"% de rayon r  :"<<endl;
	ps<<"% x y r  C "<<endl;
	ps<<"/C3{newpath  0 360 arc gsave grestore stroke}def"<<endl;
	ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
	ps<<"% Procedure Lis: ligne de largeur fixe"<<endl;
	ps<<"% du point 0 au point 1"<<endl;
	ps<<"% x0 y0 x1 y1 Lis"<<endl;
	ps<<"/Lis{gsave newpath moveto lineto stroke grestore}def"<<endl; 
	ps<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;


	//******************
	//ps<<"%%Page: 1 1"<<endl;//Particules
	ps<<"1 setlinewidth"<<endl;
	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		
		double Rmax = rsave[i];
		ps<<"1 0 0 setrgbcolor"<<endl;
		ps<<spl.body(i)->x()*zoom<<" "<<spl.body(i)->y()*zoom<<" "<<(1-eta)*Rmax*zoom<<" C3"<<endl;
		ps<<"0 0 0 setrgbcolor"<<endl;
		ps<<spl.body(i)->x()*zoom<<" "<<spl.body(i)->y()*zoom<<" "<<Rmax    *zoom<<" C3"<<endl;
	
	
		polyg * temp = dynamic_cast<polyg*>(spl.body(i));
		//cout<<temp->Rout()<<endl;
		unsigned int next;
		for( unsigned int j=0;j<temp->Vertex().size();++j)
		{
			if (j==temp->Vertex().size()-1) next =0;
			else next=j+1;
			ps<<(temp->x() + temp->Vertex(j).x())*zoom<<" "<<(temp->y() + temp->Vertex(j).y())*zoom<<" "
				<<(temp->x() + temp->Vertex(next).x())*zoom<<" "<<(temp->y() + temp->Vertex(next).y())*zoom<<" Lis"<<endl; 
		}
	//	temp->adjustCenter();
	
	}
	
	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		if(typeid(*(spl.body(i))) == typeid(polyg))
		{
		 (dynamic_cast<polyg*>(spl.body(i)))->convexify();
		cout<<" convexify OK : body "<<i<<endl;
		}
	}
	
	for(unsigned int i=0 ; i<spl.lbody().size();++i)
	{
		ps<<".5 setlinewidth"<<endl;
		
		ps<<"0 1 0 setrgbcolor"<<endl;
		
		polyg * temp = dynamic_cast<polyg*>(spl.body(i));
		//cout<<temp->Rout()<<endl;
		unsigned int next;
		for( unsigned int j=0;j<temp->Vertex().size();++j)
		{
			if (j==temp->Vertex().size()-1) next =0;
			else next=j+1;
			ps<<(temp->x() + temp->Vertex(j).x())*zoom<<" "<<(temp->y() + temp->Vertex(j).y())*zoom<<" "
				<<(temp->x() + temp->Vertex(next).x())*zoom<<" "<<(temp->y() + temp->Vertex(next).y())*zoom<<" Lis"<<endl; 
		temp->adjustCenter();
		}
	}
	ps.close();

}


void mutate_polyg(body2d *B, unsigned int nVertexMin, unsigned int nVertexMax, double Rmax, double dR, double aspectRatio,bool chord_check)
{
	if (nVertexMin < 3)  nVertexMin = 3;
	if (nVertexMax < 3)  nVertexMax = 3;

	if (aspectRatio > 1.0) aspectRatio = 1.0/aspectRatio;

	polyg * P;
	if (typeid(*B) == typeid(polyg))
		P = dynamic_cast<polyg*>(B);
	else return;

	P->Vertex().clear();
	P->Normal().clear();

// Random number of vertexes selected in range [nVertexMin nVertexMax]
	unsigned int n = (unsigned int) floor(nVertexMin + ((double)rand()/(double)RAND_MAX)*(nVertexMax-nVertexMin+1));

	double R_alea,Theta,Theta_alea;
	double Rmin=Rmax-dR;
	gdm::vertex v;
	bool badChord;
	unsigned int Ntry=0;
	double mult;
//	cout<<" new polyg"<<endl;
	for (unsigned int s=0;s<n;++s)
	{
		//cout<<"---"<<s<<"/"<<n-1<<endl;
		badChord=true;
		Ntry=0;
		mult=1;
		while(  badChord )
		{
			if( Ntry>100)
			{
				Ntry=0;
				mult+=.1;
			}
			R_alea = Rmax - ((double)rand()/(double)RAND_MAX)*dR;
			Theta = 6.28/(double) n;
			Theta_alea =  (2.0*((double)rand()/(double)RAND_MAX)-1.0)*(mult*3.0/(double)n);
			v.x() = R_alea * cos(s * Theta + Theta_alea);
			v.y() = R_alea * sin(s * Theta + Theta_alea) * aspectRatio;

			if( chord_check )
			{
				if ( s > 0)
				{
					double dx = v.x() - P->Vertex(s-1).x();
					double dy = v.y() - P->Vertex(s-1).y();
					double d  = P->Vertex(s-1).x()* v.y() - P->Vertex(s-1).y()*v.x();
					double D = Rmin*Rmin*( dx*dx+dy*dy) - d*d;
					//cout<<" D= "<<D<<endl;
					if ( D >= 0)//Intersection
						badChord=true;
					else 
						badChord=false;

					if( s == n-1 )
					{
						//cout<<" dernier v"<<endl;
						dx = v.x() - P->Vertex(0).x();
						dy = v.y() - P->Vertex(0).y();
						d  = P->Vertex(0).x()* v.y() - P->Vertex(0).y()*v.x();

						if ( Rmin*Rmin*( dx*dx+dy*dy) - d*d >= 0)
							badChord=true;
						else 
							badChord=false;
					}
				}
				else 
					badChord=false;
			}
			else 
				badChord=false;
			
			Ntry++;
			//cout<<badChord<<endl;
			//getchar();
		}
	//	cout<<" 	add vertex "<<endl;
		P->Vertex().push_back(v);
		P->Normal().push_back(v);	
	}
	//P->adjustCenter();
}

void mutate_regular_polyg(body2d *B, unsigned int nVertex, double Rmax)
{
	if (nVertex < 3)  nVertex = 3;
	polyg * P;
	if (typeid(*B) == typeid(polyg))
		P = dynamic_cast<polyg*>(B);
	else return;

	P->Vertex().clear();
	P->Normal().clear();

	gdm::vertex v;
	double offset=0. ;
	
	for (unsigned int s=0;s < nVertex;++s)
	{	
		if( s==0) 
		{
			offset =  (2.0*((double)rand()/(double)RAND_MAX)-1.0)* 6.28;
		}
			double Theta = 2.*M_PI/(double) nVertex;			
			
			v.x() = Rmax * cos(s * Theta  );
			v.y() = Rmax * sin(s * Theta  );
			//cout<<s * Theta/M_PI*180<<" "<<sqrt(v.x()*v.x() + v.y()*v.y())<<endl;
			P->rot() = offset;			
		
		P->Vertex().push_back(v);
		P->Normal().push_back(v);	
	}
	P->adjustCenter();
	P->Fill(1.);
	
	//cout<<"Area = "<<P->Area()<<" Caculus = "<<nVertex*Rmax*Rmax*pow(sin( M_PI/(double) (nVertex)),2)/tan( M_PI/(double) (nVertex))<<endl;

	
}


//! \Insertion of orientation disorder in 3 disks' cluster sample by multiplying rotation mass center by random number 
//! \authors C. Voivret and B. Saint-Cyr
vector <dof*> insert_cluster_in_disk(Sample& spl, unsigned int nDisk,double ratio)
{
	
	vector <dof*> ldof;
	double randomRot;
	unsigned int i=0;
	while(  i<spl.lbody().size()) 
	{
        
		//Random number calculus
		randomRot=rand()/(double) RAND_MAX*2.*M_PI;		
		
		if(typeid(*(spl.body(i))) == typeid(disk) && spl.body(i)->bodyDof()==NULL)
		{
			//cout<<"mutation corps "<<i<<endl;
			ldof.push_back( mutate_cluster_baptiste( spl.body(i),spl,3,ratio) );
			i=0;
			
			//Display random number generated
			cout<<" radomRot= "<<randomRot <<endl;
			
			//Mass center calculus
			ldof.back()->mcrot()=randomRot;
			
			//Filling clusters with unity mass j'ai mis un double en argument dans la fonction fill
			ldof.back()->fill(1.);
			
			//Mass center calculus
			ldof.back()->computeMassCenter();
			
			//Vertex computation -> center elementaries disks position from cluster mass center position 
			ldof.back()->computeVertex();
			
			//Applying rotation to mass center by move function
			ldof.back()->move(0.);
			
		}// test: if following body takes part to one cluster
		else if(spl.body(i)->bodyDof()!=NULL)
		{
		//cout<<"appartient deja a 1 cluster i= "<<i<<endl;
		++i;
		
		}
		
	}
	return ldof;
}


dof* mutate_cluster(body2d *B, Sample& spl,unsigned int nDisk,double ratio)
{
	if (nDisk < 2)  nDisk = 2;
	
	double r= B->sizeVerlet();
	double x= B->x();
	double y= B->y();
	spl.removeBody(B);
	
	dof * doftemp = new dof;
	for(unsigned int i=0; i< nDisk;++i)
	{
		double angle= 2.*M_PI/(double) nDisk;
		spl.lbody().push_back( new disk( x + r*(1.-ratio)*cos(i*angle),y+r*(1.-ratio)*sin(i*angle),ratio*r));
		doftemp->plugBody( spl.lbody().back());
	}
	return doftemp;
	
}

//! \Insert 3 disk's cluster in one fixed radius disk
//! \author B. Saint-Cyr 
//Mutate_cluster_baptiste function = dof type pointor
//Introduction of geometric parameter dr_r in function argues
dof* mutate_cluster_baptiste(body2d *B, Sample& spl,unsigned int nDisk,double eta)
{
	if (nDisk < 2)  nDisk = 2;
	
	cout<<" fonction baptiste "<<endl;
	//double eta;
	double alpha;
	double aire=0;
	double r= B->sizeVerlet();
	double x= B->x();
	double y= B->y();
	// effacage du corps c-a-d grand disque
	spl.removeBody(B);
	
	dof * doftemp = new dof;
	
	//cout<<" eta "<<eta<<" "<<.5*sqrt(3)<<endl;
	
	//Alpha angle calculus  
	alpha = 2*asin((sqrt(3)/2.)*(sqrt(1.-eta*eta)-eta/sqrt(3)));
	//cout<<"alpha"<<alpha<<endl;
	
	
	double angle= 2.*M_PI/(double) nDisk;
	for(unsigned int i=0; i< nDisk;++i)
	{
		
		
		//Calculus of the new coordinates of cluster's disks centers 
		spl.lbody().push_back( new disk( x + r*(1.-sin(angle/2.)/(eta+sin(angle/2.)))*cos(angle/2.+angle*i),y+r*(1.-sin(angle/2.)/(eta+sin(angle/2.)))*sin(angle/2+angle*i),r*sin(M_PI/3.)/(sin(M_PI/3.)+eta)));
		
		//Display new coordinates of cluster's disks centers and their new radius
		//cout<<" x "<< r*(1.-sin(angle/2.)/(eta+sin(angle/2.)))*cos(angle/2.+angle*i)<<" y "<<r*(1.-sin(angle/2.)/(eta+sin(angle/2.)))*sin(angle/2.+angle*i)<<" radius "<<r*sin(M_PI/3.)/(sin(M_PI/3.)+eta)<<endl;
		

		//Apply unity mass to disk's cluster
		spl.lbody().back()->Fill(1);
		doftemp->plugBody( spl.lbody().back());
		
		
	}
	
	double h = M_PI/3.;
	double a = r*sin(h)/(eta+sin(h));
	//aire =3.*M_PI*a*a;

	if ( eta>=sqrt(3)/2 && eta<=1 ) //vide interne considéré comme une phase solide
	{
		aire =  3*M_PI*a*a-6*a*a*(acos(eta)-eta*sqrt(1-eta*eta)) + a*a*(eta*eta*sqrt(3) - M_PI/2 +3*(acos(eta)-eta*sqrt(1-eta*eta)));
	}
	
	else if (eta >=0 && eta< .5*sqrt(3)) // Prise en compte de l'intersection des trois disques
	{
		aire = 3*M_PI*a*a-6*a*a*(acos(eta)-eta*sqrt(1-eta*eta)) + (3*sqrt(3)*a*a*0.5)*(sqrt(1-eta*eta)-eta/sqrt(3)) * (0.5-cos(alpha/2)) +3*a*a*alpha*0.5; 
	}
	
	//aire = 3*M_PI*a*a-6*a*a*(acos(eta)-eta*sqrt(1-eta*eta)) + (3*sqrt(3)*a*a*0.5)*(sqrt(1-eta*eta)-eta/sqrt(3)) * (0.5-cos(alpha/2)) +3*a*a*alpha*0.5;  
	
	cout<<" aire "<<aire<<endl;
	doftemp->Area()=aire;
	
	return doftemp;
	
}

