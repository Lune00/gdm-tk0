#include "polyg.hpp"

void polyg::read (istream & is)
{
  unsigned int nbVertexes = 0;
  Rout_ = 0.0;
  Rin_ = 1.0E+10;
  
  is >> grp_ >> nbVertexes;
  gdm::vertex readedVertex;
  double dist = 0.0;
  for(unsigned int i=0 ; i<nbVertexes ; ++i)
    {
    is >> readedVertex.x() >> readedVertex.y();
    vertex_.push_back(readedVertex);
    normal_.push_back(readedVertex);
    
    dist = sqrt(readedVertex.x() * readedVertex.x() + readedVertex.y() * readedVertex.y());
    if (dist > Rout_) Rout_ = dist;
    if (dist < Rin_) Rin_ = dist;
    }
  
  is >> x_  >> y_  >> rot_
     >> vx_ >> vy_ >> vrot_;
  x0_=x_;
  y0_=y_;
  rot0_=rot_;
  double nx,ny;
  double norm;
	//cout<< Vertex().size()<<endl;
  for(unsigned int i=0 ; i<nbVertexes-1 ; ++i)
    {
    nx = vertex_[i+1].y() - vertex_[i].y();
    ny = vertex_[i].x()   - vertex_[i+1].x();
    normal_[i].x() = nx * (norm = 1.0 / sqrt(nx * nx + ny * ny));
    normal_[i].y() = ny * norm;
    }
  
  nx = vertex_[0].y() - vertex_[nbVertexes-1].y();
  ny = vertex_[nbVertexes-1].x() - vertex_[0].x();
  normal_[nbVertexes-1].x() = nx * (norm = 1.0 / sqrt(nx * nx + ny * ny));
  normal_[nbVertexes-1].y() = ny * norm;  
}

void polyg::write (ostream & os)
{
  os << "polyg " << grp_ << ' ' << vertex_.size() << endl << flush;
os.precision(15);  
  for(unsigned int i=0 ; i<vertex_.size() ; ++i)
    {
    os << vertex_[i].x() << ' ' << vertex_[i].y() << endl << flush;
    }
  
  os << x_  << ' ' <<  y_  << ' ' << rot_  << ' ' 
     << vx_ << ' ' << vy_  << ' ' << vrot_
     << endl << flush;
}


void polyg::writeM(ostream &os)
{

}

void polyg::adjustCenter ()
{
  if (vertex_.empty()) return;
  
  // Center of gravity
  double AGx = 0.0;
  double AGy = 0.0;
  double inv_n = 1.0 / (vertex_.size());
  for(unsigned int i=1 ; i<vertex_.size() ; ++i)
    {
    AGx += vertex_[i].x() - vertex_[0].x();
    AGy += vertex_[i].y() - vertex_[0].y();
    }
  AGx *= inv_n;
  AGy *= inv_n;
  
  double dx = vertex_[0].x() + AGx;
  double dy = vertex_[0].y() + AGy;
  
  // Correct the vertexes and particle position
  for(unsigned int i=0 ; i<vertex_.size() ; ++i)
    {
    vertex_[i].x() -= dx;
    vertex_[i].y() -= dy;
    }
  //x_ += dx;
  //y_ += dy;
 
  //Corrige by DH NGUYEN 20/2/2012
  c=cos(rot_);
  s=sin(rot_);
  x_ += dx*c-dy*s;
  y_ += dx*s+dy*c;
  
  for(unsigned int i=0 ; i<vertex_.size() ; ++i)
 {
	double d = sqrt ( vertex_[i].x()*vertex_[i].x()+vertex_[i].y()*vertex_[i].y());
	Rout_ = Rout_ > d ? Rout_ : d;
 }
}

void polyg::writeMGP (ostream & os)
{
  os << "    <POLYG id=\"" << id_ << "\" r=\"" << Rout_ << "\">" << endl
     << "     <position x=\"" << x_ 
     << "\" y=\"" << y_ << "\" rot=\"" << rot_ << "\"/>" << endl
     << "     <velocity x=\"" << vx_ 
     << "\" y=\"" << vy_ << "\" rot=\"" << vrot_ << "\"/>" << endl; 
     
  for (unsigned int num = 0 ; num < vertex_.size() ; ++num)
    {
    os << "     <node x=\"" << vertex_[num].x() 
       << "\" y=\"" << vertex_[num].y() << "\"/>" << endl;
    }
  os << "    </POLYG>" << endl << flush;
}

void polyg::writePS (ostream & os)
{

}

polyg* polyg::duplicate()
{
  polyg* copy = new polyg;
  *copy = *this; 
  return copy;
}

void polyg::Fill(double density)
{
  double totalArea   = 0.0;
  double triangleMom = 0.0;
  
  unsigned int nVertex = vertex_.size();
  double Area[nVertex];
  double xGrav[nVertex];
  double yGrav[nVertex];
  
  // the arrays Area, xGrav and yGrav are requiered in the following
  for(unsigned int i = 1 ; i < nVertex-1 ; ++i)
    {
		Area[i] = 0.5 * ((vertex_[i].x() - vertex_[0].x()) * (vertex_[i+1].y() - vertex_[0].y())
                   -   (vertex_[i].y() - vertex_[0].y()) * (vertex_[i+1].x() - vertex_[0].x()));
		//Area[i] = 0.5 * ((vertex_[i].x() ) * (vertex_[i+1].y() )
	      //           -  (vertex_[i].y() ) * (vertex_[i+1].x() ));
		if( Area[i] <0)
    	  {
      		//gdm::fatal("Polygon have a negative Area");
      		cerr << "Aire negative !!" << endl;
      		cerr << "Area[" << i << "] = " << Area[i] << endl;
      	  }
		totalArea += Area[i];
		xGrav[i] = 0.33333333 * (vertex_[0].x() + vertex_[i].x() + vertex_[i+1].x());
		yGrav[i] = 0.33333333 * (vertex_[0].y() + vertex_[i].y() + vertex_[i+1].y());
    }
  mass_ = totalArea * density;
  Rmass_ = sqrt(totalArea/M_PI);
  
  // Central moment of inertia -  Tam quan tinh
  mom_ = 0.0;
	for(unsigned int i = 1 ; i < nVertex-1 ; ++i)
    {
		triangleMom = (vertex_[i].x() - vertex_[0].x()) * (vertex_[i].x() - vertex_[0].x())
    + (vertex_[i].y()   - vertex_[0].y()) * (vertex_[i].y()   - vertex_[0].y())
    + (vertex_[i+1].x() - vertex_[0].x()) * (vertex_[i+1].x() - vertex_[0].x())
    + (vertex_[i+1].y() - vertex_[0].y()) * (vertex_[i+1].y() - vertex_[0].y())
    - (vertex_[i].x()   - vertex_[0].x()) * (vertex_[i+1].x() - vertex_[0].x())
    - (vertex_[i].y()   - vertex_[0].y()) * (vertex_[i+1].y() - vertex_[0].y());
		triangleMom *= Area[i]/18.0;
		mom_ += triangleMom + Area[i] * (xGrav[i] * xGrav[i] + yGrav[i] * yGrav[i]);
    }
  mom_ *= density;
}

void polyg::convexify()
{
	//cout<<" deb convexify"<<endl;
	unsigned int ilower=0,ihigher=0;
	unsigned int nVertex=this->Vertex().size();
	//cout<<" Nombre de sommets"<<nVertex<<endl;
	
	for (unsigned int i=1; i< nVertex-1;++i)
	{
		if( vertex_[i].y() < vertex_[ilower].y() )
			ilower=i;
		if( vertex_[i].y() > vertex_[ihigher].y() )
				ihigher=i;
	}
	//cout<<" PLus bas sommet "<<ilower<<endl;
	//cout<<" PLus haut sommet "<<ihigher<<endl;
	polyg temp;
	temp.vertex_.push_back(Vertex(ilower));
	
	double theta,Dthetamin,thetaprec=0,thetamin;
	
	unsigned int base=ilower;
	unsigned int isuiv=ilower;
	double dx,dy,d;
	do
	{
		Dthetamin=thetamin=8;
		//cout<<" sommet de recherche "<<base<<endl;
		
		for (unsigned int j=0; j< nVertex;++j)
		{
			if( base == j ) continue;
			//Calcul de l'angle 
			dx = this->Vertex(j).x() - this->Vertex(base).x();
			dy = this->Vertex(j).y() - this->Vertex(base).y();
			d  = sqrt( dx*dx+dy*dy);
			theta= acos( dx / d );
			if( dy < 0.) theta = 2.*M_PI-theta;
			
			//if( base < ihigher && theta < M_PI ) continue;
			//cout<<"		test de "<<j<<" angle abs "<<theta<<" angle relatif  "<<theta - thetaprec;
			
			if( fabs(theta - thetaprec) < Dthetamin)
			{
				Dthetamin= fabs(theta - thetaprec);
				thetamin=theta;
				isuiv = j;
			//	cout<<"----> sommet externe "<<isuiv;
			}
			//cout<<endl;
		
		}
		if ( isuiv != ilower) temp.Vertex().push_back( this->Vertex(isuiv));
		
		thetaprec=thetamin;
		base=isuiv;
		
		//cout<<" Sommet de depart suivant "<<isuiv<<" angprec "<<thetaprec<<endl<<endl;
		
				
	}
	while ( isuiv != ilower);
	//cout<<" convexification ok "<<endl;
	
	this->Vertex().clear();
	this->Normal().clear();
	this->Vertex()=temp.Vertex();
	
	
	double nx,ny;
	  double norm;
	nVertex=this->Vertex().size();

	  for(unsigned int i=0 ; i<nVertex-1 ; ++i)
	    {
	    nx = vertex_[i+1].y() - vertex_[i].y();
	    ny = vertex_[i].x()   - vertex_[i+1].x();
	    normal_[i].x() = nx * (norm = 1.0 / sqrt(nx * nx + ny * ny));
	    normal_[i].y() = ny * norm;
	    }

	  nx = vertex_[0].y() - vertex_[nVertex-1].y();
	  ny = vertex_[nVertex-1].x() - vertex_[0].x();
	  normal_[nVertex-1].x() = nx * (norm = 1.0 / sqrt(nx * nx + ny * ny));
	  normal_[nVertex-1].y() = ny * norm;
		
}

bool polyg::Contain(gdm::vertex &v)
{
	bool logic=true;
	unsigned int i=0, start, end;
	double signv;
	double signu;
	double ux, uy, startx, starty, endx, endy;
	gdm::vertex u;
	c=cos(rot_);
	s=sin(rot_);
	double d;
	
	d=sqrt((v.x()-this->x())*(v.x()-this->x())+(v.y()-this->y())*(v.y()-this->y()));
	
	if(d>this->sizeVerlet()) return false;
	else do
	{
		start=i;
		end=(i+1)%this->Vertex().size();
		unsigned int l=(end+1)%this->Vertex().size();
		u=vertex_[l];
		ux=this->x()+u.x()*c-u.y()*s;
		uy=this->y()+u.x()*s+u.y()*c;
		startx=this->x()+vertex_[start].x()*c-vertex_[start].y()*s;
		starty=this->y()+vertex_[start].x()*s+vertex_[start].y()*c;
		endx=this->x()+vertex_[end].x()*c-vertex_[end].y()*s;
		endy=this->y()+vertex_[end].x()*s+vertex_[end].y()*c;
		
		signv=(v.x()-startx)*(endy-starty)-(v.y()-starty)*(endx-startx);
		signu=(ux-startx)*(endy-starty)-(uy-starty)*(endx-startx);
		//cout<<"i:="<<i<<"	"<<signv*signu<<endl;
		if (signv*signu<0.) logic=false;
		else i++;
	}
	while (logic &&(i<this->Vertex().size()));
	return logic; 
}

vector<gdm::vertex> polyg::SectionSegment(gdm::vertex d1,gdm::vertex d2)
{
	gdm::vertex d;
	vector<gdm::vertex> d_;
	gdm::vertex p1,p2;
	double a1, b1, c1;
	double a2, b2, c2;
	double D, Dx, Dy;
	double c=cos(rot_);
	double s=sin(rot_);
	double p1x, p1y, p2x, p2y;
	double funcp1,funcp2,funcd1,funcd2;
	
	for(unsigned int i=0;i<vertex_.size();i++)
	{
		unsigned int j=(i+1) %vertex_.size();
		p1=vertex_[i];
		p2=vertex_[j];
		p1x=this->x()+p1.x()*c-p1.y()*s;
		p1y=this->y()+p1.x()*s+p1.y()*c;
		p2x=this->x()+p2.x()*c-p2.y()*s;
		p2y=this->y()+p2.x()*s+p2.y()*c;
		
		//cout<<"p1x:= "<<p1x<<"  p1y:= "<<p1y<<"  p2x:= "<<p2x<<"  p2y:= "<<p2y<<endl;
		
		funcp1=(p1x-d1.x())*(d2.y()-d1.y())-(p1y-d1.y())*(d2.x()-d1.x());
		funcp2=(p2x-d1.x())*(d2.y()-d1.y())-(p2y-d1.y())*(d2.x()-d1.x());
		funcd1=(d1.x()-p1x)*(p2y-p1y)-(d1.y()-p1y)*(p2x-p1x);
		funcd2=(d2.x()-p1x)*(p2y-p1y)-(d2.y()-p1y)*(p2x-p1x);
		
		if( (funcp1*funcp2<0) && (funcd1*funcd2<0))
		{
			//cout<<"co giao diem!"<<endl;
			a1=d2.y()-d1.y();
			b1=d1.x()-d2.x();
			c1=d1.x()*(d2.y()-d1.y())-d1.y()*(d2.x()-d1.x());
			a2=p2y-p1y;
			b2=p1x-p2x;
			c2=p1x*(p2y-p1y)-p1y*(p2x-p1x);
			
			D=a1*b2-a2*b1;
			Dx=-b1*c2+b2*c1;
			Dy=a1*c2-a2*c1;
			d.x()=Dx/D;
			d.y()=Dy/D;
		//	cout<<"funcp1:= "<<funcp1<<"  funcp2:= "<<funcp2<<"  funcd1:= "<<funcd1<<"funcd2:= "<<funcd2<<endl;
		//	cout<<"p1x:= "<<p1x<<"  p1y:= "<<p1y<<"  p2x:= "<<p2x<<"  p2y:= "<<p2y<<endl;
		//	cout<<"d.x():="<<d.x()<<"    "<<"d.y():="<<d.y()<<endl;
			d_.push_back(d);
		}
	}
	return d_;
}



void polyg::Arrangement()
{
	gdm::vertex point;
	for(unsigned int i=0; i<vertex_.size()-2;i++)
	{
		bool logic=true;//rerurn true si on ne trouve pas le sommet i+1 
		unsigned int j=i;
		double pointl, pointk;
		do //Cherche le sommet i+1
		{
			j++;
			bool ktra=true;//tra ve gia tri dung neu dinh j thoa man
			unsigned int l=(j+1)%vertex_.size();
			//kiem tra tat ca cac dinh tru i va j xem co nam cung phia ko
			
			pointl=(vertex_[l].x()-vertex_[i].x())*(vertex_[j].y()-vertex_[i].y())- (vertex_[l].y()-vertex_[i].y())*(vertex_[j].x()-vertex_[i].x());			
			for(unsigned int k=0;k<vertex_.size();k++)
			{
				if((k!=i)&&(k!=j)&&(k!=l))
				{
					pointk=(vertex_[k].x()-vertex_[i].x())*(vertex_[j].y()-vertex_[i].y())-(vertex_[k].y()-vertex_[i].y())*(vertex_[j].x()-vertex_[i].x());	
					if (pointl*pointk<0) 
					{
					//	cout<<"pointl*pointk:="<<pointl*pointk<<endl;
						ktra=false;			
					}	
				}
			}
			
			if (ktra) 
			{	
				logic=false;
				point=vertex_[i+1];
				vertex_[i+1]=vertex_[j];
				vertex_[j]=point;
			}
		}
		while (logic);		
	}
}

//Voir NR 3th Edition ; return value positif si polygone est CCW
int polyg::ispolysimple() 
{
	int i,ii,j,jj,np,schg=0,wind=0; //Initialize sign change and winding number.
	double p0,p1,d0,d1,pp0,pp1,dd0,dd1,t,tp,t1,t2,crs,crsp=0.0;
	np = vertex_.size();
	p0 = vertex_[0].x()-vertex_[np-1].x();
	p1 = vertex_[0].y()-vertex_[np-1].y();
	for (i=0,ii=1; i<np; i++,ii++)
	{
		if (ii == np) ii = 0;
		d0 = vertex_[ii].x()-vertex_[i].x();
		d1 = vertex_[ii].y()-vertex_[i].y();
		crs = p0*d1-p1*d0;
		if (crs*crsp < 0) 
		{
			//cout<<"crs*crsp:="<<crs*crsp<<endl;
			schg = 1;
		}
		if (p1 <= 0.0) 
		{
   			if (d1 > 0.0 && crs > 0.0) wind++;
		} 
		else
		{
    		if (d1 <= 0.0 && crs < 0.0) wind--;
		}
		p0=d0;
		p1=d1;
		if (crs != 0.0) crsp = crs;
	}
	if (abs(wind) != 1) return 0;
	if (schg == 0) return (wind>0? 1 : -1);

	for (i=0,ii=1; i<np; i++,ii++) 
	{
		if (ii == np) ii=0;
		d0 = vertex_[ii].x();
		d1 = vertex_[ii].y();
		p0 = vertex_[i].x();
		p1 = vertex_[i].y();
		tp = 0.0;
		for (j=i+1,jj=i+2; j<np; j++,jj++)
		{
			if (jj == np) {if (i==0) break; jj=0;}
			dd0 = vertex_[jj].x();
			dd1 = vertex_[jj].y();
			t = (dd0-d0)*(p1-d1) - (dd1-d1)*(p0-d0);
			if (t*tp <= 0.0 && j>i+1)
			{// First loop is only to compute starting tp,
				pp0 = vertex_[j].x(); //hence test on j.
				pp1 = vertex_[j].y();
				t1 = (p0-dd0)*(pp1-dd1) - (p1-dd1)*(pp0-dd0);
				t2 = (d0-dd0)*(pp1-dd1) - (d1-dd1)*(pp0-dd0);
				if (t1*t2 <= 0.0) return 0; //Found an intersection, so done.
			}
			tp = t; 
		}
	}
	return (wind>0? 2 : -2);
}

//Polygon is CCW (+) or CW (-)
void polyg::Invert()
{
	vector<gdm::vertex> q;
	q.push_back(vertex_[0]);
	for(unsigned int j=1;j<vertex_.size();j++) q.push_back(vertex_[vertex_.size()-j]);
	vertex_=q;
	q.clear();
}
