#include "near.hpp"
// #include <typeinfo>

// nearFactory
bool near(const body2d* b1, const body2d* b2, const double dv)
{
  // const type_info &b1_type = typeid(*b1);
  // const type_info &b2_type = typeid(*b2);
  
	if ( (b1->bodyDof() != 0 ) && (b2->bodyDof() != 0 ) && ( b1->bodyDof() == b2->bodyDof() ) ) 
		return false;

  // dkdk
  if ( (typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(disk)) )
    return near(dynamic_cast<const disk*>(b1), dynamic_cast<const disk*>(b2),dv);

  // pgpg
  if ( (typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(polyg)) )
    return near(dynamic_cast<const polyg*>(b1), dynamic_cast<const polyg*>(b2),dv); 
  
  // rlrl
  //if ( (typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(rline)) )
  //  return near(dynamic_cast<const rline*>(b1), dynamic_cast<const rline*>(b2),dv);
  
  // dkrl
  if ( (typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(rline)) )
    {
    return near(dynamic_cast<const disk*>(b1),dynamic_cast<const rline*>(b2),dv);
    }
  if ( (typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(disk)) )
    {
    return near(dynamic_cast<const disk*>(b2),dynamic_cast<const rline*>(b1),dv);
    }
  
  // dkpg
  if ( (typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(polyg)) )
    {
    return near(dynamic_cast<const disk*>(b1),dynamic_cast<const polyg*>(b2),dv);
    }
  if ( (typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(disk)) )
    {
    return near(dynamic_cast<const disk*>(b2),dynamic_cast<const polyg*>(b1),dv);
    }

//pgrline
  if ( (typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(rline)) )
    {
    return near(dynamic_cast<const polyg*>(b1),dynamic_cast<const rline*>(b2),dv);
    }
  if ( (typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(polyg)) )
    {
    return near(dynamic_cast<const polyg*>(b2),dynamic_cast<const rline*>(b1),dv);
    }
	
  return false;
}

// dkdk
//dverlet : espace entre deux particules inferieur a dverlet?
bool near(const disk* d1, const disk* d2, const double dv)
{
  double lx  = d2->x() - d1->x();
  double ly  = d2->y() - d1->y();
  double sum = dv + d1->R() + d2->R();
	
  return ((lx*lx+ly*ly) < sum*sum);
}

// plpl
bool near(const polyg* p1, const polyg* p2, const double dv)
{
  double lx  = p2->x() - p1->x();
  double ly  = p2->y() - p1->y();
  double sum = dv + p1->Rout() + p2->Rout();
	
  if ((lx*lx+ly*ly) > sum*sum) return false;
  
  unsigned int v = 0;
  double c,s;
  unsigned int iCriticalVertex = 0, jCriticalVertex = 0;
  double iScalarProduct = -1.0E+10, jScalarProduct = 1.0E+10, vScalarProduct;
  double sepx, sepy, norm;
  double shadowDist;
  
  unsigned int nVertexi = p1->Vertex().size();
  unsigned int nVertexj = p2->Vertex().size();
  
  gdm::vertex * iVertex[nVertexi];
  gdm::vertex * jVertex[nVertexj];
  
  // Global vertex coordinates
  // Polygon i
  c = cos(p1->rot());
  s = sin(p1->rot());
  for(v = 0; v < nVertexi; ++v)
    {
    iVertex[v]->x() = p1->x() + c * p1->Vertex(v).x() - s * p1->Vertex(v).y();
    iVertex[v]->y() = p1->y() + s * p1->Vertex(v).x() + c * p1->Vertex(v).y();
    }
  
  // Polygon j
  c = cos(p2->rot());
  s = sin(p2->rot());	
  for(v = 0; v < nVertexj ; ++v)
    {
    jVertex[v]->x() = p2->x() + c * p2->Vertex(v).x() - s * p2->Vertex(v).y();
    jVertex[v]->y() = p2->y() + s * p2->Vertex(v).x() + c * p2->Vertex(v).y();
    }  
  
  // Separation vector (i.e. unit vector of the intercenter line)
  sepx = p2->x() - p1->x();
  sepy = p2->y() - p1->y();
  sepx *= (norm = 1.0 / sqrt(sepx * sepx + sepy * sepy));
  sepy *= norm;
  
  // Research of the critical vertexes and of the corresponding overlap 
  for(v = 0 ; v < nVertexi ; ++v)
    {
    // Max for body i
    if((vScalarProduct = sepx * iVertex[v]->x() + sepy * iVertex[v]->y()) > iScalarProduct)	
      {
      iCriticalVertex = v ;
      iScalarProduct = vScalarProduct;
      }
    }
  for(v = 0 ; v < nVertexj ; ++v)
    {
    // Min for body j
    if((vScalarProduct = sepx * jVertex[v]->x() + sepy * jVertex[v]->y()) < jScalarProduct)
      {
      jCriticalVertex = v ;
      jScalarProduct = vScalarProduct ;
      }
    }
  shadowDist = jScalarProduct - iScalarProduct;
  if(shadowDist > dv) return false;
  
  // ... 
  
  return true;
}

// rlrl
/*
bool near(const rline* d, const rline* rl, const double dv)
{  
  // TODO
  return true;
}
*/

// dkrl
bool near(const disk* d, const rline* rl, const double dv)
{
  double h   = 0.5*rl->L();
  double dxm = h*fabs(cos(rl->rot()))+rl->R()+d->R()+dv;
  double dym = h*fabs(sin(rl->rot()))+rl->R()+d->R()+dv;

  return ((fabs(d->x() - rl->x()) < dxm) && (fabs(d->y() - rl->y()) < dym));
}

// pgrl
bool near(const polyg* pg, const rline* rl, const double dv)
{
	//cout<<"near pgrl "<<endl;
  double h   = 0.5*rl->L();
  double dxm = h*fabs(cos(rl->rot()))+rl->R()+pg->Rout()+dv;
  double dym = h*fabs(sin(rl->rot()))+rl->R()+pg->Rout()+dv;

  return ((fabs(pg->x() - rl->x()) < dxm) && (fabs(pg->y() - rl->y()) < dym));
}

// dkpg
bool near(const disk* d, const polyg* p, const double dv)
{
  double lx  = p->x() - d->x();
  double ly  = p->y() - d->y();
  double sum = dv + d->R() + p->Rout();
	
  return ((lx*lx+ly*ly) < sum*sum);
}

