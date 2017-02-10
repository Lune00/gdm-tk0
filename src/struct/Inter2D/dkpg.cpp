#include "dkpg.hpp"

//vr EN TRAVAUX 

void dkpg::read(istream & is, unsigned int * Id1, unsigned int * Id2)
{
	
}

void dkpg::write(ostream & os)
{
}

double dkpg::Dist() // todo
{
  return 0.0;
}

double dkpg::Overlap() // todo
{
  return 0.0;
}

bool dkpg::Activate()
{
  if (rank_ > 0) return true;
  return false;
}

void dkpg::Frame()
{
  double cx,cy,cc;
  double dc,df,ff,dn;
  double xface,yface;
  unsigned int nbVertexes = j_->Vertex().size();
  unsigned int v,vPrev,vNext;
  
  // Global vertex coordinates and face normals for polygon j
  gdm::vertex * jVertex[nbVertexes];
  gdm::vertex * jNorm[nbVertexes];
  
  double c = cos(j_->rot());
  double s = sin(j_->rot());
  for(v = 0; v < nbVertexes ; ++v)
    {
    // Vertexes
    jVertex[v]->x() = j_->x() + c * j_->Vertex(v).x() - s * j_->Vertex(v).y();
    jVertex[v]->y() = j_->y() + s * j_->Vertex(v).x() + c * j_->Vertex(v).y();
    
    // Normales
    jNorm[v]->x() = c * j_->Normal(v).x() - s * j_->Normal(v).y();
    jNorm[v]->y() = s * j_->Normal(v).x() + c * j_->Normal(v).y();
    } 
  
  // We first check any contact with the polygon vertexes   
  for (v = 0; v < nbVertexes ; ++v)
    {
    cx = i_->x() - jVertex[v]->x();
    cy = i_->y() - jVertex[v]->y();
    cc = cx * cx + cy * cy;
    
    if (cc < i_->R()*i_->R())
      {
      vPrev = (v == 0) ? (nbVertexes - 1) : (v - 1);
      if (jNorm[vPrev]->x() * cy - jNorm[vPrev]->y() * cx < 0.0) break;
      if (jNorm[v]->x() * cy - jNorm[v]->y() * cx > 0.0) break;
      
      cc = sqrt(cc);
      dn = i_->R() - cc;
      cc = 1.0 / cc;
      nx_ = cx * cc;
      ny_ = cy * cc;
      tx_ = -ny_;
      ty_ = nx_;
      x_ = jVertex[v]->x() - 0.5 * dn * nx_;
      y_ = jVertex[v]->y() - 0.5 * dn * ny_;
	  overlap_= dn;
      	rank_ = 1;
		//rankco_=1;
      return;
      }
    
    }
  
  // There is no vertex-contact
  // So, we check contact with the polygon faces 
  for (v = 0;v<nbVertexes;++v)
    {
    cx = i_->x() - jVertex[v]->x();
    cy = i_->y() - jVertex[v]->y();
    dc = cx * jNorm[v]->x() + cy * jNorm[v]->y();
    
    if (dc > 0.0 && dc < i_->R())
      {
      vNext = (v == nbVertexes - 1) ? 0 : (v + 1);
      xface = jVertex[vNext]->x() - jVertex[v]->x();
      yface = jVertex[vNext]->y() - jVertex[v]->y();
      ff = yface * jNorm[v]->x() - xface * jNorm[v]->y();
      df = cy * jNorm[v]->x() - cx * jNorm[v]->y();
      
      if(df >= 0.0 && df <= ff)
        {
        nx_ = jNorm[v]->x();
        ny_ = jNorm[v]->y();
        tx_ = -ny_;
        ty_ = nx_;
        x_ = jVertex[v]->x() + cx - 0.5 * (dc + i_->R()) * nx_;
        y_ = jVertex[v]->y() + cy - 0.5 * (dc + i_->R()) * ny_;
        rank_ = 1;
		//rankco_=1;
        return;
        }
      }
    }
  
  rank_ = 0;
overlap_=-1;
}

void dkpg::Frame2()
{
  double cx,cy,cc;
  double dc,df,ff,dn;
  double xface,yface;
  unsigned int nbVertexes = j_->Vertex().size();
  unsigned int v,vPrev,vNext;
  
  // Global vertex coordinates and face normals for polygon j
  gdm::vertex * jVertex[nbVertexes];
  gdm::vertex * jNorm[nbVertexes];
  
  double c = cos(j_->rot());
  double s = sin(j_->rot());
  for(v = 0; v < nbVertexes ; ++v)
    {
    // Vertexes
    jVertex[v]->x() = j_->x() + c * j_->Vertex(v).x() - s * j_->Vertex(v).y();
    jVertex[v]->y() = j_->y() + s * j_->Vertex(v).x() + c * j_->Vertex(v).y();
    
    // Normales
    jNorm[v]->x() = c * j_->Normal(v).x() - s * j_->Normal(v).y();
    jNorm[v]->y() = s * j_->Normal(v).x() + c * j_->Normal(v).y();
    } 
  
  // We first check any contact with the polygon vertexes   
  for (v = 0; v < nbVertexes ; ++v)
    {
    cx = i_->x() - jVertex[v]->x();
    cy = i_->y() - jVertex[v]->y();
    cc = cx * cx + cy * cy;
    
    if (cc < i_->R()*i_->R())
      {
      vPrev = (v == 0) ? (nbVertexes - 1) : (v - 1);
      if (jNorm[vPrev]->x() * cy - jNorm[vPrev]->y() * cx < 0.0) break;
      if (jNorm[v]->x() * cy - jNorm[v]->y() * cx > 0.0) break;
      
      cc = sqrt(cc);
      dn = i_->R() - cc;
      cc = 1.0 / cc;
      nx_ = cx * cc;
      ny_ = cy * cc;
      tx_ = -ny_;
      ty_ = nx_;
      x_ = jVertex[v]->x() - 0.5 * dn * nx_;
      y_ = jVertex[v]->y() - 0.5 * dn * ny_;
	  overlap_= dn;
      	rank_ = 1;
		//rankco_=1;
      return;
      }
    
    }
  
  // There is no vertex-contact
  // So, we check contact with the polygon faces 
  for (v = 0;v<nbVertexes;++v)
    {
    cx = i_->x() - jVertex[v]->x();
    cy = i_->y() - jVertex[v]->y();
    dc = cx * jNorm[v]->x() + cy * jNorm[v]->y();
    
    if (dc > 0.0 && dc < i_->R())
      {
      vNext = (v == nbVertexes - 1) ? 0 : (v + 1);
      xface = jVertex[vNext]->x() - jVertex[v]->x();
      yface = jVertex[vNext]->y() - jVertex[v]->y();
      ff = yface * jNorm[v]->x() - xface * jNorm[v]->y();
      df = cy * jNorm[v]->x() - cx * jNorm[v]->y();
      
      if(df >= 0.0 && df <= ff)
        {
        nx_ = jNorm[v]->x();
        ny_ = jNorm[v]->y();
        tx_ = -ny_;
        ty_ = nx_;
        x_ = jVertex[v]->x() + cx - 0.5 * (dc + i_->R()) * nx_;
        y_ = jVertex[v]->y() + cy - 0.5 * (dc + i_->R()) * ny_;
        rank_ = 1;
		//rankco_=1;
        return;
        }
      }
    }
  
  rank_ = 0;
overlap_=-1;
}

void dkpg::Kin()
{
  double cix = x_ - i_->x();
  double ciy = y_ - i_->y(); 
  double cjx = x_ - j_->x();
  double cjy = y_ - j_->y();
  
  double vx = i_->vx() - j_->vx()
    -ciy * i_->vrot() + cjy * j_->vrot();
  double vy = i_->vy() - j_->vy()
    +cix * i_->vrot() - cjx * j_->vrot();
  
  vn_ = vx * nx_ + vy * ny_;
  vt_ = vx * tx_ + vy * ty_;
}

void dkpg::Res()
{
  double fx  = fn_ * nx_ + ft_ * tx_;
  double fy  = fn_ * ny_ + ft_ * ty_;
  double cix = x_ - i_->x();
  double ciy = y_ - i_->y();
  double cjx = x_ - j_->x();
  double cjy = y_ - j_->y();
  
  i_->fx()   += fx;
  i_->fy()   += fy;
  i_->frot() += (cix * fy - ciy * fx);  
  
  j_->fx()   -= fx;
  j_->fy()   -= fy;
  j_->frot() += (-cjx * fy + cjy * fx);
}

void dkpg::Res(const double dfn, const double dft, const double dfs)
{
  double dfx = dfn * nx_ + dft * tx_;
  double dfy = dfn * ny_ + dft * ty_;
  double cix = x_ - i_->x();
  double ciy = y_ - i_->y();
  double cjx = x_ - j_->x();
  double cjy = y_ - j_->y();
  
  i_->fx()   += dfx;
  i_->fy()   += dfy;
  i_->frot() += (cix * dfy - ciy * dfx);  
  
  j_->fx()   -= dfx;
  j_->fy()   -= dfy;
  j_->frot() += (-cjx * dfy + cjy * dfx);
}

void dkpg::CDcoeff() // a verifier
{
  double cix = x_ - i_->x();
  double ciy = y_ - i_->y();
  double cit = cix * tx_ + ciy * ty_;
  
  double cjx = x_ - j_->x();
  double cjy = y_ - j_->y();
  double cjt = cjx * tx_ + cjy * ty_;
  
  // temporaire
  double en = 0.0;
  //double et = 0.0;
  
  double mn = 1.0/(1.0/i_->mass()+1.0/j_->mass()
                   + (cit*cit)/i_->mom() + (cjt*cjt)/j_->mom());
  facn0_ = (1.0 + en)*mn;
  facn1_ = mn/i_->mass();
  facn2_ = mn/j_->mass();
  //facn3_ = mn*cit/i_->mom();
  facn4_ = mn*cjt/j_->mom();
}

double dkpg::An(const double dt)
{
  return (
          - facn0_ * vn_/dt + fn_
          - facn1_ * (i_->fx()*nx_ + i_->fy()*ny_)
          + facn2_ * (j_->fx()*nx_ + j_->fy()*ny_)
          //+ facn3_ * i_->frot()
          - facn4_ * j_->frot() 
          );  
}

double dkpg::At(const double dt)
{
  return 0;
}

double dkpg::As(const double dt)
{
  return (
          - facs0_ * vrot_/dt + frot_
          - facs1_ * i_->frot()
          + facs2_ * j_->frot()
          );
}


void dkpg::writeMGP (ostream & os)
{
}
