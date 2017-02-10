#ifndef _alteration_hpp
#define _alteration_hpp

#include "sample.hpp"
#include "talk.hpp"
#include "dof.hpp"
#include "NRsource.hpp"
#include <math.h>

using namespace std;

void limit_speed(Sample& spl, double vlim);
void homothetie(Sample& spl, double h, double xc = 0.0, double yc = 0.0);
void insert_polyg_in_disk(Sample& spl, unsigned int nVertexMin, unsigned int nVertexMax,double dR, double aspectRatio);
void insert_regular_polyg_in_disk(Sample& spl, unsigned int nVertex);

void insert_polyg_in_disk_CEGEO(Sample& spl, unsigned int nVertexMin, unsigned int nVertexMax,double eta);

void mutate_polyg(body2d * B, unsigned int nVertexMin, unsigned int nVertexMax,double Rmax, double dR, double aspectRatio,bool chord_check=false);
void mutate_regular_polyg(body2d *B, unsigned int nVertex, double Rmax);

vector <dof*> insert_cluster_in_disk(Sample& spl, unsigned int nDisk,double ratio);
dof* mutate_cluster(body2d * B,Sample & spl, unsigned int nDisk,double ratio);
dof* mutate_cluster_baptiste(body2d * B,Sample & spl, unsigned int nDisk,double dr_r);



#endif //_alteration_hpp
