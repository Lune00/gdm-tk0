#include "NRsource.hpp"

# define IA 16807
# define IM 2147483647
# define AM (1.0/IM)
# define IQ 127773
# define IR 2836
# define NTAB 32
# define NDIV (1+(IM-1)/NTAB)
# define EPS 1.2e-7
# define RNMX (1-EPS)

#define MAXIT 100 
#define EPS2 3.0e-7 
#define FPMIN 1.0e-30 

using namespace std;
double ran1(long *idum)
{
	long j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	
	if (*idum<=0  || !iy)
	{
		if (-(*idum)<1)  *idum=1;
		else *idum=-(*idum);
		for (j=NTAB+7;j>=0;j--)
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum<0) *idum += IM;
			if (j<NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum<0) *idum+=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum;
	if ((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

double betai(double a, double b, double x) 
//Returns the incomplete beta function Ix(a,b). 
{ 
double betacf(double a, double b, double x); 
double gammln(double xx); 
 
double bt; 
if (x <0.0 || x >1.0) printf("Bad x in routine betai");; 
if (x ==0.0 ||x == 1.0) bt=0.0; 

else 
	bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x)); 

if (x <(a+1.0)/(a+b+2.0))  
	return bt*betacf(a,b,x)/a; 
else 
	return 1.0-bt*betacf(b,a,1.0-x)/b; 
} 
 

double betacf(double a, double b, double x) 
{ 


int m,m2; 
double aa,c,d,del,h,qab,qam,qap; 

qab=a+b;
qap=a+1.0; 
qam=a-1.0; 
c=1.0;  
d=1.0-qab*x/qap; 
if (fabs(d) < FPMIN) d=FPMIN; 
d=1.0/d; 
h=d; 
for (m=1;m<=MAXIT;m++) 
{ 
m2=2*m; 
aa=m*(b-m)*x/((qam+m2)*(a+m2)); 
d=1.0+aa*d;  
if (fabs(d) < FPMIN) d=FPMIN; 
c=1.0+aa/c; 
if (fabs(c) < FPMIN) c=FPMIN; 
d=1.0/d; 
h*= d*c; 
aa =-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); 
d=1.0+aa*d; 
if (fabs(d) < FPMIN) d=FPMIN; 
c=1.0+aa/c; 
if (fabs(c) < FPMIN) c=FPMIN; 
d=1.0/d; 
del=d*c; 
h*= del; 
if (fabs(del-1.0) <EPS2) break;// Are we done? 
} 
if (m >MAXIT) printf("a or b too big, or MAXIT too small in betacf"); 
return h; 
} 


double gammln(double xx) 
{ 
double x,tmp,ser; 
static double cof[6]={76.18009172947146,-86.50532032941677, 
24.01409824083091,-1.231739572450155, 
0.1208650973866179e-2,-0.5395239384953e-5}; 
int j; 
x=xx-1.0; 
tmp=x+5.5; 
tmp -= (x+0.5)*log(tmp); 
ser=1.0; 
for (j=0;j<=5;j++) 
	{
		x+=1.0;
		ser += cof[j]/x;
	}
return -tmp+log(2.50662827465*ser); 
} 


