#include "pointSet.hpp"


void pointSet::circularise(double offset,int period)
{
	unsigned int size = pset_.size();

	for( unsigned int i=0; i< size; ++i)
	{
		dtk::point temp(pset_[i].x()+ offset,pset_[i].y());
		pset_.push_back(temp);
	}
/*	size = pset_.size()-1;


	dtk::point temp1(pset_[size].x(),pset_[size].y());
	pset_.insert(pset_.begin(),temp1);
	
	dtk::point temp(pset_[1].x(),pset_[1].y());
	pset_.push_back(temp);
	*/
	this->periodic_=true;
	/*
	for(  int i=0; i< period; ++i)
	{
		dtk::point temp(pset_[size-i].x(),pset_[size-i].y());
		pset_.insert(pset_.begin(),temp);
		size++;
	}

	for(  int i=0; i< period; ++i)
	{
		dtk::point temp(pset_[i+period].x(),pset_[i+period].y());
		pset_.push_back(temp);
	}
	*/
}

void pointSet::closeCircular()
{
	
	if(periodic_)
	{
	dtk::point temp(pset_[0].x(),pset_[0].y());
	pset_.push_back(temp);
	}
}

void pointSet::yAbs()
{
	for( unsigned int i=0; i< pset_.size() ;++i)
	{
		pset_[i].y() = fabs(pset_[i].y());
	}
}


void pointSet::upXFilter( double lim)
{
	unsigned int Nsuppr=0;
	vector<dtk::point>::iterator it;
	for(it = pset_.begin(); it != pset_.end(); )
	{
		if( (*it).x() > lim ) {pset_.erase(it);Nsuppr++;}
		else {it++;}
	}
	cout<<" upXFilter, number of erased data (> "<<lim<<") : "<<Nsuppr<<endl;
}

void pointSet::downXFilter( double lim)
{
	unsigned int Nsuppr=0;
	vector<dtk::point>::iterator it;
	for(it=pset_.begin(); it != pset_.end(); )
	{
		if( (*it).x() < lim ) {pset_.erase(it);Nsuppr++;}
		else {it++;}
	}
	cout<<" downXFilter, number of erased data (< "<<lim<<") : "<<Nsuppr<<endl;

}


void pointSet::add(double x, double y)
{	
	dtk::point temp(x,y);
//	cout<<" add : "<<temp.x()<<" "<<temp.y()<<endl;
	pset_.push_back(temp);
	//pset_.insert(pset_.end(), temp);
//	cout<<" add ok, taille "<<this->pset_.size()<<endl;
}

void pointSet::write( char * file)
{
	//cout<<" ecriture de : "<<file<<endl;
	ofstream out(file,ios::out);
	for( unsigned int i=0; i< pset_.size();++i)
	{
		out<<this->px(i)<<" "<<this->py(i)<<endl;
	}
	out.close();
}

void pointSet::read( istream & is)
{
	this->pset_.clear();
	double x,y;
	while( is )
	{
		is>> x >> y ;
		this->add(x,y);
		if( is.eof()) break;
	}
}

void pointSet::addRead( istream & is)
{
	double x,y;
	while( is )
	{
		is>> x >> y ;
		this->add(x,y);
		if( is.eof()) break;
	}
}

void pointSet::Mean()
{
	unsigned int N = pset_.size();
	if (N == 0) 
	{	
		xmean_=0.;
		ymean_=0.;
		return;
	}

	double xmean = 0.0;
	double ymean = 0.0;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		xmean += pset_[i].x();
		ymean += pset_[i].y();
	}

	xmean /= (double)N;
	ymean /= (double)N;
	xmean_=xmean;
	ymean_=ymean;

}
void pointSet::Variance()
{
	unsigned int N = pset_.size();
	if (N == 0) {xvar_=yvar_=0;return ;}

	double xvar = 0.0,yvar=0.0;
	double xs=0,ys=0;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		xs += pset_[i].x()*pset_[i].x();
		ys += pset_[i].y()*pset_[i].y();
	}

	xvar = xs/(double) (N) - xmean_*xmean_;
	yvar = ys/(double) (N) - ymean_*ymean_;


	xvar_=sqrt(xvar);
	yvar_=sqrt(yvar);
	//variance= srqrt ( < x*x > - <x>*<x> )
}
void pointSet::extractValues()
{
	Min();
	Max();
	Mean();
	Variance();
}
void pointSet::xNormalise( double ref)
{
	unsigned int N = pset_.size();
	if (N == 0) return;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		pset_[i].x() /= ref;
	}
}
void pointSet::yNormalise(double ref)
{
	unsigned int N = pset_.size();
	if (N == 0) return;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		pset_[i].y() /= ref;
	}
}
void pointSet::xoffSet( double ref)
{
	unsigned int N = pset_.size();
	if (N == 0) return;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		pset_[i].x() += ref;
	}
}
void pointSet::yoffSet(double ref)
{
	unsigned int N = pset_.size();
	if (N == 0) return;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		pset_[i].y() += ref;
	}
}
void pointSet::yCumul( )
{
	unsigned int N = pset_.size();
	if (N == 0) return;
	double cumul=0;
	for(unsigned int i=0 ; i<N ; ++i)
	{
		cumul+=pset_[i].y() ;
		pset_[i].y() = cumul;
	}
	ysum_=cumul;
}
void pointSet::Min()
{
	unsigned int N = pset_.size();
	if (N == 0)
		{xmin_=ymin_=0;return ;}

	double xmin = pset_[0].x();
	double ymin = pset_[0].y();

	for(unsigned int i=1 ; i<N ; ++i)
	{
		if(pset_[i].x() < xmin) xmin = pset_[i].x();
		if(pset_[i].y() < ymin) ymin = pset_[i].y();
	}
	xmin_=xmin;
	ymin_=ymin;
}

void pointSet::Max()
{
	unsigned int N = pset_.size();
	if (N == 0) {xmax_=ymax_=0;return ;}

	double ymax = pset_[0].y();
	double xmax = pset_[0].x();

	for(unsigned int i=1 ; i<N ; ++i)
	{
		if(pset_[i].x() > xmax) xmax = pset_[i].x();
		if(pset_[i].y() > ymax) ymax = pset_[i].y();

	}
	xmax_=xmax;
	ymax_=ymax;
}

void pointSet::decreasingXSort()
{
	stable_sort( pset_.begin(), pset_.end(),DecreasingXOrder());
}
void pointSet::decreasingYSort()
{
	stable_sort( pset_.begin(), pset_.end(),DecreasingYOrder());
}
void pointSet::increasingXSort()
{
	stable_sort( pset_.begin(), pset_.end(),IncreasingXOrder());
}
void pointSet::increasingYSort()
{
	stable_sort( pset_.begin(), pset_.end(),IncreasingYOrder());
}

//Moyenne mobile tout les "period" points avec une fenetre d'amplitude centree sur 
//le point considerer et de demi largeur "width"
pointSet pointSet::mobileMean(unsigned int period, unsigned int width)
{

	if (period == 0  )
	{
		cout<<" @ pointSet::mobileMean : bad parameters "<<endl;
		return 0;
	}
	unsigned int Np= pset_.size();
	double yval=0.;
	pointSet mob;
	if(!periodic_)
	{
		for( unsigned int i= 0; i<width ;++i)
		{
			mob.add( this->px(i), this->py(i));
		}
		
		for( unsigned int i = width ; i< Np-width ; i += period )
		{
			for ( unsigned int j = i-width; j <= i + width; ++j)
			{
				if( j <= Np)
				{
					yval+= this->py(j) ;
				}
				else
				{
					break;
				}
			}
			yval/=(double) (2*width + 1);
			mob.add( this->px(i) , yval );
			yval=0.;
		}

		for( unsigned int i= Np-width ; i<Np ;++i)
		{
			mob.add( this->px(i), this->py(i));
		}
	}
	else 
	{
		int amp=2*width+1;
		int width2 = (int) width;
		int Np2 = (int) Np;
		vector<unsigned int> rg(amp,0);
		//cout<<" periodic : "<<amp<<endl;
		
		for(  int i = 0 ; i< Np2 ; ++i )
		{
			for(int k = 0;k< amp ;++k)
			{
				if(i+k-width2<0) 
				{
					rg[k]= i-(width2)+k +Np2;
				}
				else if (i+k-width2 >Np2-1) 
				{
					
					rg[k]=i+k - width2-Np2;
				}
				else rg[k]=i+k-width2;
			}
			
			for (  int j = 0; j < amp; ++j)
			{
					yval+= this->py(rg[j]) ;
					//cout<<rg[j]<<" ";
			}
			yval/=(double) (amp);
			
			mob.add( this->px(i) , yval );
			yval=0.;
		}
		mob.periodic()=true;
		
	}

	return (mob);

}

// No periodic
pointSet pointSet::mobileMean2(unsigned int period, unsigned int width)
{

	if (period == 0  )
	{
		cout<<" @ pointSet::mobileMean : bad parameters "<<endl;
		return 0;
	}
	unsigned int Np= pset_.size();
	double yval=0.;
	double xval=0.;
	pointSet mob;

		for( unsigned int i= 0; i<width ;++i) 
		{
			xval+=this->px(i);
			yval+= this->py(i);
		}
		xval/=double(width);
		yval/=(double) (width);
		mob.add(xval, yval );
		xval=0.;
		yval=0.;
			
		for( unsigned int i = width ; i< Np-width ; i += period )
		{
			for ( unsigned int j = i-width; j <= i + width; ++j)
			{
				if( j <= Np)
				{
					yval+= this->py(j) ;
				}
				else
				{
					break;
				}
			}
			yval/=(double) (2*width + 1);
			mob.add( this->px(i) , yval );
			yval=0.;
		}

		for( unsigned int i= Np-width ; i<Np ;++i)
		{
			xval+=this->px(i);
			yval+= this->py(i);
		}
		xval/=double(width);
		yval/=(double) (width);
		mob.add(xval, yval );
		xval=0.;
		yval=0.;
	return (mob);
}

pointSet pointSet::histoBin( unsigned int Nbin  )//Nbin 
{
	this->Max();
	this->Min();
	double amp= (xmax_-xmin_) / (double) (Nbin);
	unsigned int N=this->psetSize();
	unsigned int rang=0;
	vector <unsigned int> Nev(Nbin,0);
	pointSet histo(Nbin);
	//cout<<" taille "<<histo.psetSize()<<endl;

	for (unsigned int i=0; i< Nbin; ++i)
	{
		histo.py(i) = 0.;
		histo.px(i) = xmin_ + ( i +.5)*amp;
	}

	for (unsigned int i=0; i< N; ++i)
	{
		rang = (unsigned int) floor ( (this->px(i)-xmin_)/amp);
		if( rang==Nbin ) rang--;
		histo.py(rang) += this->py(i);
		Nev[rang]++;
	}

	for (unsigned int i=0; i< Nbin; ++i)
	{
		if( Nev[i] != 0)
			histo.py(i) /= (double) (Nev[i]);
		//cout<<" i= "<<i<<" Nev "<<Nev[i]<<endl;
	}
	return histo;
}
pointSet pointSet::histoNumberBin( unsigned int Nbin  )//Nbin population not frequence
{
	this->Max();
	this->Min();
	double amp= (xmax_-xmin_) / (double) (Nbin);
	unsigned int N=this->psetSize();
	unsigned int rang=0;
	vector <unsigned int> Nev(Nbin,0);
	pointSet histo;

	for (unsigned int i=0; i< Nbin; ++i)
	{
		histo.add(0,0);
	}
	//cout<<" taille "<<histo.psetSize()<<endl;

	for (unsigned int i=0; i< Nbin; ++i)
	{
		histo.py(i) = 0.;
		histo.px(i) = xmin_ + ( i +.5)*amp;
	}

	for (unsigned int i=0; i< N; ++i)
	{
		rang = (unsigned int) floor ( (this->px(i)-xmin_)/amp);
		if( rang==Nbin ) rang--;
		histo.py(rang) += this->py(i);
		Nev[rang]++;
	}

	return histo;
}
pointSet   pointSet::histoBinVar( unsigned int Nbin)
{
	//The original set need to be sorted
	this->extractValues();
	unsigned int N=this->psetSize();
	unsigned int Nevpc = N/Nbin;
	unsigned int Current=0;
	unsigned int Nevbin;

	pointSet histo;

	for (unsigned int i=0; i< Nbin; ++i)
	{
		histo.add(0,0);
	}
	//cout<<"	histobinvar taille "<<N<<" Nevpc "<<Nevpc<<flush<<endl;

	//cout<<"avant while"<<endl;
	unsigned int i=0;
	Nevbin=0;
	while (  i< N)
	{
		if( Nevbin < Nevpc)
		{
			histo.py(Current)+= this->py(i);
			histo.px(Current)+= this->px(i);
			i++;
			Nevbin++;
		}
		else
		{
			Current++;
			//histo->px(Current)= this->px(i);
			Nevbin=0;
			//cout<<i<<" "<<histo.px( Current-1)<<" "<<histo.py( Current-1)<<endl;
		}
	}
	//cout<<" fin while"<<endl;
	for (unsigned int j=0; j< Nbin; ++j)
	{
		histo.py(j) /= (double) (Nevpc);
		histo.px(j) /= (double) (Nevpc);

	}
	//cout<<" bi var ok "<<endl;
	return  histo;

}

pointSet pointSet::slidingHisto( unsigned int Nbin, double h)
{
	//The original set need to be sorted
	this->extractValues();
	unsigned int N=this->psetSize();
	unsigned int Nev;
	double ymoy,xmoy,xi;
	if (xmax_==xmin_) return 0;
	double amp = (xmax_ - xmin_)/(double) Nbin;
	double width = (xmax_ - xmin_)*h;
	double width2 = .5*width;
	double hinit=h;

	pointSet histo;

	double dmoy=N/(xmax_-xmin_);
	double d;
	double dlim= dmoy*.05;
	unsigned int weight,sumw;
	xmoy=ymoy=sumw=0;
	//cout<<"	slidinghisto : "<<amp<<" "<<width<<endl;
	
	for (unsigned int i=0; i< Nbin; ++i)
	{

		xi = xmin_ + ( i+.5)*amp;

		//if( i>0 && histo.px(i-1) + width2 > xmax_ ) break;

		d=dlim;
		h=hinit;
		
		while( d <= dlim )
		{
			width=(xmax_ - xmin_)*h;
			width2=.5*width;
			ymoy=xmoy=0.;
			Nev=sumw=0;
			for ( unsigned int j=0;j<N;++j)
			{
				if( fabs(xi - this->px(j)) < width2 )
				{
					weight=1;//Weight( fabs( xi - this->px(j))/ width2);
					sumw+=weight;
					Nev++;
					ymoy+=weight*this->py(j);
					xmoy+=weight*this->px(j);
				}
			}

		/*	if( xi- width2 < xmin_ ) 
			{	
				d=(double) Nev/(xi-xmin_+width2);
				cout<<"depasse min "<<endl;
			}
			else if ( xi + width2 > xmax_ )
			{	
				d=(double) Nev/(xmax_-xi+width2);
				cout<<" depasse max"<<endl;
			}
			else
				{	d=(double) Nev/(width);}*/
				d=dlim+.1;		
		//if( d< dlim ) {h+=.02;cout<<"increment taille"<<endl;}
		}
		//cout<<Nev<<" "<<dmoy<<" "<< (double) Nev/width<<endl;
		
		if((double) Nev/width< dmoy ) h*=2.;

		ymoy/=(double) sumw;//Nev;
		xmoy/=(double) sumw;//Nev;

		histo.add(xmoy,ymoy);

		xmoy=ymoy=0.;
	}
	histo.increasingXSort();
	pointSet histo2=histo.mobileMean(1,4);
	return histo2;
}

double pointSet::integral ( )
{
	unsigned int N=this->pset_.size();
	double area=0.;

	for (unsigned int i=0;i<N-1;++i)
	{
		area+= .5* (max(this->px(i+1), this->px(i))-min(this->px(i+1), this->px(i)) )* ( this->py(i+1)+ this->py(i));
	}

	return( area);
}



