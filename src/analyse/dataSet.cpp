#include "dataSet.hpp"

void DataSet::write(char * file)
{
	ofstream out(file,ios::out);
	for( unsigned int i=0; i< set_.size();++i)
	{
		out<<set_[i]<<endl;
	}
	out.close();
}

DataSet::DataSet()
{
	this->Sort();
}

DataSet::DataSet(const char * name)
{
	this->Sort();
	name_ = name;
}

DataSet::DataSet(string name)
{
	this->Sort();
	name_ = name;
}

DataSet::~DataSet()
{
}

void DataSet::add(double value)
{
	this->set_.push_back(value);
}

void DataSet::write(ostream& os, unsigned int n)
{
	os << set_[n];
}

void DataSet::read(istream& is)
{
	set_.clear();
	double value;
	while( is )
	{
		is >> value;
		set_.push_back(value);
		if( is.eof()) break;
	}

}
void DataSet::addRead(istream& is)
{
	double value;
	while( is )
	{
		is >> value;
		set_.push_back(value);
		if( is.eof()) break;
	}
}

double DataSet::Mean()
{
	unsigned int N = set_.size();
	if (N == 0) 
	{	mean_=0.;
		return 0.0;
	}

	double mean = 0.0;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		mean += set_[i];
	}

	mean /= (double)N;
	mean_=mean;
	return mean;
}

/*double DataSet::Mean(unsigned int N)
{
	if (N == 0) 
	{	meandev_=0.;
		return 0.0;
	}

	double mean = 0.0;

	for(unsigned int i=0 ; i<set_.size() ; ++i)
	{
		mean += set_[i];
	}

	mean /= (double)N;
	meandev_=mean;
	return mean;

}*/

void DataSet::meanCentering()
{
	unsigned int N = set_.size();
	if (N == 0) 
	{	
		return ;
	}

	double mean = this->mean();

	for(unsigned int i=0 ; i<N ; ++i)
	{
	 	set_[i]-= mean ;
	}


	return ;
}

double DataSet::Min()
{
	unsigned int N = set_.size();
	if (N == 0)
		{return 0.0;min_=0;}

	double min = set_[0];

	for(unsigned int i=1 ; i<N ; ++i)
	{
		if(set_[i] < min) min = set_[i];
	}
	min_=min;
	return min;
}

double DataSet::Max()
{
	unsigned int N = set_.size();
	if (N == 0) {return 0.0;max_=0;}

	double max = set_[0];

	for(unsigned int i=1 ; i<N ; ++i)
	{
		if(set_[i] > max) max = set_[i];
	}
	max_=max;
	return max;
}

double DataSet::Variance()
{
	unsigned int N = set_.size();
	if (N == 0) {return 0.0;variance_=0;}

	double variance = 0.0;
	double s=0;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		s += set_[i]*set_[i];
	}

	variance = s/(double) (N) - mean_*mean_;


	variance=sqrt(variance);
	variance_= variance;//variance= srqrt ( < x*x > - <x>*<x> )
	return variance;
}
void  DataSet::extractValues()
{
	this->Max();
	this->Min();
	this->Mean();
	this->Variance();
	//this->Mean(unsigned int N);
}

void DataSet::Normalize(double ref)
{
	unsigned int N = set_.size();
	if (N == 0) return;

	for(unsigned int i=0 ; i<N ; ++i)
	{
		set_[i] /= ref;
	}
}

void DataSet::Logarithmize()
{
	unsigned int N = set_.size();
	if (N == 0) return;

	if ( set_[0] < 0.)
	{
		cout<< "negative values : can't logarithmize..."<<endl;
		return;
	}
	else
		for(unsigned int i=0 ; i<N ; ++i)
	{
		set_[i]=log (set_[i]);
	}

}

void DataSet::Sort()
{
	stable_sort( set_.begin(), set_.end() );
}

void DataSet::DecreasingSort()
{
	stable_sort( set_.begin(), set_.end(),DecreasingOrder());
}


pointSet DataSet::Rich_PDF( unsigned int nbin)// nbin classes of same number of event
{
	vector <double> pdf;
	vector <double> sigma;

	unsigned int N = set_.size();
	unsigned int Nmax= N/ nbin;
	unsigned int rest = N- Nmax*nbin;
	vector <double> borne;
	unsigned int Nev;
	unsigned int i=0;
	cout<<" number of events per classes = "<<Nmax<<endl;

	//this->DecreasingSort(); 

	borne.push_back(set_[0]);

	while (  i < N )
	{
		Nev=0;
		DataSet * data = new DataSet;
		while ( Nev<Nmax && i < N )
		{
			data->add(set_[i]);
			++i;++Nev;
		} 
		borne.push_back( set_[ i-1 ] );

		data->extractValues();
		sigma.push_back(data->variance());
		delete data;
		//cout <<" borne.size "<<borne.size()<<endl;
		//Nborne++;
		//cout <<"classe "<<Nborne-1<<": "<< borne[Nborne-1]<<" "<<borne[Nborne]<<" amp "<<borne[Nborne-1]-borne[Nborne]<<" "<<Nev<<endl;
		//cout<<borne[j]<<endl;
	}

	//Fusion de la derniere classe qui contient Nmax+rest evenements
	if (rest != 0 )
	{
		//cout <<" borne.size "<<borne.size()<<endl;
		borne[borne.size()-2]=borne[borne.size()-1];

		//borne.pop_back();
		borne.pop_back();
		//cout <<" borne.size "<<borne.size()<<endl;

		/*for( unsigned int j=0;j <borne.size()-1;++j)
		{
			cout <<"classe "<<j<<": "<< borne[j]<<" "<<borne[j+1]<<" amp "<<borne[j]-borne[j+1]<<endl;
		//cout<<borne[j]<<endl;
		}
		*/
	}
	cout<<endl;
	//correction
	double ampref=borne[borne.size()-2]-borne[borne.size()-1];
	double temp=0;

	//cout <<" ampref "<<ampref<<endl;;
	for( unsigned int j=0;j<borne.size()-1;++j)	
	{
		if( j==borne.size()-1 ) temp=rest;
		pdf.push_back( (double) (Nmax+temp)* ampref/ (borne[j]-borne[j+1]) );
		//tot+=pdf.back();
		//cout <<"classe "<<j<<": "<< borne[j]<<" "<<borne[j+1]<<" amp "<<borne[j]-borne[j+1]<<endl;
	}
	//pdf.push_back( (double) (Nmax+rest)* ampref/ (borne[borne.size()-2]-borne[borne.size()-1]) );
	//tot+=pdf.back();


	pointSet pdf2;
	for( unsigned int j=0;j< borne.size()-1;++j)	
	{
		//pdf[j]/=tot;
		//in<<.5*(borne[j]+borne[j+1])<<" "<<pdf[j]<<" "<<sigma[j]<<endl;
		pdf2.add(.5*(borne[j]+borne[j+1]),pdf[j] );

		//cout <<" j "<<j<<" pdf "<<pdf[j]<<" "<<borne[j]<<" "<<borne[j+1]<<endl;
	}
	double area = pdf2.integral();
	pdf2.yNormalise(area);
	//cout<<"	pdf "<<endl;
	return pdf2;


//return( pdf );
}

pointSet DataSet::Rich_PDF2( unsigned int Nevpc)//Nevpc events per class
{
	vector <double> pdf;
	vector <double> sigma;

	unsigned int N = set_.size();
	//unsigned int Nmax= N/ nbin;
	unsigned int rest=0;//= N- Nmax*nbin;
	vector <double> borne;
	unsigned int Nev;
	unsigned int i=0;

	//Nevpc= N/Nc;
	//cout<<" Nevpc = "<<Nevpc<<" ";

	this->DecreasingSort(); 

	borne.push_back(set_[0]);

	while (  i < N )
	{
		Nev=0;
		DataSet * data = new DataSet;
		while ( Nev<Nevpc && i < N )
		{
			data->add(set_[i]);
			++i;++Nev;
		} 
		borne.push_back( set_[ i-1 ] );

		data->extractValues();
		sigma.push_back(data->variance());
		delete data;

		if ( i==N) rest= Nev;
		//cout <<" borne.size "<<borne.size()<<endl;
		//Nborne++;
		//cout <<"classe "<<Nborne-1<<": "<< borne[Nborne-1]<<" "<<borne[Nborne]<<" amp "<<borne[Nborne-1]-borne[Nborne]<<" "<<Nev<<endl;
		//cout<<borne[j]<<endl;
	}

	//Fusion de la derniere classe qui contient Nmax+rest evenements
	if (rest != 0 )
	{
		//cout <<" borne.size "<<borne.size()<<endl;
		borne[borne.size()-2]=borne[borne.size()-1];

		//borne.pop_back();
		borne.pop_back();
		//cout <<" borne.size "<<borne.size()<<endl;

		/*for( unsigned int j=0;j <borne.size()-1;++j)
		{
			cout <<"classe "<<j<<": "<< borne[j]<<" "<<borne[j+1]<<" amp "<<borne[j]-borne[j+1]<<endl;
		//cout<<borne[j]<<endl;
		}
		*/
	}
	//correction
	double ampref=borne[borne.size()-2]-borne[borne.size()-1];

	double temp=0;

	//cout <<" ampref "<<ampref<<endl;;
	for( unsigned int j=0;j<borne.size()-1;++j)	
	{
		if( j==borne.size()-1 ) temp=rest;
		pdf.push_back( (double) (Nevpc+temp)* ampref/ (borne[j]-borne[j+1]) );
		//tot+=pdf.back();
		//cout <<"classe "<<j<<": "<< borne[j]<<" "<<borne[j+1]<<" amp "<<borne[j]-borne[j+1]<<endl;
	}
	//pdf.push_back( (double) (Nmax+rest)* ampref/ (borne[borne.size()-2]-borne[borne.size()-1]) );
	//tot+=pdf.back();


	pointSet pdf2;
	for( unsigned int j=0;j< borne.size()-1;++j)	
	{
		//pdf[j]/=tot;
		//in<<.5*(borne[j]+borne[j+1])<<" "<<pdf[j]<<" "<<sigma[j]<<endl;

		pdf2.add(.5*(borne[j]+borne[j+1]),pdf[j] );

		//cout <<" j "<<j<<" pdf "<<pdf[j]<<" "<<borne[j]<<" "<<borne[j+1]<<endl;
	}
	double area = pdf2.integral();
	pdf2.yNormalise(area);
	//cout<<"	pdf "<<endl;
	return pdf2;


//return( pdf );
}

pointSet DataSet::slidingPdf(unsigned int Nc)
{	
	unsigned int N = set_.size();
	//unsigned int Nmax= N/ nbin;
	unsigned int Nev,Nevpc,Nevpc2;
	unsigned int i=0;
	pointSet pdf;
	Nevpc= N/Nc;
	Nevpc2= Nevpc/4;
	DataSet data; 
//	cout<<" Nevpc = "<<Nevpc<<" "<<Nevpc2<<endl;
	unsigned int rest= N- Nevpc*Nc;
	//cout<<" rest "<<rest<<endl;
	//Premiere classe
	double xprec,xact=0;
	xprec= set_[0];
	pointSet var;
	while (  i < N  )
	{
		Nev=0;
		//DataSet * data = new DataSet;
		//for( unsigned int j = i - Nevpc ; j < i+ Nevpc;++j)
		while ( Nev<Nevpc && i < N )
		{
			data.add(set_[i]);
			xact=set_[i];
			++i;++Nev;
			if( i+rest == N) 
			{
			//cout<<"alerter"<<endl;
				xact=set_[N-1];
				Nev+=rest;
				i+=rest;
			}

		} 
		data.extractValues();
	//	cout<<"Nev "<<Nev<<" "<<fabs(xact-xprec)<<" "<<xact<<" "<<xprec<<endl;

		pdf.add( data.mean() , (double) (Nev) /fabs( xact-xprec));
		var.add( data.mean() , data.variance());
		xprec=xact;

		if (i<N) i-=Nevpc2;
		//cout<<" ---- i "<<i<<endl;

		// sigma.push_back(data->variance());
		data.setClear();

		//if ( i==N) rest= Nev;

	}

	double area = pdf.integral();
	pdf.yNormalise(area);

	/*	ofstream out("pdfslidvar.txt",ios::out);
	cout<<"------------------*********************************************************"<<endl;
	for( unsigned int i=0;i<pdf.psetSize();++i)
	{
		out<<pdf.px(i)<<" "<<pdf.py(i)<<" "<<var.py(i)<<endl;
	}
	*/


	return pdf;

}

pointSet DataSet::kernelPdf(unsigned int Npoint, double h)
{
	//h est la proportion de la largeur totale 

	this->extractValues();

	unsigned int N = set_.size();
	double min = this->min();
	double max = this->max();
	double amp = (max-min)/ (double) Npoint;
	double width = (max-min)*h;
	double tot=0,temp;
	double widthtemp;
	if(max<=min) return 0;
	if ( width < amp ) cout<<" Warning : kernel pdf disjoint classes......"<<endl;
	pointSet pdf;
	pointSet var;
	DataSet data;
	double xi,sum=0.;
	for (unsigned int j=0; j< Npoint ; ++j)
	{
		xi= min + (j+.5)*amp;

		if( xi- width <= min ) 
		{
			widthtemp=xi-min;
			continue;
		}
		else widthtemp=width;

		if( xi+ width >= max )
		{
			widthtemp=max - xi;
			continue;
		}
		else widthtemp=width;

		sum=0.;
		data.setClear();
		double u;
		for ( unsigned int i=0; i<N; ++i)
		{
			u=(set_[i] - xi)/ (widthtemp) ;
			sum += triangular( u );
			if( fabs(u) <1.) data.add( set_[i]);
		}
		temp = 	sum/( widthtemp*N);
		tot += temp;
		data.extractValues();
		var.add( xi , data.variance() );
		pdf.add( xi , temp);
	}
	double area = pdf.integral();
	pdf.yNormalise( area);

	return pdf;
}

void DataSet::upFilter( double lim)
{
	unsigned int Nsuppr=0;
	vector<double>::iterator it;
	for(it=set_.begin(); it != set_.end(); )
	{
		if( *it > lim ) {set_.erase(it);Nsuppr++;}
		else {it++;}
	}
	cout<<" upFilter, number of erased data (> "<<lim<<") : "<<Nsuppr<<endl;
}

void DataSet::downFilter( double lim)
{
	unsigned int Nsuppr=0;
	vector<double>::iterator it;
	for(it=set_.begin(); it != set_.end(); )
	{
		if( *it < lim ) {set_.erase(it);Nsuppr++;}
		else {it++;}
	}
	cout<<" downFilter, number of erased data (< "<<lim<<") : "<<Nsuppr<<endl;

}


DataSet * DataSet::RuningAverage(DataSet& input,unsigned int nb)
{

	DataSet * newSet = new DataSet;
	vector <unsigned int> P(nb,0);

	input.extractValues();

//double amp= ( input.max()-input.min())/(double) (nb);




	return( newSet );
}

DataSet * DataSet::Derivative(DataSet&, unsigned int n, unsigned int nb)
{
// on fit nb points avec un polynome de degres n
	DataSet * newSet = new DataSet;


	return( newSet );
}

DataSet * DataSet::Period(unsigned int nb)
{
	DataSet * newSet = new DataSet;

	for(unsigned int i=0; i<set_.size() ; i += nb)
	{
		newSet->add(set_[i]);
	}

	return newSet;
}

double DataSet::linearCorrelation(DataSet * input)
{
//Les deux set a correler doivent faire la meme taille et chaque rang doit correpondre
	this->extractValues();
	input->extractValues();
	double r=0.;
	double sumx,sumy;
	double sumx2,sumy2;
	double sumxy;
	unsigned int n = set_.size();
	sumx=sumy=sumx2=sumy2=sumxy=0.;

	for (unsigned int i=0;i< n ;++i)
	{
		//x -> this
		//y -> input
		sumxy += this->set(i)*input->set(i);
		sumx  += this->set(i);
		sumy  += input->set(i);
		sumx2 += this->set(i)*this->set(i);
		sumy2 += input->set(i)*input->set(i);

	}

	r = (n*sumxy -sumx*sumy)/ (sqrt( n*sumx2 - sumx*sumx)*sqrt( n*sumy2 - sumy*sumy));
//	r/= sqrt( n*sumx2 - sumx*sumx);
//	r/= sqrt( n*sumy2 - sumy*sumy);
	//r/=sqrt( this->variance()*(this->setSize()-1)*input->variance()*(input->setSize()-1));


	return r;

}



