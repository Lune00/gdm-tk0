#include <iostream> 
#include <string>
#include <fstream>
#include "io.hpp"
#include "alteration.hpp"

using namespace std;

int main(int argc,char **argv)
{	
	if( argv[1]==0)
	
	{
	Sample * spl;
	char SampleFile[100];
	double dR,aspectRatio;
	unsigned int nVertexMin,nVertexMax;
	string token;

	cout << "Name of the Sample file: ";
	cin >> SampleFile;
	spl = new Sample;
	sample_read(SampleFile,*spl);

	cout << "Regular polygon (yes/no) ? ";
	cin >>token;
	if ( token == "yes")
	{
		cout << "Number of vertexes for regular polygon : ";
		cin >>nVertexMin;
		insert_regular_polyg_in_disk(*spl,nVertexMin);
	}
	else
	{
	cout << "Minimum number of vertexes: ";
	cin >> nVertexMin;
	cout << "Maximum number of vertexes: ";
	cin >> nVertexMax;
	cout << "Aspect ratio: ";
	cin >> aspectRatio;
	cout << "Variation of 'radius': ";
	cin >> dR;
	insert_polyg_in_disk(*spl, nVertexMin, nVertexMax, dR, aspectRatio);
	
	}
	sample_write("polygSample.spl",*spl);

	delete spl;
	}
	else if ( argv[1]!=0)
	{
		string in=argv[1];
		if ( in=="CEGEO")
		{
			Sample * spl;
			char SampleFile[100];
			double eta;
			unsigned int nVertexMin,nVertexMax;
			
			cout << "Name of the Sample file: ";
			cin >> SampleFile;
			spl = new Sample;
			sample_read(SampleFile,*spl);
			
			cout << "Minimum number of vertexes: ";
			cin >> nVertexMin;
			cout << "Maximum number of vertexes: ";
			cin >> nVertexMax;
			cout << "Variation of 'radius': ";
			cin >> eta;

			insert_polyg_in_disk_CEGEO(*spl, nVertexMin, nVertexMax,eta);
			
			sample_write("polygSample.spl",*spl);

			delete spl;	
		}
		else
			cout<<" Unknown argument --> end"<<endl;
	}

	return 0;
}
