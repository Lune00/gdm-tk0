#include <iostream> 
#include <string>
#include <fstream>
#include "io.hpp"
#include "grid.hpp"

using namespace std;

int main(int argc,char **argv)
{
	char SampleFile[100];
	Sample * spl;
	double el;
	
	//Combien de fichier

	cout << "Name of the Sample file: ";
	cin >> SampleFile;
	cout << "Aspect ratio ? "<<endl;
	cin >> el;
	spl = new Sample;
	sample_read(SampleFile,*spl);
	
	squareGrid( *spl,el);
	
	sample_write("gridSample.spl",*spl);
	
	return 0;
}
