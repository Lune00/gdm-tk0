#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"

//! \brief Fusion between two run of the same simulation ( number of file and time file) 
//! \author C. Voivret

//fusionne deux run de la meme simulation: le dernier fichier du premier run sert de depart pour le second.
//run1 : dernier fichier à tf1
//run1 : premier fichier à tf1+k*dt

//Verif : maj liste fichier pour time mean
 
int main (int argc, char * argv[]) 
{
	
	string dir1,dir2,time1,time2;
	vector <double> t1,t2;
	double temp1,temp2;
	char  commande [100];
	cout<<" Les deux series a concatener doivent etre dans des dossier separes - au meme niveau - "<<endl;
	cout<<" Premier dossier "<<endl;
	cin>>dir1;
	cout<<" Premier fichier time "<<endl;
	cin>>time1;
	cout<<" Deuxieme dossier "<<endl;
	cin>>dir2;
	cout<<" Deuxieme fichier time "<<endl;
	cin>>time2;
	
	// lecture des fichiers time
	sprintf(commande,"%s/%s",dir1.c_str(),time1.c_str());
	
	ifstream time(commande);
	if(!time)
	{
		cerr << "Fichier de temps non trouve  " << time << endl;
		//return 0;
	}
	else
	{
		
		while(time )
		{
			
			time>> temp1 >>temp2;
			//cout<<temp2<<" ";
			t1.push_back(temp2);
			if( time.eof()) break;
			
		}
	}
	time.close();
	
	sprintf(commande,"%s/%s",dir2.c_str(),time2.c_str());
	
	ifstream timebis(commande);
	if(!timebis)
	{
		cerr << "Fichier de temps non trouve  " << timebis << endl;
		//return 0;
	}
	else
	{
		
		while(timebis )
		{
			if( timebis.eof()) break;
			timebis>> temp1 >>temp2;
			//cout<<temp2<<" ";
			t2.push_back(temp2);
			
		}
	}
	timebis.close();
	
	t1.pop_back();
	t2.pop_back();
	
/*	for (unsigned int i=0;i<t1.size();++i)
	{
		cout<<i<<" "<<t1[i]<<" "<<t2[i]<<endl;
	}
*/
	//Creation du dossier contenant la fusion
	char path[100];
	
	sprintf(path,"Fusion_%s_%s",dir1.c_str(),dir2.c_str());
	
	sprintf(commande,"mkdir %s",path);
	system( commande);
	
	sprintf(commande,"cp %s/spl_nwk_*  %s/",dir1.c_str(),path);
	system( commande);
	
	sprintf(commande,"cp %s/%s  %s/time_fus.txt",dir1.c_str(),time1.c_str(),path);
	system( commande);
	
	sprintf(commande,"%s/time_fus.txt",path);
//	system( commande);
	
	ofstream timeadd(commande,ios::app);
	unsigned int ndeb = t1.size();
	double tdeb = t1.back();
	char fich1[100],fich2[100];
		
	for( unsigned int i =0; i< t2.size();++i)
	{
		cout<<tdeb<<"  +  "<<t2[i]<<" = "<<tdeb+t2[i]<<endl;
		timeadd<<ndeb+i<<" "<<tdeb+t2[i]<<endl;
		sprintf(fich1,"%s/spl_nwk_%.4d.his",dir2.c_str(),i);
		sprintf(fich2,"%s/spl_nwk_%.4d.his",path,i+ndeb);
		cout<<fich1<<" "<<fich2<<endl;
		sprintf(commande,"cp %s %s",fich1,fich2);
		system( commande);
		
		
	}
	
	
	return 0;
}
