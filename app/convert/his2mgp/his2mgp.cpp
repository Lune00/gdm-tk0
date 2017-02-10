#include <iostream> 
#include <string>
#include <fstream>
#include "io.hpp"

using namespace std;

int main(int argc,char **argv)
{
	unsigned int first,last,nbfich,numerotation=0,per=1;
	double dt=10e-4;
	first=1;
	unsigned int Nmgp=0;
	
	//Combien de fichier
	while (!(numerotation == 3 || numerotation == 4))
	{
	cout << "Type de numerotation (3 ou 4 chiffres ) ? ";
	cin >> numerotation;
	}
	cout << "Numero du premier fichier ? ";
	cin >> first;
	cout << "Numero du dernier fichier ? ";
	cin >> last;
	cout << "Periode de traitement ? ";
	cin >> per;
	cout << "Debut de numerotation ? ";
	cin >> Nmgp;
	nbfich = last + 1 - first;
	cout << nbfich << " fichiers a traiter " << endl;

	char  nomFichier [100],nomMgp [100];

	// Chargement fichier
	Sample * spl;
	Network * nwk;
	for( unsigned int i=first;i<=last;i+=per)
	{
		spl = new Sample;
		nwk = new Network;
		if( numerotation==3 )
		sprintf(nomFichier,"spl_nwk_%.3d.his",i);
		else 
		sprintf(nomFichier,"spl_nwk_%.4d.his",i);
		
		
		sprintf(nomMgp,"mgp.out.%.3d",Nmgp);

		cout << nomFichier << " -> " << nomMgp << endl << flush;

		history_read(nomFichier,*spl,*nwk);
		write_mgpost(nomMgp,*spl,*nwk,Nmgp,Nmgp*dt);
		delete spl;
		delete nwk;
		Nmgp++;

	}
	return 0;
}
