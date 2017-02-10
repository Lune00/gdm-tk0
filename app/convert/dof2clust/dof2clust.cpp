#include <iostream> 
#include <string>
#include <fstream>
#include "io.hpp"
#include "dof.hpp"

using namespace std;

int main(int argc,char **argv)
{
	unsigned int first,last,nbfich,numerotation=0;
	first=1;
	char  includeFrom [100];
	
	cout<< " Conversion des fichiers spl avec .dof en cluster intergres "<<endl;
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
	cout << "Nom du fichier original ? ";
	cin >> includeFrom;
	nbfich = last + 1 - first;
	cout << nbfich << " fichiers a traiter " << endl;

	char  nomFichier [100];
	char  fichDof[100];
	char  comsys[100];

	// Chargement fichier
	Sample * spl;
	Network * nwk;
	GroupRelationData * grpRel;
	vector<unsigned int> d1(0),d2(0);
	unsigned int temp;
	sprintf(fichDof,"dof_0_%s.dof",includeFrom);
	ifstream dof1(fichDof);
	if ( ! dof1 )
	{
		cout<<" erreur ouverture dof1 "<<endl;
		exit(0); 
	}
	else
	{
		while (dof1)
		{
			dof1>>temp;
			d1.push_back(temp);
			cout<<temp<< ' ';
		}
		d1.pop_back();
	}
	cout<<endl;
	
	sprintf(fichDof,"dof_1_%s.dof",includeFrom);
	ifstream dof2(fichDof);
	if ( ! dof2 )
	{
		cout<<" erreur ouverture dof2 "<<endl;
		exit(0); 
	}
	else
	{
		while (dof2)
		{
			dof2>>temp;
			d2.push_back(temp);
			cout<<temp<< ' ';
			
		}
		d2.pop_back();
		
	}
	cout<<endl;
	
	dof * dtemp1, *dtemp2;
	
	spl = new Sample;
	nwk = new Network;
	dtemp1=new dof();
	dtemp2=new dof();
	cout << includeFrom << " (dof) -> (cluster) " << includeFrom << endl << flush;
	cout << includeFrom << "  ->  " << includeFrom <<"_svg "<< endl << endl << flush;
	history_read(includeFrom,*spl,*nwk);
	for (unsigned int ii=0;ii < d1.size();++ii)
	{
		//spl->body(d1[ii])->bodyDof()=dtemp1;
		dtemp1->plugBody(spl->body(d1[ii]));
	}
	for (unsigned int ii=0;ii < d2.size();++ii)
	{
		//spl->body(d2[ii])->bodyDof()=dtemp2;
		dtemp2->plugBody(spl->body(d2[ii]));
	}
	
	
	history_write(0,*spl,*nwk, *grpRel,false,true);
	sprintf(comsys,"mv %s %s_svg",includeFrom,includeFrom);
	system(comsys);
	sprintf(comsys,"mv spl_nwk_0000.his %s",includeFrom);
	system(comsys);
	
	delete spl;
	delete nwk;
	delete dtemp1;
	delete dtemp2;
	
	
	for( unsigned int i=first;i<=last;i++)
	{
		spl = new Sample;
		nwk = new Network;
		dtemp1=new dof();
		dtemp2=new dof();
		
		if( numerotation==3 )
		{
		sprintf(nomFichier,"spl_nwk_%.3d.his",i);
		
		}
		else 
		{
		sprintf(nomFichier,"spl_nwk_%.4d.his",i);
		}
		
		
		
		//sprintf(nomMgp,"mgp.out.%.3d",Nmgp);

		cout << nomFichier << " (dof) -> (cluster) " << nomFichier << endl << flush;
		cout << nomFichier << "  ->  " << nomFichier <<"_svg "<< endl << endl << flush;

		history_read(nomFichier,*spl,*nwk);
		
		sprintf(comsys,"mv %s %s_svg",nomFichier,nomFichier);
		system(comsys);
		
		for (unsigned int ii=0;ii < d1.size();++ii)
		{
			//spl->body(d1[ii])->bodyDof()=dtemp1;
			dtemp1->plugBody(spl->body(d1[ii]));
		}
		for (unsigned int ii=0;ii < d2.size();++ii)
		{
			//spl->body(d2[ii])->bodyDof()=dtemp2;
			dtemp2->plugBody(spl->body(d2[ii]));
		}
		cout<<" dof 1 = "<<dtemp1->lctrlBodies().size()<< " dof 2 = "<<dtemp2->lctrlBodies().size()<<endl;
		
		history_write(i,*spl,*nwk, *grpRel,false,true);
		
		
		delete spl;
		delete nwk;
		delete dtemp1;
		delete dtemp2;
	//	Nmgp++;

	}
	return 0;
}
