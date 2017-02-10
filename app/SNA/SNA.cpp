#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"

//! \brief Application to sequantially reload .his file and analyse 
//! \author C. Voivret


// Deux choix : 
// - monofichier : argument 1.: fichier .his( ou similaire)  argument 2.: fichier de commande d'analyse
// - multifichier: saisie au clavier de fichier spl_nwk_xxx.his ou spl_nwk_xxxx.his 
int main (int argc, char * argv[]) 
{
	
	string token,splnwk;
	
	if ( argv[1]!=0)
	{
		
		Simulation * mySimu= new Simulation();
		mySimu->read_data(argv[1]);
		
		ifstream input(argv[1]);
		if( ! input)
		{
			cerr << "Echec overture du fichier sim " << argv[1] << endl;
			
		}
		else
		{
			input >> token;
			while(input)
			{
				if (token == "includeFile")  
				{
				input >> splnwk;
				//cout<<"include trouve"<<splnwk<<endl;
				
				}
				input>>token;
			}
		}
		cout<<" splnwk extrait"<<endl;
			
		//Lecture du fichier de commande contenant les grandeurs a evaluer
		ifstream data2(argv[2]);
		if(!data2)
		{
			cerr << "Fichier de commande manquant : Analyse a partir du fichier sim " << argv[1] << endl;
			//return 0;
		}
		else
		{
			cout<<" Lecture du fichier de commande "<<endl;
			mySimu->sysA()->allFalse();
			data2 >> token;
			while(data2)
			{
				if (token == "System_A{")  
					mySimu->sysA()->read_parameters(data2);
				data2 >> token;
			}
		}
		cout<<" lecture faite "<<endl;
		
		mySimu->load_history(splnwk.c_str());
		mySimu->algo()->algoFill();
				
		cout<<" taille spl = "<<mySimu->spl()->lbody().size()<<endl;
		cout<<" taille nwk = "<<mySimu->nwk()->linter().size()<<endl;
		cout<<" nb contacts= "<<mySimu->nwk()->clist().size()<<endl;
		
		mySimu->sys()->init();
		
		mySimu->sysA()->initAnalyse();
		
		mySimu->sysA()->analyse(10.,1,1);
		
		write_mgpost("mgp.out.001",*mySimu->spl(),*mySimu->nwk(),0,0);
	//	mySimu->sysA()->writePS("visuSNA.ps");
	}
	else
	{
		
		//Premier et dernier fichier
		unsigned int Ndeb,Nfin,period,format=1;
		string fichcom,fichsim,fichtemps,include,token;
		char  nomFichier [100],commande[100],princDir[100];	
		
		cout<<"Numero du premier fichier a traiter ?"<<endl;
		cin>>Ndeb;
		cout<<"Numero du dernier fichier a traiter ?"<<endl;
		cin>>Nfin;
		cout<<"Format numerotation 3 ou 4 ? "<<endl;
		cin>>format;
		cout<<"Periode de traitement ? "<<endl;
		cin>>period;
		cout<<"Fichier sim ?"<<endl;
		cin>>fichsim;
		cout<<"Fichier de commande ?"<<endl;
		cin>>fichcom;
		cout<<"Fichier de temps ?"<<endl;
		cin>>fichtemps;
		cout<<"Nom du dossier principal a creer ?"<<endl;
		cin>>princDir;
		
		sprintf(commande,"mkdir %s",princDir);
		system( commande);
		/*cout<<"Moyenne temporelle a effectuer ? (y/n)"<<endl;
		cin>>token;
		if( token=="y")
		{
			//parametre a lire tdebut tfin
			//grandeur scalaire a evaluer : q,p,q/p,a,an,at,al,
			//grandeur vectorielle a evaluer : pdf, Xu, correlation force taille ??
			
		}
		else
		{
			cout<<"Pas de moyenne temporelle"<<endl;
		}
		*/
		
		Simulation * mySimu= new Simulation();
		
		//Lecture des parametres syteme dans fichier Sim
		//Le spl et nwk sont normalement vide... a remplir par la suite
		mySimu->read_data(fichsim.c_str());
	
		cout<<" --- lecture fichier sim ok "<<endl;
		
		//Lecture du fichier de commande contenant les grandeurs a evaluer
		ifstream com(fichcom.c_str());
		if(!com)
		{
			cerr << "Fichier de commande manquant : Analyse a partir du fichier sim " << com << endl;
			//return 0;
		}
		else
		{
			mySimu->sysA()->allFalse();
			com >> token;
			while(com )
			{
				if (token == "System_A{")  
					mySimu->sysA()->read_parameters(com);
				com >> token;
			}
		}
		cout<<" --- lecture fichier de commande ok "<<endl;
		
		//lecture du fichier time
		vector <double> t;
		double temp1,temp2;
		ifstream time(fichtemps.c_str());
		if(!time)
		{
			cerr << "Fichier de temps manquant  " << time << endl;
			//return 0;
		}
		else
		{
			
			while(time )
			{
				time>> temp1 >>temp2;
				//cout<<temp2<<" ";
				t.push_back(temp2);
			}
		}
		cout<<fichtemps<<endl;

		cout<<" --- lecture fichier temps ok "<<endl;
		
	//	ofstream analyze("Analyze.txt",ios::out);
		
		if( format==3)
		sprintf(nomFichier,"spl_nwk_%.3d.his",Ndeb);
		else if (format==4)
		sprintf(nomFichier,"spl_nwk_%.4d.his",Ndeb);
		else
		{
			cout<<" Bad format = "<<format<<endl;
			exit(0);
		}
		
		cout<<endl<<endl<<"************+++++ Chargement  :  "<<nomFichier<<endl;

		//mySimu->load_history(nomFichier );
		//mySimu->algo()->algoFill();
		
	
		
		mySimu->sys()->init();
		mySimu->sysA()->initAnalyse();		

		
		for( unsigned int i= Ndeb; i<=Nfin;i+=period)
		{
			if( format==3)
			sprintf(nomFichier,"spl_nwk/spl_nwk_%.3d.his",i);
			else if (format==4)
			sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);

			cout<<endl<<endl<<"************ Chargement  :  "<<nomFichier<<endl;
						
			//spl.includeFrom()=include;
			
			mySimu->load_history(nomFichier);
			mySimu->algo()->algoFill();
			
			//mySimu->spl()->updateBoundaries();
			//mySimu->spl()->radiusExtrema();
			//mySimu->spl()->definePeriodicityCV(true);
			
			mySimu->sysA()->plugRef();
			mySimu->sysA()->analyse(t[i],1,1);
			
			ofstream time("time.txt",ios::app);
			time<<t[i]<<endl;time.close();
			
				
			cout<<"**********************************************"<<endl;
		}
		
	}
	
	
	
	return 0;
}
