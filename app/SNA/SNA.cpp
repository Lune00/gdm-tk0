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
	cout<<"Premiere option"<<endl;	
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
		
//		write_mgpost("mgp.out.001",*mySimu->spl(),*mySimu->nwk(),0,0);
	//	mySimu->sysA()->writePS("visuSNA.ps");
	}
	else
	{
		
		//Premier et dernier fichier
		unsigned int Ndeb,Nfin,period,format=1;
		string fichcom,fichsim,fichtemps,include,token;
		char  nomFichier [100],commande[100],princDir[100];	
		
		cin>>Ndeb;
		cin>>Nfin;
		cin>>format;
		cin>>period;
		cin>>fichsim;
		cin>>fichcom;
		cin>>fichtemps;
		cin>>princDir;
		cout<<"Numeros du premier et dernier fichier a traiter :"<<Ndeb<<"-"<<Nfin<<endl;
		cout<<"Nom du dossier principal où stocker les fichiers d'analyse : "<<princDir<<endl;
		sprintf(commande,"mkdir -p %s",princDir);
		system( commande);
		cout<<"Nom du fichier temps à charger :"<<fichtemps<<endl;	

		//lecture du fichier time
		vector <double> t;
		double temp1,temp2;
		ifstream time(fichtemps.c_str());
		string line;

		if(!time.is_open())
		{
			cerr << "Fichier de temps manquant  " << time << endl;
			return 0;
		}
		else
		{
			cout<<"Lecture du fichier "<<fichtemps.c_str()<<endl;
			
			while(time)
			{
				time>> temp1 >>temp2;
				if(time.eof()) break;
				t.push_back(temp2);
			}
		}
		time.close();
		
		Simulation * mySimu= new Simulation();

		//On réécrit le fichier temps:
		ofstream timereset(fichtemps.c_str());
		unsigned int indice=0;
		if(timereset.is_open()){
			for(std::vector<double>::iterator it = t.begin();it != t.end();++it){
				timereset<<indice<<" "<<*it<<endl;
				indice++;
			}
		}
		timereset.close();


		//Lecture des parametres syteme dans fichier Sim
		//Le spl et nwk sont normalement vide... a remplir par la suite
		cout<<"--- lecture du fichier Simu.sim"<<endl;
		mySimu->read_data(fichsim.c_str());
		mySimu->sys()->init();
		//OK
	
		cout<<"--> lecture fichier sim ok "<<endl;
		
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
				cerr<<"Read parameters"<<endl;
				com >> token;
			}
		}
		cout<<"--> lecture fichier de commande ok "<<endl;
		
		cout<<"--> lecture fichier temps ok "<<endl;
		
	//	ofstream analyze("Analyze.txt",ios::out);
		
		if( format==3)
		sprintf(nomFichier,"spl_nwk_%.3d.his",Ndeb);
		else if (format==4)
		sprintf(nomFichier,"spl_nwk_%.4d.his",Ndeb);
		else if (format==5)
		sprintf(nomFichier,"spl_nwk_%.5d.his",Ndeb);
		else
		{
			cout<<" Bad format = "<<format<<endl;
			exit(0);
		}

		mySimu->sysA()->initAnalyse();		

		for( unsigned int i= Ndeb; i<=Nfin;i+=period)
		{
			if( format==3)
			sprintf(nomFichier,"spl_nwk/spl_nwk_%.3d.his",i);
			else if (format==4)
			sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);
			else if (format==5)
			sprintf(nomFichier,"spl_nwk/spl_nwk_%.5d.his",i);

			cout<<endl<<endl<<"--> Chargement  :  "<<nomFichier<<endl;
			//Load history a checker
						
			mySimu->load_history(nomFichier);
			cout<<"load_history() done"<<endl;
			mySimu->algo()->algoFill();
			cout<<"algoFill() done"<<endl;	
			mySimu->sysA()->plugRef();
			cout<<"plugRef() done"<<endl;
			cout<<"t["<<i<<"]="<<t[i]<<endl;
			mySimu->sysA()->analyse(t[i],1,1);
			cout<<"analyse done"<<endl;
			
			cout<<"*o*0ro*************0******^r******0********1**"<<endl;
		}
	//	delete mySimu ;
		
	}
	
	
	
	return 0;
}
