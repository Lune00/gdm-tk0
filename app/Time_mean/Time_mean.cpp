#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"

//! \brief Application to sequantially reload .his file and meaning scalar or/and distributions 
//! \author C. Voivret



int main (int argc, char * argv[]) 
{

	// Lire un fichier de parametre :
	// - fichiers scalaires a moyennner
	// - distributions a moyenner
	// - nom fichier listDir
	string token;

	if ( argv[1]!=0)
	{
		//lecture du fichier de parametre

	//	mySimu->sysA()->writePS("visuSNA.ps");
	}
	else
	{

		//Premier et dernier fichier
		string fichcom,fichrep,fichtemps,include,token,fichDist;
		//string anDir;
		char  nomFichier [100];	
		double tmin,tmax;
		unsigned int rgmin,rgmax;
	//	bool sf=false;

	//	cout<<"Fichier de grandeurs a evaluer ?"<<endl;
	//	cin>>fichcom;
		cout<<"Fichier contenant la liste des repertoires ?"<<endl;
		cin>>fichrep;
		//cout<<"Dossier contenant les donnees analysees ?"<<endl;
		//cin>>anDir;
		cout<<"Fichier contenant les distributions a moyenner ?"<<endl;
		cin>>fichDist;
		cout<<"Moyennes effectues entre tmin et tmax :"<<endl;
		cout<<"		tmin ?"<<endl;
		cin>>tmin;
		cout<<"		tmax ?"<<endl;
		cin>>tmax;


		ifstream inrep(fichrep.c_str());

		if(!inrep)
		{
			cerr << "Echec a l'ouverture :  " << fichrep.c_str() << endl;
			return 1;
		}
		vector <double> temps;

		inrep>>token; 
		while(inrep)
		{
			token.erase(  remove(token.begin(),token.end(),'t'),token.end());
			token.erase(  remove(token.begin(),token.end(),'='),token.end());
			istringstream iss(token);

			double temp;
			iss>>temp;
			temps.push_back(temp);

			inrep>>token;
		}
		cout<<" --- lecture fichier list rep ok "<<endl;
		inrep.close();

		//Determination des rangs
		rgmax=rgmin=1000000;

		for(unsigned int i=0;i<temps.size();++i)
		{
			if( temps[i] >= tmin && rgmin==1000000)
				rgmin=i;
			if( temps[i] >= tmax && rgmax==1000000)
				rgmax=i;
		}

		unsigned int j=0;
		while( temps[0]<=tmin)
		{	
			temps.erase( temps.begin());

		}

		j=temps.size()-1;
		while( temps[j]>=tmax)
		{	
			temps.pop_back();
			j=temps.size()-1;

			//temps.erase( temps.end());
		}

		for(unsigned int i=0;i<temps.size();++i)
		{
			cout<<temps[i]<<endl;
		}
		cout<<"rgmin = "<<rgmin<<" rgmax = "<<rgmax<<endl;


		//Lecture des valeur a moyenner
		double t1,trash;
		unsigned int rang=0;
		//DataSet z,sf,q,p,qop,t;
		/*sprintf(nomFichier,"Analyse/analyse.txt");
		ifstream in1(nomFichier);

		if(!in1)
		{
			cerr << "Echec a l'ouverture :  " << nomFichier << endl;
			return 1;
		}
		while (in1)
		{
			in1>>t1>>trash>>trash>>sftemp>>s1>>s2>>trash>>trash>>ztemp>>trash;
			if( rang>= rgmin && rang<= rgmax)
			{
				//cout<<rang<<" "<<sftemp<<endl;

				z.add(ztemp);
				sf.add(sftemp);
				p.add(s1+s2);
				q.add( max(s1,s2)-min(s1,s2));
				qop.add((max(s1,s2)-min(s1,s2))/(s1+s2));
				rang++;
			}
			else
				rang++;

			if( in1.eof()) break;
		}
		in1.close();


		sf.extractValues();
		z.extractValues();
		q.extractValues();
		p.extractValues();
		qop.extractValues();

		system("mkdir Mean");
		ofstream out1("Mean/ScalarMean.txt",ios::out);
		out1<<"sf  = "<<sf.mean()<<" "<<sf.variance()<<endl;
		out1<<"z   = "<<z.mean()<<" "<<z.variance()<<endl;
		out1<<"q   = "<<q.mean()<<" "<<q.variance()<<endl;
		out1<<"p   = "<<p.mean()<<" "<<p.variance()<<endl;
		out1<<"qop = "<<qop.mean()<<" "<<qop.variance()<<endl<<endl;
		out1.close();
		//cout<<"analyse ok"<<endl;*/



		//Anisotropy----------
		sprintf(nomFichier,"Analyse/Anisotropies/FAnisotropy.txt");
		ifstream in2(nomFichier);

		double ac1,an1,at1,al1,approx1,qop1;
		DataSet ac,an,at,al,approx,qop;
		rang=0;
		if(!in2)
		{
			cerr << "Echec a l'ouverture :  " << nomFichier << endl;
			return 1;
		}
		while (in2)
		{
			in2>>t1>>trash>>ac1>>an1>>at1>>al1>>approx1>>qop1;
			if( rang>= rgmin && rang< rgmax)
			{
				//cout<<ac1<<" "<<an1<<" "<<at1<<" "<<al1<<endl;

				ac.add(ac1);
				an.add(an1);
				at.add(at1);
				al.add(al1);
				approx.add(approx1);
				qop.add(qop1);
				rang++;
			}
			else
				rang++;

			if( in2.eof()) break;

		}
		in2.close();

		ac.extractValues();
		an.extractValues();
		at.extractValues();
		al.extractValues();
		approx.extractValues();
		qop.extractValues();
		system("mkdir Mean");

		ofstream out2("Mean/ScalarMean.txt",ios::app);
		out2<<"ac = "<<ac.mean()<<" "<<ac.variance()<<endl;
		out2<<"an = "<<an.mean()<<" "<<an.variance()<<endl;
		out2<<"at = "<<at.mean()<<" "<<at.variance()<<endl;
		out2<<"al = "<<al.mean()<<" "<<al.variance()<<endl;
		out2<<"approx = "<<approx.mean()<<" "<<approx.variance()<<endl;
		out2<<"qop = "<<qop.mean()<<" "<<qop.variance()<<endl;
		out2.close();


		//Correlations----------
		sprintf(nomFichier,"Analyse/linearCor.txt");
		ifstream in3(nomFichier);
	//	lc<<time<<" "<<epsxy<<" "<<r.mean() <<" "<<l.mean()<<" "<<lcr<<" "<<lcl<<" "<<lcld<<" "<<lclt<<" "<<lcrt<<" "<<lcRR<<" "<<lcfltemp<<endl;

		DataSet lcr,lcl,lcld,lclt,lcrt,lcRR,lcfsl,lcfwl;
		double cr,cl,cld,clt,crt,cRR,cfsl,cfwl;
		rang=0;
		if(!in3)
		{
			cerr << "Echec a l'ouverture :  " << nomFichier << endl;
			return 1;
		}
		while (in2)
		{
			in2>>t1>>trash>>trash>>trash>>cr>>cl>>cld>>clt>>crt>>cRR>>cfsl>>cfwl;
			if( rang>= rgmin && rang< rgmax)
			{
				//cout<<ac1<<" "<<an1<<" "<<at1<<" "<<al1<<endl;

				lcr .add(cr);
				lcl .add(cl);
				lcld.add(cld);
				lclt.add(clt);
				lcrt.add(crt);
				lcRR.add(cRR);
				lcfsl.add(cfsl);
				lcfwl.add(cfwl);
				
				rang++;
			}
			else
				rang++;

			if( in2.eof()) break;

		}
		in2.close();
		lcr .extractValues();
		lcl .extractValues();
		lcld.extractValues();
		lclt.extractValues();
		lcrt.extractValues();
		lcRR.extractValues();
		lcfsl.extractValues();
		lcfwl.extractValues();
		

		ofstream out3("Mean/ScalarMeanCor.txt",ios::app);
		out3<<"cr  = "<<lcr .mean()<<" "<<lcr .variance()<<endl;
		out3<<"cl  = "<<lcl .mean()<<" "<<lcl .variance()<<endl;
		out3<<"cld = "<<lcld.mean()<<" "<<lcld.variance()<<endl;
		out3<<"clt = "<<lclt.mean()<<" "<<lclt.variance()<<endl;
		out3<<"crt = "<<lcrt.mean()<<" "<<lcrt.variance()<<endl;
		out3<<"cRR = "<<lcRR.mean()<<" "<<lcRR.variance()<<endl;
		out3<<"cfsl= "<<lcfsl.mean()<<" "<<lcfsl.variance()<<endl;
		out3<<"cfwl= "<<lcfwl.mean()<<" "<<lcfwl.variance()<<endl;
		out3.close();

		cout<<" Donnees scalaires traitees "<<endl;
//--------------------------------------------------------------

		unsigned int nbinfn,perfn,wfn;
		unsigned int nbinft,perft,wft;
		unsigned int nbinsize;
		unsigned int nbinweak;
		bool pdffn,pdfft,size,fracweak;
		pdffn=pdfft=size=fracweak=false;

		if (fichDist=="Aucun") return 0;
		ifstream indist(fichDist.c_str());

		if(! indist)
		{
			cerr << "Echec a l'ouverture :  " << fichDist.c_str() << endl;
			return 1;
		}

		while(indist)
		{
			indist>>token;

			if( token == "pdffn") 
			{
				pdffn=true;
				indist>>nbinfn;
				indist>>perfn;
				indist>>wfn;
			}
			else if( token == "pdfft") 
			{
				pdfft=true;
				indist>>nbinft;
				indist>>perft;
				indist>>wft;
			}
			else if( token == "size") 
			{
				size=true;
				indist>>nbinsize;
			}
			else if( token == "fracweak") 
			{
				fracweak=true;
				indist>>nbinweak;
			}
			else
				cout<<" unknown parameter "<<endl;

			if( indist.eof()) break;

		}
		cout<<" --- lecture fichier list rep ok "<<endl;
		inrep.close();


		//Lecture des distributions a moyenner
		DataSet fn,ft;
		pointSet fdw,fds,fracw;
	//	unsigned int nbin=100,period=1,width=3;
		for( unsigned int i =0; i<temps.size();++i)
		{
			if(pdffn)
			{
				sprintf(nomFichier,"t=%.3f/pdf/fn_brut.txt",temps[i]);
				ifstream infn(nomFichier);
				if(!infn)
				{
					cerr << "Echec a l'ouverture :  " << nomFichier << endl;
				//return 1;
				}
				else
				{
					fn.addRead(infn);
					infn.close();
				//cout<<fn.setSize()<<endl;
				//fn.setClear();	
				}
			}

			if(pdfft)
			{
				sprintf(nomFichier,"t=%.3f/pdf/ft_brut.txt",temps[i]);
				ifstream inft(nomFichier);
				if(!inft)
				{
					cerr << "Echec a l'ouverture :  " << nomFichier << endl;
					pdfft=false;
				//return 1;
				}
				else
				{
					ft.addRead(inft);
					inft.close();
				//cout<<fn.setSize()<<endl;
				//fn.setClear();	
				}
			}
			if( size)
			{
				sprintf(nomFichier,"t=%.3f/forceCorr/ftotdW_brut.txt",temps[i]);
				ifstream infdw(nomFichier);
				if(!infdw)
				{
					cerr << "Echec a l'ouverture :  " << nomFichier << endl;
					size=false;
				//return 1;
				}
				else
				{
					fdw.addRead(infdw);
					infdw.close();
				//cout<<fn.setSize()<<endl;
				//fn.setClear();	
				}
				
				sprintf(nomFichier,"t=%.3f/forceCorr/ftotdS_brut.txt",temps[i]);
				ifstream infds(nomFichier);
				if(!infds)
				{
					cerr << "Echec a l'ouverture :  " << nomFichier << endl;
					size=false;
				//return 1;
				}
				else
				{
					fds.addRead(infds);
					infds.close();
				//cout<<fn.setSize()<<endl;
				//fn.setClear();	
				}
			}
			if( fracweak)
			{
				sprintf(nomFichier,"t=%.3f/Connect/fracweak_brut.txt",temps[i]);
				ifstream inw(nomFichier);
				if(!inw)
				{
					cerr << "Echec a l'ouverture :  " << nomFichier << endl;
					fracweak=false;
				//return 1;
				}
				else
				{
					fracw.addRead(inw);
					inw.close();
				//cout<<fn.setSize()<<endl;
				//fn.setClear();	
				}
				
			}
		}
		if(pdffn)
		{
			pointSet pdf = fn.kernelPdf( nbinfn,.05);
			pointSet pdfm = pdf.mobileMean(perfn,wfn);
			pointSet pdfs = (fn.slidingPdf( nbinfn) ).mobileMean(perfn,wfn);
			pdf.write("Mean/pdffn_ker.txt");
			pdfm.write("Mean/pdffn_kersli.txt");
			pdfs.write("Mean/pdffn_sli.txt");
			fn.setClear();
			pdf.psetClear();
			pdfm.psetClear();
			pdfs.psetClear();
			cout<<" pdffn ok "<<endl;
		}
		if(pdfft)
		{
			pointSet pdf = ft.kernelPdf( nbinft,.05);
			pointSet pdfm = pdf.mobileMean(perft,wft);
			pointSet pdfs = (ft.slidingPdf( nbinft) ).mobileMean(perft,wft);
			pdf.write("Mean/pdfft_ker.txt");
			pdfm.write("Mean/pdfft_kersli.txt");
			pdfs.write("Mean/pdfft_sli.txt");
			ft.setClear();
			pdf.psetClear();
			pdfm.psetClear();
			pdfs.psetClear();
			cout<<" pdfft ok "<<endl;
		}
		if(size)
		{
			pointSet S= fds.slidingHisto(nbinsize,.05);
			pointSet W= fdw.slidingHisto(nbinsize,.05);
			S.write("Mean/fds.txt");
			W.write("Mean/fdw.txt");
			fds.psetClear();
			fdw.psetClear();
			
		}
		if(fracweak)
		{
			pointSet FS= fracw.slidingHisto(nbinweak,.05);
			FS.write("Mean/fracweak.txt");
			fracw.psetClear();
		}


		cout<<" Donnees vectorielles traitees"<<endl;



	}



	return 0;
}
