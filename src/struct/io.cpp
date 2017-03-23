#include "io.hpp"

// Check whether a file exists
// Maybe this function shouldn't be fast enought
// but it is not OS dependent
bool fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
    {
    fin.close();
    return true;
    }
  fin.close();
  return false;
}


void write_mgpost(const char * fname, Sample& spl, Network& nwk, int stateId, double t)
{
  ofstream fmgpost(fname);
  unsigned int k = 0;
  unsigned int kmax = nwk.linter().size();
	
  fmgpost << "<?xml version=\"1.0\"?>" << endl
	  << " <mgpost mode=\"2D\">" << endl
	  << "  <state id=\"" << stateId 
	  << "\" time=\"" << t << "\">" << endl;
	
  for (unsigned int i=0 ; i<spl.lbody().size() ; ++i)
    {
    fmgpost << "   <body>" << endl;
    spl.body(i)->writeMGP(fmgpost);
    
    while ((k<kmax) && (nwk.inter(k)->lexifirst()->id() == spl.body(i)->id()))
      {
      nwk.inter(k)->writeMGP(fmgpost);
      ++k;
      }
		
    fmgpost << "   </body>" << endl;
    }
	
  fmgpost << "  </state>" << endl
	  << " </mgpost>" << endl;
	
}

void write_ps(const char * fname, Sample& spl, Network& nwk, int stateId, double t)
{ // EN TRAVAUX...
  ofstream fps(fname);
	
  fps << "%%!PS-Adobe-3.0 EPSF-3.0" << endl
      << "%%BoundingBox: 0 0 100 250" << endl
      << "newpath" << endl;
		
  for (unsigned int i=0;i<spl.lbody().size();i++)
    spl.body(i)->writePS(fps);		
		
}


void read_data(const char* name, 
               Algo& algo, Sample& spl, Network& nwk, System& sys, 
               GroupData& grpDat, GroupRelationData& grpRel)
{
  string token;
  ifstream datafile(name);
  if(!datafile)
    {
    cerr << "@read_data, cannot open file " << name << endl;
    return;
    }
  
  datafile >> token;
  while(datafile)
    {
     
    if(token == "includeFile")
      {
      string filename;
      datafile >> filename;
      if (fileExists(filename))
        read_data(filename.c_str(),algo,spl,nwk,sys,grpDat,grpRel);
      }
    
    if (token == "System{")            sys.read_parameters(datafile);
    if (token == "Parameters{")        algo.read_parameters(datafile);
    if (token == "GroupData{")         grpDat.read(datafile);    
    if (token == "GroupRelationData{") grpRel.read(datafile);
    if (token == "Sample{")            spl.read(datafile);
    
    datafile >> token;
    }
}


void save_data(const char* name, 
               Algo& algo, Sample& spl, Network& nwk, System& sys, 
               GroupData& grpDat, GroupRelationData& grpRel,bool saveNetwork)
{
  ofstream datafile(name);
  if(!datafile)
  {
    cerr << "@save_data, cannot open file " << name << endl;
    return;
  }
  
  datafile << "System{" << endl;
  sys.write_parameters(datafile);
  datafile << "}" << endl << endl;
  
  datafile << "Parameters{" << endl;
  algo.write_parameters(datafile);
  datafile << "}" << endl << endl;

  datafile << "GroupData{" << endl;
  grpDat.write(datafile);
  datafile << "}" << endl << endl;
  
  datafile << "GroupRelationData{" << endl;
  grpRel.write(datafile);
  datafile << "}" << endl << endl;
  
  datafile.setf(ios_base::scientific);
  
  datafile << "Sample{" << endl;
  spl.write(datafile);
  datafile << "}" << endl << endl;
  
  // Sometimes the network is very big and the user should prefer to avoid the saving
  if(saveNetwork)
    {
    datafile << "Network{" << endl;
    nwk.write(datafile);
    datafile << "}" << endl << endl;
    }
}


void history_read(const char * name, Sample& spl, Network& nwk)
{
  ifstream histFile(name);
  if(!histFile)
    {
    cerr << "@history_read, cannot open file " << name << endl;
    return;
    }  
  
  string token;
  histFile >> token;
  
  while(histFile)
    {
    if (token == "Sample{")  spl.read(histFile);
    if (token == "Network{") nwk.read(histFile);
    
    histFile >> token;
    }
  
  nwk.associate(spl);
  
  // Compute the data that have not been saved
  for (unsigned int k=0 ; k < nwk.linter().size() ; ++k)
    {	
    nwk.inter(k)->Frame();
    }
}



void history_write(unsigned int number, Sample& spl, Network& nwk, GroupRelationData& grpRel, bool TwoFiles, bool saveNetwork, bool compact)
{
  char name[50];
    //char name_spl[50];
        
    
  if(TwoFiles || !saveNetwork)
  {
      system("mkdir -p spl");
      sprintf(name,"spl/spl_%04d.his",number);
  }
    
  else
  {
      system("mkdir -p spl_nwk");
      sprintf(name,"spl_nwk/spl_nwk_%04d.his",number);
  
  }
    
  ofstream datafile(name);
  if(!datafile)
    {
    cerr << "@history_write, cannot open file " << name << endl;
    return;
    }
  
  datafile.setf(ios_base::scientific);
  
  datafile << "Sample{" << endl;
  spl.write(datafile);
  datafile << "}" << endl << endl;
  
  if(saveNetwork)
    {
    if(TwoFiles) 
      {
          
      system("mkdir nwk");
      sprintf(name,"nwk/nwk_%04d.his",number);
      ofstream nwkfile(name);
      if(!nwkfile)
        {
        cerr << "@history, cannot open file " << name << endl;
        return;
        }
      
      nwkfile << "Network{" << endl;
      nwk.write(nwkfile);
      nwkfile << "}" << endl << endl;
      }
    else
      {
	      datafile << "Network{" << endl;
   		  nwk.write(datafile);
    	  datafile << "}" << endl << endl;
      }
    }
	
    
     
     
     
    /////////////////
    
	if( compact)
	{
		char command[100];
		sprintf(command,"gzip %s",name);
		system(command);
	}
	if( grpRel.existfragmentation())
	{
      	cout<<"OK"<<endl;
      	grpRel.fragmentation()->write(datafile);
    }
	
}


void sample_read(const char * name, Sample& spl)
{
  ifstream histFile(name);
  if(!histFile)
    {
    cerr << "@sample_read, cannot open file " << name << endl;
    return;
    }  
  
  string token;
  histFile >> token;
  
  while(histFile)
    {
    if (token == "Sample{")  spl.read(histFile);
    
    histFile >> token;
    }
}

void sample_write(const char * name, Sample& spl)
{  
  ofstream datafile(name);
  if(!datafile)
    {
    cerr << "@sample_write, cannot open file " << name << endl;
    return;
    }
  
  datafile.setf(ios_base::scientific);
  
  datafile << "Sample{" << endl;
  spl.write(datafile);
  datafile << "}" << endl << endl;    
}

void network_write(const char * name, Network& nwk)
{  
  ofstream datafile(name);
  if(!datafile)
    {
    cerr << "@network_write, cannot open file " << name << endl;
    return;
    }
  
  datafile.setf(ios_base::scientific);
  
  datafile << "Network{" << endl;
  nwk.write(datafile);
  datafile << "}" << endl << endl;    
}


// Only for disks
// V. Richefeu
void sample_write_tapioK(const char * name, Sample& spl)
{  
	ofstream datafile(name);
	if(!datafile)
	{
		cerr << "@sample_write_tapioK, cannot open file " << name << endl;
		return;
	}

// Count disks and compute boundaries
	unsigned int nbdisks = 0;
	double xmin =  1.0e10;
	double xmax = -1.0e10;
	double ymin =  1.0e10;
	double ymax = -1.0e10;
	double val = 0;
	for (unsigned int i = 0 ; i < spl.lbody().size() ; ++i)
	{

		if ( (*spl.body(i)).type() == _type_disk) 
		//if (typeid(*(spl.body(i))) == typeid(disk))

		{
			++nbdisks;
			val = spl.body(i)->xmin();
			if (xmin > val) xmin = val;
			val = spl.body(i)->xmax();
			if (xmax < val) xmax = val;
			val = spl.body(i)->ymin();
			if (ymin > val) ymin = val;
			val = spl.body(i)->ymax();
			if (ymax < val) ymax = val;
		}
	}

	datafile.setf(ios_base::scientific);
	datafile << "--" << endl; // Info line
	datafile << nbdisks + 4 << " 0" <<endl;

// Four walls
	datafile << "100" << endl;
	datafile << xmin << " 0" << endl;
	
	datafile << "101" << endl;
	datafile << xmax << " 0" << endl;	

	datafile << "102" << endl;
	datafile << ymin << " 0" << endl;

	datafile << "103" << endl;
	datafile << ymax << " 0" << endl;

// Disks
	disk * D = 0;
	for (unsigned int i = 0 ; i < spl.lbody().size() ; ++i)
	{
		//if (typeid(*(spl.body(i))) == typeid(disk))
		
		if ( (*spl.body(i)).type() == _type_disk) 
		{
			D = dynamic_cast<disk*> (spl.body(i));
			datafile << "0" << endl;
			datafile << D->R() << " " << D->x()  << " " << D->y()  << " " << D->rot()  << " "
				                      << D->vx() << " " << D->vy() << " " << D->vrot() << endl;    
		}
	}  
}

