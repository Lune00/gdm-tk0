#include "simulation.hpp"

void gdm_signal_handler(int Nsig)
{
	switch(Nsig)
    {
    case SIGINT : 
      cout << endl;
      cout << "Computation has been stopped by user." << endl;
	  //quit();
      //save_dataSim("stopped.sim");
	//simulation_write("stopped.sim",this);
      cout << "Data saved in 'stopped.sim'" << endl;
      cout << "Fin simulation" << endl;
      
      exit(0);
    }
}

void Simulation::init()
{
	if (sigemptyset(&sigact_.sa_mask) == -1) 
		cout << "Signal mask can not be initialized" << endl;
		
	sigact_.sa_handler= gdm_signal_handler;
	
	if (sigaction(SIGINT, &sigact_, NULL) == -1) 
		cout << "Can not affect SIGINT" << endl;

	if (algo_ == 0)
		gdm::fatal("@Simulation::init: No algorithm!");

	if (sys_ == 0)
		gdm::fatal("@Simulation::init: No system!");
	//if (sys_->check() == 0) 
	//    gdm::fatal("@Simulation::init: Problem in system!");
	
	if (doAnalyse_ && sysA_==0) 
	 gdm::fatal("No System_A has been defined!");
  
	if (spl_ == 0)
		gdm::fatal("@Simulation::init: No sample!");

	if (nwk_ == 0)
		gdm::fatal("@Simulation::init: No network!");

	if (grpDat_ == 0)
		gdm::fatal("@Simulation::init: No Group Data!");

	if (grpRel_ == 0)
		gdm::fatal("@Simulation::init: No Group Relation Data!");
		


}

void Simulation::speak()
{
  switch(speakLevel_)
    {
    case 0:
      return;
    case 1:  
    default:
	cout <<  "Step: " << ns_ << " (" << std::setprecision(3)<< 100.0*(double)ns_/(double)(nsf_-nsi_) << " %)" << endl;
	cout << "Time: " << time_ << endl;
      return;
    }
}

void Simulation::read_data(const char* name)
{ // Probleme ou pas pour la relecture? a voir ...
 
  ifstream datafile(name);
  if(!datafile)
    {
    cerr << "@Simulation::read_data, cannot open file " << name << endl;
    return;
    }

  string token;
  string type;
  
  datafile >> token;
  while(datafile)
    {
    if(token == "includeFile")
      {
      string filename;
      datafile >> filename;
	  
      if (fileExists(filename))
		{
        Simulation::read_data(filename.c_str());
		}
      }
    
    if (token == "Simulation{")
      {
      datafile >> token;
      while(datafile)
        {
        if      (token == "speakLevel")     datafile >> speakLevel_; 
        else if (token == "ns")             datafile >> ns_;
        else if (token == "nsi")            datafile >> nsi_;
        else if (token == "nsf")            datafile >> nsf_;
        else if (token == "nSpeak")         datafile >> nSpeak_;
        else if (token == "nHist")          datafile >> nHist_;
        else if (token == "nAnalyse")       datafile >> nAnalyse_;
        else if (token == "historyNetwork") historyNetwork_ = true;
        else if (token == "numFileHist")    datafile >> numFileHist_;
        else if (token == "twoFilesHist")   twoFilesHist_ = true;
        else if (token == "compactHist")    compactHist_ = true;
        else if (token == "}")              break;
        else cerr << "@Simulation::read_data, Unknown parameter: " << token << endl;
        
        datafile >> token;
        }
      }
    
    if (token == "System{")
      {
      datafile >> type; 
      if (sys_ == 0) sys_ = System::factory(type);
      sys_->read_parameters(datafile);

      }

	if (token == "System_A{")
      {
	
	  doAnalyse_= true;
      datafile >> type; 
      if (sysA_ == 0) sysA_ = System_A::factory(type);
      sysA_->read_parameters(datafile);

      }

    if (token == "Parameters{" || token == "Algo{")
      {
      datafile >> type;
      if (algo_ == 0) algo_ = Algo::factory(type);
      algo_->read_parameters(datafile);

      }
	
    if (token == "GroupData{")
      {
      if (grpDat_ == 0) grpDat_ = new GroupData();
      grpDat_->read(datafile);
    
      }
	
    if (token == "GroupRelationData{")
      {
      if (grpRel_ == 0) grpRel_ = new GroupRelationData();
      grpRel_->read(datafile);
      }

	
    if (token == "Sample{")
      {
      if (spl_ == 0) spl_ = new Sample();
      sys_->ldof() = spl_->read(datafile);
      }
	

    if (token == "Network{")
      {
      if (nwk_ == 0) nwk_ = new Network();
      nwk_->read(datafile);
      //cout << "nwk_->read(datafile)" << " " << nwk_->linter().size() << endl;
      //if (spl_ == 0) nwk_->read(datafile);
      //else nwk_->read(datafile,*spl_);
      }
	
    datafile >> token;
    }

  // a deplacer dans init()
  
  // Check for existence of needed entities
  if (!sys_)     gdm::fatal("No System has been defined!");
  if (!algo_)    gdm::fatal("No Algorithm has been defined!");
  if (!grpDat_)  gdm::fatal("No GroupData has been defined!");
  if (!grpRel_)  gdm::fatal("No GroupRelationData has been defined!");
  if (!spl_)  
    {
    spl_ = new Sample();
    gdm::warning("No Sample defined! an empty sample has been build");
    }
  
  if (!nwk_) nwk_ = new Network();
  else nwk_->associate(*spl_);

  if (doAnalyse_ && !sysA_) gdm::fatal("No System_A has been defined!");
    //else gdm::warning("No System_A defined! No analysis will be made ");
  

  // Here we plug all this entities togheter with !!!
  assert(spl_);
  assert(nwk_);
  
  sys_->plugSample(spl_);
  sys_->plugNetwork(nwk_);
  sys_->plugGroupRelationData(grpRel_);
  algo_->plugSample(spl_);
  algo_->plugNetwork(nwk_);
  algo_->plugSystem(sys_);
  algo_->plugGroupData(grpDat_);
  algo_->plugGroupRelationData(grpRel_);

  if (doAnalyse_)  sysA_->plugSystem(sys_);

//cout<<" taille dof simu "<<sys_->ldof().size()<<endl;

}


void Simulation::save_data(const char* name)
{
  /*
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
  if(historyNetwork_)
    {
    datafile << "Network{" << endl;
    nwk.write(datafile);
    datafile << "}" << endl << endl;
    }
   */
}

void Simulation::run()
{
	gdm::fun_GDM();
	char name[15];
	
	algo_->algoFill();
	sys_->init();
	algo_->algoFill();

	if (doAnalyse_) sysA_->initAnalyse();
	algo_->stand();

	algo_->look();
	cout<<" nHist "<<nHist_<<endl;
	
	for (ns_ = nsi_+1 ; ns_ <= nsf_ ; ++ns_)
	{
		time_ += algo_->dt();
		algo_->step();
		algo_->look();
		algo_->hand(ns_);
		
		if ( doAnalyse_ && ns_%nAnalyse_  == 0 )
		{
			sysA_->analyse(time_,ns_,nsf_);
		}

		if (ns_%nSpeak_ == 0) 
		{
			this->speak();
			algo_->speak();
			nwk_->speak();
		}

		if (ns_%nHist_  == 0)
		{
			cout<<"Ecriture "<<numFileHist_<<endl;
			
			history_write(numFileHist_, *spl_, *nwk_, *grpRel_, twoFilesHist_, historyNetwork_, compactHist_);
			ofstream time("time.txt",ios::app);
			time<<numFileHist_<<" "<<time_<<endl;
			time.close();
			
			sprintf((char *) name, "mgp.out.%04d",numFileHist_);
			//write_mgpost(name,*spl_,*nwk_,numFileHist_,time_);
			numFileHist_ += 1;
            
        
		}
	}
	cerr<<"La simulation est terminee."<<endl;
}

void Simulation::load_history(const char * fname )
{
  ifstream histFile(fname);
  if(!histFile)
    {
    cerr << "@history_read, cannot open file " << fname << endl;
    return;
    }  
  
	purge(this->spl()->lbody());
	purge(this->nwk()->linter());

	this->spl()->lbody().clear();
	this->nwk()->linter().clear();
	this->nwk()->clist().clear();
  string token;
  histFile >> token;
  
  while(histFile)
    {
    if (token == "Sample{")  this->sys()->ldof() = this->spl()->read(histFile);
    if (token == "Network{") this->nwk()->read(histFile);
    if (token == "Fragmentation{") this->sys()->grpRel()->fragmentation()->load_history(histFile);
    
    histFile >> token;
    }
  
  this->nwk()->associate(*this->spl());
//cout<<"taille bande spl = "<<this->spl()->leftband().size()<<endl;
  cout<<"taille nwk = "<<this->nwk()->linter().size()<<endl;
// Compute the data that have not been saved
  for (unsigned int k=0 ; k < this->nwk()->linter().size() ; ++k)
 {	
    this->nwk()->inter(k)->Frame2();
 }
}
