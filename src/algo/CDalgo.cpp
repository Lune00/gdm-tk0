
//    step --> update step num., update positions, drive & share functions of the system,
//             update system velocities
//    look --> contact,fres,mass matrix & array of coefficients,
//             iter, velocities update
//    hand --> increment time, update verlet list 
//
//
#include "CDalgo.hpp"

using namespace std;

void CDalgo::read_parameters(istream & is)
{
	string token;

	is >> token;
	while(is)
	{	
		if      (token == "speakLevel") is >> speakLevel_;
		else if (token == "dt")         is >> dt_;
		else if (token == "epsf")       is >> epsf_;
		else if (token == "niterconv")  is >> niterconv_;
		else if (token == "nitermx")    is >> nitermx_;
		else if (token == "nitermn")    is >> nitermn_;
		else if (token == "nver")       is >> nver_;
		else if (token == "nsuper")     is >> nsuper_;
		else if (token == "}")          break;
		else cerr << "@CDalgo::read_parameters, Unknown parameter: " << token << endl;

		is >> token;
	}
}

void CDalgo::write_parameters(ostream & os)
{
	os << "CD" << endl;
	os << "speakLevel " << speakLevel_ << endl;
	os << "dt      "    << dt_         << endl;
	os << "epsf    "    << epsf_       << endl;
	os << "nitermn "    << nitermn_    << endl;
	os << "nitermx "    << nitermx_    << endl;
	os << "nver    "    << nver_       << endl;
}

void CDalgo::speak()
{
	switch(speakLevel_)
	{
		case 0:
		return;

		default:
		cout << "Number of iterations: " << niter_ << endl;
		cout << flush;
		return;
	}  
}


// Take care of contacts and forces
void CDalgo::look()
{

	//cout<<"look"<<endl;
// Contact list, local frames and relative velocities
	contact();

// Compute force resultants
	fres();

// Mass matrix, arrays of coefficients
//Seulement pour les disk appel des bons parametres de restitution
	for(unsigned int k = 0 ; k < nwk_->clist().size() ; ++k)
	{
	  //cout<<"test dkdk? 1"<<endl;
	  if(nwk_->inter(nwk_->clist(k))->type() == 0 )
	  {
	    nwk_->inter(nwk_->clist(k))->CDcoeff(grpRel_);
	  }
	  else
	  {
	    nwk_->inter(nwk_->clist(k))->CDcoeff();
	  }
	}


	// The main loop: compute forces
	niter_ = (this->*iterFunc)();
	//cout<<"niter:=	"<<niter_<<endl;

	// Mise à jours les fractures
	if (grpRel_->existfragmentation()) resistance();

	// .......... New velocities

	// Imposed forces
	double invmassi,momi;
	for (unsigned int i = 0 ; i < sys_->lctrl().size() ; ++i)	
	{
	  invmassi = 1.0 / spl_->body(i)->mass();
	  momi = spl_->body(i)->mom();

	  if(sys_->ctrl(i).x()   == _FORCE)
	  {
	    //cout<<" body "<<i<<endl;
	    spl_->body(i)->vx() += dt_ * spl_->body(i)->fx()   * invmassi;
	  }
	  if(sys_->ctrl(i).y()   == _FORCE) 
	    spl_->body(i)->vy() += dt_ * spl_->body(i)->fy()   * invmassi;
	  if(sys_->ctrl(i).rot() == _FORCE)
	    spl_->body(i)->vrot() += dt_ * spl_->body(i)->frot() / momi;

	}

	//sys_->share(); //vr en test...

	// Free bodies
	for (unsigned int i = sys_->lctrl().size() ; i < spl_->lbody().size() ; ++i)
	{
	  if( spl_->body(i)->bodyDof()==NULL)
	  {

	    invmassi = 1./spl_->body(i)->mass();
	    momi  = spl_->body(i)->mom();

	    spl_->body(i)->vx()   += dt_ * spl_->body(i)->fx()   * invmassi;
	    spl_->body(i)->vy()   += dt_ * spl_->body(i)->fy()   * invmassi;
	    spl_->body(i)->vrot() += dt_ * spl_->body(i)->frot() / momi;
	  }

	}	
	// bodies controled by a dof (different of CONTROL)
	for( unsigned int i=0;i < sys_->ldof().size();++i)
	{
	  sys_->ldof(i)->imposeVelocityOfForce(dt_);
	  //cout<<" dof "<<i<<endl;
	}
}

// Update the contact list, local frames and relative velocities
void CDalgo::contact()
{
  // No-contact case
  if (nwk_->linter().empty()) return;

  nwk_->clist().clear();
  /*if (grpRel_->existfragmentation()) 
    for (unsigned int k=0 ; k < nwk_->linter().size() ; ++k)
    {	
    nwk_->inter(k)->dAct()= (grpRel_->getLaw(nwk_->inter(k)->first()->grp(),nwk_->inter(k)->second()->grp()))->dAct(nwk_->inter(k));
    nwk_->inter(k)->Frame();

    unsigned int id1,id2,cluster_id;
    bool stop=false;
    bool samecluster=false;

    id1=nwk_->inter(k)->first()->id()-4;
    id2=nwk_->inter(k)->second()->id()-4;

  //Cherche cluster que id1 apparait
  unsigned int i=0;
  while (!stop)
  {
  unsigned int j=0;
  while((!stop) &&(j<grpRel_->fragmentation()->Cluster()[i].size() ))
  {
  if (grpRel_->fragmentation()->Cluster()[i][j]==id1) {stop=true;}
  if (!stop) j++;
  }
  if (!stop) i++;
  }

  cluster_id=i;

  for (unsigned int j=0; j<grpRel_->fragmentation()->Cluster()[cluster_id].size(); j++)
  if (grpRel_->fragmentation()->Cluster()[cluster_id][j]==id2) samecluster=true;

  //cout<<"same?:="<<samecluster<<endl;

  if( (!samecluster || (samecluster && (nwk_->inter(k)->rang()>1))) && (nwk_->inter(k)->Activate())  )
  {   
  //	cout<<"dAct:="<<nwk_->inter(k)->dAct()<<endl;
  nwk_->clist().push_back(k);
  nwk_->inter(k)->CDcoeff();
  nwk_->inter(k)->Kin();
  }
  else nwk_->inter(k)->clearForceAndMoment();
  }
  else*/
  for (unsigned int k=0 ; k < nwk_->linter().size() ; ++k)
  {	
    nwk_->inter(k)->dAct()= (grpRel_->getLaw(nwk_->inter(k)->first()->grp(),nwk_->inter(k)->second()->grp()))->dAct(nwk_->inter(k));
    nwk_->inter(k)->Frame();

    if (nwk_->inter(k)->Activate())
    {   
      nwk_->clist().push_back(k);
      if(nwk_->inter(k)->type() == 0 )
      {
//	cout<<"DKDK -> appel des bons coeffs ? 2"<<endl;
	nwk_->inter(k)->CDcoeff(grpRel_);
      }
      else
      {
	nwk_->inter(k)->CDcoeff();
      } 

      nwk_->inter(k)->Kin();
    }
    else nwk_->inter(k)->clearForceAndMoment();
  }
}

void CDalgo::fres()
{	  
  // External forces
  // Force-controled bodies are NOT subjected to acceleration field
  for (unsigned int i=0 ; i < sys_->lctrl().size() ; ++i)
  {	
    if(sys_->ctrl(i).x() == _FORCE)
      spl_->body(i)->fx()   = sys_->ctrl(i).xval();

    if(sys_->ctrl(i).y() == _FORCE)
      spl_->body(i)->fy()   = sys_->ctrl(i).yval();

    if(sys_->ctrl(i).rot() == _FORCE)
      spl_->body(i)->frot() = sys_->ctrl(i).rotval(); 
  } 

  sys_->computeAccField();//ne fait rien pour le moment
  // rq: il ne s'agit pas vraiment de n'importe quel 'champs d'accélération'.
  //     Par exemple, on ne peut pas avoir un champs convergent vers un point
  //     (au mieux on peut faire gx = cos(theta) et gy = sin(theta) sur chacunes des particules)
  // External forces of free bodies 
  for (unsigned int i=sys_->lctrl().size() ; i < spl_->lbody().size() ; ++i)
  {	
    if(spl_->body(i)->bodyDof()==NULL)
    {
      spl_->body(i)->fx()   =  sys_->gx() * spl_->body(i)->mass();//GRAVITY
      spl_->body(i)->fy()   =  sys_->gy() * spl_->body(i)->mass();//GRAVITY
      spl_->body(i)->frot() =  0.;
    }
  }

  for( unsigned int i=0;i < sys_->ldof().size();++i)
  {
    sys_->ldof(i)->imposeForce();
  }

  // Contact forces
  for (unsigned int c=0 ; c < nwk_->clist().size() ; ++c)
    nwk_->inter(nwk_->clist(c))->Res();

  // Force resultants on bloqued degrees of freedom
  for(unsigned int i=0;i<sys_->lctrl().size();++i)
  {
    if (sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
    if (sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
    if (sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;

  }
  for( unsigned int i=0;i < sys_->ldof().size();++i)
  {
    sys_->ldof(i)->imposeForceOfVelocity();
  }
}

// CETTE PROCEDURE QUI EST APPELEE DANS LE CODE!
// Gauss-Seidel iteration, original procedure
unsigned int CDalgo::iter_0()
{
	//cout<<"iter_0"<<endl;
	// No-contact case
	if (nwk_->clist().empty()) return 0;
	unsigned int c;

	// Variables
	double fnij0,anij;
	double ftij0,atij;
	double frotij0,asij;

	double muij;
	double ahij;
	double rsij; 

	double dfnij = 0.0, dftij = 0.0, dfrotij = 0.0;
	double fz1 = 0.0, fz2 = 0.0, dfz = 0.0;
	unsigned int nstop = 0;

	inter2d* oxo;

	unsigned int muId = grpRel_->getId("mu");
	//	unsigned int ahId = grpRel_->getId("ah");
	unsigned int rsId = grpRel_->getId("rs");
	//unsigned int clId = grpRel_->getIdLaw("");

	// Main CD loop
	for(unsigned int kiter=1 ; kiter <= nitermx_ ; ++kiter)
	{
		c = 0;
		while(c < nwk_->clist().size())
		{
			oxo = nwk_->inter(nwk_->clist(c));
			rsij = grpRel_->getParameterQuickly(rsId,oxo->first()->grp(),oxo->second()->grp());

			// Si on veut modéliser la fragmentation
			if (grpRel_->existfragmentation()) 
			{
				double fco,fri,ft;
				fco=grpRel_->fragmentation()->fco(oxo); // force cohésion
				fri=grpRel_->fragmentation()->frict(oxo); // coefficient de friction
				//frot=grpRel_->fragmentation()->frot(oxo);

				//le cas en collérant
				if ((fco!=0) && (fri!=0))
				{
					if( oxo->rang()==0 )//distance interaction
					{
						fnij0 = oxo ->fn(); 
						oxo->fn() = -fco;
						dfnij = oxo ->fn() - fnij0;
						// .......... On peut parfaitement rajouter le frottement

						// .......... Normal force sum
						fz2 += fabs(oxo->fn());

						// .......... Write new values of force resultants
						oxo->Res(dfnij,dftij,dfrotij);

						for(unsigned int i=0;i<sys_->lctrl().size();++i)
						{
							if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
							if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
							if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
						}

						if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
						if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();
					}
					else//contact
						for(unsigned int r=0 ; r < oxo->rang() ; ++r)
						{
							oxo->current() = r;

							// .......... Normal forces
							fnij0 = oxo->fn();
							anij = oxo->An(dt_);
							if(anij > -fco)
								oxo->fn() = anij;
							else
								oxo->fn() = -fco;

							dfnij = oxo->fn() - fnij0;

							// .......... Forces tangentiel

							ftij0 = oxo->ft();
							atij = oxo->At(dt_);
							ft=fri*oxo->fn()+grpRel_->fragmentation()->res_cisaillement(oxo);

							if (atij >= ft)
								oxo->ft() = ft;
							else if(atij <= -ft) 
								oxo->ft() = -ft;
							else if((atij > -ft) && (atij < ft)) 
								oxo->ft() = atij;

							dftij =  oxo->ft() - ftij0;

							// .......... Normal force sum
							fz2 += fabs(oxo->fn());

							// .......... Write new values of force resultants
							oxo->Res(dfnij,dftij,dfrotij);

							for(unsigned int i=0;i<sys_->lctrl().size();++i)
							{
								if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
								if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
								if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
							}

							if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
							if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();
						} //end of for r
				}
				else
				{ 
					ahij=(grpRel_->getLaw(oxo->first()->grp(),oxo->second()->grp()))->fco(oxo);//Cohésion
					muij = grpRel_->getParameterQuickly(muId,oxo->first()->grp(),oxo->second()->grp());//coefficience de friction
					cerr<<"Appel de mu : "<< oxo->first()->grp()<< " "<< oxo->second()->grp()<< " : "<<muij<<endl;

					if( oxo->rang()==0 )//distance interaction
					{
						fnij0 = oxo ->fn(); 
						oxo->fn() = -ahij;// - ahij---> - fco( oxo)
						dfnij = oxo ->fn() - fnij0;
						// .......... On peut parfaitement rajouter le frottement

						// .......... Normal force sum
						fz2 += fabs(oxo->fn());

						// .......... Write new values of force resultants
						oxo->Res(dfnij,dftij,dfrotij);

						for(unsigned int i=0;i<sys_->lctrl().size();++i)
						{
							if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
							if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
							if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
						}

						if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
						if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();
					}
					else//contact
						for(unsigned int r=0 ; r < oxo->rang() ; ++r)
						{
							oxo->current() = r;

							// .......... Normal forces
							fnij0 = oxo->fn();
							anij = oxo->An(dt_);
							if(anij > -ahij)
								oxo->fn() = anij;
							else
								oxo->fn() = -ahij;

							dfnij = oxo->fn() - fnij0;

							// .......... Friction forces

							if(muij != 0.0)
							{
								double mufn = muij * (oxo->fn() + ahij);  
								ftij0 = oxo->ft();

								atij = oxo->At(dt_);

								if (atij >= mufn)
									oxo->ft() = mufn;
								else if(atij <= -mufn) 
									oxo->ft() = -mufn;
								else if((atij > -mufn) && (atij < mufn)) // teste util ? a tester...
									oxo->ft() = atij;

								dftij =  oxo->ft() - ftij0;
							}

							// .......... Torques

							if(rsij != 0.0)
							{
								frotij0 = oxo->frot();

								asij = oxo->As(dt_);

								double rsfn = rsij * (oxo->fn() + ahij);

								if(asij >= rsfn)
								{
									oxo->frot() = rsfn;
								}

								if(asij <= -rsfn)
								{
									oxo->frot() = -rsfn;
								}

								if(asij>-rsfn && asij<rsfn)
								{ 
									oxo->frot() = asij;
								}

								dfrotij = oxo->frot() - frotij0;
							}

							// .......... Normal force sum
							fz2 += fabs(oxo->fn());

							// .......... Write new values of force resultants
							oxo->Res(dfnij,dftij,dfrotij);

							for(unsigned int i=0;i<sys_->lctrl().size();++i)
							{
								if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
								if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
								if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
							}

							if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
							if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();
						} //end of for r
				}
			}
			// Si on ne considere pas la fragmentation
			else
			{
				rsij = grpRel_->getParameterQuickly(rsId,oxo->first()->grp(),oxo->second()->grp());
				muij = grpRel_->getParameterQuickly(muId,oxo->first()->grp(),oxo->second()->grp());//coefficience de friction
				//cout<<"muij "<<oxo->first()->grp()<<" "<<oxo->second()->grp()<<" = "<<muij<<endl;
				ahij=(grpRel_->getLaw(oxo->first()->grp(),oxo->second()->grp()))->fco(oxo);//Cohésion

				//cerr<<"Appel de mu : "<< oxo->first()->grp()<< " "<< oxo->second()->grp()<< " : "<<muij<<endl;
				if( oxo->rang()==0 )//distance interaction
				{
					cout<<"rang = 0"<<endl;
					fnij0 = oxo ->fn(); 
					oxo->fn() = -ahij;// - ahij---> - fco( oxo)
					dfnij = oxo ->fn() - fnij0;
					//if(oxo->fn()!=0) cout<<"Cohésion à distance:="<<oxo->fn()<<endl;
					// .......... On peut parfaitement rajouter le frottement

					// .......... Normal force sum
					fz2 += fabs(oxo->fn());

					// .......... Write new values of force resultants
					oxo->Res(dfnij,dftij,dfrotij);

					for(unsigned int i=0;i<sys_->lctrl().size();++i)
					{
						if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
						if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
						if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
					}

					if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
					if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();
				}
				else//contact
					for(unsigned int r=0 ; r < oxo->rang() ; ++r)
					{
						oxo->current() = r;
						// .......... Normal forces

						fnij0 = oxo->fn();
						anij = oxo->An(dt_);
						if(anij > -ahij)
							oxo->fn() = anij;
						else
							oxo->fn() = -ahij;

						dfnij = oxo->fn() - fnij0;

						// .......... Friction forces

						if(muij != 0.0)
						{
							//cout<<"ahij="<<ahij<<endl;
							double mufn = muij * (oxo->fn() + ahij);  
							ftij0 = oxo->ft();

							atij = oxo->At(dt_);

							if (atij >= mufn)
								oxo->ft() = mufn;
							else if(atij <= -mufn) 
								oxo->ft() = -mufn;
							else if((atij > -mufn) && (atij < mufn)) // teste util ? a tester...
								oxo->ft() = atij;

							dftij =  oxo->ft() - ftij0;
							//On a verifie que si musij=0 alors oxo->ft() = 0 donc c'est ok.
						}		

						// .......... Torques
						if(rsij != 0.0)
						{
							frotij0 = oxo->frot();
							asij = oxo->As(dt_);
							//cout<<"asij:="<<asij<<endl;
							double rsfn = rsij * (oxo->fn() + ahij);

							if(asij >= rsfn)
							{
								oxo->frot() = rsfn;
							}

							if(asij <= -rsfn)
							{
								oxo->frot() = -rsfn;
							}

							if(asij>-rsfn && asij<rsfn)
							{ 
								oxo->frot() = asij;
							}
							dfrotij = oxo->frot() - frotij0;
						}

						// .......... Normal force sum
						fz2 += fabs(oxo->fn());

						// .......... Write new values of force resultants
						oxo->Res(dfnij,dftij,dfrotij);

						for(unsigned int i=0;i<sys_->lctrl().size();++i)
						{
							if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
							if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
							if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
						}

						if ( oxo->first()->bodyDof() != NULL ) oxo->first()->bodyDof()->imposeForceOfVelocity();
						if ( oxo->second()->bodyDof() != NULL ) oxo->second()->bodyDof()->imposeForceOfVelocity();

					}
			}
			++c;
		}

		// .......... Check convergence

		if( kiter >= nitermn_ )
		{
			if(fz2 != 0.0)
				dfz = fabs(fz2-fz1)/fz2;
			else
				dfz = 0.0;

			if(dfz < epsf_)
			{
				++nstop;
				if(nstop == niterconv_)
					return kiter;
			}
			else
				nstop = 0;
		}

		fz1 = fz2;
		fz2 = 0.0;
	}
	return nitermx_;
}

// Gauss-Seidel iteration, J J Moreau 
unsigned int CDalgo::iter_1()
{
	cout<<"Ham iter_1"<<endl;

	// No-contact case
	if (nwk_->clist().empty()) return 0;

	// Variables
	unsigned int verif,c;
	unsigned int totalcon = nwk_->clist().size();

	unsigned int it,k;
	double fnij0,anij;
	double ftij0,atij;

	double muij;
	double ahij;
	//double rsij; 

	double dfnij = 0.0, dftij = 0.0; 

	unsigned int muId = grpRel_->getId("mu");
	unsigned int ahId = grpRel_->getId("ah");

	it = 0;
	verif = 0;
	c = 0;

	while ( verif++ < 2 * totalcon )
	{
		if ( it++ >  nitermx_ ) break;
		if ( c    >= totalcon ) c = 0;

		k = nwk_->clist(c++);

		muij = grpRel_->getParameterQuickly(muId,nwk_->inter(k)->first()->grp(),nwk_->inter(k)->second()->grp());
		ahij = grpRel_->getParameterQuickly(ahId,nwk_->inter(k)->first()->grp(),nwk_->inter(k)->second()->grp());
		//rsij = 0.0; 


		for(unsigned int r=0 ; r < nwk_->inter(k)->rang() ; ++r)
		{
			nwk_->inter(k)->current() = r;  

			// .......... Normal forces

			fnij0 = nwk_->inter(k)->fn();

			anij = nwk_->inter(k)->An(dt_);

			if(anij > -ahij)
				nwk_->inter(k)->fn() = anij;
			else
				nwk_->inter(k)->fn() = -ahij; 

			dfnij = nwk_->inter(k)->fn() - fnij0;

			if( fabs(dfnij/fnij0)  > epsf_) verif = 0;

			// .......... Friction forces

			if(muij != 0.0)
			{
				cout<<"THERE"<<endl;
				double mufn = muij * (nwk_->inter(k)->fn() + ahij);
				ftij0 = nwk_->inter(k)->ft();

				atij = nwk_->inter(k)->At(dt_);

				if(atij >= mufn)
					nwk_->inter(k)->ft() = mufn;
				else if(atij <= -mufn) 
					nwk_->inter(k)->ft() = -mufn;
				else if((atij > -mufn) && (atij < mufn))
					nwk_->inter(k)->ft() = atij;

				dftij =  nwk_->inter(k)->ft() - ftij0;
			}

			nwk_->inter(k)->Res(dfnij,dftij,0.0);

			// .......... Bloqued degrees of freedom
			for(unsigned int i=0;i<sys_->lctrl().size();++i)
			{
				if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
				if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
				if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;        
			}
			if ( nwk_->inter(k)->first()->bodyDof() != NULL || nwk_->inter(k)->second()->bodyDof() != NULL )
				nwk_->inter(k)->first()->bodyDof()->imposeForceOfVelocity();

			/*for( unsigned int i=0;i < sys_->ldof().size();++i)
			  {
			  sys_->ldof(i)->imposeForceOfVelocity();
			  }*/
		} // r (current contact point)
	} // verif

	return it;
}

// Pour etude speciale
unsigned int CDalgo::iter_2()
{
	cout<<"Ham iter_2"<<endl;

	// No-contact case
	if (nwk_->clist().empty()) return 0;
	unsigned int c;

	// on efface toutes les forces normales
	inter2d* oxo;
	c = 0;
	while(c < nwk_->clist().size())
	{
		oxo = nwk_->inter(nwk_->clist(c));
		oxo->fn() = 0.0;
		++c;
	}

	// Variables
	double fnij0,anij;
	double ftij0,atij;
	double frotij0,asij;

	double muij;
	double ahij;
	double rsij; 

	double dfnij = 0.0, dftij = 0.0, dfrotij = 0.0;
	double fz1 = 0.0, fz2 = 0.0, dfz = 0.0;
	unsigned int nstop = 0;

	//char name[50];

	unsigned int muId = grpRel_->getId("mu");
	unsigned int ahId = grpRel_->getId("ah");
	unsigned int rsId = grpRel_->getId("rs");


	// Main CD loop
	for(unsigned int kiter=1 ; kiter <= nitermx_ ; ++kiter)
	{
		c = 0;

		while(c < nwk_->clist().size())
		{
			oxo = nwk_->inter(nwk_->clist(c));

			muij = grpRel_->getParameterQuickly(muId,oxo->first()->grp(),oxo->second()->grp());
			ahij = grpRel_->getParameterQuickly(ahId,oxo->first()->grp(),oxo->second()->grp());
			rsij = grpRel_->getParameterQuickly(rsId,oxo->first()->grp(),oxo->second()->grp());

			for(unsigned int r=0 ; r < oxo->rang() ; ++r)
			{
				oxo->current() = r;

				// .......... Normal forces

				fnij0 = oxo->fn();

				anij = oxo->An(dt_);

				if(anij > -ahij)
					oxo->fn() = anij;
				else
					oxo->fn() = -ahij;

				dfnij = oxo->fn() - fnij0;

				// .......... Friction forces

				if(muij != 0.0)
				{
					double mufn = muij * (oxo->fn() + ahij);
					ftij0 = oxo->ft();

					atij = oxo->At(dt_);

					if (atij >= mufn)
						oxo->ft() = mufn;
					else if(atij <= -mufn) 
						oxo->ft() = -mufn;
					else if((atij > -mufn) && (atij < mufn)) // teste util ? a tester...
						oxo->ft() = atij;

					dftij =  oxo->ft() - ftij0;
				}


				// .......... Torques

				if(rsij != 0.0)
				{
					frotij0 = oxo->frot();

					asij = oxo->As(dt_);

					double rsfn = rsij * (oxo->fn() + ahij);

					if(asij >= rsfn)
					{
						oxo->frot() = rsfn;
					}

					if(asij <= -rsfn)
					{
						oxo->frot() = -rsfn;
					}

					if(asij>-rsfn && asij<rsfn)
					{ 
						oxo->frot() = asij;
					}


					dfrotij = oxo->frot() - frotij0;
				}

				// .......... Bloqued degrees of freedom
				for(unsigned int i=0;i<sys_->lctrl().size();++i)
				{
					if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->fx()   = 0.0;
					if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->fy()   = 0.0;
					if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->frot() = 0.0;
				}

				// .......... Write new values of force resultants
				oxo->Res(dfnij,dftij,dfrotij);

				// .......... Normal force sum
				fz2 += fabs(oxo->fn());     

			}

			++c;
		}

		// .......... Sauvegarde des forces

		//sprintf(name,"nwk_%03d.his",kiter);
		// history_write(kiter, *spl_, *nwk_, false, true);
		//network_write(name,*nwk_);

		// .......... Check convergence
		if(fz2 != 0.0)
			dfz = fabs(fz2-fz1)/fz2;
		else
			dfz = 0.0;

		if(dfz < epsf_)
		{
			++nstop;
			if(nstop == nitermn_)
				return kiter;
		}
		else
			nstop = 0;

		fz1 = fz2;
		fz2 = 0.0;
	}

	return nitermx_;
}

// Update positions and velocities for one step
void CDalgo::step()
{
	//cout<<"step"<<endl;
	unsigned int N = spl_->lbody().size();

	// Update externally imposed velocities and forces (enslaved)
	sys_->drive(); 

	// Update enslaved velocities
	sys_->share(); 

	// Update positions	of all bodies
	for (unsigned int i=0 ; i<N ; ++i)
	{
		if( spl_->body(i)->bodyDof()==NULL)
		{
			spl_->body(i)->x()   += spl_->body(i)->vx()   * dt_;
			spl_->body(i)->y()   += spl_->body(i)->vy()   * dt_;
			spl_->body(i)->rot() += spl_->body(i)->vrot() * dt_;
		}
	}
	for( unsigned int i=0;i < sys_->ldof().size();++i)
	{
		sys_->ldof(i)->move(dt_);
	}

	// Set velocities of controlled bodies
	for (unsigned int i=0 ; i<sys_->lctrl().size() ; i++)
	{
		if(sys_->ctrl(i).x()   == _VELOCITY) spl_->body(i)->vx()   = sys_->ctrl(i).xval();
		if(sys_->ctrl(i).y()   == _VELOCITY) spl_->body(i)->vy()   = sys_->ctrl(i).yval();
		if(sys_->ctrl(i).rot() == _VELOCITY) spl_->body(i)->vrot() = sys_->ctrl(i).rotval();
	}
	for( unsigned int i=0;i < sys_->ldof().size();++i)
	{
		sys_->ldof(i)->imposeVelocity();
	}

}

// Increment time, Update Verlet list (if necessary)
void CDalgo::hand(unsigned int ns)
{  
	//	cout<<"hand"<<endl;
	if( nwk_->useSuperList() )//Si superliste
	{

		if (ns % nsuper_ == 0  ) //si maj super liste
		{ 
			cout<<" ------ Update SuperList ------ "<<endl;
			//	Nmajsuper_=0;
			if( spl_->isMonoPeriodic() ) 
			{
				//cout<<"		Nmajsuper = "<<Nmajsuper_<<" "<<tm_<<"---------------------"<<endl;
				spl_->translate();
				spl_->updateBands();
				nwk_->buildSuperList(spl_, grpRel_);
				nwk_->buildSuperListP(spl_, grpRel_);
			}
			else
			{
				nwk_->buildSuperList(spl_, grpRel_);
			}

			nwk_->stock(forcesSafe_);

			nwk_->verlet(spl_, grpRel_);

			if (spl_->isMonoPeriodic()) 
			{    
				nwk_->verletP(spl_, grpRel_);
			}
			contact();
			nwk_->retrieve(forcesSafe_);



		}
		else //Pas maj super mais maj verlet
		{
			if ( ns % nver_ == 0  )
			{
				//cout<<" ------ Use of superList : Update Verlet = "<<ns<<" ------ "<<endl;
				//Nmaj_=0;

				nwk_->stock(forcesSafe_);

				nwk_->verlet(spl_, grpRel_);

				if (spl_->isMonoPeriodic()) 
				{    
					nwk_->verletP(spl_, grpRel_);
				}

				nwk_->retrieve(forcesSafe_);
				contact();


			}


		}
	}

	else 
	{
		if ( ns % nver_ == 0  )
		{
			cout<<" ------ No use of superList : Update Verlet = "<<ns<<" ------ "<<endl;
			//Nmaj_=0;


			nwk_->stock(forcesSafe_);

			if (spl_->isMonoPeriodic() ) 
			{ 
				spl_->translate();
				spl_->updateBands();
				//	cout<<"av verlet"<<endl;
			}

			nwk_->verlet(spl_, grpRel_);
			//	cout<<"ap verlet"<<endl;

			if (spl_->isMonoPeriodic()) 
			{    
				nwk_->verletP(spl_, grpRel_);
				//	cout<<"ap verletP"<<endl;

			}

			contact();
			nwk_->retrieve(forcesSafe_);
			//	cout<<"apres cont"<<endl;
		}
	}
}

// Initialize the computation
void CDalgo::stand()
{
	assert(sys_ != 0);
	assert(grpRel_ != 0);

	// Check existence of the required parameters
	if(!grpRel_->exist("mu"))
		cerr << "@CDalgo::stand, required parameter 'mu' does not exist" << endl << flush;
	//	if(!grpRel_->exist("ah"))
	//		cerr << "@CDalgo::stand, required parameter 'ah' does not exist" << endl << flush;
	//if(!grpRel_->exist("en"))
	//cerr << "@CDalgo::stand, required parameter 'en' does not exist" << endl << flush;
	//if(!grpRel_->exist("et"))
	//cerr << "@CDalgo::stand, required parameter 'et' does not exist" << endl << flush;
	if(!grpRel_->exist("rs"))
		cerr << "@CDalgo::stand, required parameter 'rs' does not exist" << endl << flush;
	//	if(!grpRel_->exist("dact"))
	//		cerr << "@CDalgo::stand, required parameter 'dact' does not exist" << endl << flush;


	// Take into account system data  
	//	sys_->init();
	sys_->trans();
	sys_->share();

	// Compute mass of bodies
	double density = 2500.0;

	if(grpDat_->exist("density"))
	{
		for (unsigned int i=0;i<spl_->lbody().size();++i)
		{
			density = grpDat_->getParameter("density",spl_->body(i)->grp());
			spl_->body(i)->Fill(density);
		}


	}
	else
	{
		spl_->fill(2500.0);
		gdm::warning("@CDalgo::stand, GroupData density not defined (so it is set to 2500)");
	}

	for(unsigned int i=0;i<sys_->ldof().size();++i)
	{
		sys_->ldof(i)->fill(density);
		sys_->ldof(i)->setGravity(sys_->gx(),sys_->gy() );
		sys_->ldof(i)->computeVertex();
		sys_->ldof(i)->computeMassCenter();
		sys_->ldof(i)->computeMoment();
		//sys_->ldof(i)->print();
	}


	// Set up the initial contact list and the force list
	if (nwk_->useSuperList()) 
		nwk_->buildSuperList(spl_, grpRel_);

	if( nwk_->linter().empty()) 
		nwk_->verlet(spl_, grpRel_);
	else
	{
		nwk_->stock(forcesSafe_);
		nwk_->verlet(spl_, grpRel_);
		contact();
		nwk_->retrieve(forcesSafe_);		
	}

	if (spl_->isMonoPeriodic()) 
	{ 
		if (nwk_->useSuperList()) nwk_->buildSuperListP(spl_, grpRel_);
		nwk_->verletP(spl_, grpRel_);
	}



	//contact();	
}

unsigned int CDalgo::contactStatut( inter2d * inter)
{

	//retourne le statut d'un contact
	//0 contact frottant non glissant
	//1 contact frottant glissant
	unsigned int muId = grpRel_->getId("mu");
	//unsigned int ahId = grpRel_->getId("ah");
	//unsigned int rsId = grpRel_->getId("rs");

	double muij = grpRel_->getParameterQuickly(muId,inter->first()->grp(),inter->second()->grp());
	//double ahij = grpRel_->getParameterQuickly(ahId,inter->first()->grp(),inter->second()->grp());
	//double rsij = grpRel_->getParameterQuickly(rsId,inter->first()->grp(),inter->second()->grp());

	if (muij != 0.)
	{
		if( fabs(inter->ft()) >= muij*inter->fn()*(.99))
			return (1);
		else
			return (0);
	}


	return (0);
}

//Mise à jour le table cohésion

void CDalgo::resistance()
{
	// No-contact case
	if (nwk_->clist().empty()) return;

	//Vérifier le résistance et mise à jour le table cohésion
	inter2d* oxo;
	for(unsigned int i=0; i<nwk_->clist().size(); i++)
	{
		oxo = nwk_->inter(nwk_->clist(i));

		unsigned int first=oxo->first()->id()-4;
		unsigned int second=oxo->second()->id()-4;

		//Si fracture
		if( (grpRel_->fragmentation()->fco(oxo)!=0) && grpRel_->fragmentation()->fracture(oxo))
		{
			unsigned int index =first/grpRel_->fragmentation()->nsite();
			grpRel_->fragmentation()->numfissure()++;
			grpRel_->fragmentation()->lgrain()[index]=true;
			cout<<"first:="<<first<<"	second:="<<second<<endl;
			for(unsigned int r=0;r < oxo->rang();++r)
			{
				oxo->current() = r;
				cout<<"current:="<<oxo->current()<<"	fn("<<r<<"):="<<oxo->fn()<<"		ft("<<r<<"):="<<oxo->ft()<<endl;
			}

			unsigned int j=0;
			bool logic=true;		
			do
			{
				if (grpRel_->fragmentation()->ActiveCohesion()[first][j]==second) 
				{
					grpRel_->fragmentation()->ActiveCohesion()[first].erase(grpRel_->fragmentation()->ActiveCohesion()[first].begin()+j);
					logic=false;
				}
				j++;
			}
			while (logic);

			j=0;
			logic=true;					
			do
			{
				if (grpRel_->fragmentation()->ActiveCohesion()[second][j]==first) 
				{			
					grpRel_->fragmentation()->ActiveCohesion()[second].erase(grpRel_->fragmentation()->ActiveCohesion()[second].begin()+j);
					logic=false;
				}
				j++;
			}
			while (logic);	
		}		
	}

	//Mise à jour le list Cluster
	unsigned int Ncluster=grpRel_->fragmentation()->Cluster().size();
	unsigned int index2;
	vector< vector<unsigned int> > NewCluster;
	vector<unsigned int> ListPolygone;	//Cluster temporaire
	vector<unsigned int> cluster;		//Cluster initial

	NewCluster.clear();
	for (unsigned int i=0; i<Ncluster; i++)
	{
		cluster=grpRel_->fragmentation()->Cluster()[i];

		while(!cluster.empty())
		{
			unsigned int first,second;
			inter2d * oxo;

			ListPolygone.clear();
			ListPolygone.push_back(cluster[0]);
			cluster.erase(cluster.begin());

			for (unsigned int j=0; j<ListPolygone.size(); j++)
			{
				first=ListPolygone[j]+4;

				for(unsigned int k=0; k<cluster.size(); k++)
				{
					second=cluster[k]+4;
					for(unsigned int l=0;l<nwk_->clist().size();l++)
					{
						oxo = nwk_->inter(nwk_->clist(l));
						if ( ((oxo->first()->id()==first && oxo->second()->id()==second) || (oxo->first()->id()==second && oxo->second()->id()==first)) && (oxo->rang()==2) )
							ListPolygone.push_back(cluster[k]);
					}
				}

				//Supprimer les particules dans cluster initial que on a ajouté à ListPolygone

				for (unsigned int k=0; k<ListPolygone.size(); k++)
				{
					unsigned int l=0;
					bool stop=false;

					while((l<cluster.size()) && (!stop))
					{
						if (cluster[l]==ListPolygone[k]) {cluster.erase(cluster.begin()+l); stop=true;}
						else l++;
					}	
				}	
			}
			NewCluster.push_back(ListPolygone);	
		}		
	}
	grpRel_->fragmentation()->Cluster()=NewCluster;

	for(unsigned int i=0;i<NewCluster.size();i++)
		if (NewCluster[i].size()<grpRel_->fragmentation()->nsite())
		{
			index2=NewCluster[i][0]/grpRel_->fragmentation()->nsite();
			grpRel_->fragmentation()->Casse()[index2]=true;
		}
}

