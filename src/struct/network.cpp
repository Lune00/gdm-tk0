#include "network.hpp"

void Network::associate(Sample & spl)
{
	if (iNum_.empty() || jNum_.empty()) return;

	unsigned int N = iNum_.size();
	for (unsigned int n = 0 ; n < N ; ++n)
	{
		linter_[n]->plug(spl.body(iNum_[n]),spl.body(jNum_[n]));
	}
	clearTmpBodyNum();
}
// attention spl doit deja exister
void Network::read(istream & is/*, Sample & spl*/)
{

	cout<<"Read network : "<<endl;

	string token;
	unsigned int n = 0;
	unsigned int Id1, Id2;
	inter2d * o_o = 0;
	clearTmpBodyNum();
	is >> token;	
	while(is)
	{
		if      (token == "dverlet") {is >> dverlet_;}
		else if (token == "dsuperList") 
		{
			is >> dsuperList_;
			//is >> trash;cout<<"read of dsuperlist deactivated by CV....."<<endl;
			useSuperList_ = true;
		}
		else if (token == "speakLevel") is >> speakLevel_;
		else if (token == "}") break;		
		else
		{


			o_o = inter2d::factory(token);
			if (o_o != 0)
			{
				linter_.push_back(o_o);
				linter_[n]->read(is,&Id1,&Id2);
				iNum_.push_back(Id1);
				jNum_.push_back(Id2);
				//cout<<o_o->fn()<<" "<<Id2<<endl;
				clist_.push_back(n);
				++n;
			}
		}

		is >> token;
	}


	cout<<"Nombre de contacts chargé = "<<n<<endl;
	cout<<"Taille de clist_ = "<<clist_.size()<<endl;
}
void Network::write(const char * fname)
{
	ofstream fnetwork(fname);

	fnetwork << "Network{" << endl;

	for (unsigned int i=0;i<linter_.size();++i)
		if (linter_[i]->fn() != 0.0) linter_[i]->write(fnetwork);
	fnetwork << "}" << endl;
}
void Network::write(ostream & os)
{	
	os << "dverlet " << dverlet_ << endl;
	if (useSuperList_)
		os << "dsuperList " << dsuperList_ << endl;
	double fn;
	for (unsigned int i=0;i<linter_.size();++i)
	{
		fn=0.;
		for(unsigned int r=0;r<linter_[i]->rang();r++)
		{
			linter_[i]->current()=r;
			fn+=linter_[i]->fn();
		}
		linter_[i]->current()=0;
		if (fn!=0.0) linter_[i]->write(os);
	}
}
void Network::speak()
{
	switch(speakLevel_)
	{
		case 0:
			return;

		default:
			cout << "Number of potential interactions: " << linter_.size() << endl;
			cout << flush;
			return;
	}
}
// Stock all the forces (in memory)
void Network::stock(vector<fsafe>& fs)
{
	fs.clear();

	for (unsigned int k=0 ; k < linter_.size() ; ++k)
	{
		linter_[k]->current()=0;
		fs.push_back(fsafe(linter_[k]));
	}
	/*
	   for (unsigned int k=0 ; k < clist_.size() ; ++k)
	   {
	   linter_[clist_[k]]->current()=0;
	   fs.push_back(fsafe(linter_[clist_[k]]));
	   }
	 */

	//	cout<<"-----------fin stockage: "<<fs.size()<<" interactions stockee------------"<<endl;
}

// Stock the forces belonging to a special list
void Network::stock(vector<fsafe>& fs, vector<unsigned int>& sublist)
{
	fs.clear();

	for (unsigned int k=0 ; k < sublist.size() ; ++k)
	{

		fs.push_back(fsafe(linter_[sublist[k]]));
	}
}


// Retrieve the forces previously stocked
void Network::retrieve(vector<fsafe> & fs)
{
	if (fs.empty() || linter_.empty()) return;

	unsigned int i=0;
	unsigned int Nnew = linter_.size();
	unsigned int Nold = fs.size();

	unsigned int NdkdkOld=0;//indice de la premiere interaction periodique dans fs
	unsigned int NdkdkNew=0;//indice de la premiere interaction periodique dans linter

	while( NdkdkNew < Nnew && linter_[NdkdkNew]->type() != _type_dkdkP   )
	{ 
		//cout<<linter_[NdkdkNew]->type()<<endl;
		++NdkdkNew;
	}
	// cout<<" NdkdkNew "<<NdkdkNew<<endl;
	while( NdkdkOld < Nold && fs[NdkdkOld].type() != _type_dkdkP) 
	{
		//cout<<fs[NdkdkOld].type()<<endl;
		++NdkdkOld;
	}
	++NdkdkNew;
	++NdkdkOld;
	//cout<<" Nold "<<Nold<<" Nnew= "<<Nnew<<endl;
	//cout<<" NdkdkOld = "<<NdkdkOld<<" NdkdkNew "<<NdkdkNew<<endl; //Plus un pour la taille du tableau
	/*for (unsigned j=min(NdkdkNew,NdkdkOld)-10;j<max(Nnew,Nold);++j)
	  {
	  if(j<Nnew)
	  cout<<" New : first "<<linter_[j]->lexifirst()->id()<<"  second "<<linter_[j]->lexisecond()->id();

	  if(j<Nold)
	  cout<<"---- Old : first "<<fs[j].first()->id()<<"  second "<<fs[j].second()->id();
	  cout<<" j = "<<j;
	  if (j==NdkdkNew) cout<<" == NdkdkNew ";
	  if (j==NdkdkOld) cout<<" == NdkdkOld ";
	  cout<<endl;
	  }*/

	for (unsigned int k=0 ; k < Nnew ; ++k)
	{	

		if( i == NdkdkOld && k<NdkdkNew ) k = NdkdkNew;
		if( k == NdkdkNew && i<NdkdkOld ) i = NdkdkOld;

		while( i < Nold && fs[i].first()->id() < linter_[k]->lexifirst()->id()) ++i;
		if( i == NdkdkOld && k<NdkdkNew) continue;
		if (i == Nold) break;

		while( i < Nold
				&& fs[i].first()->id()  == linter_[k]->lexifirst()->id()  
				&& fs[i].second()->id() <  linter_[k]->lexisecond()->id()) ++i;
		if( i == NdkdkOld && k<NdkdkNew) continue;
		if (i == Nold) break;

		/*	if (fs[i].type() != linter_[k]->type() ) 
			{
			cout<<"type different "<< linter_[k]->lexifirst()->id() <<" et "<< linter_[k]->lexisecond()->id()<<endl;
			cout<<" 	type safe"<<fs[i].type()<<" linter.type "<<linter_[k]->type()<<endl;
			}*/

		if(    fs[i].first()->id()  == linter_[k]->lexifirst()->id() 
				&& fs[i].second()->id() == linter_[k]->lexisecond()->id()	)
		{
			//linter_[k]->rang() = fs[i].rang();
			linter_[k]->current()=0;
			linter_[k]->fn() = fs[i].fn();
			linter_[k]->ft() = fs[i].ft();


			if( fs[i].rang()==2)
			{
				//   cout<<"retrieve double"<<endl;
				linter_[k]->rang()=2;
				linter_[k]->current()=1;
				linter_[k]->fn()=fs[i].fn2();
				linter_[k]->ft()=fs[i].ft2();
				linter_[k]->current()=0;

				//	cout<<" interaction double retrouvee entre "<< linter_[k]->lexifirst()->id() <<" et "<< linter_[k]->lexisecond()->id()<<" ";
				//	cout<<" fn1= "<<linter_[k]->fn()<<" ft1= "<<linter_[k]->ft()<<endl;
				//	linter_[k]->current()=1;
				//	cout<<" fn2= "<<linter_[k]->fn()<<" ft2= "<<linter_[k]->ft()<<endl;
				//	linter_[k]->current()=0;
			}
			/*	if(linter_[k]->lexifirst()->id()<100 )//&& linter_[k]->lexifirst()->id()>2
				{
				cout<<" interaction simple retrouvee entre "<< linter_[k]->lexifirst()->id() <<" et "<< linter_[k]->lexisecond()->id()<<" ";
				cout<<" type safe "<<fs[i].type()<<" type inter "<<linter_[k]->type();
				cout<<" fn1= "<<linter_[k]->fn()<<" ft1= "<<linter_[k]->ft()<<endl;
				}
			 */
			//if(k>min(NdkdkNew,NdkdkOld)-10)
			++i;
		}
		//else cout<<" interaction non retrouvee"<<endl;

	}
}


void Network::buildSuperList(Sample * spl,  GroupRelationData * grpRel)
{
	//cerr<<"---Building SuperList"<<endl;
	superList_.clear();//Reset le list des pairs particules 

	unsigned int N=spl->lbody().size();
	unsigned int i,j;

	struct particle_pair O___O;
	for (i=0 ; i<N ; ++i)
	{
		for (j=i+1 ; j<N ; ++j)
		{
			if (grpRel->act(spl->body(i)->grp(), spl->body(j)->grp()))
			{
				if (near(spl->body(i), spl->body(j), dsuperList_))
				{
					O___O.i = spl->body(i);
					O___O.j = spl->body(j);
					superList_.push_back(O___O);
				}
			}
		}
	}
}

// Build the list of potential periodics interactions
void Network::buildSuperListP(Sample* spl, GroupRelationData* grpRel)
{
	//cerr<<"---Building SuperlistP"<<endl;
	//Il faut s'assurer que la largeur de bande est au moins egale a la dsuperlistP
	superListP_.clear();
	unsigned int NL =spl->leftband().size();
	unsigned int NR =spl->rightband().size();
	double P = spl->rightBoundary() - spl->leftBoundary();
	unsigned int i,j;
	unsigned int r,l;

	body2d* ghost=NULL;
	struct particle_pair O___O;
	O___O.i = NULL ;
	O___O.j = NULL ;

	for (i = 0 ; i < NR ; ++i)
	{
		for (j = 0 ; j < NL ; ++j)
		{
			r = spl->rightband(i);
			l = spl->leftband(j);

			if (grpRel->act(spl->body(r)->grp(), spl->body(l)->grp()))
			{
				ghost = spl->body(l)->duplicate();
				ghost->x() += P;

				//if (near(spl->body(r), ghost, dVerletVar(min( ghost->sizeVerlet(),spl->body(r)->sizeVerlet())) ))
				if (near(spl->body(r), ghost, dsuperListP_ ))
				{
					O___O.i = spl->body(r);
					O___O.j = spl->body(l);
					superListP_.push_back(O___O);
				}

				delete ghost;
			}
		}
	}
}


// Build the list of potential interactions
void Network::verlet(Sample * spl, GroupRelationData * grpRel)
{
	purge(linter_);
	linter_.clear();// linter_ declaree en tant que donnee membre
//cerr<<"---Building verlet list"<<endl;
	if (useSuperList_)
	{
		unsigned int k;
		unsigned int N = superList_.size();
		inter2d* o_o = NULL;

		for (k=0 ; k<N ; ++k)
		{
			if (grpRel->act(superList_[k].i->grp(), superList_[k].j->grp()))
			{
				if (near(superList_[k].i, superList_[k].j, dverlet_))
				{
					o_o = inter2d::factory(superList_[k].i, superList_[k].j);
					if(o_o != 0) 
					{
						linter_.push_back(o_o);
					}
				}
			}
		}
	}
	else
	{
		unsigned int N=spl->lbody().size();
		unsigned int i,j;
		inter2d* o_o = NULL;

		for (i=0 ; i<N ; ++i)
		{
			for (j=i+1 ; j<N ; ++j)
			{
				if (grpRel->act(spl->body(i)->grp(), spl->body(j)->grp()))
				{
					if (near(spl->body(i), spl->body(j), dverlet_))
					{
						o_o = inter2d::factory(spl->body(i), spl->body(j));
						if(o_o != 0) 
						{
							linter_.push_back(o_o);
						}
					}
				}
			}
		}

	}

}

// Build the list of potential periodics interactions
void Network::verletP(Sample* spl, GroupRelationData* grpRel)
{
	//cerr<<"Building verletP"<<endl;
	unsigned int NL =spl->leftband().size();
	unsigned int NR =spl->rightband().size();
	double P = spl->rightBoundary() - spl->leftBoundary();
	unsigned int i,j;
	unsigned int r,l;

	inter2d* o_o = NULL;
	body2d* ghost = NULL;

	if (useSuperList_)
	{
		unsigned int k;
		unsigned int N = superListP_.size();

		for (k=0 ; k<N ; ++k)
		{
			if (grpRel->act(superListP_[k].i->grp(), superListP_[k].j->grp()))
			{
				ghost = superListP_[k].j->duplicate();
				ghost->x() += P;

				if (near(superListP_[k].i, ghost, dverlet_))
				{
					o_o = inter2d::factoryP(superListP_[k].i, superListP_[k].j,P);
					if(o_o != 0) 
					{
						linter_.push_back(o_o);
					}
				}
				delete ghost;
			}
		}

	}
	else
	{
		for (i = 0 ; i < NR ; ++i)
		{
			for (j = 0 ; j < NL ; ++j)
			{
				r = spl->rightband(i);
				l = spl->leftband(j);

				if (grpRel->act(spl->body(r)->grp(), spl->body(l)->grp()))
				{
					ghost = spl->body(l)->duplicate();
					ghost->x() += P;

					if (near(spl->body(r), ghost, dverlet_ ))
					{
						o_o = inter2d::factoryP(spl->body(r), spl->body(l), P);
						if(o_o != 0) 
						{
							//link->Frame(); // indispendsable ???
							linter_.push_back(o_o);
						}
					}

					delete ghost;
				}
			}
		}
	}
}

