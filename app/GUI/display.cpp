#include "display.h"
#include <math.h>



const float _pi     = 3.14159265358979323846;
const float _pi_2   = 1.57079632679489661923;
const float _pi_8   = 0.39269908169872415480;
const float _15pi_8 = 5.89048622548086232200;
const float _3pi_2  = 4.71238898038468985769;
const float _2pi    = 6.28318530717958647692;


Display::Display(int x,int y,int w,int h,const char *l)
	: Fl_Gl_Window(x,y,w,h,l)
{
	_vAng = 0.0;
	_hAng = 0.0;

//spl.updateBoundaries();

	_size   = 1.0;
	_xshift = 0.0;
	_yshift = 0.0;

	_scalf = 1.0;
	_scalv = 1.0;

	_disp_shapes    = true;
	_disp_vel       = false;
	_disp_orient    = false;
	_disp_cluster   = true;

	_disp_fn        = true;
	_disp_fn_pos    = true;
	_disp_fn_neg    = true;

//modif le 09/03/09 
	_disp_fresNormal	= true;
	//_disp_fresNormal_pos	= true;
//---------------------------------	
	_disp_ft        = false;
	_disp_verlet    = false;
	_disp_locFrame  = false; 

	ifstream fconfig("./display.opt");
	if(fconfig != NULL) read_parameters(fconfig);
}

void Display::read_parameters(istream & is)
{
	string token;

	is >> token;	
	while(is)
	{	
		if      (token == "shapes")   is >> _disp_shapes;
		else if (token == "vel")      is >> _disp_vel;
		else if (token == "orient")   is >> _disp_orient;
		else if (token == "cluster")  is >> _disp_cluster;
		else if (token == "fn")       is >> _disp_fn;
		else if (token == "ft")       is >> _disp_ft;
		else if (token == "verlet")   is >> _disp_verlet;
		else if (token == "locFrame") is >> _disp_locFrame;

		else if (token == "size")    is >> _size;
		else if (token == "xshift")  is >> _xshift;
		else if (token == "yshift")  is >> _yshift;

		else if (token == "scalf") is >> _scalf;
		else if (token == "scalv") is >> _scalv;
		
//modif le 09/03/09
		else if (token == "fresNormal")		is >> _disp_fresNormal;
//-----------------------------------------------------------------		
		//else if (token == "nwiter") is>> _nwiter;

		else if (token == "}")        break;
		else cerr << "@Display::read_parameters, Unknown parameter: " << token << endl;

		is >> token;
	}
}

void Display::write_parameters(ostream & os)
{
	os << "shapes   " << _disp_shapes   << endl;
	os << "vel      " << _disp_vel      << endl;
	os << "orient   " << _disp_orient   << endl;
	os << "cluster  " << _disp_cluster  << endl;
	os << "fn       " << _disp_fn       << endl;
	os << "ft       " << _disp_ft       << endl;
	os << "verlet   " << _disp_verlet   << endl;
	os << "locFrame " << _disp_locFrame << endl;
//modif le 09/03/09
	os << "fresNormal	" << _disp_fresNormal << endl;
//-----------------------------------------------------
	os << endl;
	os << "size     " << _size          << endl;
	os << "xshift   " << _xshift        << endl;
	os << "yshift   " << _yshift        << endl;
	os << endl;
	os << "scalf    " << _scalf         << endl;
	os << "scalv    " << _scalv         << endl;
}

void Display::drawSystem()
{
	double x,y;
	if( _disp_cluster)
	{
	for(unsigned int i=0;i<sim()->sys()->ldof().size();++i)
	{
		x=sim()->sys()->ldof(i)->mcx();
		y=sim()->sys()->ldof(i)->mcy();
		glLineWidth(3.0f);
		glBegin(GL_LINES);
		glColor3f (0.0f, .5f, .5f);
		for(unsigned int j=0;j<sim()->sys()->ldof(i)->lctrlBodies().size();++j)
		{

			glVertex3f(x , y, 0.0f);
			glVertex3f(sim()->sys()->ldof(i)->ctrlBody(j)->x(),
				sim()->sys()->ldof(i)->ctrlBody(j)->y(), 0.0f);


		}
		glEnd();
		glBegin(GL_POLYGON);
		double L = sim()->sys()->ldof(i)->ctrlBody(0)->sizeVerlet();
		for (float t=0.;t<_2pi;t+=_pi_8)
		{
			glVertex3f (x+L*cos(t),
				y+L*sin(t), 0.0f);
		}
		glEnd();
		
		
	}
	}

	if (typeid(*sim()->sys()) == typeid( shearP_CD ))
	{
		if ( sim()->doAnalyse() && typeid(*sim()->sysA()) == typeid(shearP_CD_A))
		{
			body2d * bodyi;
			double L;
			shearP_CD_A * shearSysA = dynamic_cast<shearP_CD_A*> ( sim()->sysA());

			glLineWidth(1.0f);
			glBegin(GL_LINES);
			glColor3f (0.5f, 0.0f, 1.0f);
			glVertex3f (sim()->spl()->xmin(),shearSysA->totalProbe().h1(), 0.0f);
			glVertex3f (sim()->spl()->xmax(),shearSysA->totalProbe().h1(), 0.0f);
			glEnd();
			glLineWidth(1.0f);
			glBegin(GL_LINES);
			glColor3f (0.5f, 0.0f, 1.0f);
			glVertex3f (sim()->spl()->xmin(),shearSysA->totalProbe().h2(), 0.0f);
			glVertex3f (sim()->spl()->xmax(),shearSysA->totalProbe().h2(), 0.0f);
			glEnd();

			bodyi=shearSysA->Partref();
			L = (bodyi->xmax() - bodyi->x()) *0.8;
			glColor4f (1.0f, 1.0f, 0.0f, 0.2f);
			glBegin(GL_POLYGON);
			for (float t=0.;t<_2pi;t+=_pi_8)
			{
				glVertex3f (bodyi->x()+L*cos(t),
					bodyi->y()+L*sin(t), 0.0f);
			}
			glEnd();

			if( shearSysA->ngap1() != NULL && shearSysA->ngap2() != NULL)
			{
				bodyi=shearSysA->ngap1();
				L = (bodyi->xmax() - bodyi->x())*.7 ;
				glColor4f (1.0f, 0.0f, 0.0f, 1.f);
				glBegin(GL_POLYGON);
				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+L*cos(t),
						bodyi->y()+L*sin(t), 0.0f);
				}
				glEnd();


				bodyi=shearSysA->ngap2();
				L = (bodyi->xmax() - bodyi->x())*.7 ;
				glColor4f (1.0f, 0.0f, 0.0f, 1.0f);
				glBegin(GL_POLYGON);
				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+L*cos(t),
						bodyi->y()+L*sin(t), 0.0f);
				}
				glEnd();
			}

		}//sys_A	
	}//sys
	
	else if (typeid(*sim()->sys()) == typeid(Biaxial))
	{

		if ( sim()->doAnalyse() && typeid(*sim()->sysA()) == typeid(Biaxial_A))
		{
			body2d * bodyi;

			Biaxial_A * biaxA = dynamic_cast<Biaxial_A*> ( sim()->sysA());

			circularProbe * cProbe= dynamic_cast<circularProbe* > (& biaxA->totalProbe());
			double  L = cProbe->R();
			glColor3f (1.0f, 0.0f, 0.0f);
			glBegin(GL_LINE_LOOP);
			for (float t=0.;t<_2pi;t+=_pi_8)
	          {
	          glVertex3f (cProbe->x()+L*cos(t),
	                      cProbe->y()+L*sin(t), 0.0f);
	          }
	        glEnd();
	
			rectangularProbe * rProbe= dynamic_cast<rectangularProbe* > (& biaxA->rProbe());
			glColor3f (1.0f, 0.0f, 0.0f);
			glBegin(GL_LINE_LOOP);
			
	        glVertex3f (rProbe->x()-rProbe->hl(),rProbe->y()-rProbe->hh(), 0.0f);
	        glVertex3f (rProbe->x()-rProbe->hl(),rProbe->y()+rProbe->hh(), 0.0f);
	        glVertex3f (rProbe->x()+rProbe->hl(),rProbe->y()+rProbe->hh(), 0.0f);
	        glVertex3f (rProbe->x()+rProbe->hl(),rProbe->y()-rProbe->hh(), 0.0f);
	          
	        glEnd();
			
			
			
			if( biaxA->ngap1() != NULL && biaxA->ngap2() != NULL)
			{
				bodyi=biaxA->ngap1();
				L = bodyi->sizeVerlet()*.5 ;
				glColor4f (1.0f, 0.0f, 0.0f, 1.f);
				glBegin(GL_POLYGON);
				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+L*cos(t),
						bodyi->y()+L*sin(t), 0.0f);
				}
				glEnd();


				bodyi=biaxA->ngap2();
				L = bodyi->sizeVerlet()*.5 ;
				glColor4f (1.0f, 0.0f, 0.0f, 1.0f);
				glBegin(GL_POLYGON);
				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+L*cos(t),
						bodyi->y()+L*sin(t), 0.0f);
				}
				glEnd();
			}
			
				//interdof * interdofk ;
			
				/*for (unsigned int j=0 ; j<biaxA->clusterA()->linterdof().size() ; ++j)
				{
					interdofk = biaxA->clusterA()->catchInterdof(j);
					if(_disp_fresNormal)//
					{
						glLineWidth(1.0f);
					
							if (interdofk->rank() == 1)
							{
								glColor3f (1.0f, 0.0f, 0.0f);
								
								glBegin(GL_POLYGON);
								glVertex3f (interdofk->first()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,//la projection de fresNormal_ dans le repere globale 
									interdofk->first()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);//donne l'épaisseur du trait.
								glVertex3f (interdofk->second()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f); 
								glVertex3f (interdofk->second()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glVertex3f (interdofk->first()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->first()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glEnd();
							}
							else if (interdofk->rank() == 2)
							{
								glColor3f (0.0f, 1.0f, 0.0f);//

								glBegin(GL_POLYGON);
								glVertex3f (interdofk->first()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->first()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);
								glVertex3f (interdofk->second()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f); 
								glVertex3f (interdofk->second()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glVertex3f (interdofk->first()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->first()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glEnd();

							}
							else if (interdofk->rank() == 3)
							{
							    glColor3f (0.0f, 0.0f, 1.0f);
								
								glBegin(GL_POLYGON);
								glVertex3f (interdofk->first()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->first()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);
								glVertex3f (interdofk->second()->mcx()+interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()+interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f); 
								glVertex3f (interdofk->second()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->second()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glVertex3f (interdofk->first()->mcx()-interdofk->tx()*interdofk->fresNormal()*_scalf,
									interdofk->first()->mcy()-interdofk->ty()*interdofk->fresNormal()*_scalf, 0.0f);    
								glEnd();

							}
					}

			}//fin for*/
		}//sys_A	
	}//sys
	

}


void Display::drawSample( ) 
{
	glShadeModel(GL_FLAT);

	glLineWidth(1.0f);

	unsigned int N = sim()->spl()->lbody().size();
	body2d* bodyi;

	for (unsigned int i = 0 ; i < N ; ++i)
	{
	// idee1 : fonction surcharge -> display(dynamic_cast<disk*>(spl.body(i))) 

		bodyi=sim()->spl()->body(i);
		
		if( _fill_bodies )
		{
			if (typeid(*bodyi) == typeid(disk)) // ----------
			{
				float R=dynamic_cast<disk*>(bodyi)->R();

				if (_disp_shapes)
				{
					glColor3f (0.5f, 0.5f, 0.5f);
					
					//glColor4f (1.0f, 1.0f, 0.0f, 0.2f);
					glBegin(GL_POLYGON);

					for (float t=0.;t<_2pi;t+=_pi_8)
					{
						glVertex3f (bodyi->x()+R*cos(t),
							bodyi->y()+R*sin(t), 0.0f);
					}
					glEnd();
				}

				if(_disp_orient)
				{
					glColor3f (0.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f (bodyi->x()+R*cos(bodyi->rot()), 
						bodyi->y()+R*sin(bodyi->rot()), 0.0f);
					glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
					glEnd();
				}

				if(_disp_vel)
				{
					glColor3f (0.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f (bodyi->x()+bodyi->vx()*_scalv, 
						bodyi->y()+bodyi->vy()*_scalv, 0.0f);
					glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
					glEnd();
				}

			} //disk

			if (typeid(*bodyi) == typeid(polyg)) // ----------
			{
				polyg* Polyg = dynamic_cast<polyg*>(bodyi);
				float rot = bodyi->rot();
				float c   = cos(rot);
				float s   = sin(rot);

				if (_disp_shapes)
				{
					glColor3f (0.5f, 0.5f, 0.5f);
					glBegin(GL_POLYGON);

					for (unsigned int i=0 ; i<Polyg->Vertex().size() ; ++i)
					{
						float x = Polyg->Vertex(i).x();
						float y = Polyg->Vertex(i).y();

						glVertex3f (bodyi->x() + (x*c - y*s),
							bodyi->y() + (x*s + y*c), 0.0f);
					}
					glEnd();
				}
				
				if(_disp_orient)
				{
					glColor3f (0.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f (bodyi->x()+.7*bodyi->sizeVerlet()*cos(bodyi->rot()), 
						bodyi->y()+.7*bodyi->sizeVerlet()*sin(bodyi->rot()), 0.0f);
					glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
					glEnd();
				}
				
				if(_disp_vel)
				{
					glColor3f (0.0f, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex3f (bodyi->x()+bodyi->vx()*_scalv, 
						bodyi->y()+bodyi->vy()*_scalv, 0.0f);
					glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
					glEnd();
				}

			} // polyg

			if (typeid(*bodyi) == typeid(rline))
			{
				float L_2 = dynamic_cast<rline*>(bodyi)->L()*0.5;
				float R   = dynamic_cast<rline*>(bodyi)->R();
				float rot = bodyi->rot();
				float c   = cos(rot);
				float s   = sin(rot);

				if(_disp_shapes)
				{
					glLineWidth(1.0f);
					glColor3f (0.5f, 0.5f, 0.5f);

					glBegin(GL_POLYGON);
					for (float t=rot-_pi_2 ; t< rot+_pi_2 ; t += _pi_8)
						glVertex3f (bodyi->x()+L_2*c+R*cos(t),
						bodyi->y()+L_2*s+R*sin(t), 0.0f);

					for (float t=rot+_pi_2 ; t<rot+_3pi_2 ; t += _pi_8)
						glVertex3f (bodyi->x()-L_2*c+R*cos(t),
						bodyi->y()-L_2*s+R*sin(t), 0.0f);
					glEnd();
				}

			}  //rline
			
		}
		else//void bodies
		{
		if (typeid(*bodyi) == typeid(disk)) // ----------
		{
			float R=dynamic_cast<disk*>(bodyi)->R();

			if (_disp_shapes)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINE_LOOP);

				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+R*cos(t),
						bodyi->y()+R*sin(t), 0.0f);
				}
				glEnd();
			}

			if(_disp_orient)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINES);
				glVertex3f (bodyi->x()+R*cos(bodyi->rot()), 
					bodyi->y()+R*sin(bodyi->rot()), 0.0f);
				glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
				glEnd();
			}

			if(_disp_vel)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINES);
				glVertex3f (bodyi->x()+bodyi->vx()*_scalv, 
					bodyi->y()+bodyi->vy()*_scalv, 0.0f);
				glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
				glEnd();
			}

		} //disk

		if (typeid(*bodyi) == typeid(polyg)) // ----------
		{
			polyg* Polyg = dynamic_cast<polyg*>(bodyi);
			float rot = bodyi->rot();
			float c   = cos(rot);
			float s   = sin(rot);

			if (_disp_shapes)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINE_LOOP);

				for (unsigned int i=0 ; i<Polyg->Vertex().size() ; ++i)
				{
					float x = Polyg->Vertex(i).x();
					float y = Polyg->Vertex(i).y();

					glVertex3f (bodyi->x() + (x*c - y*s),
						bodyi->y() + (x*s + y*c), 0.0f);
				}
				glEnd();



		// >>> temporaire
		/*
				float R = dynamic_cast<polyg*>(bodyi)->Rout();
				glBegin(GL_LINE_LOOP);

				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+R*cos(t),
						bodyi->y()+R*sin(t), 0.0f);
				}
				glEnd();
		*/
		// <<<


			}
			
			if(_disp_orient)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINES);
				glVertex3f (bodyi->x()+.7*bodyi->sizeVerlet()*cos(bodyi->rot()), 
					bodyi->y()+.7*bodyi->sizeVerlet()*sin(bodyi->rot()), 0.0f);
				glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
				glEnd();
			}
			
			if(_disp_vel)
			{
				glColor3f (0.0f, 0.0f, 0.0f);
				glBegin(GL_LINES);
				glVertex3f (bodyi->x()+bodyi->vx()*_scalv, 
					bodyi->y()+bodyi->vy()*_scalv, 0.0f);
				glVertex3f (bodyi->x(), bodyi->y(), 0.0f);
				glEnd();
			}

		} // polyg

		if (typeid(*bodyi) == typeid(rline))
		{
			float L_2 = dynamic_cast<rline*>(bodyi)->L()*0.5;
			float R   = dynamic_cast<rline*>(bodyi)->R();
			float rot = bodyi->rot();
			float c   = cos(rot);
			float s   = sin(rot);

			if(_disp_shapes)
			{
				glLineWidth(1.0f);
				glColor3f (0.0f, 0.0f, 0.0f);

				glBegin(GL_LINE_LOOP);
				for (float t=rot-_pi_2 ; t< rot+_pi_2 ; t += _pi_8)
					glVertex3f (bodyi->x()+L_2*c+R*cos(t),
					bodyi->y()+L_2*s+R*sin(t), 0.0f);

				for (float t=rot+_pi_2 ; t<rot+_3pi_2 ; t += _pi_8)
					glVertex3f (bodyi->x()-L_2*c+R*cos(t),
					bodyi->y()-L_2*s+R*sin(t), 0.0f);
				glEnd();
			}

		}  //rline
	}
	} //for

			//Periodicity
	if( sim()->spl()->isMonoPeriodic())
	{
		unsigned int N=sim()->spl()->leftband().size();
		double P=sim()->spl()->boundWidth();
		double L;
		for( unsigned int i=0;i<N;++i)
		{
			bodyi=sim()->spl()->body( sim()->spl()->leftband(i));
			glColor3f (0.0f, 1.0f, 0.0f);
			glBegin(GL_LINE_LOOP);
			L = (bodyi->xmax() - bodyi->x()) ;
			for (float t=0.;t<_2pi;t+=_pi_8)
			{
				glVertex3f (bodyi->x()+P+L*cos(t),
					bodyi->y()+L*sin(t), 0.0f);
			}
			glEnd();
		}

			//Grains to translate
		N=sim()->spl()->lbody().size();
		for( unsigned int i=0;i<N;++i)
		{
			bodyi=sim()->spl()->body( i);
			if( bodyi->x() < sim()->spl()->leftBoundary() ||  bodyi->x() > sim()->spl()->rightBoundary() )
			{
				glColor4f (0.0f, 1.0f, 0.0f , 0.5f);
				L = (bodyi->xmax() - bodyi->x()) *0.8;
				glBegin(GL_POLYGON);
				for (float t=0.;t<_2pi;t+=_pi_8)
				{
					glVertex3f (bodyi->x()+L*cos(t),
						bodyi->y()+L*sin(t), 0.0f);
				}
				glEnd();
			}
		}
	}
	
	
	if (_disp_tracking_body && trackingBody_ != NULL)
	{
		glColor4f (1.0f, 0.5f, 0.0f , 0.9f);
		double L = trackingBody_->sizeVerlet() *0.8;
		glBegin(GL_POLYGON);
		for (float t=0.;t<_2pi;t+=_pi_8)
		{
			glVertex3f (trackingBody_->x()+L*cos(t),
				trackingBody_->y()+L*sin(t), 0.0f);
		}
		glEnd();
		
		
		glColor3f (0.0f, 0.0f, 1.f );
		L = trackingBody_->sizeVerlet() + sim()->nwk()->dverlet() ;
		glBegin(GL_LINE_LOOP);
		for (float t=0.;t<_2pi;t+=_pi_8)
		{
			glVertex3f (trackingBody_->x()+L*cos(t),
				trackingBody_->y()+L*sin(t), 0.0f);
		}
		glEnd();
		
		if (sim()->nwk()->useSuperList())
		{
			glColor3f (0.0f, 0.0f, 1.f );
			L = trackingBody_->sizeVerlet() + sim()->nwk()->dsuperList() ;
			glBegin(GL_LINE_LOOP);
			for (float t=0.;t<_2pi;t+=_pi_8)
			{
				glVertex3f (trackingBody_->x()+L*cos(t),
					trackingBody_->y()+L*sin(t), 0.0f);
			}
			glEnd();
		}
	}
	
	if (_disp_tracking_cluster && trackingCluster_ != NULL)
	{
		for( unsigned int i=0;i< trackingCluster_->lctrlBodies().size();++i)
		{
		glColor4f (1.0f, 1.0f, .0f , 0.4f);
		double L = trackingCluster_->ctrlBody(i)->sizeVerlet() *0.8;
		glBegin(GL_POLYGON);
		for (float t=0.;t<_2pi;t+=_pi_8)
		{
			glVertex3f (trackingCluster_->ctrlBody(i)->x()+L*cos(t),
				trackingCluster_->ctrlBody(i)->y()+L*sin(t), 0.0f);
		}
		glEnd();
		}
	}
}

void Display::drawNetwork( ) 
{
	glShadeModel(GL_FLAT);

	unsigned int N = sim()->nwk()->linter().size();
	inter2d* interk;

	for (unsigned int k = 0 ; k < N ; ++k)
	{
		interk=sim()->nwk()->inter(k);

		if (typeid(*interk) == typeid(dkdk))
		{
			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 0.0f, 1.0f);
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->second()->x(), 
					interk->second()->y(), 0.0f);
				glEnd();
			}


			if(_disp_fn)
			{
				glLineWidth(1.0f);
		// scalf = spl.Rmoy() / nwk.Maxfn()
				if(_disp_fn_pos && interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
							// mais le pb c'est que cette list est dans CDalgo au lieu de Network
							// il faudrait corriger ca...
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+interk->tx()*interk->fn()*_scalf,
						interk->second()->y()+interk->ty()*interk->fn()*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-interk->tx()*interk->fn()*_scalf,
						interk->second()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);    
					glEnd();
				}
				else if (_disp_fn_neg && interk->fn()<0.0)
				{
					//cout<<"force negative"<<interk->fn()<<endl;
					double fn = fabs(interk->fn() );
					glColor3f (0.0f, 0.0f, 1.0f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*fn*_scalf,
						interk->first()->y()+interk->ty()*fn*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+interk->tx()*fn*_scalf,
						interk->second()->y()+interk->ty()*fn*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-interk->tx()*fn*_scalf,
						interk->second()->y()-interk->ty()*fn*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*fn*_scalf,
						interk->first()->y()-interk->ty()*fn*_scalf, 0.0f);    
					glEnd();
				}
			}

			if(_disp_ft)
			{
				glLineWidth(1.0f);
		// scalf = spl.Rmoy() / nwk.Maxfn()
				if(interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
							// mais le pb c'est que cette list est dans CDalgo au lieu de Network
							// il faudrait corriger ca...
				{
					glColor3f (1.0f, 0.5f, 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+interk->tx()*interk->ft()*_scalf,
						interk->second()->y()+interk->ty()*interk->ft()*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-interk->tx()*interk->ft()*_scalf,
						interk->second()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);    
					glEnd();
				}
			}     
	/*
			if(_disp_ft)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (0.0f, 0.0f, 1.0f);

					float R = 0.5*( dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<disk*>(interk->second())->R());
					glBegin(GL_POLYGON);
					glVertex3f (interk->x()+interk->nx()*interk->ft()*_scalf+R*interk->tx(),
						interk->y()+interk->ny()*interk->ft()*_scalf+R*interk->ty(), 0.0f);
					glVertex3f (interk->x()-interk->nx()*interk->ft()*_scalf+R*interk->tx(),
						interk->y()-interk->ny()*interk->ft()*_scalf+R*interk->ty(), 0.0f); 
					glVertex3f (interk->x()-interk->nx()*interk->ft()*_scalf-R*interk->tx(),
						interk->y()-interk->ny()*interk->ft()*_scalf-R*interk->ty(), 0.0f);
					glVertex3f (interk->x()+interk->nx()*interk->ft()*_scalf-R*interk->tx(),
						interk->y()+interk->ny()*interk->ft()*_scalf-R*interk->ty(), 0.0f);
					glEnd();

				}
			}
	*/

			if(_disp_locFrame)
			{	if( interk->activate())
				{
				if( interk->rang()>0)
				{
					float R = 0.25*( dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<disk*>(interk->second())->R());
					glLineWidth(2.0f);
					glBegin(GL_LINES);
					glColor3f (1.0f, 0.0f, 0.0f);
					glVertex3f(interk->x(), interk->y(), 0.0f);
					glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
					glColor3f (0.0f, 0.0f, 1.0f);
					glVertex3f(interk->x(), interk->y(), 0.0f);
					glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
					glEnd();
				}
				else
				{
					glLineWidth(1.0f);
					glBegin(GL_LINES);
					glColor3f (0.0f, 0.0f, 1.0f);
					glVertex3f (interk->first()->x(), 
						interk->first()->y(), 0.0f);
					glVertex3f (interk->second()->x(), 
						interk->second()->y(), 0.0f);
					glEnd();
				}
				}
			}

		} // dkdk
		
	


		if (typeid(*interk) == typeid(dkdkP))
		{
			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 1.0f, 0.0f);
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->second()->x()+ sim()->spl()->boundWidth(), 
					interk->second()->y(), 0.0f);
				glEnd();
			}


			if(_disp_fn)
			{
				glLineWidth(1.0f);
		// scalf = spl.Rmoy() / nwk.Maxfn()
				if(_disp_fn_pos && interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
							// mais le pb c'est que cette list est dans CDalgo au lieu de Network
							// il faudrait corriger ca...
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()+interk->tx()*interk->fn()*_scalf,
						interk->second()->y()+interk->ty()*interk->fn()*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()-interk->tx()*interk->fn()*_scalf,
						interk->second()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);   
					glEnd();
				}
				else if (_disp_fn_neg && interk->fn()<0.0)
				{
					//cout<<"force negative"<<interk->fn()<<endl;
					double fn = fabs(interk->fn() );
					glColor3f (0.0f, 0.0f, 1.0f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*fn*_scalf,
						interk->first()->y()+interk->ty()*fn*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()+interk->tx()*fn*_scalf,
						interk->second()->y()+interk->ty()*fn*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()-interk->tx()*fn*_scalf,
						interk->second()->y()-interk->ty()*fn*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*fn*_scalf,
						interk->first()->y()-interk->ty()*fn*_scalf, 0.0f);    
					glEnd();
				}
			}

			if(_disp_ft)
			{
				glLineWidth(1.0f);
		// scalf = spl.Rmoy() / nwk.Maxfn()
				if(interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
							// mais le pb c'est que cette list est dans CDalgo au lieu de Network
							// il faudrait corriger ca...
				{
					glColor3f (1.0f, 0.5f, 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()+interk->tx()*interk->ft()*_scalf,
						interk->second()->y()+interk->ty()*interk->ft()*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()+sim()->spl()->boundWidth()-interk->tx()*interk->ft()*_scalf,
						interk->second()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);   
					glEnd();
				}
			}     
	/*
			if(_disp_ft)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (0.0f, 0.0f, 1.0f);

					float R = 0.5*( dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<disk*>(interk->second())->R());
					glBegin(GL_POLYGON);
					glVertex3f (interk->x()+interk->nx()*interk->ft()*_scalf+R*interk->tx(),
						interk->y()+interk->ny()*interk->ft()*_scalf+R*interk->ty(), 0.0f);
					glVertex3f (interk->x()-interk->nx()*interk->ft()*_scalf+R*interk->tx(),
						interk->y()-interk->ny()*interk->ft()*_scalf+R*interk->ty(), 0.0f); 
					glVertex3f (interk->x()-interk->nx()*interk->ft()*_scalf-R*interk->tx(),
						interk->y()-interk->ny()*interk->ft()*_scalf-R*interk->ty(), 0.0f);
					glVertex3f (interk->x()+interk->nx()*interk->ft()*_scalf-R*interk->tx(),
						interk->y()+interk->ny()*interk->ft()*_scalf-R*interk->ty(), 0.0f);
					glEnd();

				}
			}
	*/

			if(_disp_locFrame)
			{
				if( interk->activate())
					{
					if( interk->rang()>0)
					{
						float R = 0.25*( dynamic_cast<disk*>(interk->first())->R()
							+dynamic_cast<disk*>(interk->second())->R());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
					}
					else
					{
						glLineWidth(1.0f);
						glBegin(GL_LINES);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f (interk->first()->x(), 
							interk->first()->y(), 0.0f);
						glVertex3f (interk->second()->x()+ sim()->spl()->boundWidth(), 
							interk->second()->y(), 0.0f);
						glEnd();
					}
					}
				
			}

		} // dkdkP

		if (typeid(*interk) == typeid(pgpg))
		{
			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 0.0f, 1.0f);
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->second()->x(), 
					interk->second()->y(), 0.0f);
				glEnd();        
			}
			
			if(_disp_fn)
			{
				glLineWidth(1.0f);
				double norm,tx,ty;
		// scalf = spl.Rmoy() / nwk.Maxfn()
				if(interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
							// mais le pb c'est que cette list est dans CDalgo au lieu de Network
							// il faudrait corriger ca...
				{
					double fn = interk->fn();
					if ( interk->rang()==2)
					{
						interk->current()=1;
						fn+= interk->fn();
						interk->current()=0;
						
					}
				 	norm = sqrt( interk->Vbranchx()*interk->Vbranchx() + interk->Vbranchy()*interk->Vbranchy());
					tx = -(interk->Vbranchy())/norm;
					ty = interk->Vbranchx()/norm;
					glColor3f (1.0f, 0.0f, 0.0f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+tx*fn*_scalf,
						interk->first()->y()+ty*fn*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+tx*fn*_scalf,
						interk->second()->y()+ty*fn*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-tx*fn*_scalf,
						interk->second()->y()-ty*fn*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-tx*fn*_scalf,
						interk->first()->y()-ty*fn*_scalf, 0.0f); 
					/*glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+interk->tx()*interk->fn()*_scalf,
						interk->second()->y()+interk->ty()*interk->fn()*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-interk->tx()*interk->fn()*_scalf,
						interk->second()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);  */
					glEnd();
				}
			}
			if(_disp_ft)
			{
				glLineWidth(1.0f);
				double norm,tx,ty;
				
			// scalf = spl.Rmoy() / nwk.Maxfn()
				if(interk->fn()>0.0) // ne devrait pas etre la !!!! IL FAUDRAIT PARCOURIR LA _clist
								// mais le pb c'est que cette list est dans CDalgo au lieu de Network
								// il faudrait corriger ca...
				{
						double ft = interk->ft();
						if ( interk->rang()==2)
						{
							interk->current()=1;
							ft+= interk->ft();
							interk->current()=0;

						}
					 norm = sqrt( interk->Vbranchx()*interk->Vbranchx() + interk->Vbranchy()*interk->Vbranchy());
					 tx = -(interk->Vbranchy())/norm;
					 ty = interk->Vbranchx()/norm;
					
					glColor3f (1.0f, 0.5f, 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+tx*ft*_scalf,
						interk->first()->y()+ty*ft*_scalf, 0.0f);
					glVertex3f (interk->second()->x()+tx*ft*_scalf,
						interk->second()->y()+ty*ft*_scalf, 0.0f); 
					glVertex3f (interk->second()->x()-tx*ft*_scalf,
						interk->second()->y()-ty*ft*_scalf, 0.0f);    
					glVertex3f (interk->first()->x()-tx*ft*_scalf,
						interk->first()->y()-ty*ft*_scalf, 0.0f);   
					glEnd();
				}
			}
			if(_disp_locFrame)
			{
				if(interk->fn() > 0.)
				{
					if( interk->rang()==1)

					{
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
					}
					else if ( interk->rang()==2)
					{
						interk->current()=0;
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();

						interk->current()=1;
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (.0f, 1.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
						interk->current()=0;
					}
				}
			}
	/*
			if(_disp_locFrame)
			{
				for(unsigned int r=0;r<interk->rang();++r)
				{
					interk->current() = r;
		//if(interk->fn() > 0.)
					{
			//float Rout = 2.*( dynamic_cast<polyg*>(interk->first())->Rout()
			//                 +dynamic_cast<polyg*>(interk->second())->Rout());
						float Rout = 0.001;
			// Rin
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x()-interk->nx()*Rout, interk->y()-interk->ny()*Rout, 0.0f);
						glVertex3f(interk->x()+interk->nx()*Rout, interk->y()+interk->ny()*Rout, 0.0f);
						glColor3f (0.0f, 0.5f, 1.0f);
						glVertex3f(interk->x()-interk->tx()*Rout, interk->y()-interk->ty()*Rout, 0.0f);
						glVertex3f(interk->x()+interk->tx()*Rout, interk->y()+interk->ty()*Rout, 0.0f);
						glEnd();
					}
				}
			}     
	*/
		} // pgpg

		if (typeid(*interk) == typeid(dkpg))
		{

			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 0.0f, 1.0f);
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->second()->x(), 
					interk->second()->y(), 0.0f);
				glEnd();        
			}

			if(_disp_locFrame && interk->rang() == 1)
			{

			//float Rout = 2.*( dynamic_cast<polyg*>(interk->first())->Rout()
			//                 +dynamic_cast<polyg*>(interk->second())->Rout());
				float Rout = 0.001;
			// Rin

				glLineWidth(2.0f);
				glBegin(GL_LINES);
				glColor3f (1.0f, 0.0f, 0.0f);
				glVertex3f(interk->x()-interk->nx()*Rout, interk->y()-interk->ny()*Rout, 0.0f);
				glVertex3f(interk->x()+interk->nx()*Rout, interk->y()+interk->ny()*Rout, 0.0f);
				glColor3f (0.0f, 0.5f, 1.0f);
				glVertex3f(interk->x()-interk->tx()*Rout, interk->y()-interk->ty()*Rout, 0.0f);
				glVertex3f(interk->x()+interk->tx()*Rout, interk->y()+interk->ty()*Rout, 0.0f);
				glEnd();


			}      
		} // dkpg

		if (typeid(*interk) == typeid(dkrl))
		{
			if(_disp_fn)
			{
				glLineWidth(1.0f);
				if(_disp_fn_pos && interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					float R = dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);   
					glEnd();
				}
				else if (_disp_fn_neg && interk->fn()<0.0)
				{
					double fn = fabs( interk->fn());
					glColor3f (0.0f, 0.0f, 1.0f);
					float R = dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*fn*_scalf,
						interk->first()->y()+interk->ty()*fn*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*fn*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*fn*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*fn*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*fn*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*fn*_scalf,
						interk->first()->y()-interk->ty()*fn*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_ft)
			{
				glLineWidth(1.0f);
		//if(interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					float R = dynamic_cast<disk*>(interk->first())->R()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_locFrame)
			{
				if(interk->fn() > 0.)
				{
					float R = 0.5*( dynamic_cast<disk*>(interk->first())->R());
					glLineWidth(3.0f);
					glBegin(GL_LINES);
					glColor3f (0.0f, 1.0f, 0.0f);
					glVertex3f(interk->x()-interk->nx()*R, interk->y()-interk->ny()*R, 0.0f);
					glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
					glColor3f (0.0f, 0.0f, 1.0f);
					glVertex3f(interk->x()-interk->tx()*R, interk->y()-interk->ty()*R, 0.0f);
					glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
					glEnd();
				}
			}

			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 1.0f, 0.0f);
				float R = dynamic_cast<disk*>(interk->first())->R()
					+ dynamic_cast<rline*>(interk->second())->R();
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->first()->x()-R*interk->nx(), 
					interk->first()->y()-R*interk->ny(), 0.0f);
				glEnd();        
			}

		}

		if (typeid(*interk) == typeid(pgrl))
		{
			//pgrl * inter= dynamic_cast<pgrl*>(interk);
			if(_disp_fn)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					float R = dynamic_cast<polyg*>(interk->first())->Rout()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_ft)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.5f, 0.5f);							
					float R = dynamic_cast<polyg*>(interk->first())->Rout()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_locFrame)
			{
				if(interk->fn() > 0.)
				{
					if( interk->rang()==1)

					{
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
					}
					else if ( interk->rang()==2)
					{
						interk->current()=0;
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();

						interk->current()=1;
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 1.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
						interk->current()=0;
					}
				}
			}

			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 1.0f, 0.0f);
				float R = dynamic_cast<polyg*>(interk->first())->Rout()
					+ dynamic_cast<rline*>(interk->second())->R();
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->first()->x()-R*interk->nx(), 
					interk->first()->y()-R*interk->ny(), 0.0f);
				glEnd();        
			}


		}
		if (typeid(*interk) == typeid(pgrl))
		{
			//pgrl * inter= dynamic_cast<pgrl*>(interk);
			if(_disp_fn)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.0f, 0.0f);
					float R = dynamic_cast<polyg*>(interk->first())->Rout()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->fn()*_scalf,
						interk->first()->y()-interk->ty()*interk->fn()*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_ft)
			{
				glLineWidth(1.0f);
				if(interk->fn()>0.0)
				{
					glColor3f (1.0f, 0.5f, 0.5f);							
					float R = dynamic_cast<polyg*>(interk->first())->Rout()
						+dynamic_cast<rline*>(interk->second())->R();
					glBegin(GL_POLYGON);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf, 0.0f);
					glVertex3f (interk->first()->x()+interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()+interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f); 
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf-interk->nx()*R,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf-interk->ny()*R, 0.0f);
					glVertex3f (interk->first()->x()-interk->tx()*interk->ft()*_scalf,
						interk->first()->y()-interk->ty()*interk->ft()*_scalf, 0.0f);   
					glEnd();
				}
			}

			if(_disp_locFrame)
			{
				if(interk->fn() > 0.)
				{
					if( interk->rang()==1)
					
					{
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
					glLineWidth(2.0f);
					glBegin(GL_LINES);
					glColor3f (1.0f, 0.0f, 0.0f);
					glVertex3f(interk->x(), interk->y(), 0.0f);
					glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
					glColor3f (0.0f, 0.0f, 1.0f);
					glVertex3f(interk->x(), interk->y(), 0.0f);
					glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
					glEnd();
					}
					else if ( interk->rang()==2)
					{
						interk->current()=0;
						float R = 0.5*( dynamic_cast<polyg*>(interk->first())->Rout());
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 0.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
						
						interk->current()=1;
						glLineWidth(2.0f);
						glBegin(GL_LINES);
						glColor3f (1.0f, 0.0f, 0.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->nx()*R, interk->y()+interk->ny()*R, 0.0f);
						glColor3f (0.0f, 1.0f, 1.0f);
						glVertex3f(interk->x(), interk->y(), 0.0f);
						glVertex3f(interk->x()+interk->tx()*R, interk->y()+interk->ty()*R, 0.0f);
						glEnd();
						interk->current()=0;
					}
				}
			}

			if(_disp_verlet)
			{
				glLineWidth(1.0f);
				glBegin(GL_LINES);
				glColor3f (0.0f, 1.0f, 0.0f);
				float R = dynamic_cast<polyg*>(interk->first())->Rout()
					+ dynamic_cast<rline*>(interk->second())->R();
				glVertex3f (interk->first()->x(), 
					interk->first()->y(), 0.0f);
				glVertex3f (interk->first()->x()-R*interk->nx(), 
					interk->first()->y()-R*interk->ny(), 0.0f);
				glEnd();        
			}
			

		}

	} //for
}

void Display::draw() 
{
	if (!valid()) 
	{
		glLoadIdentity();
		glViewport(0,0,w(),h());
		glOrtho(-10,10,-10,10,-20020,10000);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	glClearColor(1.0f,1.0f,1.0f,0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	glTranslatef(_xshift, _yshift, 0);
	glRotatef(_hAng,0,1,0); glRotatef(_vAng,1,0,0);
	glScalef(float(_size),float(_size),float(_size));

//glRotatef(_hAng,0,1,0); glRotatef(_vAng,1,0,0);
//glScalef(float(_size),float(_size),float(_size));
//glTranslatef(_xshift, _yshift, 0);

//	glutKeyboardFunc( display_keyboard_wrapper(this) );	
	drawSample();
	drawNetwork();
	drawSystem();
	glPopMatrix();
}

