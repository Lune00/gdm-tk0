#ifndef _display_h
#define _display_h

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.H>
#include <FL/glut.H>
#include "simulation.hpp"

#include "shearP_CD.hpp"
#include "shearP_CD_A.hpp"

#include "biaxial.hpp"
#include "biaxial_A.hpp"

#include "dkdkP.hpp"
#include "dkdk.hpp" // a corriger
#include "pgpg.hpp"
#include "dkrl.hpp" // a corriger
#include "pgrl.hpp"
#include "dkpg.hpp"
#include "interdof.hpp"

class Display : public Fl_Gl_Window 
{

public:


    
    Display(int x,int y,int w,int h,const char *l=0);
    
    void read_parameters  (istream&);
    void write_parameters (ostream&);
	
	Simulation*  sim()  const { return sim_;}
	Simulation* &sim()        { return sim_;}

    void draw(); 
    
    void size(float val)         { _size = val;          }
    float size()                 { return _size;         }
    void vAng(float val)         { _vAng = val;          }
    float vAng()                 { return _vAng;         }
    void hAng(float val)         { _hAng = val;          }
    float hAng()                 { return _hAng;         }
    void xshift(float val)       { _xshift = val;        }
    float xshift()               { return _xshift;       }
    void yshift(float val)       { _yshift = val;        }
    float yshift()               { return _yshift;       }

    void addsize(float val)         { _size += val;          }
	void addxshift(float val)       { _xshift += val;        }
	void addyshift(float val)       { _yshift += val;        }
    
    void disp_shapes(bool opt)   { _disp_shapes = opt;   }
    bool disp_shapes()           { return _disp_shapes;  }
    void disp_orient(bool opt)   { _disp_orient = opt;   }
    bool disp_orient()           { return _disp_orient;  }
    
    void disp_vel(bool opt)      { _disp_vel = opt;      }
    bool disp_vel()              { return _disp_vel;     }
    
    void disp_locFrame(bool opt) { _disp_locFrame = opt; }
    bool disp_locFrame()         { return _disp_locFrame;}

    void disp_fn(bool opt)       { _disp_fn = opt;       }
    bool disp_fn()               { return _disp_fn;      }
	void disp_fn_pos(bool opt)   { _disp_fn_pos = opt;   }
    bool disp_fn_pos()           { return _disp_fn_pos;  }
	void disp_fn_neg(bool opt)   { _disp_fn_neg = opt;   }
    bool disp_fn_neg()           { return _disp_fn_neg;  }

    void disp_ft(bool opt)       { _disp_ft = opt;       }
    bool disp_ft()               { return _disp_ft;      }
    void disp_verlet(bool opt)   { _disp_verlet = opt;   }
    bool disp_verlet()           { return _disp_verlet;  }
    
    void fixScalf(bool opt)      { _fixScalf = opt;      }
    bool fixScalf()              { return _fixScalf;     }
    void scalf(float val)        { _scalf = val;         }
    float scalf()                { return _scalf;        }
    void fixScalv(bool opt)      { _fixScalv = opt;      }
    bool fixScalv()              { return _fixScalv;     }
    void scalv(float val)        { _scalv = val;         }
    float scalv()                { return _scalv;        }

	void fill_bodies(bool opt)   { _fill_bodies = opt;   }
    bool fill_bodies()           { return _fill_bodies;  }


	void disp_trackingBody(bool opt)   { _disp_tracking_body = opt;   }
    bool disp_trackingBody()           { return _disp_tracking_body;  }

	void disp_trackingCluster(bool opt)   { _disp_tracking_cluster = opt;   }
    bool disp_trackingCluster()           { return _disp_tracking_cluster;  }

	void disp_cluster(bool opt)   { _disp_cluster = opt;   }
    bool disp_cluster()           { return _disp_cluster;  }
	
	//le 10/03/09 modif
	void disp_fresNormal(bool opt)		{_disp_fresNormal = opt;}
	bool disp_fresNormal()				{return _disp_fresNormal;}
	
	//void disp_fresNormal_pos(bool opt)	{_disp_fresNormal_pos = opt;}
	//bool disp_fresNormal_pos()			{return _disp_fresNormal;}
//-----------------------------------------------------------------------------
    void trackingBody( unsigned int i ) 
		{
			if( i < sim()->spl()->lbody().size() )
			trackingBody_ = sim()->spl()->body(i);
			else 
				cout<<" trackingBody over ... max = "<<sim()->spl()->lbody().size()<<endl;
			//return tracking_;
		}
	void trackingCluster( unsigned int i ) 
			{
				if( i < sim()->sys()->ldof().size() )
				trackingCluster_ = sim()->sys()->ldof(i);
				else 
					cout<<" trackingCluster over ... max = "<<sim()->sys()->ldof().size()<<endl;
				//return tracking_;
			}
    
private:

    void drawSample();
    void drawNetwork();
	void drawSystem();
		
	void keyboard(unsigned char key,int x,int y) {cout<<"apuiage de "<<key<<endl;}

//	Sample * spl_;
//	Network * nwk_;
//	System * sys_;
	Simulation * sim_; 
 
    float _size;
    float _vAng,_hAng;
    float _xshift,_yshift;
    
    float _scalf; bool _fixScalf;
    float _scalv; bool _fixScalv;

	bool    _disp_tracking_body;
	body2d* trackingBody_;
	
	bool    _disp_tracking_cluster;
	dof* trackingCluster_;
	
	bool _disp_cluster;
	
	bool _fill_bodies;
    
    bool _disp_shapes;
    bool _disp_orient;
    
    bool _disp_vel;
    
    bool _disp_locFrame;
    bool _disp_fn;
	bool _disp_fn_pos;
	bool _disp_fn_neg;
	
    bool _disp_ft;    
    bool _disp_verlet;
	
//modif du 09/03/09

	bool _disp_fresNormal;
	//bool _disp_fresNormal_pos;
};



#endif // _display_h

