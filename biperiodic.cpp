#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<vector>

using namespace std;

class Tensor2x2{
	private:
		double xx_;
		double xy_;
		double yx_;
		double yy_;
	public:
		Tensor2x2(): xx_(0.0), xy_(0.0),yx_(0.0), yy_(0.0){};
		Tensor2x2(double xx, double xy, double yx, double yy): xx_(xx), xy_(xy), yx_(yx), yy_(yy){};
		~Tensor2x2(){};
		
		void set(double xx, double xy, double yx, double yy);
		double getDet(){ return (xx_*yy_ - xy_*yx_) ;}

		//Accesseurs:
		double getxx() const { return xx_;} 
		double getxy() const { return xy_;} 
		double getyx() const { return yx_;} 
		double getyy() const { return yy_;} 
};

void Tensor2x2::set(double xx, double xy, double yx, double yy){
	xx_ = xx;
	xy_ = xy;
	yx_ = yx;
	yy_ = yy;
}

class Vecteur{
	private:
		double x_;
		double y_;
	public:
		Vecteur(): x_(0.), y_(0.){};
		Vecteur(double x, double y): x_(x), y_(y){};
		~Vecteur(){};

		Vecteur productTensor(Tensor2x2&);
		void print(){cout<<x_<<" "<<y_<<endl;}
		void update(double x,double y){ x_=x; y_=y;}
		void add(double dx, double dy) { x_ += dx ; y_ += dy;}

		//Accessors:
		double getNorme() const {return sqrt( x_*x_ + y_ * y_);}
		double getx() const {return x_;}
		double gety() const {return y_;}
		//Mutators:
		void setx(double x) {x_ = x;}
		void sety(double y) {y_ = y;}
};

Vecteur Vecteur::productTensor(Tensor2x2& T){
	Vecteur r;
		r.x_ = T.getxx() * this->x_ + T.getxy() * this->y_;
		r.y_ = T.getyx() * this->x_ + T.getyy() * this->y_;
	return r;
}

class Particle{
	private:
		//Position absolue
		Vecteur r_;
		//Reduced coordinate
		Vecteur s_;
		//Vitesse:
		Vecteur v_;
		//Acceleration:
		Vecteur a_;
	public:
		Particle(): r_(), v_() {}; 
		Particle(Vecteur r, double L, Vecteur v) {r_ = r; v_ = v; calculReduced(L);}
		~Particle(){};
		void updateSimple(double dt);
		Vecteur& getR() { return r_;}
		void write(ofstream&);
		void calculReduced(double L);
		void affiche();
};

//sx et sy doivent etre compris entre 0 et 1 (initialiement dans la cellule)
void Particle::calculReduced(double L){
	s_.update(r_.getx() / L,r_.gety() / L);
}

void Particle::write(ofstream& of){
	of<< r_.getx()<<" "<<r_.gety()<<endl;
}

void Particle::affiche(){
	cout<<"Coordonnees absolues: "<<r_.getx()<<" "<<r_.gety()<<endl;
	cout<<"Coordonnees reduites: "<<s_.getx()<<" "<<s_.gety()<<endl;
}	


//Euler methode: update acceleration, vitesse and position
//Kinematics: update juste position a partir de la vitesse
void Particle::updateSimple(double dt){
	double x = r_.getx() + v_.getx() * dt;
	double y = r_.gety() + v_.gety() * dt;
	r_.update(x,y);
}

//Cellule periodique
class Cell{

	private:
		//Dimension:
		double L_;
		//Coordonnees centre:
		double xc_;
		double yc_;
		//Metrics: collective degrees of freedom
		Tensor2x2 h;
		Tensor2x2 hdot;
	public:
		Cell();
		Cell(double L): L_(L),xc_(L/2.),yc_(L/2.){initAffine();}
		~Cell();

		//Met a jours la periodicite 
		void PeriodicBoundaries(Particle&);
		//Distorsion cellule sur un pas de temps(calcul champ homogene)
		void distorsion(double dt);
		void initAffine();
		//Renvoie le volume de la cellule
		double getVolume();
		double getL() const { return L_;}
		double getL2() const { return L_/2.;}
};

Cell::Cell(){
	L_ = 1. ;
	xc_ = 0.5 * L_ ;
	yc_ = 0.5 * L_ ;
}


//Init h et hdot
void Cell::initAffine(){
	//Cellule initiale:
	h.set(L_,0.,0.,L_);
	//Taux de deformation de la cellue:
	hdot.set(1.,0.,0.,1.);
}

//Le volume est donné par det(h) (deux vecteurs de base de la cellule)
double Cell::getVolume(){
	return h.getDet();
}
	

Cell::~Cell(){};

//Applique le deplacement affine
void Cell::distorsion(double dt){




}


//PeriodicBoundary Conditions
void Cell::PeriodicBoundaries(Particle& p){

	//Periodic position
	//if( p.getR().getx() > getL2() ) p.getR().add(-getL(),0.);
	//if( p.getR().getx() < -getL2() ) p.getR().add(getL(),0.);
	//if( p.getR().gety() < -getL2() ) p.getR().add(0.,getL());
	//if( p.getR().gety() > getL2() ) p.getR().add(0.,-getL());

	//Other method without IF
	p.getR().add( - getL() *  floor( p.getR().getx()/ getL()) ,0. );
	p.getR().add(0., - getL() * floor( p.getR().gety()/ getL()) );
}


void verletalgo(Cell& cell,std::vector<Particle> ps,double dt){


}

int main(){

	//Initialisation:
	double const L = 1.;
	double dt = .05 ;
	Cell cell(L);

	Vecteur r0(L/3,L/3);
	Vecteur v0(0.,0.);
	//Initialise coordonnees reduites, absolues, vitesse fluctuante
	Particle p(r0,L,v0);
	p.affiche();

	ofstream tmp("particle.txt");
	//Kinetics Time integration
	for(int t = 0; t < 10; t++){
		//Check for boundary:
		cell.PeriodicBoundaries(p);
		//Update position:
		//Calcul du champ homogene: ie def cellule
		cell.distorsion(dt);
		//Application champ aux particules:

		//Print:
		cout<<cell.getVolume()<<endl;
		p.write(tmp);
	}
	tmp.close();
}

