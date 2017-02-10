#include "interdof.hpp"

//le 15/01/09 B.SAINT-CYR
void interdof::Frame()
{
	//cout<<i_->mcx()<<" "<<i_->mcy()<<endl;
	vbranchx_=i_->mcx()-j_->mcx();
	vbranchy_=i_->mcy()-j_->mcy();
	vbranch_ = sqrt(vbranchx_*vbranchx_+vbranchy_*vbranchy_);
	double invVbranch = 1.0/vbranch_;
	
	nx_ = vbranchx_ * invVbranch; 
	ny_ = vbranchy_ * invVbranch;
	tx_ = -ny_;
	ty_ = nx_;
}

void interdof::calcFres()//projection de la résultante des forces normales sur les vecteurs intercentres
{
	for(unsigned int i=0;i<linter_.size();++i)
	{
		fresNormal_ += linter_[i]->fx() * nx_ + linter_[i]->fy() * ny_;
		fresOrthoNorm_ += linter_[i]->fx() * tx_ + linter_[i]->fy() * ty_;
		fresNormalx_ = fresNormal_ * nx_ + fresOrthoNorm_ * tx_;
		fresNormaly_ = fresNormal_ * ny_ + fresOrthoNorm_ * ty_; 
		//cout<<i<<" "<<linter_[i]->fx()<<"---------- "<<linter_[i]->fy()<<endl;

	}
}

