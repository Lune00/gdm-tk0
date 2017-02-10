#ifndef _anisotropy_hpp
#define _anisotropy_hpp

#include "sample.hpp"
#include "network.hpp"
#include "tensor.hpp"
#include "probe.hpp"
#include "heightProbe.hpp"
#include "dkdkP.hpp"
#include "pointSet.hpp"
#include "dataSet.hpp"

#include <iomanip>

using namespace std;

gdm::Tensor2x2 * Fabric(Sample& spl, Network& nwk);

gdm::Tensor2x2 * FabricInProbe(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * FabricInProbe_ss(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * FabricInProbe_sv(Probe& prb, Sample& spl, Network& nwk);

gdm::Tensor2x2 * Fabric2InProbe(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * Fabric2InProbe_ss(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * Fabric2InProbe_sv(Probe& prb, Sample& spl, Network& nwk);

gdm::Tensor2x2 * SimpleContactFabricInProbe(Probe& prb, Sample& spl, Network& nwk);//need to be normalised

gdm::Tensor2x2 * DoubleContactFabricInProbe(Probe& prb, Sample& spl, Network& nwk);//need to be normalised

gdm::Tensor2x2 * fnAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * fnAnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * fnAnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * fnAnisoInProbe2(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * fn2AnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * fn2AnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * fn2AnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * ftAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ftAnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ftAnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * ftAnisoInProbe2(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * ft2AnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ft2AnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ft2AnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * lengthAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * length2AnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * length2AnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * length2AnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * lnAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * lnAnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * lnAnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * ltAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ltAnisoInProbe_ss(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ltAnisoInProbe_sv(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * fnlAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);
gdm::Tensor2x2 * ftlAnisoInProbe(Probe & prb, Sample& spl, Network & nwk, unsigned int Nbin);

gdm::Tensor2x2 * fnlcAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);
gdm::Tensor2x2 * ftlcAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);



pointSet  PthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);

pointSet  CorrelationthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);
gdm::Tensor2x2 * CorrelationAnisoInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);

/*     
pointSet  FNthetaInProbe( Probe & prb, Sample &spl, Network& nwk);
         
pointSet  FThetaInProbe( Probe & prb, Sample &spl, Network& nwk);
 */
pointSet  SimpleContactPthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);
                                                                           
pointSet  DoubleContactPthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nc);

gdm::Tensor2x2 * fnlthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nt,unsigned int nl);

gdm::Tensor2x2 * ftlthetaInProbe( Probe & prb, Sample &spl, Network& nwk,unsigned int nt,unsigned int nl);

#endif // _anisotropy_hpp


