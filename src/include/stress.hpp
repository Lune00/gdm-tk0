#ifndef _stress_hpp
#define _stress_hpp

#include "sample.hpp"
#include "network.hpp"
#include "tensor.hpp"
#include "probe.hpp"
#include "heightProbe.hpp"
#include "dkdkP.hpp"
#include "dkdk.hpp"
#include "disk.hpp"
#include "interdof.hpp"

using namespace std;

//! \brief 
//! \author C. Voivret

gdm::Tensor2x2 * Stress(Sample& spl, Network& nwk);
gdm::Tensor2x2 * StressInProbe(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * StressInProbe_sv(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * StressInProbe_ss(Probe& prb, Sample& spl, Network& nwk);
gdm::Tensor2x2 * PartialLengthStressInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Npoint,double lmoy,double rmax,double pressure);
gdm::Tensor2x2 * PartialNormalForceStressInProbe(Probe& prb, Sample& spl, Network& nwk,unsigned int Npoint,double fnmoy,double pressure);

void IMT_Body(Sample &spl, Network &nwk, vector <gdm::Tensor2x2 > &bodyStress,double rho=1.);

#endif // _stress_hpp
