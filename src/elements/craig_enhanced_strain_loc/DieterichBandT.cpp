#include "DieterichBandT.h"

using namespace Tahoe;

/* constructor */
DieterichBandT::DieterichBandT(dArrayT normal, dArrayT slipDir, dArrayT
perpSlipDir, dArrayT coords,
				   double fH_delta_0, double
				   residCohesion, ArrayT<dSymMatrixT> stressList,
				   SSEnhLocDieterichT* element, double
			       theta_0):
BandT(normal, slipDir, perpSlipDir, coords,
				   fH_delta_0, residCohesion, stressList,
      element),
fTheta(theta_0),
fLastTheta(theta_0),
fSlipRate(0.0),
fLastSlipRate(0.0),
fLastThetaRate(0.0) //what should this be?
//fLastJumpIncrement(0.0),
//fDeltaTheta(0.0)
{}


double DieterichBandT::Theta()
{
  return fTheta;
}

/*
double DieterichBandT::DeltaTheta()
{
  return fDeltaTheta;
}
*/

void DieterichBandT::StoreTheta(double theta)
{
  fTheta = theta;
}

double DieterichBandT::SlipRateLast()
{
  return fLastSlipRate;
}

void DieterichBandT::CloseStep()
{
  /* inherited */
  BandT::CloseStep();

  /* update ISV theta */
  fLastTheta = fTheta;
  /* store jump increment for next time step */
  fLastSlipRate = fSlipRate;
  //fLastThetaRate = 1.0 -  //get d_c
}
