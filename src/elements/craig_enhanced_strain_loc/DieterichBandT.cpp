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
fLastJumpIncrement(0.0),
fDeltaTheta(0.0)
{}


double DieterichBandT::Theta()
{
  return fTheta;
}

double DieterichBandT::DeltaTheta()
{
  return fDeltaTheta;
}

void DieterichBandT::StoreDeltaTheta(double deltaTheta)
{
  fDeltaTheta = deltaTheta;
}

double DieterichBandT::JumpIncrLast()
{
  return fLastJumpIncrement;
}

void DieterichBandT::CloseStep()
{
  /* inherited */
  BandT::CloseStep();

  /* update ISV theta */
  fTheta += fDeltaTheta;
  /* store jump increment for next time step */
  fLastJumpIncrement = JumpIncrement();
}
