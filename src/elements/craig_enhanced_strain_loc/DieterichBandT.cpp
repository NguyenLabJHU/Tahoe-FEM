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
fTheta(theta_0)
{}
