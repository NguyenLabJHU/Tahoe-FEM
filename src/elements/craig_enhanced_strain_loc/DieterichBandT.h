/* base class */
#include "BandT.h"

#include "SSEnhLocDieterichT.h"

#ifndef _DIETERICH_BAND_T_H_
#define _DIETERICH_BAND_T_H_

namespace Tahoe {

  /* forward delaration */
  class SSEnhLocDieterichT;

  class DieterichBandT: public BandT
    { 

    public:

      DieterichBandT(dArrayT normal, dArrayT slipDir, dArrayT
		     perpSlipDir, dArrayT coords,
				   double fH_delta_0, double
				   residCohesion, ArrayT<dSymMatrixT> stressList,
				   SSEnhLocDieterichT* element, double
	       theta_0);

      virtual double Theta();
      virtual double ThetaLast() {return fLastTheta;};
      virtual void StoreTheta(double theta);
      virtual double SlipRate() {return fSlipRate;};
      virtual double SlipRateLast();
      virtual void StoreSlipRate(double slipRate) {fSlipRate = slipRate;};
      virtual double ThetaRateLast() {return fLastThetaRate;};
      virtual void CloseStep();
      virtual void UpdateThetaRate(double d_c);

    private:

      double fLastTheta;
      double fTheta;
      double fSlipRate;
      double fLastSlipRate;
      double fLastThetaRate;

    };

} //end namespace Tahoe

#endif //_DIETERICH_BAND_T_H_
