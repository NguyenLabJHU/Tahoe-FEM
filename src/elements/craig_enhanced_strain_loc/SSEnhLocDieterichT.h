#ifndef _SS_ENH_LOC_DIETERICH_T_H
#define _SS_ENH_LOC_DIETERICH_T_H

/* base class */
#include "SSEnhLocCraigT.h"

#include "DieterichBandT.h"

namespace Tahoe {

  /* forward declarations */
  class SSMatSupportT;
  
  class SSEnhLocDieterichT: public SSEnhLocCraigT
    {
    public:

      /* constructor */
      SSEnhLocDieterichT(const ElementSupportT& support);

      /** describe the parameters needed by the interface */
      virtual void DefineParameters(ParameterListT& list) const;

      /** accept parameter list */
      virtual void TakeParameterList(const ParameterListT& list);

    protected:
      virtual void FormStiffness(double constK);
      virtual double CalculateJumpIncrement();
      virtual bool IsBandActive();

      virtual BandT* FormNewBand(dArrayT normal, dArrayT slipDir,
				 dArrayT perpSlipDir, dArrayT coords, double residCohesion, ArrayT<dSymMatrixT>
				 stressList);
      
      /* to update Theta */
      virtual void CloseStep(void);

    private:

      /*band parameters*/
      double fMu_star;
      double fTheta_star;
      double fV_star;
      double fA;
      double fB;
      double fD_c;
      double fTheta_0; //perhaps same as fTheta_star?
    }; //end class declaration

}//end namespace Tahoe

#endif // _SS_ENH_LOC_DIETERICH_T_H_


