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
      
      //virtual void CloseStep(void);
      virtual void LoadBand(int elementNumber);



      /* math functions for jump increment */
      virtual double DeltaTheta(double jumpIncrement, double deltaTheta);
      virtual double DeltaG(double jumpIncr, double deltaTheta);
      virtual double DdeltaGdJump(double jumpIncr, double deltaTheta);
      virtual double DdeltaGdJumpGlobal(double jumpIncr, double deltaTheta);
     virtual dSymMatrixT StressIncrOnBand(double jumpIncrement);
      virtual dSymMatrixT LastStressOnBand();
      virtual dSymMatrixT AvgStrainRelaxation(double jumpIncrement);
      virtual double BigConstant(double jumpIncrement, double deltaTheta);

      virtual double DdeltaGdJumpAtConstTheta(double jumpIncrement,
					      double deltaTheta); 
      virtual double DdeltaGdTheta(double jumpIncrement, double
				   deltaTheta);
      virtual double DdeltaGdThetaGlobal(double jumpIncrement, double
				   deltaTheta);
      virtual double DThetaDJump(double jumpIncrement, double deltaTheta);
      virtual dSymMatrixT FormdGdSigma(int ndof, double
					     jumpIncrement, double deltaTheta);

    private:

      //DieterichBandT* &fDieterichBand;
      DieterichBandT* fDieterichBand;

      /*band parameters*/
      double fMu_star;
      double fTheta_star;
      double fV_star;
      double fFrictionA;
      double fFrictionB;
      double fD_c;
      double fTheta_0; //perhaps same as fTheta_star?
    }; //end class declaration

}//end namespace Tahoe

#endif // _SS_ENH_LOC_DIETERICH_T_H_


