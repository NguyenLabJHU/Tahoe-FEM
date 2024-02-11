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
				 dArrayT perpSlipDir, dArrayT coords, double area);
      
      virtual void CloseStep(void);


      virtual void LoadBand(int elementNumber);

      /* math functions for jump increment */
	  virtual double NewCohesion(double slipRate, double jumpIncrement, double thetaNew);
      virtual dSymMatrixT StressIncrOnBand(double jumpIncrement);
      virtual dSymMatrixT LastStressOnBand();
      virtual dSymMatrixT AvgStrainRelaxation(double jumpIncrement);

      virtual dSymMatrixT FormdGdSigma(int ndof, double
					     slipRate, double thetaNew);



      virtual double JumpIncrement(double slipRate);
      virtual double ThetaNew(double slipRate);
      virtual double Phi(double slipRate, double jumpIncrement, double thetaNew);
      virtual double DPhidSlipRate(double slipRate, double jumpIncr, double thetaNew);
	  virtual double DCohesiondSlipRate(double slipRate, double jumpIncr, double thetaNew);
	  
	  virtual double DSecondOrderSofteningDSlipRate(double slipRate, double jumpIncr, double thetaNew);
	  virtual double D2MuDSlipRate2(double slipRate, double thetaNew);
	  virtual double D2ThetaDSlipRate2(double slipRate);
	  
	  virtual double DNormalStressdSlipRate(double jumpIncr);
      virtual dSymMatrixT DSigmadSlipRate(double JumpIncrement);
      virtual double DjumpdSlipRate();
      virtual double DmudSlipRate(double slipRate, double thetaNew);
      virtual double ArcSinhArg(double slipRate, double theta);
      virtual double DthetadSlipRate(double slipRate);
      virtual double NormalStress(double jumpIncr);
	  virtual double ShearStress(double jumpIncr);
      virtual double FrictionCoeff(double slipRate, double theta);
      virtual double arcsinh(double arg);

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
      double fBeta_zeta;
      double fBeta_theta;
	  bool fSimpleSoftening;
	  bool fNoFrictionInTension;
	  double fInitialSlipRate;
    }; //end class declaration

}//end namespace Tahoe

#endif // _SS_ENH_LOC_DIETERICH_T_H_


