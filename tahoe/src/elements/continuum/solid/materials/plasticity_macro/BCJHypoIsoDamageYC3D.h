/* $Id: BCJHypoIsoDamageYC3D.h,v 1.4 2003-01-29 07:35:06 paklein Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_YC_3D_H_
#define _BCJ_HYPO_ISO_DAMAGE_YC_3D_H_

#include "BCJHypo3D.h"
#include "VoidGrowthModelImp.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "SpectralDecompT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;

class BCJHypoIsoDamageYC3D : public BCJHypo3D
{
 public:
  // constructor
  BCJHypoIsoDamageYC3D(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~BCJHypoIsoDamageYC3D();

  // initialize arrays and void growth model
  virtual void Initialize();

  // Cauchy stress
  virtual const dSymMatrixT& s_ij(); 
 
  // form residual
  virtual void FormRHS(const dArrayT& variab, dArrayT& rhs);

  // form Jacobian
  virtual void FormLHS(const dArrayT& variab, dMatrixT& lhs);

  // update/reset state
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // output related methods
  virtual int NumOutputVariables() const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:
  // initial value of variables
  virtual void InitializeVariables(ElementCardT& element);

  // initialize void growth model (implicit VGM)
  virtual void SetVoidGrowthModel();

  // backward integration of constitutive equations
  virtual void IntegrateConstitutiveEqns(bool& converged, int subIncr, int totSubIncrs);

  // trial stresses
  virtual void ElasticTrialStress();

  // Cauchy stress and backstress
  virtual void UpdateStresses();

  // tangent moduli
  virtual void TangentModuli();

  // forward gradient estimate for primary unknowns (DEQPe, DEQPh, ALPH, KAPP, DAMG)
  virtual void ForwardGradientEstimate();

  // cheack for negative values of solution variables
  virtual bool IsSolnVariableNegative();

 private:
  // indexes to access internal variable (scalars) array
  enum InternalVariables { kDEQPe = 0,   // deviatoric plastic strain increment
			   kDEQPh = 1,   // volumetric plastic strain increment
			   kALPH  = 2,   // norm of back stress tensor
			   kKAPP  = 3,   // isotropic hardening stress
                           kDAMG  = 4 }; // void volume fraction

  enum EQValues { kEQPe_n  = 0,      // deviatoric equivalent plastic strain
                  kEQPe    = 1,
		  kEQPh_n  = 2,      // volumetric equivalent plastic strain
		  kEQPh    = 3,
		  kEQP_n   = 4,      // total equivalent plastic strain
		  kEQP     = 5,
		  kEQXie_n = 6,      // deviatoric equivalent stress
		  kEQXie   = 7,
		  kEQXih_n = 8,      // hydrostatic stress
		  kEQXih   = 9 };

  // internal quantities
  void ComputeInternalQntsRHS(const dArrayT& array);
  void ComputeInternalQntsLHS(const dArrayT& array);

 protected:
  // code for void growth model
  int fVGMCode;

  // initial damage and strain rate sensitivity exponent (use in VGModel)
  double fDamg0;
  double fm;

  // equivalent stress & equivalent plastic strain of porous material
  double fEQXi;
  double fDEQP;

  // pointer to void growth model class
  VoidGrowthModelImp* fVoidGrowthModel;

  // trial equivalent deviatoric/hydrostatic stresses
  double fEQXieTr;       // deviatoric
  double fEQXihTr;       // hydrostatic

  // coefficients for isotropic void growth (damage) model
  double fA1Dmg;         // A1
  double fA2Dmg;         // A2
  double fA1iDmg;        // 1/A1
  double fA2iDmg;        // 1/A2
  double fdA1Dmg;        // dA1/dDmg
  double fdA2Dmg;        // dA2/dDmg

  // derivatives of equivalent stresses
  dArrayT fdEQXie;
  dArrayT fdEQXih;
};

} // namespace Tahoe 
#endif /* _BCJ_HYPO_ISO_DAMAGE_YC_3D_ */
