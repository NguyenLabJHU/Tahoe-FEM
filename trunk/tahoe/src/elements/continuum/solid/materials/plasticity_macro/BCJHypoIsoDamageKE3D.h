/* $Id: BCJHypoIsoDamageKE3D.h,v 1.5 2004-07-15 08:29:14 paklein Exp $ */
#ifndef _BCJ_HYPO_ISO_DAMAGE_KE_3D_H_
#define _BCJ_HYPO_ISO_DAMAGE_KE_3D_H_

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

class BCJHypoIsoDamageKE3D : public BCJHypo3D
{
 public:
  // constructor
  BCJHypoIsoDamageKE3D(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~BCJHypoIsoDamageKE3D();

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

  // forward gradient estimate for primary unknowns (EQXie, EQXih, ALPH, KAPP, DAMG)
  virtual void ForwardGradientEstimate();

  // cheack for negative values of solution variables
  virtual bool IsSolnVariableNegative();

 private:
  // indexes to access internal variable (scalars) array
  enum InternalVariables { kEQXie = 0,   // deviatoric effective stress
			   kEQXih = 1,   // hydrostatis stress
			   kALPH  = 2,   // norm of back stress tensor
			   kKAPP  = 3,   // isotropic hardening stress
                           kDAMG  = 4 }; // void volume fraction

  enum EQValues { kEQPe_n    = 0,      // deviatoric effective plastic strain
                  kEQPe      = 1,
		  kEQPh_n    = 2,      // volumetric plastic strain
		  kEQPh      = 3,
		  kEQPeDot_n = 4,      // effective deviatoric plastic strain rate
		  kEQPeDot   = 5 };

  // internal quantities
  void ComputeInternalQntsRHS(const dArrayT& array);
  void ComputeInternalQntsLHS(const dArrayT& array);

 protected:
  // code for void growth model
  int fVGMCode;

  // initial damage and strain rate sensitivity exponent (use in VGModel)
  double fDamg0;
  double fm;

  // equivalent stress & equivalent plastic strain rates of porous material
  double fEQXi;
  double fEQPDot;
  double fEQPeDot;
  double fEQPhDot;

  // pointer to void growth model class
  VoidGrowthModelImp* fVoidGrowthModel;

  // trial equivalent stresses
  double fEQXieTr;       // deviatoric
  double fEQXihTr;       // hydrostatic

  // coefficients for isotropic void growth (damage) model
  double fA1Dmg;         // A1
  double fA2Dmg;         // A2
  double fdA1Dmg;        // dA1/dDmg
  double fdA2Dmg;        // dA2/dDmg

  // derivatives of: (i) strain rates and (ii) factor eta
  dArrayT fdEQPeDot;
  dArrayT fdEQPhDot;
  dArrayT fdEta;
};

} // namespace Tahoe 
#endif /* _BCJ_HYPO_ISO_DAMAGE_KE_3D_ */
