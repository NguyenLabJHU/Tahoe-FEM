/*
  File: VoidGrowthModelImp.h
*/

#ifndef _VOID_GROWTH_MODEL_IMP_H_
#define _VOID_GROWTH_MODEL_IMP_H_

#include "ios_fwd_decl.h"
#include "EVPFDBaseT.h"

namespace Tahoe {

class VoidGrowthModelImp
{
 public:
  // enumeration to select type of void growth model
  enum VoidGrowthModel { kCocks    = 1,   // Cocks' model (1989)
                         kDuvaCrow = 2,    // Duva & Crow's model (1992)
                         kSofronis = 3 }; // Sofronis & McMeeking model (1992)

  // constructor/virtual destructor
  VoidGrowthModelImp();
  virtual ~VoidGrowthModelImp();

  // coefficients A1 and A2 of pressure dependent viscoplastic potential
  virtual void ACoefficients(const double& vvf, double& A1, double& A2) const = 0;

  // derivatives of A1 and A2 
  virtual void ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const = 0;
  
  // set rate sensitivity exponent
  void SetRateSensitivity(const double m);
  
 protected:
  // strain rate sensitivity exponent
  double fm;
};

inline void VoidGrowthModelImp::SetRateSensitivity(const double m) { fm = m; }

/* DERIVED CLASSES */

class CocksVGModel: public VoidGrowthModelImp
{
 public:
  // constructor/virtual destructor
  CocksVGModel(EVPFDBaseT& macro);
  ~CocksVGModel();

  // coefficients A1 and A2 of pressure dependent viscoplastic potential
  virtual void ACoefficients(const double& vvf, double& A1, double& A2) const;

  // derivatives of A1 and A2 
  virtual void ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const;
  
 private:
};

class DuvaCrowVGModel: public VoidGrowthModelImp
{
 public:
  // constructor/virtual destructor
  DuvaCrowVGModel(EVPFDBaseT& macro);
  ~DuvaCrowVGModel();
  
  // coefficients A1 and A2 of pressure dependent viscoplastic potential
  virtual void ACoefficients(const double& vvf, double& A1, double& A2) const;

  // derivatives of A1 and A2 
  virtual void ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const;
  
 private:
};

class SofronisVGModel: public VoidGrowthModelImp
{
 public:
  // constructor/virtual destructor
  SofronisVGModel(EVPFDBaseT& macro);
  ~SofronisVGModel();
  
  // coefficients A1 and A2 of pressure dependent viscoplastic potential
  virtual void ACoefficients(const double& vvf, double& A1, double& A2) const;

  // derivatives of A1 and A2 
  virtual void ADerivCoefficients(const double& vvf, double& dA1, double& dA2) const;
  
 private:
};

} // namespace Tahoe 
#endif /*  _VOID_GROWTH_MODEL_IMP_H_ */
