/* $Id: CrystalElasticity.h,v 1.5 2004-07-15 08:29:06 paklein Exp $ */
#ifndef _CRYSTAL_ELASTICTIY_H_
#define _CRYSTAL_ELASTICTIY_H_

#include "ios_fwd_decl.h"

namespace Tahoe {

class PolyCrystalMatT;
class dArrayT;
class dMatrixT;

class CrystalElasticity
{
 public:
  // enumeration to select type of crystal elasticity
  enum CrysElastType { kIsoCrysElast = 1,   // isotropic
                       kCubCrysElast = 2 }; // orthotropic (cubic)

  // constructor/virtual destructor
  CrystalElasticity();
  virtual ~CrystalElasticity();

  // compute elastic material constants (e.g. Lame's constants)
  virtual void ElasticityProps(dArrayT& matprop) const;

  // compute elasticity matrix
  virtual void ComputeModuli(dMatrixT& moduli) const;
  
  // query for isotropic/anisotropic elasticity (default: false)
  virtual bool IsIsotropic() const;
  
 protected:
  // general stiffness coefficients
  double fC11;
  double fC12;
  double fC44;    // fC44=0.5*(fC11-fC12) for isotropic elasticity
};

class IsotropicCrystalElast: public CrystalElasticity
{
 public:
  // constructor/virtual destructor
  IsotropicCrystalElast(PolyCrystalMatT& poly);
  ~IsotropicCrystalElast();

  // query for isotropic elasticity
  virtual bool IsIsotropic() const;

 private:
  double fYoung;
  double fPoisson;
};

class CubicCrystalElast: public CrystalElasticity
{
 public:
  // constructor/virtual destructor
  CubicCrystalElast(PolyCrystalMatT& poly);
  ~CubicCrystalElast();
  
 private:
};

} // namespace Tahoe 
#endif /* _CRYSTAL_ELASTICTIY_H_ */
