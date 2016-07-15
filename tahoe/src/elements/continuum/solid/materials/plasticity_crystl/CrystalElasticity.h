/* $Id: CrystalElasticity.h,v 1.6 2016-07-15 13:07:10 tdnguye Exp $ */
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
                        kCubCrysElast = 2,  // orthotropic (cubic)
                        kPECrysElast = 3 }; // Polyethylene crystal

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
  // general stiffness coefficients for isotropic and cubic crystal
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

class PECrystalElast: public CrystalElasticity
{
public:
    // constructor/virtual destructor
    PECrystalElast(PolyCrystalMatT& poly);
    ~PECrystalElast();
    
    // compute elastic material constants (e.g. Lame's constants)
    virtual void ElasticityProps(dArrayT& matprop) const;
    
    // compute elasticity matrix
    virtual void ComputeModuli(dMatrixT& moduli) const;
    
    // query for isotropic/anisotropic elasticity (default: false)
    virtual bool IsIsotropic() const;
protected:
    //additional constants for orthotropic crystal
    double fC22, fC33, fC55,fC66, fC13, fC23;

private:
};

} // namespace Tahoe 
#endif /* _CRYSTAL_ELASTICTIY_H_ */
