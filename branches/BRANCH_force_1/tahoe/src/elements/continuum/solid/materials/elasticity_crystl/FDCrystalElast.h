/* $Id: FDCrystalElast.h,v 1.6 2004-07-15 08:28:16 paklein Exp $ */
#ifndef _FD_CRYSTAL_ELAST_H_
#define _FD_CRYSTAL_ELAST_H_

#include "CrystalElast.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "LAdMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;
class UpLagr_ExternalFieldT;

class FDCrystalElast : public CrystalElast
{
 public:
  // constructor
  FDCrystalElast(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~FDCrystalElast();

  // number of crystal variables to be stored
  virtual int NumVariablesPerElement();

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // output related methods
  virtual int NumOutputVariables() const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

 protected:
  // initial value of crystal variables
  virtual void InitializeCrystalVariables();

  // recover crystal variables
  virtual void LoadCrystalData(ElementCardT& element, int intpt);

  // Cauchy stress
  void CrystalS_ij();

  // stress due to thermal strains 
  void InverseFthermal();

  // elastic moduli
  void CrystalC_ijkl();

  // 4th order tensor: c_ijkl=0.5*(b_ik b_jl + b_il b_jk)
  void Set_I_b_Tensor(const dSymMatrixT& b, dMatrixT& c);

  // 4th order tensor transformation: Co_ijkl = F_iI F_jJ F_kK f_lL Ci_IJKL
  void FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F);

  // vector multiplications
  void a_i_b_i(double& C, dArrayT a, dArrayT b);

 protected:

  //element type with external field
  const UpLagr_ExternalFieldT* fExFieldElement;

  // temperature
  double fTemp_DegC;
  dArrayT array1;

  // normals to slip plane
  dArrayT fVecNorm;           //specimen orientation
  dArrayT fVecNormC;          //crystal orientation

  // stress normal to slip plane
  double fsnorm;

  // deformation gradients
  dMatrixT fF;   
  dMatrixT fFe;   
  dMatrixT fFthi;   

  // tensors in material configuration
  dSymMatrixT fCe;
  dSymMatrixT fBe;
  dSymMatrixT fEe;
  dSymMatrixT fS;

  // tensors in polar decomposition of Fe
  dSymMatrixT fUe;
  dMatrixT fRe;

  // rotation matrix from Euler angles
  dMatrixT fRotMat;

  // crystal Cauchy stress
  dSymMatrixT fs_ij;
  
  // anisotropic crystal elasticity matrix
  dMatrixT fc_ijkl;

  // anisotropic contribution to crystal elasticity matrix
  dMatrixT fCanisoLat;    // lattice frame
  dMatrixT fCanisoSpl;    // sample configuration

  // second order identity tensor
  dSymMatrixT fISym;

  // principal values of Cauchy stress
  dArrayT fsEigs;

  // tensor of thermal expansion
  dMatrixT falpha; 

  // general workspaces
  dMatrixT fmatx1;
  dMatrixT fmatx2;
  dMatrixT fmatx3;
  dMatrixT fRank4;
  dSymMatrixT fsymmatx1;
  dArrayT fvector1;
};

} // namespace Tahoe 
#endif /* _FD_CRYSTAL_ELAST_H_ */

