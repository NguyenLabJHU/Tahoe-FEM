/* $Id: LocalCrystalPlast_C.h,v 1.6 2003-01-29 07:35:05 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_C_H_
#define _LOCAL_CRYSTAL_PLAST_C_H_

#include "LocalCrystalPlast.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;

class LocalCrystalPlast_C : public LocalCrystalPlast
{
 public:
  // constructor
  LocalCrystalPlast_C(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~LocalCrystalPlast_C();

  // initialize arrays
  virtual void Initialize();

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // update/reset crystal state
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
 
  // deformation gradient
	const dMatrixT& DeformationGradient(const LocalArrayT& disp);

  // deformation gradient at center of element
  const dMatrixT& DefGradientAtCenter(const LocalArrayT& disp);

 private:

  // number of crystal variables to be stored
  virtual int NumVariablesPerElement();
 
   // initial value of crystal variables
  virtual void InitializeCrystalVariables(ElementCardT&);

 protected:
  // number of nodes/element: 4-node Quad (2D) and 8-node Hexa (3D)
  const int fNNodes;

  // references to initial coords
  const LocalArrayT& fLocInitX;

  // arrays for shape funtion derivatives at center
  dArray2DT fLNa;
  dArray2DT fLDNa;
  dArray2DT fGDNa;

  // displacement gradient
  dMatrixT fGradU;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_C_H_ */
