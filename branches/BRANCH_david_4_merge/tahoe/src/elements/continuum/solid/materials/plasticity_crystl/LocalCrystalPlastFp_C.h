/* $Id: LocalCrystalPlastFp_C.h,v 1.6 2005-01-21 16:51:22 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_FP_C_H_
#define _LOCAL_CRYSTAL_PLAST_FP_C_H_

#include "LocalCrystalPlastFp.h"

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

class LocalCrystalPlastFp_C : public LocalCrystalPlastFp
{
 public:
  // constructor
  LocalCrystalPlastFp_C(void);

  // destructor
  ~LocalCrystalPlastFp_C();

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

 protected:
 
  // deformation gradient at center of element
  void DeformationGradient_C(const LocalArrayT& disp, dMatrixT& F_3D);

 private:

  // number of crystal variables to be stored
  virtual int NumVariablesPerElement();
 
   // initial value of crystal variables
  virtual void InitializeCrystalVariables(ElementCardT& element);

 protected:

  // number of nodes/element: 4-node Quad (2D) and 8-node Hexa (3D)
  //const int fNNodes;
  int fNNodes;

  // references to initial coords
  const LocalArrayT* fLocInitX;

  // arrays for shape funtion derivatives at center
  dArray2DT fLNa;
  dArray2DT fLDNa;
  dArray2DT fGDNa;

  // displacement gradient
  dMatrixT fGradU;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_FP_C_H_ */
