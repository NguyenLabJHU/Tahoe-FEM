/*
  File GradientTools_C.h
*/

#ifndef _GRADIENT_TOOLS_C_H_
#define _GRADIENT_TOOLS_C_H_

#include <iostream.h>
#include <ctype.h>

#include "Utils.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"


namespace Tahoe {

class GradientTools_C 
{
 public:
   // constructor
   GradientTools_C(int numip, int numsd);

   // dNa/dE (fLDNa) -> dNa/dX (fGDNa): coords = Initial_X
   // dNa/dE (fLDNa) -> dNa/dx (fGDNa): coords = Current_X
   void ComputeGDNa(const LocalArrayT& coords);

   // interpolate variables at center of element
   void InterpolateTensorAtCenter(const ArrayT<dMatrixT>& ipvalues, dMatrixT& cvalue);

   // spatial gradient of a tensor (dA/dx) at center of element
   void GradientOfTensorAtCenter(const ArrayT<dMatrixT>& nodal, ArrayT<dMatrixT>& tensorGrad);

   // curl of a tensor at center of element
   void CurlOfTensorAtCenter(const ArrayT<dMatrixT>& Tnodal, dMatrixT& CurlT);

 private: 
   // set shape functions and their derivatives at center for 2D
   void SetLocalShape2D();

   // set shape functions and their derivatives at center for 3D
   void SetLocalShape3D();

   // jacobian for change of derivatives
   void Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNa);
   
   // dNa/dE -> dNa/dX : tensor = Jac^(-1) at Center
   // dNa/dX -> dNa/dx : tensor = F^(-1) at Center
   void ChangeOfVariables(const dMatrixT& tensor, const dArray2DT& DNa_known,
                          dArray2DT& DNa_unknown);

 private:
   // some dimensions
   int fNumIP;
   int fNumSD;

   // matrix of shape functions (parent domain)
   dArray2DT fNa;             // (1 x #IP)

   // array of matrices for shape funtion derivatives
   dArray2DT fLDNa;   // (#sd x #ip) (local, parent domain)
   dArray2DT fGDNa;   // (#sd x #ip) (global, physical domain)

   //workspace arrays
   dArrayT fArray1;           // (#sd) 
   dArrayT fArray2;           // (#sd)
 
   // workspace matrices
   dMatrixT fJac;             // (#sd x #sd)
   dMatrixT fMatx1;           // (#sd x #sd)

   // workspace for spatial gradient of a tensor T
   ArrayT<dMatrixT> fGradT;   // [knsd x (knsd x knsd)], knsd=3
};

} // namespace Tahoe 
#endif /* _GRADIENT_TOOLS_C_H_ */
