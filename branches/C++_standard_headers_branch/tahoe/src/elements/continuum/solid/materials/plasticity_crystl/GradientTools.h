/*
  File GradientTools.h
*/

#ifndef _GRADIENT_TOOLS_H_
#define _GRADIENT_TOOLS_H_

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

class GradientTools 
{
 public:
   // constructor
   GradientTools(int numip, int numnodes, int numsd);

   // dNa/dE (fLDNa) -> dNa/dX (fGDNa): coords = Initial_X
   // dNa/dE (fLDNa) -> dNa/dx (fGDNa): coords = Current_X
   void ComputeGDNa(const LocalArrayT& coords);

   // interpolate variables at IPs
   void InterpolateTensor(const ArrayT<dMatrixT>& npvalues, ArrayT<dMatrixT>& ipvalues);

   // extrapolate variables to nodes
   void ExtrapolateTensor(const ArrayT<dMatrixT>& ipvalues, ArrayT<dMatrixT>& npvalues);

   // spatial gradient of a tensor (dA/dx) at integration point
   void GradientOfTensorAtIP(const ArrayT<dMatrixT>& nodal, ArrayT<dMatrixT>& tensorGrad, 
                             int ip);

   // curl of a tensor at integration point
   void CurlOfTensorAtIP(const ArrayT<dMatrixT>& Tnodal, dMatrixT& CurlT, int ip);

 private: 
   // set shape functions and their derivatives at IPs for 2D
   void SetLocalShape2D();

   // set shape functions and their derivatives at IPs for 3D
   void SetLocalShape3D();

   // set extrapolation matrix of ipvariables to npvariables
   void SetNodalExtrapolation();

   // jacobian for change of derivatives
   void Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNa);
   
   // dNa/dE -> dNa/dX : tensor = Jac^(-1) at IP
   // dNa/dX -> dNa/dx : tensor = F^(-1) at IP
   void ChangeOfVariables(const dMatrixT& tensor, const dArray2DT& DNa_known,
                          dArray2DT& DNa_unknown);

   void ChangeOfVariables();

 private:
   // some dimensions
   int fNumIP;
   int fNumNodes;
   int fNumSD;

   // matrix of shape functions (parent domain)
   dArray2DT fNa;             // (#IP x #nodes)

   // array of matrices for shape funtion derivatives
   ArrayT<dArray2DT> fLDNa;   // [#IP x (#sd x #nodes)] (local, parent domain)
   ArrayT<dArray2DT> fGDNa;   // [#IP x (#sd x #nodes)] (global, physical domain)

   // matrix for nodal extrapolation (least square smoothing)
   dArray2DT fNodalExtrap;    // (#nodes x #IP) 

   // weights at IP
   dArrayT fWeights;

   //workspace arrays
   dArrayT fArray1;           // (#sd) 
   dArrayT fArray2;           // (#sd)
 
   // workspace matrices
   dMatrixT fJac;             // (#sd x #sd)
   dMatrixT fMatx1;           // (#sd x #sd)

   // workspace for spatial gradient of a tensor T
   ArrayT<dMatrixT> fGradT;   // (knsd x knsd), knsd=3
};

} // namespace Tahoe 
#endif /* _GRADIENT_TOOLS_H_ */
