/* $Id: Utils.h,v 1.3 2003-01-27 07:00:29 paklein Exp $ */
#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream.h>
#include "ArrayT.h"

namespace Tahoe {

class ifstreamT;
class StringT;
class LocalArrayT;
class dArrayT;
class dArray2DT;
class dMatrixT;
class iArrayT;

// auxiliar functions for input/output
ifstreamT& OpenExternal(ifstreamT& in, ifstreamT& in2, const char* name);
void OpenExternal(ifstreamT& input, StringT& filename, const char* name);
void SetStreamPrefs(ostream& stream);

// auxiliar function to compute deformation gradient at center of element
void SetLocalShape_C(dArray2DT& Na, dArray2DT& DNa); 
void ComputeGDNa_C(const LocalArrayT& coords, const dArray2DT& LDNa, dArray2DT& GDNa);

// auxiliar functions to compute spatial gradient terms at IPs
void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x, dArray2DT& nodal_extrap);
void ComputeGDNa(const LocalArrayT& coords, const ArrayT<dArray2DT>& LDNa, ArrayT<dArray2DT>& GDNa);
void Jacobian(const LocalArrayT& nodal, const dArray2DT& DNa, dMatrixT& jac);
void ChangeOfVariables(const dMatrixT& tensor, const dArray2DT& DNa_known, dArray2DT& DNa_unknown);
void Interpolate(const ArrayT<dMatrixT>& nodal, ArrayT<dMatrixT>& interp, const dArray2DT& Na);
void Extrapolate(const ArrayT<dMatrixT>& ipvalues, ArrayT<dMatrixT>& extrap, const dArray2DT& Ma);

// auxiliar functions for matrix inversion (from Numerical Recipes)
void LUDecomposition(dMatrixT& a, const int& n, iArrayT& indx, double& d);
void LUBackSubstitution(dMatrixT& a, const int& n, const iArrayT& indx, dArrayT& b);
dMatrixT MatrixInversion(dMatrixT& a);

// auxiliar functions to print general messages, warnings, and error messages
void writeMessage(const char* msg);
void writeWarning(const char* msg);
void throwRunTimeError(const char* msg);
void throwMemoryError(const char* msg);

// some macros
#ifndef max
inline int max(int i1, int i2) {return i1 >= i2 ? i1 : i2;}
inline double max(double d1, double d2) {return d1 >= d2 ? d1 : d2;}
#endif

#ifndef min
inline int min(int i1, int i2) {return i1 <= i2 ? i1 : i2;}
inline double min(double d1, double d2) {return d1 <= d2 ? d1 : d2;}
#endif

#ifndef sign
inline int sign(int x) {return x > 0 ? 1 : (x == 0 ? 0 : -1);}
#endif

} // namespace Tahoe 
#endif /* _UTILS_H_ */

