#ifndef _PARADYN_T_H_
#define _PARADYN_T_H_

/* direct members */
#include "GeometryT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class ifstreamT;
class StringT;
class iArrayT;
class dArrayT;
class iArray2DT;
class dArray2DT;
template <class TYPE> class AutoArrayT;

class ParaDynT
{
public:
enum VariableTypeT {kScalarElemental,
		    kVectorElemental,
		    kScalarNodal,
		    kVectorNodal};

 ParaDynT (ostream& fMainOut);

 // Compute
 dArrayT ComputeLatticeParameter(dArray2DT local);
 dArray2DT ComputeMinMax(dArray2DT local);
 dArray2DT ComputeBounds(const iArrayT& periodic_cond,
			 dArray2DT local,
			 const dArray2DT& BoundsCoord,
			 const dArrayT& lattice_parameter);

 // geometry
 void WriteHeader (ostream& fgeo, ArrayT<StringT>& header) const;
 
 void WriteBoundHeader (ostream& fgeo) const;
 void WriteBounds (ostream& fgeo, const dArray2DT& bounds) const;

 void WriteCoordinateHeader (ostream& fgeo) const;
 void WriteCoordinates (ostream& fgeo, 
			const dArray2DT& coords,
			const iArrayT& types,
			const iArrayT& parts) const;
 
 void WriteTime (ostream& fvar) const;
 
 private:
 ostream& fOut;    // error messaging
};

} // namespace Tahoe 
#endif
