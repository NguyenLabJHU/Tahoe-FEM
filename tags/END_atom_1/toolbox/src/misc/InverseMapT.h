/* $Id: InverseMapT.h,v 1.1.2.1 2002-12-16 09:13:47 paklein Exp $ */
#ifndef _INVERSE_MAP_T_H_
#define _INVERSE_MAP_T_H_

/* base class */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class nArrayT;

/** class to construct and access an inverse map. Given a forward map
 * \f$ j = g(i) \f$, the class constructs and handles access to the inverse
 * map \f$ i = g^{-1}(j) \f$. Inverse of a value not in the forward map
 * throws an exception.
 */
class InverseMapT: private AutoArrayT<int>
{
public:

	/** constructor */
	InverseMapT(void);

	/** construct the inverse map */
	void SetMap(const nArrayT<int>& forward);
	
	/** clear the map. Sets InverseMapT::Length to zero without
	 * necessarily freeing any memory. Use InverseMapT::Free to
	 * release allocated memory */
	void ClearMap(void) { Dimension(0); };

	/** map the global index to the local index */
	int Map(int global) const;
	
	/** release memory */
	void Free(void);
	
	/** return the logical size of the map */
	int Length(void) { return AutoArrayT<int>::Length(); };
	
private:

	/** minimum value in the forward map. The index shift allows some
	 * saving in memory since global tags less than the shift are not stored
	 * in the map. */	
	int fShift;	
};

/* inlines */

/* constructor */
inline InverseMapT::InverseMapT(void): fShift(0) {}

/* map the global index to the local index */
inline int InverseMapT::Map(int global) const
{
	int map = (*this)[global - fShift];
	if (map == -1) 
		ExceptionT::GeneralFail("InverseMapT::Map", "%d was not in the forward map", global);
	return map;
}

/* release memory */
inline void InverseMapT::Free(void) { AutoArrayT<int>::Free(); }

} /* namespace Tahoe */
 
#endif /* _INVERSE_MAP_T_H_ */
