/* $Id: InverseMapT.h,v 1.10 2005-02-24 17:34:04 paklein Exp $ */
#ifndef _INVERSE_MAP_T_H_
#define _INVERSE_MAP_T_H_

/* base class */
#include "AutoArrayT.h"

/* forward declarations */
#include "nArrayT.h"

namespace Tahoe {

/** class to construct and access an inverse map. Given a forward map
 * \f$ j = g(i) \f$, the class constructs and handles access to the inverse
 * map \f$ i = g^{-1}(j) \f$. By default, InverseMapT::Map of a value not in 
 * the forward map throws an exception. This behavior can be changed by
 * setting the 
 */
class InverseMapT: private AutoArrayT<int>
{
public:

	/** enum for what to do about trying to map values that 
	 * are out of range */
	enum SettingT {
		Throw,   /**< throw ExceptionT::OutOfRange */
		MinusOne /**< return -1 */
	};

	/** constructor */
	InverseMapT(void);

	/** \name construct the inverse map */
	/*@{*/
	void SetMap(const nArrayT<int>& forward);
	void SetMap(const ArrayT<int>& forward);
	/*@}*/
	
	/** recover the forward map */
	void Forward(ArrayT<int>& forward) const;
	
	/** set the flag for handling calls to InverseMapT::Map that
	 * are out of range */
	void SetOutOfRange(SettingT setting) { fOutOfRange = setting; };
	
	/** return the out of range flag */
	SettingT OutOfRange(void) const { return fOutOfRange; };
	
	/** clear the map. Sets InverseMapT::Length to zero without
	 * necessarily freeing any memory. Use InverseMapT::Free to
	 * release allocated memory */
	void ClearMap(void) { Dimension(0); };

	/** map the global index to the local index */
	int Map(int global) const;
	
	/** release memory */
	void Free(void);
	
	/** return the logical size of the map */
	int Length(void) const { return AutoArrayT<int>::Length(); };

	/** return the number of entrees in the map */
	int Entrees(void) const { return fEntrees; };
	
private:

	/** minimum value in the forward map. The index shift allows some
	 * saving in memory since global tags less than the shift are not stored
	 * in the map. */	
	int fShift;	
	
	/** how to handle out of range */
	SettingT fOutOfRange;
	
	/** number of mapped entries */
	int fEntrees;
};

/* inlines */

/* constructor */
inline InverseMapT::InverseMapT(void): 
	fShift(0),
	fEntrees(0),
	fOutOfRange(Throw)
{

}

inline void InverseMapT::SetMap(const ArrayT<int>& forward) {
	nArrayT<int> fwd;
	fwd.Alias((ArrayT<int>&) forward);
	SetMap(fwd);
}

/* map the global index to the local index */
inline int InverseMapT::Map(int global) const
{
	int dex = global - fShift;
	int map = -1;
	if (dex > -1 && dex < fLength) map = (*this)[dex];

	if (map == -1 && fOutOfRange == Throw) 
		ExceptionT::OutOfRange("InverseMapT::Map", "%d was not in the forward map", global);
	return map;
}

/* release memory */
inline void InverseMapT::Free(void) { 
	AutoArrayT<int>::Free(); 
	fShift = 0;
	fEntrees = 0;
}

} /* namespace Tahoe */
 
#endif /* _INVERSE_MAP_T_H_ */
