/* $Id: SolidMatSupportT.h,v 1.2.2.1 2002-12-10 17:06:03 paklein Exp $ */
#ifndef _STRUCT_MAT_SUPPORT_T_H_
#define _STRUCT_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

namespace Tahoe {

/* forward declarations */
class ElasticT;

/** support for the small strain Tahoe materials classes */
class SolidMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	SolidMatSupportT(int nsd, int ndof, int nip);

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const ElasticT* Elastic(void) const { return fElastic; };

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;

	/** nodal temperatures. Returns NULL if not available */
	const LocalArrayT* Temperatures(void) const { return fTemperatures; };

	/** nodal temperatures from the last time step. Returns NULL if 
	 * not available */
	const LocalArrayT* LastTemperatures(void) const { return fLastTemperatures; };
	/*@}*/

	/** \name set host code information */
	/*@{*/
	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);

	/** set pointer local array */
	virtual void SetLocalArray(const LocalArrayT& array);

	void SetTemperatures(const LocalArrayT& temperatures);
	void SetLastTemperatures(const LocalArrayT& last_temperatures);
	/*@}*/

  private:

  	/** pointer to the solid element */
	const ElasticT* fElastic;
	
	/** \name pointers to local arrays */
	/*@{*/
	const LocalArrayT* fLastDisp;
	const LocalArrayT* fVel;
	const LocalArrayT* fAcc;

	const LocalArrayT* fTemperatures;
	const LocalArrayT* fLastTemperatures;
	/*@}*/
};

/* inlines */
inline void SolidMatSupportT::SetTemperatures(const LocalArrayT& temperatures)
{
	fTemperatures = &temperatures;
}

inline void SolidMatSupportT::SetLastTemperatures(const LocalArrayT& last_temperatures)
{
	fLastTemperatures = &last_temperatures;
}

} /* namespace Tahoe */
#endif /* _STRUCT_MAT_SUPPORT_T_H_ */