/* $Id: SolidMatSupportT.h,v 1.2 2002-11-14 17:06:21 paklein Exp $ */
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

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;

	/** nodal temperatures. Returns NULL if not available */
	const LocalArrayT* Temperatures(void) const;

	/** nodal temperatures from the last time step. Returns NULL if 
	 * not available */
	const LocalArrayT* LastTemperatures(void) const;
	/*@}*/

  private:

  	/** pointer to the solid element */
	const ElasticT* fElastic;	
};

} /* namespace Tahoe */
#endif /* _STRUCT_MAT_SUPPORT_T_H_ */
