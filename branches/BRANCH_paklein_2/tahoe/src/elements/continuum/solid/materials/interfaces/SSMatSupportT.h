/* $Id: SSMatSupportT.h,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#ifndef _FD_MAT_SUPPORT_T_H_
#define _FD_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class LocalArrayT;
class SmallStrainT;

/** support for the small strain Tahoe materials classes */
class SSMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	SSMatSupportT(int nsd, int ndof, int nip);

	/** destructor */
	~SSMatSupportT(void);

	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/

	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;
	/*@}*/

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const SmallStrainT* SmallStrain(void) const { return fSmallStrain; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);

	/** return a pointer the specified local array, or NULL if the array is not
	 * available. During calls the materials routines these will contain the
	 * values for the current element. */
	virtual const LocalArrayT* LocalArray(LocalArrayT::TypeT t) const;
	/*@}*/

  private:

  	/** \name return values */
	/*@{*/
  	ArrayT<dSymMatrixT> fStrain_List;
  	ArrayT<dSymMatrixT> fStrain_last_List;
	/*@}*/

  	/** pointer to the small strain element */
	const SmallStrainT* fSmallStrain;	
};

} /* namespace Tahoe */
#endif /* _FD_MAT_SUPPORT_T_H_ */
