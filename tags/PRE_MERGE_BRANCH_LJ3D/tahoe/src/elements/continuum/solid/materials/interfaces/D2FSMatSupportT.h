/* $Id: D2FSMatSupportT.h,v 1.4 2003-01-29 07:34:57 paklein Exp $ */
#ifndef _D2_FD_MAT_SUPPORT_T_H_
#define _D2_FD_MAT_SUPPORT_T_H_

/* base class */
#include "FSMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class D2MeshFreeFSSolidT;

/** support for the 2nd gradient, finite strain Tahoe materials classes */
class D2FSMatSupportT: public FSMatSupportT
{
  public:

	/** constructor */
	D2FSMatSupportT(int nsd, int ndof, int nip);

	/** \name host code information */
	/*@{*/
	/** return a pointer to the host element. Returns NULL if no
	 * no element information in available. The ContinuumElementT
	 * pointer is set using MaterialSupportT::SetContinuumElement. */
	const D2MeshFreeFSSolidT* D2MeshFreeFDElastic(void) const { return fD2MeshFreeFDElastic; };

	/** set the element group pointer */
	virtual void SetContinuumElement(const ContinuumElementT* p);
	/*@}*/

  private:
  
	/** pointer to the host element */
	const D2MeshFreeFSSolidT* fD2MeshFreeFDElastic;
};

} /* namespace Tahoe */
#endif /* _D2_FD_MAT_SUPPORT_T_H_ */
