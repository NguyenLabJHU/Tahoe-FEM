/* $Id: ElementListT.h,v 1.5.12.1 2003-06-14 18:46:24 paklein Exp $ */
/* created: paklein (04/20/1998) */
#ifndef _ELEMENTLIST_T_H_
#define _ELEMENTLIST_T_H_

/* base class */
#include "pArrayT.h"
#include "iArrayT.h"

/* direct members */
#include "ElementSupportT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ElementBaseT;
class eIntegratorT;
class StringT;
class FEManagerT;
class ElementSupportT;

/** list of elements. Constructs list of element objects and
 * provides some attributes. */
class ElementListT: public pArrayT<ElementBaseT*>
{
public:

	/** constructor */
	ElementListT(void);

	/** destructor */
	~ElementListT(void);

	/** echo data from the I/O streams */
	void EchoElementData(ifstreamT& in, ostream& out, FEManagerT& fe);
	
	/** returns true of ALL element groups have interpolant DOF's */
	bool InterpolantDOFs(void) const;

	/** returns true if contact group present */
	bool HasContact(void) const;

	/** change the active element groups.
	 * \param mask list with length of the \e total number of element
	 *        groups with true|false determining whether the element
	 *        group is active. */
	void SetActiveElementGroupMask(const ArrayT<bool>& mask);

private:

	/** data needed for element contructors */
	ElementSupportT fSupport;

	/** cached pointers to element groups */
	ArrayT<ElementBaseT*> fAllElementGroups;
};

} /* namespace Tahoe */

#endif /* _ELEMENTLIST_T_H_ */
