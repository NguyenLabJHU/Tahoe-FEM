/* $Id: ElementListT.h,v 1.1.1.1.8.2 2002-04-29 02:45:09 paklein Exp $ */
/* created: paklein (04/20/1998) */

#ifndef _ELEMENTLIST_T_H_
#define _ELEMENTLIST_T_H_

/* base class */
#include "pArrayT.h"
#include "iArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ElementBaseT;
class eControllerT;
class StringT;
class ElementSupportT;

/** list of elements. Constructs list of element objects and
 * provides some attributes. */
class ElementListT: public pArrayT<ElementBaseT*>
{
public:

	/* constructor */
	ElementListT(const ElementSupportT& support);
	
	/* echo data from the I/O streams */
	void EchoElementData(ifstreamT& in, ostream& out,
		eControllerT* e_controller);
	
	/* returns true of ALL element groups have interpolant DOF's */
	bool InterpolantDOFs(void) const;

	/* returns true if contact group present */
	bool HasContact(void) const;

private:

	/* data needed for element contructors */
	const ElementSupportT& fSupport;
};

#endif /* _ELEMENTLIST_T_H_ */
