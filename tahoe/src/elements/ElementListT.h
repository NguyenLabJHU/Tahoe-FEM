/* $Id: ElementListT.h,v 1.1.1.1.8.4 2002-04-30 08:21:59 paklein Exp $ */
/* created: paklein (04/20/1998) */

#ifndef _ELEMENTLIST_T_H_
#define _ELEMENTLIST_T_H_

/* base class */
#include "pArrayT.h"
#include "iArrayT.h"

/* direct members */
#include "ElementSupportT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class ElementBaseT;
class eControllerT;
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

	/* echo data from the I/O streams */
	void EchoElementData(ifstreamT& in, ostream& out, FEManagerT& fe);
	
	/* returns true of ALL element groups have interpolant DOF's */
	bool InterpolantDOFs(void) const;

	/* returns true if contact group present */
	bool HasContact(void) const;

private:

	/* data needed for element contructors */
	ElementSupportT fSupport;
};

#endif /* _ELEMENTLIST_T_H_ */
