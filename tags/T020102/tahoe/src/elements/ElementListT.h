/* $Id: ElementListT.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (04/20/1998)                                          */

#ifndef _ELEMENTLIST_T_H_
#define _ELEMENTLIST_T_H_

/* base class */
#include "pArrayT.h"
#include "iArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class FEManagerT;
class ElementBaseT;
class eControllerT;
class StringT;

class ElementListT: public pArrayT<ElementBaseT*>
{
public:

	/* constructor */
	ElementListT(FEManagerT& fe_manager);
	
	/* echo data from the I/O streams */
	void EchoElementData(ifstreamT& in, ostream& out,
		eControllerT* e_controller);
	
	/* returns true of ALL element groups have interpolant DOF's */
	bool InterpolantDOFs(void) const;

	/* returns true if contact group present */
	bool HasContact(void) const;

private:

	/* data needed for element contructors */
	FEManagerT& fFEManager;
};

#endif /* _ELEMENTLIST_T_H_ */
