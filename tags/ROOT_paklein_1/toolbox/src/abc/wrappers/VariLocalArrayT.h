/* $Id: VariLocalArrayT.h,v 1.3 2002-07-05 22:26:22 paklein Exp $ */
/* created: paklein (04/26/1999)                                          */
/* Wrapper for dynamically re-sizing the number of nodes in               */
/* a LocalArrayT's.                                                       */

#ifndef _VARI_LOCALARRAY_T_H_
#define _VARI_LOCALARRAY_T_H_

/* base class */
#include "VariBaseT.h"

namespace Tahoe {

/* forward declarations */
class LocalArrayT;

class VariLocalArrayT: public VariBaseT<double>
{
public:

	/* constructors */
	VariLocalArrayT(void);
	VariLocalArrayT(int headroom, LocalArrayT& ward, int minordim);
	
	/* set the managed array - can only be set ONCE */
	void SetWard(int headroom, LocalArrayT& ward, int minordim);

	/* set number of nodes */
	void SetNumberOfNodes(int numnodes);

	/* return the ward */
	const LocalArrayT& TheWard(void) const;
	
private:

	/* dimensions */
	int fMinorDim;	
	
	/* the managed array */
	LocalArrayT* fWard;
};

/* inlines */

/* return the ward */
inline const LocalArrayT& VariLocalArrayT::TheWard(void) const
{
	return *fWard;
}

} // namespace Tahoe 
#endif /* _VARI_LOCALARRAY_T_H_ */