/* $Id: VariLocalArrayT.h,v 1.4 2003-01-27 06:42:46 paklein Exp $ */
/* created: paklein (04/26/1999) */
#ifndef _VARI_LOCALARRAY_T_H_
#define _VARI_LOCALARRAY_T_H_

/* base class */
#include "VariBaseT.h"

namespace Tahoe {

/* forward declarations */
class LocalArrayT;

/** Wrapper for dynamically re-sizing the number of nodes in
 * a LocalArrayT's */
class VariLocalArrayT: public VariBaseT<double>
{
public:

	/** \name constructors */
	/*@{*/
	VariLocalArrayT(void);
	VariLocalArrayT(int headroom, LocalArrayT& ward, int minordim);
	/*@}*/
	
	/** set the managed array - can only be set ONCE */
	void SetWard(int headroom, LocalArrayT& ward, int minordim);

	/** return true if the ward is already set */
	bool HasWard(void) const { return fWard != NULL; };

	/** set number of nodes */
	void SetNumberOfNodes(int numnodes);

	/** return the ward */
	const LocalArrayT& TheWard(void) const;
	
private:

	/** dimensions */
	int fMinorDim;	
	
	/** the managed array */
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
