/* $Id: VIB.h,v 1.4 2002-10-20 22:48:54 paklein Exp $ */
/* created: paklein (10/30/1997)                                          */
/* Base class for isotropic VIB solvers.                                  */

#ifndef _VIB_H_
#define _VIB_H_

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class C1FunctionT;
class dSymMatrixT;

class VIB
{
public:

	/* constructor */
	VIB(ifstreamT& in, int nsd, int numstress, int nummoduli);

	/* destructor */
	virtual ~VIB(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	

protected:

	/* allocate memory for all the tables */
	void Dimension(int numbonds);

protected:

	/* number of spatial dimensions */
	int fNumSD;

	/* potential function */
	C1FunctionT* fPotential;

	/* length table */
	dArrayT	fLengths;

	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	/* jacobian table */
	dArrayT	fjacobian;

	/* STRESS angle tables - by associated stress component */
	int fNumStress;
	dArray2DT fStressTable;
	  	
	/* MODULI angle tables */
	int fNumModuli; 	
	dArray2DT fModuliTable;	
};

} // namespace Tahoe 
#endif /* _VIB_H_ */
