/* $Id: VIB.h,v 1.5 2004-07-15 08:27:40 paklein Exp $ */
/* created: paklein (10/30/1997) */
#ifndef _VIB_H_
#define _VIB_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class C1FunctionT;
class dSymMatrixT;

/** base class for Virtual Internal Bond calculations */
class VIB: virtual public ParameterInterfaceT
{
public:

	/** constructor */
	VIB(int nsd, int numstress, int nummoduli);

	/* destructor */
	virtual ~VIB(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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
