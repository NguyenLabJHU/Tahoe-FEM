/* $Id: SWDataT.h,v 1.4 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (03/22/1997) */
/* Container class for Stillinger-Weber potential parameters */
#ifndef _SW_DATA_T_H_
#define _SW_DATA_T_H_

/* base class */
#include "ParameterInterfaceT.h"

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

class SWDataT: public ParameterInterfaceT
{
/* classes that use this container */
friend class MixedSWDiamondT;
friend class ModCB2D_SWT;
friend class SW2BodyT;
friend class SW3BodyT;

public:

	/** \name constructors
	 * construct either from stream of using the ParameterInterfaceT interface */
	/*@{*/
	SWDataT(void);
	SWDataT(ifstreamT& in);
	/*@}*/

	/* I/O functions */
	void Read(ifstreamT& in);
	void Write(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
 
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* unit scaling */
	double	feps;

	/* 2 body potential */
	double	fA;
	double	fdelta;
	 	
	/* 3 body potential */
	double	fgamma;
	double	flambda;
	
	double	frcut;		/* cut-off distance */
	double	fa;			/* lattice spacing parameter */
	
	/* derived values */
	double	fB;
};

} /* namespace Tahoe */
#endif /* _SW_DATA_T_H_ */
