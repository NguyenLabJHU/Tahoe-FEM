/* $Id: SWDataT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/22/1997)                                          */
/* Container class for Stillinger-Weber potential parameters              */

#ifndef _SW_DATA_T_H_
#define _SW_DATA_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;

class SWDataT
{
/* classes that use this container */
friend class MixedSWDiamondT;
friend class ModCB2D_SWT;
friend class SW2BodyT;
friend class SW3BodyT;

public:

	/* Constructors */
	SWDataT(void);
	SWDataT(ifstreamT& in);

	/* I/O functions */
	void Read(ifstreamT& in);
	void Write(ostream& out) const;

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

#endif /* _SW_DATA_T_H_ */
