/* $Id: ErcolessiAdamsAl.h,v 1.1.1.1 2001-01-29 08:20:24 paklein Exp $ */
/* created: paklein (12/04/1996)                                          */
/* ErcolessiAdamsAl.h                                                     */

#ifndef _ERCOLESSIADAMS_AL_H_
#define _ERCOLESSIADAMS_AL_H_

/* base class */
#include "EAM.h"

class ErcolessiAdamsAl: public EAM
{
public:

	/* Constructor */
	ErcolessiAdamsAl(CBLatticeT& lattice);

	/*
	 * Unstressed lattice parameter.
	 */
	 virtual double LatticeParameter(void) const;

private:

	/*
	 * Set the spline data - called by the constructor
	 */
	virtual void SetPairPotential(void);
	virtual void SetEmbeddingEnergy(void);
	virtual void SetElectronDensity(void); 	
	
};

#endif /* _ERCOLESSIADAMS_AL_H_ */
