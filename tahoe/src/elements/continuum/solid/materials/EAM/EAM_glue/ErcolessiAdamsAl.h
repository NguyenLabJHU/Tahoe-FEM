/* $Id: ErcolessiAdamsAl.h,v 1.2.56.1 2004-06-16 00:31:53 paklein Exp $ */
/* created: paklein (12/04/1996) */
#ifndef _ERCOLESSIADAMS_AL_H_
#define _ERCOLESSIADAMS_AL_H_

/* base class */
#include "EAM.h"

namespace Tahoe {

/** Ercolessi and Adams EAM aluminum potentials */
class ErcolessiAdamsAl: public EAM
{
public:

	/* Constructor */
	ErcolessiAdamsAl(CBLatticeT& lattice, int nsd);

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

} // namespace Tahoe 
#endif /* _ERCOLESSIADAMS_AL_H_ */
