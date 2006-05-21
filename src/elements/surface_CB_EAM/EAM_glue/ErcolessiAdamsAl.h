/* $Id: ErcolessiAdamsAl.h,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
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

} // namespace Tahoe 
#endif /* _ERCOLESSIADAMS_AL_H_ */
