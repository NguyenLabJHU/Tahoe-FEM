/* $Id: ErcolessiAdamsAl.h,v 1.4 2007-07-05 00:03:19 paklein Exp $ */
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

	/** unstressed lattice parameter */
	 virtual double LatticeParameter(void) const;

	/** atomic mass */
	 virtual double Mass(void) const;

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
