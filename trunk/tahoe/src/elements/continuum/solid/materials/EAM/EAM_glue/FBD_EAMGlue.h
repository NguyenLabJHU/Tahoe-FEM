/* $Id: FBD_EAMGlue.h,v 1.5 2007-07-05 00:03:19 paklein Exp $ */
/* created: paklein (01/30/2000) */
#ifndef _FBD_EAM_GLUE_H_
#define _FBD_EAM_GLUE_H_

/* base class */
#include "EAM.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** EAM glue functions in the form Paradyn form */
class FBD_EAMGlue: public EAM
{
public:

	/* constructor */
	FBD_EAMGlue(CBLatticeT& lattice, ifstreamT& in);

	/* ustressed lattice parameter */
	 virtual double LatticeParameter(void) const { return fLatticeParameter; }

	/* atomic mass */
	 virtual double Mass(void) const { return fMass; };

private:

	/* set glue */
	virtual void SetPairPotential(void);
	virtual void SetEmbeddingEnergy(void);
	virtual void SetElectronDensity(void); 	

private:

	/* lattice parameters */
	double fLatticeParameter;
	
	/* atomic mass (amu) */
	double fMass;
};

} // namespace Tahoe 
#endif /* _FBD_EAM_GLUE_H_ */
