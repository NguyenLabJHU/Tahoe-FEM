/* $Id: FBD_EAMGlue.h,v 1.1.1.1.10.1 2002-06-27 18:03:08 cjkimme Exp $ */
/* created: paklein (01/30/2000)                                          */
/* FBD_EAMGlue.h                                                          */
/* EAM glue functions in the form used in PRL v##, n##, 1986.             */

#ifndef _FBD_EAM_GLUE_H_
#define _FBD_EAM_GLUE_H_

/* base class */
#include "EAM.h"

/* forward declarations */

namespace Tahoe {

class ifstreamT;

class FBD_EAMGlue: public EAM
{
public:

	/* constructor */
	FBD_EAMGlue(CBLatticeT& lattice, ifstreamT& in);

	/* ustressed lattice parameter */
	 virtual double LatticeParameter(void) const;

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
