/* $Id: EAMFCC3D.h,v 1.4 2004-04-09 02:02:58 hspark Exp $ */
/* created: paklein (12/02/1996)                                          */
/* EAMFCC3D.h                                                             */

#ifndef _EAMFCC3D_H_
#define _EAMFCC3D_H_

/* base class */
#include "CBLatticeT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dMatrixT;
class dSymMatrixT;
class EAM;
class EAM_particle;

/* bond parameters */
const int kEAMFCC3DNumBonds			= 54;
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;

class EAMFCC3D: public CBLatticeT
{
public:

	/* EAM glue functions */
	enum GlueTypeT {kErcolessiAdamsAl = 0,
                         kVoterChenAl = 1,
                         kVoterChenCu = 2,
                     kFoilesBaskesDaw = 3,
					     kEAMParticle = 4};

	/* constructor */
	EAMFCC3D(ifstreamT& in, int EAMcode, int numspatialdim,
		int numbonds = kEAMFCC3DNumBonds);
	EAMFCC3D(ifstreamT& in, const dMatrixT& Q, int EAMcode, int numspatialdim,
		int numbonds = kEAMFCC3DNumBonds);

	/* destructor */
	virtual ~EAMFCC3D(void);

	/* strain energy density */
	double EnergyDensity(const dSymMatrixT& strain);

	/* return the material tangent moduli in Cij */
	void Moduli(dMatrixT& Cij, const dSymMatrixT& strain);

	/* return the symmetric 2nd PK stress tensor */
	void SetStress(const dSymMatrixT& strain, dSymMatrixT& stress);

	/* I/O functions */
	virtual void Print(ostream& out) const;

	/* calculate electron density at ghost atom */
	void ElectronDensity(const dSymMatrixT& strain, double& edensity, double& embforce);

	/* initialize bond tables */
	void InitBondTables(void);

protected:

	/* initialize bond table values */
	virtual void LoadBondTable(void);

private:

	/* Set glue functions */
	void SetGlueFunctions(ifstreamT& in);
	 	   	    	
protected:   	    	

	int		fEAMcode;
	double	fLatticeParameter;
	double	fCellVolume;
	
	/* embedded atom solvers */
	EAM*	fEAM;
	EAM_particle* fEAM_particle;
};

} // namespace Tahoe 
#endif /* _EAMFCC3D_H_ */
