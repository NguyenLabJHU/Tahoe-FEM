/* $Id: EAMFCC3D.h,v 1.3.56.1 2004-06-16 00:31:52 paklein Exp $ */
/* created: paklein (12/02/1996) */
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

class EAMFCC3D: public CBLatticeT
{
public:

	/* EAM glue functions */
	enum GlueTypeT {kErcolessiAdamsAl = 0,
                         kVoterChenAl = 1,
                         kVoterChenCu = 2,
                     kFoilesBaskesDaw = 3};

	/* constructor */
	EAMFCC3D(ifstreamT& in, int EAMcode, int nsd);
//	EAMFCC3D(ifstreamT& in, const dMatrixT& Q, int EAMcode, int nsd);

	/* destructor */
	virtual ~EAMFCC3D(void);

	/* strain energy density */
	double EnergyDensity(const dSymMatrixT& strain);

	/* return the material tangent moduli in Cij */
	void Moduli(dMatrixT& Cij, const dSymMatrixT& strain);

	/* return the symmetric 2nd PK stress tensor */
	void SetStress(const dSymMatrixT& strain, dSymMatrixT& stress);

protected:

	/* initialize bond table values */
	virtual void LoadBondTable(void);

private:

	/* Set glue functions */
	void SetGlueFunctions(ifstreamT& in, int EAMcode, int nsd);
	 	   	    	
protected:   	    	

	double	fLatticeParameter;
	double	fCellVolume;
	
	/* embedded atom solver */
	EAM* fEAM;	
};

} // namespace Tahoe 

#endif /* _EAMFCC3D_H_ */
