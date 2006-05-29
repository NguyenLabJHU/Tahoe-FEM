/* $Id: EAMFCC3D_surf.h,v 1.3 2006-05-29 17:22:56 paklein Exp $ */
/* created: paklein (12/02/1996) */
#ifndef _EAMFCC3D_SURF_H_
#define _EAMFCC3D_SURF_H_

/* base class */
#include "FCCLatticeT_Surf.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dMatrixT;
class dSymMatrixT;
class EAM;
class EAM_particle;

class EAMFCC3D_surf: public FCCLatticeT_Surf
{
public:

	/** EAM glue functions */
	enum GlueTypeT {kErcolessiAdamsAl = 0,
                         kVoterChenAl = 1,
                         kVoterChenCu = 2,
                     kFoilesBaskesDaw = 3,
					     kEAMParticle = 4};

	/** constructor */
	EAMFCC3D_surf(int nshells, int normal);

	/* destructor */
	virtual ~EAMFCC3D_surf(void);

	/* strain energy density */
	double EnergyDensity(const dSymMatrixT& strain);

	/* return the material tangent moduli in Cij */
	void Moduli(dMatrixT& Cij, const dSymMatrixT& strain);

	/* return the symmetric 2nd PK stress tensor */
	void SetStress(const dSymMatrixT& strain, dSymMatrixT& stress);

	/* calculate electron density at ghost atom */
	void ElectronDensity(const dSymMatrixT& strain, double& edensity, double& embforce);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* initialize bond table values */
	virtual void LoadBondTable(void);

protected:   	    	

	double	fLatticeParameter;
	double	fCellVolume;

	/** \name embedded atom solver */
	/*@{*/
	EAM* fEAM;
	EAM_particle* fEAM_particle;	
	/*@}*/
};

} // namespace Tahoe 

#endif /* _EAMFCC3D_SURF_H_ */
