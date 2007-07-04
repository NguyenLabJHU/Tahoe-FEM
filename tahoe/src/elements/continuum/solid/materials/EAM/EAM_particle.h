/* $Id: EAM_particle.h,v 1.4 2007-07-04 17:51:54 paklein Exp $ */
/* created: hspark(02/25/2004) */
#ifndef _EAM_PARTICLE_H_
#define _EAM_PARTICLE_H_

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "EAMPropertyT.h"

namespace Tahoe {

/* forward declarations */
class CBLatticeT;
class iArrayT;
class dSymMatrixT;

/** EAM calculations for Cauchy-Born constitutive models using the EAMPropertyT
 * potential functions */
class EAM_particle
{
public:

	/* constructor */
	EAM_particle(CBLatticeT& lattice, const StringT& param_file);

	/* destructor */
	virtual ~EAM_particle(void);

	/** set "glue" functions and dimension work space */
	void Initialize(int nsd, int numbonds);

	/* compute unit strain energy density:
	 *
	 *     unit strain energy = energy/atom
	 */
	double ComputeUnitEnergy(void);
	
	/* compute unit 2nd PK stress:
	 *
	 *     unit 2nd PK stress = SIJ*(volume per cell/atoms per cell)
	 */
	void ComputeUnitStress(dSymMatrixT& stress);
	   	    	
	/* compute unit material tangent moduli:
	 *
	 *   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
	 */
	void ComputeUnitModuli(dMatrixT& moduli); 	    	

	/* unstressed lattice parameter */
	double LatticeParameter(void) const { return fLatticeParameter; };

	/** atomic mass */
	double Mass(void) const { return fMass; }

	/* compute the total electron density - moved public by HSP 3/5/04 */
	double TotalElectronDensity(void);

	/* return the embedding force for a given electron density */
	double ReturnEmbeddingForce(double rho);

private:

	/* form matrix of mixed pair potential and embedding
	 * energy derivatives.  NOTE: computes the UPPER triangle
	 * ONLY */
	void FormMixedDerivatives(double rho);	

	/* Moduli tensor contributions */
	void FormSingleBondContribution(double rho, dMatrixT& moduli);
	void FormMixedBondContribution(double rho, dMatrixT& moduli);
	
private:   	

	/** EAM property and function pointers */
	/*@{*/
	EAMPropertyT* fEAMProperty;
	EAMPropertyT::PairEnergyFunction    fPairEnergy;
	EAMPropertyT::PairForceFunction     fPairForce;
	EAMPropertyT::PairStiffnessFunction fPairStiffness;
	EAMPropertyT::EmbedStiffnessFunction fEmbedStiffness;	
	EAMPropertyT::EmbedEnergyFunction   fEmbedEnergy;
	EAMPropertyT::EmbedForceFunction fEmbedForce;
	EAMPropertyT::EDEnergyFunction fEDEnergy;
	EAMPropertyT::EDForceFunction fEDForce;
	EAMPropertyT::EDStiffnessFunction fEDStiffness;
	/*@{*/

	CBLatticeT&	fLattice;
//	const iArrayT&	fCounts;
//	const dArrayT&	fBonds;

	/* parameters */
//	int		fNumSpatialDim;
//	int		fNumBonds;
//	int		fModuliDim;
	double  fLatticeParameter;
	double  fMass;
	
	dMatrixT	fBondTensor4;
	dMatrixT	fAmn; /* mixed derivative matrix */

	dArrayT		fBondTensor2;
//	dArrayT		fBondTensor2b;

	/* 2nk rank bond component tensor */
	dArray2DT	fTensor2Table;	

	/* for batch evaluation of bond data */
	dArrayT	fBond1;
	dArrayT	fBond2;
	dArrayT	fBond3;
};

} // namespace Tahoe 
#endif /* _EAM_PARTICLE_H_ */
