/* $Id: EAM.h,v 1.9 2007-07-05 00:10:54 hspark Exp $ */
/* created: paklein (12/02/1996) */
#ifndef _EAM_H_
#define _EAM_H_

/* direct members */
#include "dMatrixT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class CBLatticeT;
class C1FunctionT;
class iArrayT;
class dSymMatrixT;

/** Cauchy-Born calculations for EAM potentials */
class EAM
{
public:

	/** constructor.
	* \param lattice crystal used to compute Cauchy-Born geometry information
	* \param nsd number of spatial dimensions needed for the stress and modulus calculations.
	*        This can be different from the dimensions of lattice, i.e., for 2D plane strain
	*        calculations nsd = 2, but the lattice will have 3 dimensions. */
	EAM(CBLatticeT& lattice);

	/* destructor */
	virtual ~EAM(void);

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

	/* compute unit 2nd PK surface stress:
	 *
	 *     unit 2nd PK surface stress = SIJ*(area per cell/atoms per cell)
	 */
	void ComputeUnitSurfaceStress(dSymMatrixT& stress);
	   	    	
	/* compute unit material tangent moduli:
	 *
	 *   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
	 */
	void ComputeUnitModuli(dMatrixT& moduli); 	    	

	/* unstressed lattice parameter */
	virtual double LatticeParameter(void) const = 0;

	/* atomic mass */
	virtual double Mass(void) const = 0;

	/** \name access glue functions */
	/*@{*/
	const C1FunctionT* PairPotential(void) { return fPairPotential; };
	const C1FunctionT* EmbeddingEnergy(void) { return fEmbeddingEnergy; };
	const C1FunctionT* ElectronDensity(void) { return fElectronDensity; };
	/*@}*/
	
	/** compute the total electron density */
	double TotalElectronDensity(void);
	
	/** calculate representative electron densities for bulk, surface1 and surface2 atoms */
	void ComputeElectronDensity(void);
	
private:

	/* form matrix of mixed pair potential and embedding
	 * energy derivatives.  NOTE: computes the UPPER triangle
	 * ONLY */
	void FormMixedDerivatives(double rho);	

	/* Moduli tensor contributions */
	void FormSingleBondContribution(double rho, dMatrixT& moduli);
	void FormMixedBondContribution(double rho, dMatrixT& moduli);

	/* set the glue function pointers - called by Initialize() */
	virtual void SetPairPotential(void) = 0;
	virtual void SetEmbeddingEnergy(void) = 0;
	virtual void SetElectronDensity(void) = 0; 	

protected:

	/** \name glue functions */
	/*@{*/
	C1FunctionT* fPairPotential;
	C1FunctionT* fEmbeddingEnergy;
	C1FunctionT* fElectronDensity;
	/*@}*/
	
private:   	

	CBLatticeT&	fLattice;
//	const iArrayT& fCounts;		
//	const dArrayT& fBonds;

	/* parameters */
//	int	fNumSpatialDim;
//	int	fNumBonds;
//	int	fModuliDim;
	
	dMatrixT fBondTensor4;
	dMatrixT fAmn; /* mixed derivative matrix */

	dArrayT fBondTensor2;
	dArrayT fBondTensor2b;

	/* 2nk rank bond component tensor */
	dArray2DT fTensor2Table;	

	/* interactiont table for surface clusters */
	dArray2DT fIntType;

	/* for batch evaluation of bond data */
	dArrayT	fBond1;
	dArrayT	fBond2;
	dArrayT	fBond3;
	dArrayT fBond4;
	dArrayT fBond5;
	dArrayT fBond6;
	dArrayT fRepRho;
};

} // namespace Tahoe 
#endif /* _EAM_H_ */
