/* $Id: EAM.h,v 1.3.42.2 2004-03-06 01:22:47 hspark Exp $ */
/* created: paklein (12/02/1996)                                          */
/* EAM.h                                                                  */

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

class EAM
{
public:

	/* constructor */
	EAM(CBLatticeT& lattice);

	/* destructor */
	virtual ~EAM(void);

	/* set "glue" functions */
	void SetGlueFunctions(void);

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
	 virtual double LatticeParameter(void) const = 0;

	/** \name access glue functions */
	/*@{*/
	const C1FunctionT* PairPotential(void) { return fPairPotential; };
	const C1FunctionT* EmbeddingEnergy(void) { return fEmbeddingEnergy; };
	const C1FunctionT* ElectronDensity(void) { return fElectronDensity; };
	/*@}*/
	
	/* compute the total electron density - moved public by HSP 3/5/04 */
	double TotalElectronDensity(void);
	
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

	/* glue functions */
	C1FunctionT*	fPairPotential;
	C1FunctionT*	fEmbeddingEnergy;
	C1FunctionT*	fElectronDensity;
	
private:   	

	CBLatticeT&	fLattice;
	const iArrayT&	fCounts;		
	const dArrayT&	fBonds;

	/* parameters */
	int		fNumSpatialDim;
	int		fNumBonds;
	int		fModuliDim;
	
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
#endif /* _EAM_H_ */