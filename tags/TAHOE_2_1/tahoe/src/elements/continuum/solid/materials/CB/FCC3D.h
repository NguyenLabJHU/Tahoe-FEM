/* $Id: FCC3D.h,v 1.5 2004-07-15 08:26:42 paklein Exp $ */
#ifndef _FCC_3D_H_
#define _FCC_3D_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class FCCLatticeT;
class PairPropertyT;
class BondLatticeT;

/** 3D Cauchy-Born material for FCC crystals with pair potential interactions. */
class FCC3D: public NL_E_MatT
{
public:

	/** constructor */
	FCC3D(void);	
	
	/** destructor */
	~FCC3D(void);

	/** \name Cauchy-Born parameters */
	/*@{*/
	/** return a reference to the bond lattice */
	const BondLatticeT& BondLattice(void) const;

	/** reference volume */
	double CellVolume(void) const { return fAtomicVolume; };

	/** nearest neighbor distance */
	double NearestNeighbor(void) const { return fNearestNeighbor; };
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/** symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/** strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

	/** return the equi-axed stretch at which the stress is zero. This method
	 * assumes the material is isotropic when subject to equi-axed stretch. */
	double ZeroStressStretch(void);

private:

	/** nearest neighbor distance */
	double fNearestNeighbor;

	/** bond information */
	FCCLatticeT* fFCCLattice;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** atomic volume */
	double fAtomicVolume;

	/** dummy full bond density array */
	dArrayT fFullDensity;
	
	/** flag to indicate whether stress calculation for output should include
	 * the full bond density */
	bool fFullDensityForStressOutput;
};

} /* namespace Tahoe */

#endif /* _FCC_3D_H_ */
