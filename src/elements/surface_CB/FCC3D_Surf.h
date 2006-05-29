/* $Id: FCC3D_Surf.h,v 1.7 2006-05-29 17:22:15 paklein Exp $ */
#ifndef _FCC_3D_SURF_H_
#define _FCC_3D_SURF_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class FCCLatticeT_Surf;
class PairPropertyT;
class BondLatticeT;

/** 3D Cauchy-Born material for FCC crystals with pair potential interactions. */
class FCC3D_Surf: public NL_E_MatT
{
public:

	/** constructor */
	FCC3D_Surf(void);	
	
	/** destructor */
	~FCC3D_Surf(void);

	/** \name Cauchy-Born parameters */
	/*@{*/
	/** return a reference to the bond lattice */
	const BondLatticeT& BondLattice(void) const;

	/** reference volume */
	double CellVolume(void) const { return fAtomicVolume; };

	/** nearest neighbor distance */
	double NearestNeighbor(void) const { return fNearestNeighbor; };
	/*@}*/

	/** thickness of surface layer to subtract off of bulk */
	double SurfaceThickness(void) const { return fSurfaceThickness; };

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

//TEMP
/** spatial tangent modulus */
//virtual const dMatrixT& c_ijkl(void) { return FSSolidMatT::c_ijkl(); };

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

	/** surface layer thickness */
	double fSurfaceThickness;

	/** bond information */
	FCCLatticeT_Surf* fFCCLattice_Surf;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** atomic volume */
	double fAtomicVolume;

	/** atomic area for surface cauchy-born */
	double fAtomicArea;

	/** dummy full bond density array */
	/* THIS IS THE REFERENCE VOLUME/AREA */
	dArrayT fFullDensity;
	
	/** flag to indicate whether stress calculation for output should include
	 * the full bond density */
	bool fFullDensityForStressOutput;
};

} /* namespace Tahoe */

#endif /* _FCC_3D_SURF_H_ */
