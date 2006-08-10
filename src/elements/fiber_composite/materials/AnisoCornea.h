/* $Id: AnisoCornea.h,v 1.2 2006-08-10 01:46:53 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ANISO_CORNEA_2D_H_
#define _ANISO_CORNEA_2D_H_

/* base classes */
#include "FSFiberMatT.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;
class C1FunctionT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class AnisoCornea: public FSFiberMatT
{
public:

	/* constructor */
	AnisoCornea(void);

	/* destructor */
	~AnisoCornea(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

//	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, SubListT& sub_lists) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:
	/*calculates  matrix contribution to 2PK stress*/
	virtual void ComputeMatrixStress(const dSymMatrixT& C, dSymMatrixT& Stress);

	/*calculates matrix contribution to modulus*/
	virtual void ComputeMatrixMod(const dSymMatrixT& C, dSymMatrixT& Stress, dMatrixT& Mod);
	
	/*computes integrated fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress);
	
	/*computes integrated moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress,
					dSymMatrixT& FiberMod);

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dSymMatrixT& stretch);
	
private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/*constitutive values for matrix*/
	double fMu;
	double fGamma;
	
	/* integration point generator */
	CirclePointsT*	fCircle;
	
	/* potential function */
	C1FunctionT* fPotential;

	/* fibril distribution function */
	C1FunctionT* fDistribution;

	/* length table */
	/*I4*/
	dArrayT	fI4;

	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	/* jacobian table */
	dArrayT	fjacobian;

	/* STRESS angle tables for fiber stress - by associated stress component */
	dArray2DT fStressTable;
	  	
	/* MODULI angle tables for fiber moduli */
	dArray2DT fModuliTable;	

};

} // namespace Tahoe 
#endif /* _ISO_VIB_2D_H_ */
