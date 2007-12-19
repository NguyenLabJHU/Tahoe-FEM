/* $Id: AnisoFiber3D.h,v 1.1 2007-12-19 23:35:54 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ANISO_FIBER_3D_H_
#define _ANISO_FIBER_3D_H_

/* base classes */
#include "FSFiberMatT.h"
#include "SolidMaterialsConfig.h"

#if defined(VIB_MATERIAL)

namespace Tahoe {

/* forward declarations */
class SpherePointsT;
class C1FunctionT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class AnisoFiber3D: public FSFiberMatT
{
public:

	/* constructor */
	AnisoFiber3D(void);

	/* destructor */
	~AnisoFiber3D(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/*compute output variables*/
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

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
					dMatrixT& FiberMod);
	
	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dSymMatrixT& FiberStretch);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/*constitutive values for matrix*/
	double fMu;
	double fGamma;

	/* potential function */
	C1FunctionT* fPotential;

	/* fibril distribution function */
//	C1FunctionT* fDistribution;

	/* length table */
	/*I4*/
	dArrayT	fI4;

	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	/* jacobian table */
	dArrayT	fjacobian;

	/* STRESS angle tables - by associated stress component */
	int fNumStress;
	dArray2DT fStressTable;
	  	
	/* MODULI angle tables */
	int fNumModuli; 	
	dArray2DT fModuliTable;	

	/* integration point generator */
	SpherePointsT* fSphere;

};

} // namespace Tahoe 
#endif /* _ISO_VIB_2D_H_ */
#endif /*VIB_MATERIAL*/
