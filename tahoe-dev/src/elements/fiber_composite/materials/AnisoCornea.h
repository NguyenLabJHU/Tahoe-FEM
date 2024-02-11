/* $Id: AnisoCornea.h,v 1.11 2010-06-24 13:49:16 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ANISO_CORNEA_2D_H_
#define _ANISO_CORNEA_2D_H_

/* base classes */
#include "SolidMaterialsConfig.h"
 
/* base class */
#include "FSFiberMatSplitT.h"

#if defined(VIB_MATERIAL)
namespace Tahoe {

/* forward declarations */
class CirclePointsT;
class C1FunctionT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class AnisoCornea: public FSFiberMatSplitT
{
public:

	enum DistributionT {kCircumferential = 0,
							   kMeridional,
							   kMixed,
							   kFile};
	/* constructor */
	AnisoCornea(void);

	/* destructor */
	~AnisoCornea(void);
	
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

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/*@}*/

	/*calculates  matrix contribution to 2PK stress*/
	virtual void ComputeDevMatrixStress(const dSymMatrixT& Cbar,  dSymMatrixT& Stress);

	/*calculates matrix contribution to modulus*/
	virtual void ComputeDevMatrixMod(const dSymMatrixT& Cbar, dSymMatrixT& Stress, dMatrixT& Mod);
	
	/*calculates  matrix contribution to 2PK stress*/
	virtual double ComputeVolMatrixStress(const double I3);

	/*calculates matrix contribution to modulus*/
	virtual double ComputeVolMatrixMod(const double I3);

	/*computes integrated fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);
	
	/*computes integrated moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod);


	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dSymMatrixT& stretch);
	
private:
	/* initialize angle tables */
	void Construct(void);
	
protected:
	
	/*constitutive values for matrix*/
	double fMu;
	double fKappa;
	
	/* integration point generator */
	CirclePointsT*	fCircle;

	/* potential function */
	C1FunctionT* fPotential;

	
	/* length table */
	/*I4 */
	dArrayT	fI4; /*C:M*/
	
	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	int fDType; // flag

	/* for inhomogeneous material */
	ArrayT<dArrayT> fjacobians; // for an inhomogeneous material
	double fk, fphi, fm; // paramters of von mises distribution function
	

	/* STRESS angle tables for fiber stress - by associated stress component */
	dArray2DT fStressTable;
	  	
	/* MODULI angle tables for fiber moduli */
	dArray2DT fModuliTable;	

};

}
#endif 
#endif /*VIB_MATERIAL*/
