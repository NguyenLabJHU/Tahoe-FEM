/* $Id: AnisoCornea.h,v 1.1 2006-08-03 01:10:41 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ANISO_CORNEA_2D_H_
#define _ANISO_CORNEA_2D_H_

/* base classes */
#include "FSFiberMatT.h"

/* direct members */
#include "C1FunctionT.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class AnisoCornea: public FSFiberMatT
{
public:

	/* constructor */
	AnisoCornea(void);

	/* destructor */
	~AnisoCornea(void);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Sum()/3.0; };
	/*@}*/

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

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

protected:
	/* allocate memory for all the tables */
	void Dimension(int numbonds);

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dSymMatrixT& stretch);

private:

	/* initialize angle tables */
	void Construct(void);
	
	void AssembleFiberStress(const dArrayT& fib_stress, dSymMatrixT& global_stress);
	void AssembleFiberModuli(const dSymMatrixT& fib_mod, dMatrixT& global_mod);

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
	dArrayT	fLengths;

	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	/* jacobian table */
	dArrayT	fjacobian;
	dArrayT fangles;

	int fNumSD;
			
	/* STRESS angle tables for fiber stress - by associated stress component */
	int fNumFibStress;
	dArray2DT fStressTable;
	  	
	/*rotation matrix from cartesian to fibril coordinates*/
//	dMatrixT fQ;
//	dArrayT fNT; /*NT-orientation*/
//	dArrayT fIS; /*IS-orientation*/
//	dArrayT fOP;  /*normal orientation*/

	/* MODULI angle tables for fiber moduli */
	int fNumFibModuli; 	
	dArray2DT fModuliTable;	

	/* return values */
	dSymMatrixT fMatStress;
	dMatrixT	fMatMod;
	dMatrixT    fModulus;
	dSymMatrixT fStress;

private:

	/* stretch */
	dSymMatrixT fC;
	dArrayT fFiberStress;
	dSymMatrixT fFiberMod;
};

} // namespace Tahoe 
#endif /* _ISO_VIB_2D_H_ */
