/* $Id: IsoCorneaModel.h,v 1.1 2005-04-15 22:49:44 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ISO_CORNEA_H_
#define _ISO_CORNEA_H_

/* base classes */
#include "ParameterInterfaceT.h"
#include "OgdenIsotropicT.h"


namespace Tahoe {

/* forward declarations */
class C1FunctionT;

/* forward declarations */
class SpherePointsT;

/** 3D Isotropic VIB using Ogden's spectral formulation */
 class IsoCorneaModel: public OgdenIsotropicT, virtual public ParameterInterfaceT
{
public:

	/* constructor */
	IsoCorneaModel(void);

	/* destructor */
	~IsoCorneaModel(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);

        /** \name implementation of the ParameterInterfaceT interface */
        /*@{*/
        /** information about subordinate parameter lists */
        virtual void DefineSubs(SubListT& sub_list) const;

        /** a pointer to the ParameterInterfaceT of the given subordinate */
        virtual ParameterInterfaceT* NewSub(const StringT& name) const;

        /** accept parameter list */
        virtual void TakeParameterList(const ParameterListT& list);
        /*@}*/

protected:

	/* principal values given principal values of the stretch tensors,
	 * i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod);

	/* allocate memory for all the tables */
	void Dimension(int numbonds);

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:

	/*material parameters*/
	double fR0; /*undeformed RMS distance*/

	/*Bulk Function*/
	double fC0;
	
	/* number of spatial dimensions */
	int fNumSD;

	/* potential function */
	C1FunctionT* fPotential;

	/* length table */
	dArrayT	fLengths;

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
#endif /* _ISO_CORNEA_H_ */
