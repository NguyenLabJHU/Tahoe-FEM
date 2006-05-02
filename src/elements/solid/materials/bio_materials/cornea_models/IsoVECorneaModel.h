/* $Id: IsoVECorneaModel.h,v 1.1 2006-05-02 00:58:53 thao Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _IsoVECornea_
#define _IsoVECornea_

/* base classes */
#include "ParameterInterfaceT.h"
#include "RGViscoelasticityT.h"
#include "SpectralDecompT.h"


namespace Tahoe {

/* forward declarations */
class C1FunctionT;
class SpherePointsT;

/** 3D Isotropic VIB using Ogden's spectral formulation */
 class IsoVECorneaModel: public RGViscoelasticityT, virtual public ParameterInterfaceT
{
public:

	/* constructor */
	IsoVECorneaModel(void);

	/* destructor */
	~IsoVECorneaModel(void);
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);

	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

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
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress, 
		const C1FunctionT* potential, const double R0, const double C0);
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod, const C1FunctionT* potential, const double R0, const double C0);

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs, const double R0);

  	/* allocate memory for all the tables */
	void Dimension(int numbonds);

private:

	void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus,
					   const C1FunctionT* potential, const double R0, const double C0,
						const double ietaS, const double ietaB);
	void ComputeiKAB(const dSymMatrixT& eigenmodulus, const double ietaS, const double ietaB);

	/* initialize angle tables */
	void Construct(void);  

   protected:
	/* return values */
	dMatrixT	fModulus;
	dSymMatrixT     fStress;

	/* free energy potential */
	C1FunctionT* fPot_EQ;
    C1FunctionT* fPot_NEQ;

	double fR0_EQ; /*undeformed RMS distance*/
	double fR0_NEQ;
	double fC0_EQ;	/*Bulk Function*/
	double fC0_NEQ;	

	/* length table */
	dArrayT	fLengths;

	/* potential tables */
	dArrayT	fU;
	dArrayT	fdU;
	dArrayT	fddU;

	/* jacobian table */
	dArrayT	fjacobian;

	int fNumSD;
	/* STRESS angle tables - by associated stress component */
	int fNumStress;
	dArray2DT fStressTable;
	  	
	/* MODULI angle tables */
	int fNumModuli; 	
	dArray2DT fModuliTable;	

	/* integration point generator */
	SpherePointsT* fSphere;

	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;
	
	/*work space*/
	dSymMatrixT fb;
	dSymMatrixT fbe;
	dSymMatrixT fiC;
	
	dArrayT     fEigs;
	dArrayT		fEigs_e;
	dArrayT	    fdWdE_EQ;
	dArrayT     fdWdE_NEQ;
	dSymMatrixT fddWddE_EQ;
	dSymMatrixT fddWddE_NEQ;

	dMatrixT fCalg;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
  	/*viscosities*/
	double fietaS;
	double fietaB;
};

} // namespace Tahoe 
#endif /* _IsoVECornea_ */
