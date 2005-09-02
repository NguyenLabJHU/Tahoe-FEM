/* $Id: ABAQUS_UMAT_SS_BaseT.h,v 1.5 2004-08-01 20:41:53 paklein Exp $ */
#ifndef _ABAQUS_UMAT_SS_BASE_T_H_
#define _ABAQUS_UMAT_SS_BASE_T_H_

/* base classes */
#include "ABAQUS_BaseT.h"
#include "SSIsotropicMatT.h"

/* library support options */
#ifdef __F2C__

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

/* f2c */
#include "f2c.h"

namespace Tahoe {

/** wrapper for using ABAQUS UMAT's with small strain elements */
class ABAQUS_UMAT_SS_BaseT: protected ABAQUS_BaseT, public SSIsotropicMatT 
{
public:

	/** constructor */
	ABAQUS_UMAT_SS_BaseT(void);

	/** required parameter flags */
	virtual bool Need_Strain_last(void) const { return true; };

	/** material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* materials initialization */
	virtual bool NeedsPointInitialization(void) const; // false by default
	virtual void PointInitialize(void);                // per ip

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // per element
	virtual void ResetHistory(void);  // per element

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fPressure; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, i.e., state variables. Returns 0
	 * by default */
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/* set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels) = 0;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* I/O functions */
	virtual void PrintProperties(ostream& out) const;

private:

	/* load element data for the specified integration point */
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/* make call to the UMAT */
	void Call_UMAT(double t, double dt, int step, int iter);
	void Reset_UMAT_Increment(void); // set back to last converged
	void Set_UMAT_Arguments(void);   // compute strains, rotated stresses, etc.
	void Store_UMAT_Modulus(void);   // write modulus to storage

	/* UMAT function wrapper */
	virtual void UMAT(doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, char*,
		integer*, integer*, integer*, integer*,
		doublereal*, integer*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*,
		integer*, integer*, integer*, integer*, integer*,
		integer*, ftnlen) = 0;

protected:
	
	GlobalT::SystemTypeT fTangentType;

	/** properties array */
	nArrayT<doublereal> fProperties;

private:

	/** if true, writes debugging info to output */
	bool fDebug;

	/** if true, uses the modulus computed by the UMAT; otherwise, returns
	 * the finite difference approximation computed by FSSolidMatT */
	bool fUseUMATModulus;

	/** number of "elastic" iterations */
	int fNumElasticIterations;

	/* material name */
	StringT fUMAT_name;
	//other options:
	//  strain type
	//  orientation (*ORIENTATION)
	//  expansion   (*EXPANSION)

	/* work space */
	dSymMatrixT fStress; /**< integration point coordinates */
	dArrayT fIPCoordinates; /**< integration point coordinates */
	double fPressure; /**< pressure for the most recent calculation of the stress */
	
	/* material output data */
	iArrayT fOutputIndex;
	ArrayT<StringT> fOutputLabels;
	
	/* dimensions */
	int fModulusDim; // dimension of modulus storage --- need this???
	int fBlockSize;  // storage block size (per ip)

	/* UMAT dimensions */
	integer ndi;    // number of direct stress components (always 3)
	integer nshr;   // number of engineering shear stress components (2D: 1, 3D: 3)
	integer ntens;  // stress/strain array dimension: ndi + nshr
	integer nstatv; // number of state variables
	
	/* UMAT array arguments */
	nMatrixT<doublereal> fddsdde;
	nArrayT<doublereal>  fdstran;
	nMatrixT<doublereal> fdrot;
	nMatrixT<doublereal> fdfgrd0;
	nMatrixT<doublereal> fdfgrd1;
	nArrayT<doublereal>  fcoords;
	
	/* UMAT stored array arguments */
	nArrayT<doublereal> fstress;
	nArrayT<doublereal> fstrain;
	nArrayT<doublereal> fsse_pd_cd;
	nArrayT<doublereal> fstatv;

	/* stored modulus */
	nArrayT<doublereal> fmodulus;

	/* reset-able history */
	nArrayT<doublereal> fstress_last;
	nArrayT<doublereal> fstrain_last;
	nArrayT<doublereal> fsse_pd_cd_last;
	nArrayT<doublereal> fstatv_last;
	
	/* UMAT argument array storage */
	nArrayT<doublereal> fArgsArray;
	
	/* polar decomposition work space */
	dMatrixT    fmat_nsd;
	dSymMatrixT fsym_mat_nsd;
};

/* inlines */
inline GlobalT::SystemTypeT ABAQUS_UMAT_SS_BaseT::TangentType(void) const
{
	return fTangentType;
}

} /* namespace Tahoe */

#else /* __F2C__ */

#ifndef __MWERKS__
#error "ABAQUS_UMAT_SS_BaseT requires __F2C__"
#endif

#endif /* __F2C__ */


#endif /* _ABAQUS_UMAT_SS_BASE_T_H_ */
