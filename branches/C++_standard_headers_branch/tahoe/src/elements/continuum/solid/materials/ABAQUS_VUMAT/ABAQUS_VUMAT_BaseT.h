/* $Id: ABAQUS_VUMAT_BaseT.h,v 1.19 2004-08-01 20:41:58 paklein Exp $ */
#ifndef _ABAQUS_VUMAT_BASE_T_H_
#define _ABAQUS_VUMAT_BASE_T_H_

/* base class */
#include "ABAQUS_BaseT.h"
#include "FSIsotropicMatT.h"

/* library support options */
#ifdef __F2C__

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

/* f2c */
#include "f2c.h"

namespace Tahoe {

/* forward declarations */
class SpectralDecompT;

/** interface for ABAQUS/Explicit VUMAT's. The class is derived
 * from IsotropicT because the VUMAT interface assumes elastic
 * response is approximately isotropic */
class ABAQUS_VUMAT_BaseT: protected ABAQUS_BaseT, public FSIsotropicMatT
{
public:

	/* constructor */
	ABAQUS_VUMAT_BaseT(void);

	/* destructor */
	~ABAQUS_VUMAT_BaseT(void);

	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

	/** material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* materials initialization */
	virtual bool NeedsPointInitialization(void) const; // false by default
	virtual void PointInitialize(void);                // per ip

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // per element
	virtual void ResetHistory(void);  // per element

	/** \name spatial description */
	/*@{*/
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
	void Load(const ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/* make call to the VUMAT */
	void Call_VUMAT(double t, double dt, int step, int iter);
	void Reset_VUMAT_Increment(void); // set back to last converged
	void Set_VUMAT_Arguments(void);   // compute strains, rotated stresses, etc.

	/* VUMAT function wrapper */
	virtual void VUMAT(integer*, integer*, integer*, integer*, integer*, integer*, integer*, doublereal*,
		doublereal*, doublereal*, char*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
		doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*, doublereal*, doublereal*, doublereal*, doublereal*, doublereal*,
                doublereal*) = 0;
protected:

	GlobalT::SystemTypeT fTangentType;

	/** properties array */
	nArrayT<doublereal> fProperties;

private:

	/** if true, writes debugging info to output */
	bool fDebug;

	/* material name */
	StringT fVUMAT_name;
	//other options:
	//  strain type
	//  orientation (*ORIENTATION)
	//  expansion   (*EXPANSION)

	/* work space */
	dArrayT fIPCoordinates; /**< integration point coordinates */
	double fPressure; /**< pressure for the most recent calculation of the stress */

	/* material output data */
	iArrayT fOutputIndex;
	ArrayT<StringT> fOutputLabels;
	
	/* dimensions */
	int fModulusDim; // dimension of modulus storage --- need this???
	int fBlockSize;  // storage block size (per ip)

	/* VUMAT dimensions */
	integer ndi;    // number of direct stress components (always 3)
	integer nshr;   // number of engineering shear stress components (2D: 1, 3D: 3)
	integer ntens;  // stress/strain array dimension: ndi + nshr
	integer nstatv; // number of state variables
	
	/* VUMAT array arguments */
	//nMatrixT<doublereal> fddsdde;
	nArrayT<doublereal>  fdstran;
	nMatrixT<doublereal> fdrot;
	nMatrixT<doublereal> fdfgrd0;
	nMatrixT<doublereal> fdfgrd1;
	nArrayT<doublereal>  fcoords;
	nMatrixT<doublereal> fRelSpin;
	nArrayT<doublereal>  fUOld2;
	nArrayT<doublereal>  fUNew2;

	/* VUMAT stored array arguments */
	nArrayT<doublereal> fstress;
	nArrayT<doublereal> fstrain;
	//nArrayT<doublereal> fsse_pd_cd;
	nArrayT<doublereal> fstatv;

	/* stored modulus */
	//nArrayT<doublereal> fmodulus;

	/* reset-able history */
	nArrayT<doublereal> fstress_last;
	nArrayT<doublereal> fstrain_last;
	//nArrayT<doublereal> fsse_pd_cd_last;
	nArrayT<doublereal> fstatv_last;
	
	/* VUMAT argument array storage */
	nArrayT<doublereal> fArgsArray;
	
	/* polar decomposition work space */
	SpectralDecompT* fDecomp;
	dMatrixT fF_rel;
	dMatrixT fA_nsd, fROld, fRNew;
	dSymMatrixT fU1, fU2, fU1U2, fUOld, fUNew;
};

/* inlines */
inline GlobalT::SystemTypeT ABAQUS_VUMAT_BaseT::TangentType(void) const
{
	return fTangentType;
}

} /* namespace Tahoe */

#else /* __F2C__ */

#ifndef __MWERKS__
#error "ABAQUS_VUMAT_BaseT requires __F2C__"
#endif

#endif /* __F2C__ */

#endif /* _ABAQUS_VUMAT_BASE_T_H_ */
