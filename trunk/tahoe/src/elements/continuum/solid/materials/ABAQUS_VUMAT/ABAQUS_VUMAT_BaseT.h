/* $Id: ABAQUS_VUMAT_BaseT.h,v 1.17 2004-07-15 08:26:37 paklein Exp $ */
#ifndef _ABAQUS_VUMAT_BASE_T_H_
#define _ABAQUS_VUMAT_BASE_T_H_

/* base class */
#include "ABAQUS_BaseT.h"
#include "FSSolidMatT.h"
#include "IsotropicT.h"

/* library support options */
#ifdef __F2C__

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

//TEMP - debugging
#include <fstream.h>

/* f2c */
#include "f2c.h"

namespace Tahoe {

/* forward declarations */
class SpectralDecompT;

/** interface for ABAQUS/Explicit VUMAT's. The class is derived
 * from IsotropicT because the VUMAT interface assumes elastic
 * response is approximately isotropic */
class ABAQUS_VUMAT_BaseT: protected ABAQUS_BaseT, public FSSolidMatT, public IsotropicT
{
public:

	/* constructor */
	ABAQUS_VUMAT_BaseT(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	~ABAQUS_VUMAT_BaseT(void);

	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialization */
	virtual void Initialize(void);

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

	/* material description */
	virtual const dMatrixT& C_IJKL(void);  // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

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
	
private:

	//debugging
	ofstream flog;
	
	/* material name */
	StringT fVUMAT_name;
	GlobalT::SystemTypeT fTangentType;
	//other options:
	//  strain type
	//  orientation (*ORIENTATION)
	//  expansion   (*EXPANSION)

	/* work space */
	dMatrixT    fModulus;            // return value
	dSymMatrixT fStress;             // return value
	dArrayT fIPCoordinates;          // integration point coordinates
	double fPressure; /**< pressure for the most recent calculation of the stress */

	/* properties array */
	nArrayT<doublereal> fProperties;
	
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
