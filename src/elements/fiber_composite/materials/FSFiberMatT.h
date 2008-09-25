/* $Id: FSFiberMatT.h,v 1.6 2008-09-25 13:30:13 thao Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FD_FIB_MAT_T_H_
#define _FD_FIB_MAT_T_H_

/* base class */
#include "FSSolidMatT.h"
#include "SpectralDecompT.h"

/* direct members */
#include "FSFiberMatSupportT.h"

namespace Tahoe {

/* forward declarations */
class UpLagFiberCompT;

/** base class for finite deformation fiber composite constitutive models. The interface *
 * provides access to the element-computed fiber orientation vectors in the global (lab) *
 * cartesian coordinates.                                                                */
class FSFiberMatT: public FSSolidMatT
{
public:

	/** constructor */
	FSFiberMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetFSFiberMatSupport(const FSFiberMatSupportT* support);

	/** fiber materials support */
	const FSFiberMatSupportT& FiberMatSupportT(void) const;
	/*@}*/

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

	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

protected:
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;

	/*retrieves fiber rotation matrix*/
	virtual const dMatrixT& GetRotation(void);
	
	/**rotate fiber stress and fiber moduli from local fiber coords to global cartesian coords*/
	void ComputeFiberStretch(const dSymMatrixT& global_stretch, dSymMatrixT& fib_stretch);
	void AssembleFiberStress(const dSymMatrixT& fib_stress, dSymMatrixT& global_stress, 
				const int fillmode = dSymMatrixT::kAccumulate);
	void AssembleFiberModuli(const dMatrixT& fib_mod, dMatrixT& global_mod,
				const int fillmode = dSymMatrixT::kAccumulate);
	/*calculates  matrix contribution to 2PK stress*/
	virtual void ComputeMatrixStress(const dSymMatrixT& C, dSymMatrixT& Stress) = 0;

	/*calculates matrix contribution to modulus*/
	virtual void ComputeMatrixMod(const dSymMatrixT& C, dSymMatrixT& Stress, dMatrixT& Mod) = 0;
	
	/*computes integrated fiber stress in local frame*/
	virtual void ComputeFiberStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress) = 0;
	
	/*computes integrated moduli in local frame*/
	virtual void ComputeFiberMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod) = 0;

protected:

	/** support for finite strain materials */
	const FSFiberMatSupportT* fFSFiberMatSupport;

	/* stretch */
	dSymMatrixT fC;
	
	/* return values */
	dMatrixT    fModulus;
	dSymMatrixT fStress;

	/*dimension*/
	int fNumSD;
	int fNumStress;
	int fNumModuli;
	int fNumFibStress;
	int fNumFibModuli; 	
	
	/*Rotation Tensor*/
	dMatrixT fQ;
	
	/*stretch, stress, and moduli defined in local fiber coords*/
	dSymMatrixT fFiberStretch;
	dSymMatrixT fFiberStress;
	dMatrixT fFiberMod;
};

/* fiber element materials support */
inline const FSFiberMatSupportT& FSFiberMatT::FiberMatSupportT(void) const
{ 
#if __option(extended_errorcheck)
	if (!fFSFiberMatSupport) 
		ExceptionT::GeneralFail("FSFiberMatT::FSFiberMatSupport", "pointer not set");
#endif

	return *fFSFiberMatSupport; 
}

} /* namespace Tahoe */

#endif /* _FD_STRUCT_MAT_T_H_ */
