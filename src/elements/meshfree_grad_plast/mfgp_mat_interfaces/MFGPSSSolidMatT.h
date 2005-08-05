/* $Id: MFGPSSSolidMatT.h,v 1.3 2005-08-05 07:18:17 kyonten Exp $  */
#ifndef _MFGP_SS_STRUCT_MAT_T_H_
#define _MFGP_SS_STRUCT_MAT_T_H_

/* base class */
#include "MFGPMaterialT.h"

/* direct members */
#include "dSymMatrixT.h"
#include "DetCheckT.h"

namespace Tahoe {

/* forward declarations */
class MFGPMatSupportT;

/** defines the interface for small strain continuum materials */
class MFGPSSSolidMatT: public MFGPMaterialT
{
public:

	/** constructor */
	MFGPSSSolidMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetMFGPMatSupport(const MFGPMatSupportT* support);

	/* required parameter flags */
	virtual bool Need_Strain(void) const { return true; };
	virtual bool Need_Strain_last(void) const { return false; };
	virtual bool Need_Lap_Strain(void) const { return true; };
	virtual bool Need_Lap_Strain_last(void) const { return false; };
	virtual bool Need_Lambda(void) const { return true; };
	virtual bool Need_Lambda_last(void) const { return false; };
	virtual bool Need_Lap_Lambda(void) const { return true; };
	virtual bool Need_Lap_Lambda_last(void) const { return false; };

	/*@{*/
	/** elastic strain */
	const dSymMatrixT& e(void);

	/** elastic strain at the given integration point */
	const dSymMatrixT& e(int ip);

	/** elastic strain */
	const dSymMatrixT& e_last(void);

	/** elastic strain at the given integration point */
	const dSymMatrixT& e_last(int ip);
	
	/** laplacian of elastic strain */
	const dSymMatrixT& lap_e(void);

	/** laplacian of elastic strain at the given integration point */
	const dSymMatrixT& lap_e(int ip);

	/** laplacian of elastic strain */
	const dSymMatrixT& lap_e_last(void);

	/** laplacian of elastic strain at the given integration point */
	const dSymMatrixT& lap_e_last(int ip);
	/*@}*/
	
	/*@{*/
	/** plastic multiplier */
	const dArrayT& pm(void);

	/** plastic multiplier at the given integration point */
	const dArrayT& pm(int ip);

	/** plastic multiplier */
	const dArrayT& pm_last(void);

	/** plastic multiplier at the given integration point */
	const dArrayT& pm_last(int ip);
	
	/** laplacian of plastic multiplier */
	const dArrayT& lap_pm(void);

	/** laplacian plastic multiplier at the given integration point */
	const dArrayT& lap_pm(int ip);

	/** laplacian of plastic multiplier */
	const dArrayT& lap_pm_last(void);

	/** laplacian of plastic multiplier at the given integration point */
	const dArrayT& lap_pm_last(int ip);
	/*@}*/
	

	/* material description */
	virtual const dMatrixT& C_IJKL(void);  // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
		// since infinitesimal strain materials don't make this
		// distinction, these functions just call the respective
		// c_ijkl and s_ij. this means 2 virtual function calls
		// per call. for efficiency, each derived SSSolidMatT should
		// overload these.

	/** return modulus. This default implementation computes the material
	 * modulus using a finite difference approach */
	virtual const dMatrixT& c_ijkl(void);	 
	
	/** spatial elastic modulus */
	virtual const dMatrixT& ce_ijkl(void);

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

	/** return the strain in the material at the current integration point. 
	 * Returns the small strain tensor. */
	virtual void Strain(dSymMatrixT& strain) { strain = e(); };
	
	/** return the laplacian of strain in the material at the current integration point. 
	 * Returns the small strain tensor. */
	virtual void LapStrain(dSymMatrixT& strain) { strain = lap_e(); };
	
	/** return the plastic multiplier in the material at the current integration point.*/
	virtual void PlasticMultiplier(dArrayT& lambda) { lambda = pm(); };
	
	/** return the laplacian of plastic multiplier in the material at the current integration point.*/ 
	virtual void LapPlasticMultiplier(dArrayT& lap_lambda) { lap_lambda = lap_pm(); };

	/*inquire if dissipation variables used in material force calculation are needed*/
	virtual bool HasDissipVar(void) const {return false;};
	
	/** test for localization assuming small strains. check for bifurcation using current
	 * Cauchy stress and the spatial tangent moduli.
	 * \param normal orientation of the localization if localized
	 * \return true if the determinant of the acoustical tensor is negative
	 * or fals if the determinant is positive. */
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
							AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact);
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double &detA);
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);

protected:

	/** small strain material support */
	const MFGPMatSupportT* fMFGPMatSupport;

	/** return value for the modulus */
	dMatrixT fModulus;
};

/* inlines */
/* strain - returns the elastic strain */
inline const dSymMatrixT& MFGPSSSolidMatT::e(void)
{
	return fMFGPMatSupport->LinearStrain();
}

/* elastic strain at the given integration point */
inline const dSymMatrixT& MFGPSSSolidMatT::e(int ip)
{
	return fMFGPMatSupport->LinearStrain(ip);
}

/* strain - returns the elastic strain */
inline const dSymMatrixT& MFGPSSSolidMatT::e_last(void)
{
	return fMFGPMatSupport->LinearStrain_last();
}

/* elastic strain at the given integration point */
inline const dSymMatrixT& MFGPSSSolidMatT::e_last(int ip)
{
	return fMFGPMatSupport->LinearStrain_last(ip);
}

/* strain - returns the laplacian of elastic strain */
inline const dSymMatrixT& MFGPSSSolidMatT::lap_e(void)
{
	return fMFGPMatSupport->LapLinearStrain();
}

/* laplacian of elastic strain at the given integration point */
inline const dSymMatrixT& MFGPSSSolidMatT::lap_e(int ip)
{
	return fMFGPMatSupport->LapLinearStrain(ip);
}

/* strain - returns the laplacian of elastic strain */
inline const dSymMatrixT& MFGPSSSolidMatT::lap_e_last(void)
{
	return fMFGPMatSupport->LapLinearStrain_last();
}

/* laplacian of elastic strain at the given integration point */
inline const dSymMatrixT& MFGPSSSolidMatT::lap_e_last(int ip)
{
	return fMFGPMatSupport->LapLinearStrain_last(ip);
}

/* strain - returns the plastic multiplier */
inline const dArrayT& MFGPSSSolidMatT::pm(void)
{
	return fMFGPMatSupport->LambdaPM();
}

/* plastic multiplier at the given integration point */
inline const dArrayT& MFGPSSSolidMatT::pm(int ip)
{
	return fMFGPMatSupport->LambdaPM(ip);
}

/* strain - returns the plastic multiplier */
inline const dArrayT& MFGPSSSolidMatT::pm_last(void)
{
	return fMFGPMatSupport->LambdaPM_last();
}

/* plastic multiplier at the given integration point */
inline const dArrayT& MFGPSSSolidMatT::pm_last(int ip)
{
	return fMFGPMatSupport->LambdaPM_last(ip);
}

/* strain - returns the laplacian of plastic multiplier */
inline const dArrayT& MFGPSSSolidMatT::lap_pm(void)
{
	return fMFGPMatSupport->LapLambdaPM();
}

/* laplacian of plastic multiplier at the given integration point */
inline const dArrayT& MFGPSSSolidMatT::lap_pm(int ip)
{
	return fMFGPMatSupport->LapLambdaPM(ip);
}

/* strain - returns the laplacian of plastic multiplier */
inline const dArrayT& MFGPSSSolidMatT::lap_pm_last(void)
{
	return fMFGPMatSupport->LapLambdaPM_last();
}

/* laplacian of plastic multiplier at the given integration point */
inline const dArrayT& MFGPSSSolidMatT::lap_pm_last(int ip)
{
	return fMFGPMatSupport->LapLambdaPM_last(ip);
}

} // namespace Tahoe 
#endif /* _MFGP_SS_STRUCT_MAT_T_H_ */
