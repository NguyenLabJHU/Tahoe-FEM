/* $Id: SSSolidMatT.h,v 1.15 2004-09-22 06:12:57 paklein Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _SS_STRUCT_MAT_T_H_
#define _SS_STRUCT_MAT_T_H_

/* base class */
#include "SolidMaterialT.h"

/* direct members */
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** defines the interface for small strain continuum materials */
class SSSolidMatT: public SolidMaterialT
{
public:

	/** constructor */
	SSSolidMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetSSMatSupport(const SSMatSupportT* support);

	/* required parameter flags */
	virtual bool Need_Strain(void) const { return true; };
	virtual bool Need_Strain_last(void) const { return false; };

	/** elastic strain */
	const dSymMatrixT& e(void);

	/** elastic strain at the given integration point */
	const dSymMatrixT& e(int ip);

	/** elastic strain */
	const dSymMatrixT& e_last(void);

	/** elastic strain at the given integration point */
	const dSymMatrixT& e_last(int ip);

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

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

	/** return the strain in the material at the current integration point. 
	 * Returns the small strain tensor. */
	virtual void Strain(dSymMatrixT& strain) { strain = e(); };

	/*inquire if dissipation variables used in material force calculation are needed*/
	virtual bool HasDissipVar(void) const {return false;};

	/** small strain material support additive split of thermal strains by default */
	virtual bool SupportsThermalStrain(void) const { return true; };

protected:

	/* return the acoustical tensor and wave speeds */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);

	/* acoustical tensor routines */
	void Q_2D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const;
	void Q_3D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const;

protected:

	/** small strain material support */
	const SSMatSupportT* fSSMatSupport;

	/** return value for the modulus */
	dMatrixT fModulus;

private:
	
	/* work space */
	dSymMatrixT	fStrainTemp; // elastic strain (w/o thermal)
	dSymMatrixT fQ;          // return value

	/* thermal strain: e_elastic = e_total - e_thermal */
	bool        fHasThermalStrain;
	dSymMatrixT fThermalStrain;
};

} // namespace Tahoe 
#endif /* _SS_STRUCT_MAT_T_H_ */
