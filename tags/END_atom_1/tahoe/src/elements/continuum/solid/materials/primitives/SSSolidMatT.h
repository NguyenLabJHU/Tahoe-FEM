/* $Id: SSSolidMatT.h,v 1.5 2002-11-14 17:06:39 paklein Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _SS_STRUCT_MAT_T_H_
#define _SS_STRUCT_MAT_T_H_

/* base class */
#include "StructuralMaterialT.h"

/* direct members */
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SSMatSupportT;

/** defines the interface for small strain continuum materials */
class SSSolidMatT: public StructuralMaterialT
{
public:

	/** constructor */
	SSSolidMatT(ifstreamT& in, const SSMatSupportT& support);

	/* I/O functions */
	virtual void PrintName(ostream& out) const;

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

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

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
	const SSMatSupportT& fSSMatSupport;

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
