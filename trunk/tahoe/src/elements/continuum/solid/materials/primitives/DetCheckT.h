/* $Id: DetCheckT.h,v 1.20 2005-01-25 23:13:22 raregue Exp $ */
/* created: paklein (09/11/1997) */

#ifndef _DETCHECK_T_H_
#define _DETCHECK_T_H_

/* 5 degree increment */
#define sweepIncr 5         // should be an integral divisor of 360
#define numThetaChecks 72    // should be 360/sweepIncr and should be even
#define numPhiChecks  19      // should be 90/sweepIncr+1

/* 1 degree increment
#define sweepIncr 1         // should be an integral divisor of 360
#define numThetaChecks 360    // should be 360/sweepIncr and should be even
#define numPhiChecks  91      // should be 90/sweepIncr+1
*/


#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;
class dMatrixT;
class dMatrixEXT;
class dArrayT;
class dTensor4DT;
//class ContinuumElementT;
class SolidMatSupportT;

/** class to support checks of loss of ellipticity.  \note this class does 
  * not dynamically allocate memory on construction */

class DetCheckT
{
public:

	/** constructor
	 * \param s_jl Cauchy stress
	 * \param c_ijkl spatial tangent modulus
	 * \param ce_ijkl spatial elastic modulus */
	DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl, const dMatrixT& ce_ijkl);

	/** check ellipticity of tangent modulus.
	 * \return 1 if acoustic tensor isn't positive definite,
	 * and returns the normal to the surface of localization.
	 * returns 0, otherwise */
	int IsLocalized(dArrayT& normal);

	/** check ellipticity of tangent modulus using closed form algorithm
	 * taken from R.A.Regueiro's SPINLOC.
	 * \return 1 if acoustic tensor isn't positive definite,
	 * and returns the normal to the surface of localization.
	 * returns 0, otherwise */
	bool IsLocalized_SS(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);
//	bool IsLocalized_SS(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double detA);
	
	/** set pointer to the calling element group */
	//void SetElementGroup(const ContinuumElementT* element) { fElement = *element; };

	/** set pointer to the support of the calling material */
	void SetfStructuralMatSupport(const SolidMatSupportT& support) { fStructuralMatSupport = &support; };

private:

	/** closed-form check for localization, assuming plane strain conditions.
	 * Taken from R.A.Regueiro's SPINLOC.
	 * 1 is 11
	 * 2 is 22
	 * 3 is 12
	 * angle theta subtends from the x1 axis to the band normal */
	bool SPINLOC_localize(const double *c__, double *thetan, bool *loccheck);

	/*3D Small Strain check for localization */
//	bool DetCheck3D_SS(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double detA);
	bool DetCheck3D_SS(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);

	/* auxiliary functions to DetCheck3D_SS */
	void FindApproxLocalMins(double detA [numThetaChecks] [numPhiChecks],
				 int localmin [numThetaChecks] [numPhiChecks], dTensor4DT& C);
	dArrayT ChooseNewNormal(dArrayT& prevnormal, dMatrixEXT& J);
	dArrayT ChooseNormalFromNormalSet(AutoArrayT <dArrayT> &normalSet, dTensor4DT &C);

	/* 2D determinant check function */
	int DetCheck2D(dArrayT& normal);

	/* compute coefficients of det(theta) function */
	void ComputeCoefficients(void);

	/* determinant function and derivatives */
	double det(double t) const;
	double ddet(double t) const;
	double dddet(double t) const;

private:

	const dSymMatrixT& fs_jl;	/* Cauchy stress          */
	const dMatrixT& fc_ijkl;	/* spatial tangent moduli */
	const dMatrixT& fce_ijkl;	/* spatial elastic moduli */

	double phi2, phi4;	/* phase shifts */
	double A0, A2, A4;	/* amplitudes   */

	/** pointer to the support of the calling material */
	const SolidMatSupportT* fStructuralMatSupport;
	
	/** pointer to calling element group */
//	const ContinuumElementT* fElement;
	
	/** flag to indicate first pass */
	static bool fFirstPass;
};

} // namespace Tahoe 
#endif /* _DETCHECK_T_H_ */
