/* $Id: DetCheckT.h,v 1.6 2002-02-26 01:47:43 raregue Exp $ */
/* created: paklein (09/11/1997) */

#ifndef _DETCHECK_T_H_
#define _DETCHECK_T_H_

/* forward declarations */
class dSymMatrixT;
class dMatrixT;
class dArrayT;
class dTensor4DT;

/** class to support checks of loss of ellipticity.  \note this class does 
 * not dynamically allocate memory on construction */
class DetCheckT
{
public:

	/** constructor
	 * \param s_jl Cauchy stress
	 * \param c_ijkl spatial tangent modulus */
	DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl);

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
	int IsLocalized_SS(dArrayT& normal);

private:

	/** closed-form check for localization, assuming plane strain conditions.
	 * Taken from R.A.Regueiro's SPINLOC.
	 * 1 is 11
	 * 2 is 22
	 * 3 is 12
	 * angle theta subtends from the x1 axis to the band normal */
	int SPINLOC_localize(double *c__, double *thetan, int *loccheck);

	/*3D Small Strain check for localization */
	int DetCheck3D_SS(dArrayT& normal);

	/* 2D determinant check function */
	int DetCheck2D(dArrayT& normal);

	/* compute coefficients of det(theta) function */
	void ComputeCoefficients(void);

	/* determinant function and derivatives */
	double det(double t) const;
	double ddet(double t) const;
	double dddet(double t) const;

private:

	const dSymMatrixT& 	fs_jl;	/* Cauchy stress          */
	const dMatrixT&		fc_ijkl;/* spatial tangent moduli */

	double phi2, phi4;	/* phase shifts */
	double A0, A2, A4;	/* amplitudes   */
};

#endif /* _DETCHECK_T_H_ */
