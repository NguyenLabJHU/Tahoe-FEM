/* $Id: DetCheckT.h,v 1.2 2001-07-03 01:35:41 paklein Exp $ */
/* created: paklein (09/11/1997)                                          */
/* Note: Class does not dynamically allocate memory on construction       */

#ifndef _DETCHECK_T_H_
#define _DETCHECK_T_H_

/* forward declarations */
class dSymMatrixT;
class dMatrixT;
class dArrayT;

class DetCheckT
{
public:

	/* Constructor */
	DetCheckT(const dSymMatrixT& s_jl, const dMatrixT& c_ijkl);

	/*
	 * Returns 1 if acoustic tensor isn't positive definite,
	 * and returns the normal to the surface of localization.
	 * Returns 0, otherwise.
	 */
	int IsLocalized(dArrayT& normal);

private:

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
