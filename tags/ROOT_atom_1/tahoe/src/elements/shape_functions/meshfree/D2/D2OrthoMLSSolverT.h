/* $Id: D2OrthoMLSSolverT.h,v 1.4 2002-07-02 19:57:02 cjkimme Exp $ */
/* created: paklein (10/10/1999)                                          */
/* base class for orthogonal basis, moving least squares, interpolants    */
/* (from Lu et al, Comp Meth App Mech Eng, 126, 1995, 131-153)            */

#ifndef _D2_ORTHO_MLS_SOLVER_T_H_
#define _D2_ORTHO_MLS_SOLVER_T_H_

/* base class */
#include "OrthoMLSSolverT.h"

/* direct members */
#include "dSymMatrixT.h"


namespace Tahoe {

class D2OrthoMLSSolverT: public OrthoMLSSolverT
{
public:

	/* constructor */
	D2OrthoMLSSolverT(int nsd, int complete);
	
	/* class dependent initializations */
	void Initialize(void);
	
	/* set MLS at coords given sampling points and influence of each, returns 1
	 * if successful and 0 if not */
	int SetField(const dArray2DT& nodalcoords, const ArrayT<double>& dmax,
		const dArrayT& samplept);
	
	/* return field value and derivatives - valid AFTER SetField() */
	const dArray2DT& DDphi(void) const;	

	/* return field value and derivatives - valid AFTER SetField() */
	const dArray2DT& DDb(void) const;	

	/* the weight function */
	const dArray2DT& DDw(void) const;	

	/* the weight function */
	const dArray2DT& DDq(void) const;	

private:

	/* return monomials evaluated at coords */
	virtual void _SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp,
		dArray2DT& DDp) = 0;

	/* set weight functions and derivatives - returns the number
	 * of active neighbors */
	int SetWeight(const dArray2DT& localcoords, const ArrayT<double>& dmax);

	/* set basis functions and derivatives */
	void SetBasis(const dArray2DT& localcoords, const dArrayT& samplept);
	
	/* set b and derivatives */
	void Setb(int pterm);
	
	/* set Schmidt orthogonalization */
	void Setalpha(int row, int maxcol);

	/* configure solver for current number of neighbors */
	void Dimension(int numneighbors);

protected:	

	/* weight function derivatives */
	dArray2DT fDDw;
	
	/* b (2.11 b) derivatives */
	dArray2DT fDDb;
	
	/* (orthogonal) basis function derivatives */
	dArray2DT fDDq;

	/* monomial derivatives */
	dArray2DT fDDp;
	
	/* nodal q derivatives */
	ArrayT<dArray2DT> fDDqJ;

	/* rows of alpha derivatives */
	dArray2DT fDDa;
	
	/* return values {DDphi} of all nodes at int pt */
	dArray2DT fDDphi;
	
	/* variable memory managers */
	nArray2DGroupT<double> fArray2DGroup3;
	
private:

	/* work space */
	dSymMatrixT fNSDsym;  	
};

/* inlines */

/* return field value and derivatives */
inline const dArray2DT& D2OrthoMLSSolverT::DDphi(void) const { return fDDphi; }

inline const dArray2DT& D2OrthoMLSSolverT::DDw(void) const { return fDDw; }
inline const dArray2DT& D2OrthoMLSSolverT::DDb(void) const { return fDDb; }
inline const dArray2DT& D2OrthoMLSSolverT::DDq(void) const { return fDDq; }

} // namespace Tahoe 
#endif /* _D2_ORTHO_MLS_SOLVER_T_H_ */
