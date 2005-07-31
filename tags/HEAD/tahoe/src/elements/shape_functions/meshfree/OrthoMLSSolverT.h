/* $Id: OrthoMLSSolverT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (07/03/1998)                                          */
/* base class for orthogonal basis, moving least squares, interpolants    */
/* (from Lu et al, Comp Meth App Mech Eng, 126, 1995, 131-153)            */

#ifndef _ORTHO_MLS_SOLVER_T_H_
#define _ORTHO_MLS_SOLVER_T_H_

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "nArrayGroupT.h"
#include "nArray2DGroupT.h"
#include "nVariArray2DT.h"

class OrthoMLSSolverT
{
public:

	/* constructor */
	OrthoMLSSolverT(int nsd, int complete);
	
	/* destructor */
	virtual ~OrthoMLSSolverT(void);
	
	/* class dependent initializations */
	void Initialize(void);
	
	/* set MLS at coords given sampling points and influence of each, returns 1
	 * if successful and 0 if not */
	int SetField(const dArray2DT& nodalcoords, const dArrayT& dmax,
		const dArrayT& samplept);
	
	/* return field value and derivatives - valid AFTER SetField() */
	const dArrayT& phi(void) const;	
	const dArray2DT& Dphi(void) const;	

	/* debugging functions */
	
	/* return field value and derivatives - valid AFTER SetField() */
	const dArrayT& b(void) const;	
	const dArray2DT& Db(void) const;	

	/* orthogonality checking dimensions mat */
	void ComputeOrtho(dMatrixT& mat) const;
	void ComputeDOrtho(int deriv, dMatrixT& mat) const; //deriv ortho check

	/* the weight function */
	const dArrayT& w(void) const;	
	const dArray2DT& Dw(void) const;	

	/* the weight function */
	const dArrayT& q(void) const;	
	const dArray2DT& Dq(void) const;	

	/* return the number of monomial terms */
	int NumberOfMonomials(void) const;

protected:

	/* configure solver for current number of neighbors */
	void Dimension(int numneighbors);

	/* return the number of monomial terms for the given completeness */
	virtual int NumberOfMonomials(int completeness) const = 0;

	/* return monomials evaluated at coords */
	virtual void SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp) = 0;

	/* set weight functions and derivatives - returns the number
	 * of active neighbors */
	int SetWeight(const dArray2DT& localcoords, const ArrayT<double>& dmax);

	/* set basis functions and derivatives */
	void SetBasis(const dArray2DT& localcoords, const dArrayT& samplept);
	
	/* set b and derivatives */
	void Setb(int pterm);
	
	/* set Schmidt orthogonalization */
	void Setalpha(int row, int maxcol);

protected:	

	const int fNumSD;        //spatial dimension
	const int fComplete;     //order of completeness in basis
	
	int fNumNeighbors; //(current) number of neighbors
	
	/* local nodal coordinates (centered at current int pt) */
	dArray2DT fLocCoords;

	/* weight functions and derivatives */
	dArrayT   fw;
	dArray2DT fDw;
	
	/* b (2.11 b) and derivatives */
	dArrayT   fb;
	dArray2DT fDb;
	
	/* (orthogonal) basis functions (2.4) and derivatives */
	dArrayT   fq;
	dArray2DT fDq;

	/* monomials and derivatives */
	dArrayT   fp;
	dArray2DT fDp;
	
	/* nodal p */
	dArray2DT fpJ;
	
	/* nodal q and derivatives */
	dArray2DT         fqJ;
	ArrayT<dArray2DT> fDqJ;

	/* rows of alpha and derivatives */
	dArrayT	  fa;
	dArray2DT fDa;
	
	/* return values {phi,Dphi} of all nodes at int pt */
	dArrayT   fphi;
	dArray2DT fDphi;
	
	/* variable memory managers */
	nArrayGroupT<double>   fArrayGroup;
	nArray2DGroupT<double> fArray2DGroup1;
	nArray2DGroupT<double> fArray2DGroup2;
	nVariArray2DT<double>  fLocCoords_man;
};

/* inlines */

/* return field value and derivatives */
inline const dArrayT& OrthoMLSSolverT::phi(void) const { return fphi; }
inline const dArray2DT& OrthoMLSSolverT::Dphi(void) const { return fDphi; }

inline const dArrayT&   OrthoMLSSolverT::w(void) const { return fw; }
inline const dArray2DT& OrthoMLSSolverT::Dw(void) const { return fDw; }

inline const dArrayT&   OrthoMLSSolverT::b(void) const { return fb; }
inline const dArray2DT& OrthoMLSSolverT::Db(void) const { return fDb; }

inline const dArrayT&   OrthoMLSSolverT::q(void) const { return fq; }
inline const dArray2DT& OrthoMLSSolverT::Dq(void) const { return fDq; }

/* return the number of monomial terms */
inline int OrthoMLSSolverT::NumberOfMonomials(void) const
{
	return NumberOfMonomials(fComplete);
}

#endif /* _ORTHO_MLS_SOLVER_T_H_ */
