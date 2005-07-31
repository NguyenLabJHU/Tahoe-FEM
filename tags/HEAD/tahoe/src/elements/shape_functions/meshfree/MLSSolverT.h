/* $Id: MLSSolverT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */
/* base class for moving least squares, interpolants                      */

#ifndef _MLS_SOLVER_T_H_
#define _MLS_SOLVER_T_H_

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "nArrayGroupT.h"
#include "nArray2DGroupT.h"
#include "nMatrixGroupT.h"
#include "nVariArray2DT.h"

/* forward declarations */
class BasisT;
class C1FunctionT;

class MLSSolverT
{
public:

	/* constructor */
	MLSSolverT(int nsd, int complete);
	
	/* destructor */
	virtual ~MLSSolverT(void);
	
	/* class dependent initializations */
	void Initialize(void);
	
	/* set MLS at fieldpt given sampling points and influence of each, returns 1
	 * if successful and 0 if not */
	int SetField(const dArray2DT& coords, const dArrayT& dmax,
		const dArrayT& volume, const dArrayT& fieldpt, int order);
	
	/* return field value and derivatives - valid AFTER SetField() */
	const dArrayT& phi(void) const;	
	const dArray2DT& Dphi(void) const;	
	const dArray2DT& DDphi(void) const;	

	/* basis dimension */
	int BasisDimension(void) const;

	//TEMP: debugging functions
	
	/* return field value and derivatives - valid AFTER SetField() */
	/* the weight function */
	const dArrayT& w(void) const;	
	const dArray2DT& Dw(void) const;	
	const dArray2DT& DDw(void) const;	

	/* correction function coefficients */
	const dArrayT& b(void) const;	
	const dArrayT& Db(int component) const;	
	const dArrayT& DDb(int component) const;	

	/* correction function */
	const dArrayT& C(void) const;
	const dArray2DT& DC(void) const;	
	const dArray2DT& DDC(void) const;	

private:

	/* configure solver for current number of neighbors */
	void Dimension(void);

	/* set window functions and derivatives - returns the number
	 * of active neighbors */
	int SetWindow(const dArrayT& dmax);

	/* set moment matrix, inverse, and derivatives */
	int SetMomentMartrix(const dArrayT& volume);
	void ComputeM(const dArrayT& volume);
	void ComputeDM(const dArrayT& volume);
	void ComputeDDM(const dArrayT& volume);

	/* set correction function coefficients */
	void SetCorrectionCoefficient(void);

	/* set correction function */
	void SetCorrection(void);

	/* set shape functions */
	void SetShapeFunctions(const dArrayT& volume);

	/* replace M with its inverse */
	int SymmetricInverse3x3(dMatrixT& M);
	int SymmetricInverse4x4(dMatrixT& M);
	
protected:	

	const int fNumSD;    // spatial dimension
	const int fComplete; // order of completeness in basis

	/* runtime parameters */
	int fOrder;        // requested order of derivatives
	int fNumNeighbors; //(current) number of neighbors
	
	/* basis functions */
	BasisT* fBasis;
	
	/* window function */
	C1FunctionT* fWindow;
	
	/* local nodal coordinates (centered at current field pt) */
	dArray2DT fLocCoords;

	/* window function */
	dArrayT   fw;   // [nnd]
	dArray2DT fDw;  // [nsd] x [nnd]
	dArray2DT fDDw; // [nstr] x [nnd]

	/* correction function coefficients */
	dArrayT fb;           // [nbasis]
	ArrayT<dArrayT> fDb;  // [nsd] x [nbasis]
	ArrayT<dArrayT> fDDb; // [nstr] x [nbasis]

	/* inverse of moment matrix */
	dMatrixT fMinv;        // [nbasis] x [nbasis]
	ArrayT<dMatrixT> fDM;  // [nsd] x [nbasis] x [nbasis]
	ArrayT<dMatrixT> fDDM; // [nstr] x [nbasis] x [nbasis]
	
	/* correction function */
	dArrayT   fC;   // [nnd]
	dArray2DT fDC;  // [nsd] x [nnd]
	dArray2DT fDDC; // [nstr] x [nnd]
	
	/* return values of all nodes at field pt */
	dArrayT   fphi;   // [nnd]
	dArray2DT fDphi;  // [nsd] x [nnd]
	dArray2DT fDDphi; // [nstr] x [nnd]
	
	/* variable memory managers */
	nArrayGroupT<double>   fArrayGroup;    // [nnd]
	nArray2DGroupT<double> fArray2DGroup2; // [nsd] x [nnd]
	nArray2DGroupT<double> fArray2DGroup3; // [nstr] x [nnd]	
	nVariArray2DT<double>  fLocCoords_man;

private:

	/* work space */
	dSymMatrixT fNSDsym;
	dMatrixT    fMtemp;
	dArrayT     fbtemp1, fbtemp2, fbtemp3, fbtemp4;
};

/* inlines */

/* return field value and derivatives */
inline const dArrayT& MLSSolverT::phi(void) const { return fphi; }
inline const dArray2DT& MLSSolverT::Dphi(void) const { return fDphi; }
inline const dArray2DT& MLSSolverT::DDphi(void) const { return fDDphi; }

inline const dArrayT&   MLSSolverT::w(void) const { return fw; }
inline const dArray2DT& MLSSolverT::Dw(void) const { return fDw; }
inline const dArray2DT& MLSSolverT::DDw(void) const { return fDDw; }

inline const dArrayT& MLSSolverT::C(void) const { return fC; }	
inline const dArray2DT& MLSSolverT::DC(void) const { return fDC; }
inline const dArray2DT& MLSSolverT::DDC(void) const { return fDDC; }

/* correction function coefficients */
inline const dArrayT& MLSSolverT::b(void) const { return fb; }
inline const dArrayT& MLSSolverT::Db(int component) const
{
	return fDb[component];
}
inline const dArrayT& MLSSolverT::DDb(int component) const
{
	return fDDb[component];
}

#endif /* _MLS_SOLVER_T_H_ */
