/* $Id: MLSSolverT.h,v 1.2 2001-06-19 23:22:03 paklein Exp $ */
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
#include "WindowT.h"
#include "MeshFreeT.h"

/* forward declarations */
class BasisT;

/** class to calculate MLS shape functions and derivatives */
class MLSSolverT
{
public:

	/** constructor.
	 * \param nsd number of spatial dimensions
	 * \param complete order of completeness for the basis functions
	 * \param window_type window function specifier
	 * \param window_params array of window function parameters */
	MLSSolverT(int nsd, int complete, MeshFreeT::WindowTypeT window_type, 
		const dArrayT& window_params);
	
	/** destructor */
	virtual ~MLSSolverT(void);
	
	/** write parameters */
	virtual void WriteParameters(ostream& out) const;
	
	/** class dependent initializations */
	void Initialize(void);
	
	/** computevshape function and derivatives. 
	 * \param coords coordinates of the neighborhood nodes: [nnd] x [nsd]
	 * \param nodal_params support parameters for each node: [nnd] x [nparam] 
	 * \param volume array of nodal volumes 
	 * \param fieldpt point of field evaluation
	 * \order highest order field derivative to compute
	 * \return 1 if successful, or 0 otherwise */
	int SetField(const dArray2DT& coords, const dArray2DT& nodal_param,
		const dArrayT& volume, const dArrayT& fieldpt, int order);

	/* return field value and derivatives from previous SetField */

	/** shape function values.
	 * \return array of nodal shape functions: [nnd] */
	const dArrayT& phi(void) const;
	
	/** shape function derivatives.
	 * \return array of shape functions derivatives: [nsd] x [nnd] */
	const dArray2DT& Dphi(void) const;	

	/** shape function second derivatives.
	 * \return array of shape functions second derivatives: [nstr] x [nnd] */
	const dArray2DT& DDphi(void) const;	
		
	/** neighbor search type needed by the window function */
	WindowT::SearchTypeT SearchType(void) const;

	/** coverage test */
	bool Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const;

	/** basis dimension.
	 * \return the number of basis functions */
	int BasisDimension(void) const;
	
	/** number of nodal field parameters */
	int NumberOfSupportParameters(void) const;
	
	/** "synchronization" of nodal field parameters */
	void SynchronizeSupportParameters(dArray2DT& params_1, dArray2DT& params_2);

	/** modify nodal shape function parameters */
	void ModifySupportParameters(dArray2DT& nodal_params) const;

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
	int SetWindow(const dArray2DT& support_params);

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
	MeshFreeT::WindowTypeT fWindowType;
	WindowT* fWindow;
	
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

/* number of nodal field parameters */
inline int MLSSolverT::NumberOfSupportParameters(void) const
{
#if __option(extended_errorcheck)
	if (!fWindow) throw eGeneralFail;
#endif
	return fWindow->NumberOfSupportParameters();
}

/* coverage test */
inline bool MLSSolverT::Covers(const dArrayT& x_n, const dArrayT& x, 
	const dArrayT& param_n) const
{
#if __option(extended_errorcheck)
	if (!fWindow) throw eGeneralFail;
#endif
	return fWindow->Covers(x_n, x, param_n);
}

/* neighbor search type needed by the window function */
inline WindowT::SearchTypeT MLSSolverT::SearchType(void) const
{
#if __option(extended_errorcheck)
	if (!fWindow) throw eGeneralFail;
#endif
	return fWindow->SearchType();
}

/* modify nodal shape function parameters */
inline void MLSSolverT::ModifySupportParameters(dArray2DT& nodal_params) const
{
	fWindow->ModifySupportParameters(nodal_params); 
};

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
