/* $Id: MLSSolverT.h,v 1.4 2004-09-13 23:33:44 raregue Exp $ */
/* created: paklein (12/08/1999) */
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

namespace Tahoe {

/* forward declarations */
class BasisT;

/** class to calculate MLS shape functions and derivatives. Before
 * MLSSolverT::SetField can used to calculate MLS functions, MLSSolverT::Initialize
 * must be called to initialize some internal work space */
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
	
	/** compute shape function and derivatives. 
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
	
	/** shape function third derivatives.
	 * \return array of shape functions third derivatives: [nstr] x [nnd] */
	const dArray2DT& DDDphi(void) const;	//added  [nstr] should change?? KY
		
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

	/** \name translate support parameters to support sizes */
	/*@{*/
	/** compute spherical support size */
	virtual double SphericalSupportSize(const dArrayT& param_n) const;

	/** compute spherical support size in batch */
	virtual void SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const;

	/** compute rectangular support size in batch */
	virtual void RectangularSupportSize(const dArray2DT& param_n, dArray2DT& support_size) const;
	/*@}*/

	/* debugging functions */
	
	/* return field value and derivatives - valid AFTER SetField() */
	/* the weight function */
	const dArrayT& w(void) const;	
	const dArray2DT& Dw(void) const;	
	const dArray2DT& DDw(void) const;
	const dArray2DT& DDDw(void) const;	// kyonten

	/* correction function coefficients */
	const dArrayT& b(void) const;	
	const dArrayT& Db(int component) const;	
	const dArrayT& DDb(int component) const;
	const dArrayT& DDDb(int component) const; // kyonten

	/* correction function */
	const dArrayT& C(void) const;
	const dArray2DT& DC(void) const;	
	const dArray2DT& DDC(void) const;
	const dArray2DT& DDDC(void) const; // kyonten

private:

	/** error checking accessor to the window function */
	WindowT& Window(void) const;
	
	/* configure solver for current number of neighbors */
	void Dimension(void);

	/* set moment matrix, inverse, and derivatives */
	int SetMomentMatrix(const dArrayT& volume);
	void ComputeM(const dArrayT& volume);
	void ComputeDM(const dArrayT& volume);
	void ComputeDDM(const dArrayT& volume);
	void ComputeDDDM(const dArrayT& volume); // kyonten

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

	/* parameters */
	const int fNumSD;    /**< number of spatial dimensions   */
	const int fComplete; /**< order of completeness in basis */

	/* runtime parameters */
	int fOrder;        /**< number of derivatives to calculate at the current field point */
	int fNumNeighbors; /**< number of neighbors at the current field point */
	
	/** basis functions */
	BasisT* fBasis;
	
	/** \name window function */
	/*@{*/
	MeshFreeT::WindowTypeT fWindowType; /**< window function type   */
	WindowT* fWindow;                   /**< window function object */
	dArrayT fOrigin;
	/*@}*/
	
	/** local nodal coordinates centered at current field point */
	dArray2DT fLocCoords;

	/* window function */
	dArrayT   fw;   /**< values of window function at the current field point: [nnd] */
	dArray2DT fDw;  /**< values of window function gradient at the current field point: [nsd] x [nnd] */
	dArray2DT fDDw; /**< second gradient of window functions at the current field point: [nstr] x [nnd] */
	dArray2DT fDDDw; /**< third gradient of window functions at the current field point: [nstr] x [nnd] */
						//[nstr] should change for third derivative?? KY
	
	/* correction function coefficients */
	dArrayT fb;           /**< correction function coefficients at the current field point: [nbasis] */
	ArrayT<dArrayT> fDb;  /**< gradient of correction function coefficients: [nsd] x [nbasis] */
	ArrayT<dArrayT> fDDb; /**< second gradient of correction function coefficient: [nstr] x [nbasis] */
	ArrayT<dArrayT> fDDDb; /**< third gradient of correction function coefficient: [nstr] x [nbasis] */
							//[nstr] should change for third derivative?? KY
	
	/* inverse of moment matrix */
	dMatrixT fMinv;        /**< moment matrix at the current field point: [nbasis] x [nbasis] */
	ArrayT<dMatrixT> fDM;  /**< gradient of moment matrix: [nsd] x [nbasis] x [nbasis] */
	ArrayT<dMatrixT> fDDM; /**< second gradient of moment matrix: [nstr] x [nbasis] x [nbasis] */
	ArrayT<dMatrixT> fDDDM; /**< third gradient of moment matrix: [nstr] x [nbasis] x [nbasis] */
								//[nstr] should change for third derivative?? KY
	/* correction function */
	dArrayT   fC;   /**< correction function at the current field point: [nnd] */
	dArray2DT fDC;  /**< gradient of the correction function: [nsd] x [nnd] */
	dArray2DT fDDC; /**< second gradient of the correction function: [nstr] x [nnd] */
	dArray2DT fDDDC; /**< third gradient of the correction function: [nstr] x [nnd] */
						//[nstr] should change for third derivative?? KY
	
	/* return values of all nodes at field pt */
	dArrayT   fphi;   /**< nodal shape functions at the current field point: [nnd] */
	dArray2DT fDphi;  /**< nodal shape function gradients at the current field point: [nsd] x [nnd] */
	dArray2DT fDDphi; /**< second gradient of nodal shape functions: [nstr] x [nnd] */
	dArray2DT fDDDphi; /**< third gradient of nodal shape functions: [nstr] x [nnd] */
						//[nstr] should change for third derivative?? KY
	
	/* variable memory managers */
	nArrayGroupT<double>   fArrayGroup;    /**< variable memory manager for arrays length [nnd] */
	nArray2DGroupT<double> fArray2DGroup2; /**< variable memory manager for 2D arrays length [nsd] x [nnd] */
	nArray2DGroupT<double> fArray2DGroup3; /**< variable memory manager for 2D arrays length [nstr] x [nnd]	*/
	nArray2DGroupT<double> fArray2DGroup4; /**< variable memory manager for 2D arrays length [nstr] x [nnd]	*/
	nVariArray2DT<double>  fLocCoords_man; /**< variable memory manager for local coordinates array */

private:

	/* work space */
	dSymMatrixT fNSDsym;
	dMatrixT    fMtemp;
	dArrayT     fbtemp1, fbtemp2, fbtemp3, fbtemp4;
	dArrayT     fbtemp5, fbtemp6, fbtemp7, fbtemp8, fbtemp9, fbtemp10;// kyonten (DDDb)
};

/* inlines */

/* error checking accessor to the window function */
inline WindowT& MLSSolverT::Window(void) const {
#if __option(extended_errorcheck)
	if (!fWindow) ExceptionT::GeneralFail("MLSSolverT::Window", "window not set");
#endif
	return *fWindow;
}

/* number of nodal field parameters */
inline int MLSSolverT::NumberOfSupportParameters(void) const
{
	return Window().NumberOfSupportParameters();
}

/* coverage test */
inline bool MLSSolverT::Covers(const dArrayT& x_n, const dArrayT& x, 
	const dArrayT& param_n) const
{
	return Window().Covers(x_n, x, param_n);
}

/* neighbor search type needed by the window function */
inline WindowT::SearchTypeT MLSSolverT::SearchType(void) const
{
	return Window().SearchType();
}

/* compute spherical support size */
inline double MLSSolverT::SphericalSupportSize(const dArrayT& param_n) const {
	return Window().SphericalSupportSize(param_n);
}

/* compute spherical support size in batch */
inline void MLSSolverT::SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const {
	Window().SphericalSupportSize(param_n, support_size);
}

/* compute rectangular support size in batch */
inline void MLSSolverT::RectangularSupportSize(const dArray2DT& param_n, dArray2DT& support_size) const {
	Window().RectangularSupportSize(param_n, support_size);
}

/* return field value and derivatives */
inline const dArrayT&   MLSSolverT::phi(void) const { return fphi; }
inline const dArray2DT& MLSSolverT::Dphi(void) const { return fDphi; }
inline const dArray2DT& MLSSolverT::DDphi(void) const { return fDDphi; }
inline const dArray2DT& MLSSolverT::DDDphi(void) const { return fDDDphi; }// kyonten

inline const dArrayT&   MLSSolverT::w(void) const { return fw; }
inline const dArray2DT& MLSSolverT::Dw(void) const { return fDw; }
inline const dArray2DT& MLSSolverT::DDw(void) const { return fDDw; }
inline const dArray2DT& MLSSolverT::DDDw(void) const { return fDDDw; }// kyonten

inline const dArrayT&   MLSSolverT::C(void) const { return fC; }	
inline const dArray2DT& MLSSolverT::DC(void) const { return fDC; }
inline const dArray2DT& MLSSolverT::DDC(void) const { return fDDC; }
inline const dArray2DT& MLSSolverT::DDDC(void) const { return fDDDC; }// kyonten

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
inline const dArrayT& MLSSolverT::DDDb(int component) const // kyonten
{
	return fDDDb[component];
}

} // namespace Tahoe 
#endif /* _MLS_SOLVER_T_H_ */
