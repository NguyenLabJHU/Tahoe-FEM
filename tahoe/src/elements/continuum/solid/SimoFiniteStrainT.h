/* $Id: SimoFiniteStrainT.h,v 1.1 2001-07-11 01:02:16 paklein Exp $ */

#ifndef _SIMO_FINITE_STRAIN_T_H_
#define _SIMO_FINITE_STRAIN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "dMatrixT.h"

/* forward declarations */
class SimoShapeFunctionT;

/** enhanced strain, finite deformation quad/hex element.
 * formulation due to Simo, Armero, and Taylor, CMAME \b 110, 359-386, 1993 */
class SimoFiniteStrainT: public FiniteStrainT
{
public:

	/** constructor */
	SimoFiniteStrainT(FEManagerT& fe_manager);

	/** data initialization */
	virtual void Initialize(void);
		
protected:

	/** write element parameter to out */
	virtual void PrintControlData(ostream& out) const;

	/** construct shape function */
	virtual void SetShape(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

private:

	/** compute modified, enhanced deformation gradient including the
	 * additional incompressible mode which appears in 3D only.
     * \note this function takes the place of the base class implementation
     * of SetGlobalShape because the enhancement to the deformation gradient 
     * isn't simply additive. A modified differential operator is applied to 
     * the standard part of the deformation gradient as well, see (2.14) and 
     * Section 3.4 */
	void ModifiedEnhancedDeformation(void);

protected:

	/* user-defined parameters */
	bool fIncompressibleMode; /**< flag to include incompressible mode (3D only) */
	int  fLocalIterationMax;  /**< sub-iterations to solve for the element modes */
	
	/* derived parameters */
	int  fNumModeShapes; /**< number of mode shapes per element */
	
	/* element degrees of freedom */
	dArray2DT fElementModes;     /**< all element modes */
	dArray2DT fCurrElementModes; /**< modes for current element */

	/** enhanced shape functions */
	SimoShapeFunctionT* fEnhancedShapes;

//need:
//(1) element, enhanced DOF's
//    - include incompressible mode for 3D
//(2) shape functions with enhanced gradients
//    - include both standard and higher accuracy integration schemes
//(3) calculate standard and enhanced parts of element force and stiffness
//(4) internal iteration for enhanced modes
//(5) element output of enhanced mode information??

	/* workspace */
	dMatrixT fStressMat;   /**< space for a stress tensor */
	dMatrixT fStressStiff; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;      /**< shape function gradients matrix */
	
	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x;
};

#endif /* _SIMO_FINITE_STRAIN_T_H_ */
