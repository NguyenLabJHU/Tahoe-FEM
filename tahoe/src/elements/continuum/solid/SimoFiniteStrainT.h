/* $Id: SimoFiniteStrainT.h,v 1.3 2001-07-20 00:58:01 paklein Exp $ */

#ifndef _SIMO_FINITE_STRAIN_T_H_
#define _SIMO_FINITE_STRAIN_T_H_

/* base class */
#include "FiniteStrainT.h"

/* direct members */
#include "LAdMatrixT.h"

/* forward declarations */
class SimoShapeFunctionT;

/** enhanced strain, finite deformation quad/hex element.
 * formulation due to Simo, Armero, and Taylor, CMAME \b 110, 359-386, 1993 */
class SimoFiniteStrainT: public FiniteStrainT
{
public:

	/** constructor */
	SimoFiniteStrainT(FEManagerT& fe_manager);

	/** destructor */
	~SimoFiniteStrainT(void);

	/** data initialization */
	virtual void Initialize(void);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual void ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;
		
protected:

	/** write element parameter to out */
	virtual void PrintControlData(ostream& out) const;

	/** construct shape function */
	virtual void SetShape(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d"). Uses the
	 * stored values of the stress tensor calculated during solution
	 * of the internal modes in order to spare evaluations of the
	 * stress. */
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

	/** compute enhanced part of F and total F */
	void ComputeEnhancedDeformation(bool need_F, bool need_F_last);

	/** calculate the residual from the internal force */
	void FormKd_enhanced(ArrayT<dMatrixT>& PK1_list, dArrayT& RHS_enh);

	/** form the stiffness associated with the enhanced modes */
	void FormStiffness_enhanced(dMatrixT& K_22);

protected:

	/* user-defined parameters */
	bool fIncompressibleMode; /**< flag to include incompressible mode (3D only) */
	int  fLocalIterationMax;  /**< sub-iterations to solve for the element modes */
	double fAbsTol; /**< absolute tolerance on residual of enhanced modes */
	double fRelTol; /**< relative tolerance on residual of enhanced modes */
	
	/* derived parameters */
	int  fNumModeShapes; /**< number of mode shapes per element */
	
	/* element degrees of freedom */
	dArray2DT   fElementModes;     /**< all element modes */
	LocalArrayT fCurrElementModes; /**< modes for current element */

	/* element degrees of freedom from last time step */
	dArray2DT   fElementModes_last;     /**< all element modes */
	LocalArrayT fCurrElementModes_last; /**< modes for current element */

	/** enhanced shape functions */
	SimoShapeFunctionT* fEnhancedShapes;

  	/* return values */
  	ArrayT<dMatrixT> fF_enh_List;
  	dArrayT          fF_enh_all;
  	ArrayT<dMatrixT> fF_enh_last_List;
  	dArrayT          fF_enh_last_all;

	/* Galerkin part of the deformation gradient */
  	ArrayT<dMatrixT> fF_Galerkin_List;
  	dArrayT          fF_Galerkin_all;
  	ArrayT<dMatrixT> fF_Galerkin_last_List;
  	dArrayT          fF_Galerkin_last_all;

//need:
//(1) element, enhanced DOF's
//    - include incompressible mode for 3D
//(2) shape functions with enhanced gradients
//    - include both standard and higher accuracy integration schemes
//(3) calculate standard and enhanced parts of element force and stiffness
//(4) internal iteration for enhanced modes
//(5) element output of enhanced mode information??

	/* list of IP stresses used my FormKd */
	ArrayT<dMatrixT> fPK1_list;

	/* workspace */
	dMatrixT fStressMat;   /**< space for a stress tensor */
	dMatrixT fStressStiff; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;      /**< shape function gradients matrix */
	
	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x;
	
	/* work space for enhanced modes */
	dMatrixT fWP_enh;
	dMatrixT fGradNa_enh;
	dArrayT  fRHS_enh; 
	LAdMatrixT fK22;
};

#endif /* _SIMO_FINITE_STRAIN_T_H_ */
