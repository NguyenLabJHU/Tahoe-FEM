/* $Id: SimoFiniteStrainT.h,v 1.4 2001-08-20 06:47:16 paklein Exp $ */

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

	/** increment current element */
	virtual bool NextElement(void);	

	/** construct shape function */
	virtual void SetShape(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix. \note This function is very
	 * similar to TotalLagrangianT::FormStiffness, except that the stress
	 * and material tangent are retrieved from fPK1_list and fc_ijkl_list,
	 * respectively instead of being calculated in place. */
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

	/** form the stiffness associated with the enhanced modes
	 * \param K_22 destination for the 2,2 block of the stiffness matrix
	 * \param K_12 destination for the 1,2 (or transposed 2,1) block of 
	 *             the stiffness matrix. Passing NULL skips calculation */
	void FormStiffness_enhanced(dMatrixT& K_22, dMatrixT* K_12);
	
protected:

	/* user-defined parameters */
	bool fIncompressibleMode; /**< flag to include incompressible mode (3D only) */
	int  fLocalIterationMax;  /**< sub-iterations to solve for the element modes */
	double fAbsTol; /**< absolute tolerance on residual of enhanced modes */
	double fRelTol; /**< relative tolerance on residual of enhanced modes */
	
	/* derived parameters */
	int fNumModeShapes; /**< number of mode shapes per element */
	
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

	/** storage for the 1st Piola-Kirchhoff stresses at all integration points
	 * of all elements. These are "loaded" fPK1_list for element calculations
	 * during SimoFiniteStrainT::NextElement. */
	dArray2DT fPK1_storage;

	/** 1st Piola-Kirchhoff stresses at the integration points of the current 
	 * element. These are computed during SimoFiniteStrainT::SetGlobalShape in
	 * SimoFiniteStrainT::FormKd_enhanced while solving for the
	 * enhanced modes */
	ArrayT<dMatrixT> fPK1_list;

	/** storage for the 1st Piola-Kirchhoff stresses at all integration points
	 * of all elements. These are "loaded" fc_ijkl_list for element calculations
	 * during SimoFiniteStrainT::NextElement. */
	dArray2DT fc_ijkl_storage;

	/** material tangent moduli at the integration points of the current 
	 * element. These are computed during SimoFiniteStrainT::SetGlobalShape in
	 * SimoFiniteStrainT::FormStiffness_enhanced while solving for the
	 * enhanced modes */
	ArrayT<dMatrixT> fc_ijkl_list;

	/* workspace */
	dMatrixT fStressMat;      /**< space for a stress tensor */
	dMatrixT fStressStiff_11; /**< compact stress stiffness contribution */
	dMatrixT fStressStiff_12; /**< compact stress stiffness contribution */
	dMatrixT fStressStiff_22; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;         /**< shape function gradients matrix */
	
//	dArrayT   fTemp2;
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x, fDNa_x_enh;
	
	/* work space for enhanced modes */
	dMatrixT fWP_enh;
	dMatrixT fGradNa_enh;
	dArrayT  fRHS_enh;
	dMatrixT fB_enh;
	LAdMatrixT fK22;	
	dMatrixT   fK12;
};

#endif /* _SIMO_FINITE_STRAIN_T_H_ */
