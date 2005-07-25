/* $Id: SimoFiniteStrainT.h,v 1.16.30.1 2005-07-25 02:36:55 paklein Exp $ */
#ifndef _SIMO_FINITE_STRAIN_T_H_
#define _SIMO_FINITE_STRAIN_T_H_

/* base classes */
#include "FiniteStrainT.h"
#include "DOFElementT.h"

/* direct members */
#include "LAdMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SimoShapeFunctionT;

/** enhanced strain, finite deformation quad/hex element.
 * formulation due to Simo, Armero, and Taylor, CMAME \b 110, 359-386, 1993. 
 * The enhanced element modes can be solved by three methods, either as part
 * of the global system of equations, using a static condensation, or the 
 * local iteration procedure described in the paper above. Solution of
 * the enhanced modes as part of the global system of equations is supported
 * by inheritance of the DOFElementT interface. */
class SimoFiniteStrainT: public FiniteStrainT, public DOFElementT
{
public:

	/** solution method for the enhanced modes */
	enum SolutionMethodT {
	            kMonolithic = 0, /**< solve internal mode as part of global equations */
		kStaticCondensation = 1, /**< solve internal mode by static condensation */
            kLocalIteration = 2  /**< solve internal mode with a staggered, local iteration */
		};                

	/** constructor */
	SimoFiniteStrainT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~SimoFiniteStrainT(void);

	/** data initialization */
	virtual void Initialize(void);

	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(ifstreamT& in);

	/** write restart information from stream */
	virtual void WriteRestart(ofstreamT& out) const;

	/** return field connectivities. Returns connectivities including
	 * the tags for the enhanced element modes when using the monolithic
	 * solver method. */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/** append element equations numbers to the list. Returns equations 
	 * including the tag equation numbers for the enhanced element modes 
	 * when using the monolithic solver method. */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** determine number of tags needed. See DOFElementT. */
	virtual void SetDOFTags(void);
	
	/** return the array tag numbers in the specified set currently 
	 * used/need by the group. See DOFElementT. */
	virtual iArrayT& DOFTags(int tag_set);

	/** generate nodal connectivities. See DOFElementT. */
	virtual void GenerateElementData(void);

	/** return the connectivities associated with the element generated
	 * degrees of freedom. See DOFElementT. */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/** restore the element degrees of freedom. See DOFElementT. */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/** check element group for tag reconfiguration. See DOFElementT. 
	 * This function returns 0. No reconfiguration is needed unless
	 * the number of elements is changing. */
	virtual int Reconfigure(void) { return 0; };		

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;

protected:

	/** increment current element */
	virtual bool NextElement(void);	

	/** construct shape function */
	virtual void SetShape(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** write all current element information to the stream. used to generate
	 * debugging information after runtime errors */
	virtual void CurrElementInfo(ostream& out) const;

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

	/** form the contribution to the stiffness associated with the Galerkin
	 * part of the deformation gradient */
	void FormStiffness_staggered(double constK);

	/** form the stiffness associated with the enhanced modes
	 * \param K_22 destination for the 2,2 block of the stiffness matrix
	 * \param K_12 destination for the 1,2 (or transposed 2,1) block of 
	 *        the stiffness matrix. Passing NULL skips calculation */
	void FormStiffness_enhanced(dMatrixT& K_22, dMatrixT* K_12);

	/** compute and assemble the element stiffness for the monolithic
	 * solution scheme */
	void FormStiffness_monolithic(double constK);

	/** form the contribution to the the residual force associated with the 
	 * Galerkin part of the deformation gradient */
	void FormKd_staggered(double constK);

	/** compute and assemble the residual force for the monolithic
	 * solution scheme */
	void FormKd_monolithic(double constK);

	/** calculate the residual from the internal force */
	void FormKd_enhanced(ArrayT<dMatrixT>& PK1_list, dArrayT& RHS_enh);
	
protected:

	/* user-defined parameters */
	bool fIncompressibleMode; /**< flag to include incompressible mode (3D only) */
	SolutionMethodT fModeSolveMethod; /**< approach to solving for the enhanced modes */
	
	/* parameters needed for the kLocalIteration mode solver method */
	int fLocalIterationMax;  /**< sub-iterations to solve for the element modes */
	double fAbsTol; /**< absolute tolerance on residual of enhanced modes */
	double fRelTol; /**< relative tolerance on residual of enhanced modes */
	bool fModesConverged; /**< true if local iteration satisfies convergence tolerances */
	
	/* derived parameters */
	int fNumModeShapes; /**< number of mode shapes per element */
	
	/** tags for the element modes. This array is used only if the
	 * monolithic solution method is used. There is one tag per
	 * element. */
	iArrayT fEnhancedModeTags;

	/** connectivities to link the tags for the enhanced modes to the
	 * node number of the elements. */
	iArray2DT fEnhancedConnectivities;
	
	/** all element modes stored in \a local \a ordering */	
	dArray2DT fElementModes;

	/** modes for current element */
	LocalArrayT fCurrElementModes; 

	/* element degrees of freedom from last time step */
	dArray2DT   fElementModes_last;     /**< all element modes stored in \a local \a ordering */
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
	dMatrixT fStressStiff_21; /**< compact stress stiffness contribution */
	dMatrixT fStressStiff_22; /**< compact stress stiffness contribution */
	dMatrixT fGradNa;         /**< shape function gradients matrix */
	
	dMatrixT  fTempMat1, fTempMat2;
	dArray2DT fDNa_x, fDNa_x_enh;
	
	/* work space for enhanced modes */
	dMatrixT fWP_enh;
	dMatrixT fGradNa_enh;
	dArrayT  fRHS_enh;
	dMatrixT fB_enh;
	LAdMatrixT fK22;	
	dMatrixT fK12, fK11, fK21;
};

} // namespace Tahoe 
#endif /* _SIMO_FINITE_STRAIN_T_H_ */
