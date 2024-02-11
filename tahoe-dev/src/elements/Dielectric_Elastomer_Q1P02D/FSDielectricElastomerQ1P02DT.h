
#if !defined(_FSDielectricElastomerQ1P02DT_)
#define _FSDielectricElastomerQ1P02DT_

#include <cassert>

#include "FiniteStrainT.h"
#include "FSDEMatQ1P02DT.h"
#include "dMatrixT.h"

namespace Tahoe {

  /* Forward declarations */
  class FSDEMatSupportQ1P02DT;

  // interface for finite deformation dielectric elastomers 
  // based on 2008 JMPS paper of Suo et al.
  
  class FSDielectricElastomerQ1P02DT: public FiniteStrainT {

  public:

    // constructor
    FSDielectricElastomerQ1P02DT(const ElementSupportT& support);

    // destructor
    virtual ~FSDielectricElastomerQ1P02DT();

    // specify parameters needed by the interface
    virtual void DefineParameters(ParameterListT& list) const;

    // accept parameter list
    virtual void TakeParameterList(const ParameterListT& list);

	/* define total # of DOFs/node, i.e. 4 (3 mech, 1 electric) */
	virtual int TotalNumDOF() const;

    virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
        AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;

    // Electric field at current integration point
    const dArrayT& ElectricField() const;

    // Electric field at given integration point
    const dArrayT& ElectricField(int ip) const;

    // increment current element
    virtual bool NextElement();

    // element stiffness matrix
    virtual void FormStiffness(double constK);

    // internal force
    virtual void FormKd(double constK);

	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/* Q1P0 STUFF */
	/** finalize current step - step is solved */
	virtual void CloseStep(void);
	
	/** restore last converged state */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/** read restart information from stream */
	virtual void ReadRestart(istream& in);

	/** write restart information from stream */
	virtual void WriteRestart(ostream& out) const;

  protected:

    // \param p an existing MaterialSupportT to be initialized. If
    // 0, allocate a new MaterialSupportT and initialize it.
    virtual MaterialSupportT*
    NewMaterialSupport(MaterialSupportT* p = 0) const;

    // Return a pointer to a new material list. Recipient is
    // responsible for freeing the pointer.
    // \param name list identifier
    // \param size length of the list
    virtual MaterialListT* NewMaterialList(const StringT& name, int size);

    // form shape functions and derivatives
    virtual void SetGlobalShape(void);

    // write all current element information to the stream. used to
    // generate debugging information after runtime errors
    virtual void CurrElementInfo(ostream& out) const;

    // Initialize local arrays
    virtual void SetLocalArrays();
	virtual void SetShape();

    // driver for calculating output values
    virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
        const iArrayT& e_codes, dArray2DT& e_values);

	/** accumulate the element mass matrix
	 * \param ip_weight array of weights per integration point or NULL
	 *        if no additional weighting is needed beyond those defined by
	 *        the integration scheme */
	virtual void FormMass(MassTypeT mass_type, double constM, bool axisymmetric,
		const double* ip_weight);

	/** element body force contribution 
	 * \param mass_type mass matrix type of ContinuumElementT::MassTypeT
	 * \param constM pre-factor for the element integral
	 * \param nodal nodal values. Pass NULL for no nodal values: [nen] x [ndof]
	 * \param ip_values integration point source terms. Pass NULL for no integration
	 *        point values : [nip] x [ndof]
	 * \param ip_weight array of weights per integration point or NULL
	 *        if no additional weighting is needed beyond those defined by
	 *        the integration scheme */
	virtual void FormMa(MassTypeT mass_type, double constM, bool axisymmetric,
		const LocalArrayT* nodal_values,
		const dArray2DT* ip_values,
		const double* ip_weight);

  private:

    void Workspace();
	void MassMatrix();

    //void Set_G(const dArray2DT& DNaX, dMatrixT& B_C);

	/** compute mean shape function gradient, H (reference volume), and
	 * current element volume, equation (2.20) */
	void SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const;

	/** special mixed index term in the tangent. Needed to compute
	 * the term in the tangent resulting from
	 * \f$ \nabla \mathbf{u} \textrm{:} \left( \nabla \boldsymbol{\eta} \right)^T \f$.
	 */
	void bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const;

  protected:

    // electric field
    ArrayT<dArrayT> fE_List;
    dArrayT fE_all;

    // The material support used to construct materials lists. This
    // pointer is only set the first time
    // FSDielectricElastomerQ1P02DT::NewMaterialList is called.
    FSDEMatSupportQ1P02DT* fFSDEMatSupportQ1P02D;

	/** current coords with local ordering */
	LocalArrayT fLocCurrCoords;

	/* BELOW IS Q1P0 STUFF COPIED FROM SimoQ1P0.h */
	/** \name element volume */
	/*@{*/
	/** deformed element volume */
	dArrayT fElementVolume;

	/** deformed element volume from the last time step */
	dArrayT fElementVolume_last;
	/*@}*/
	
	/** element pressure. Calculated during SimoQ1P0::FormKd. */
	dArrayT fPressure;

	/** determinant of the deformation gradient for the current element */
	dArrayT fJacobian;

	/** \name work space */
	/*@{*/
	dArray2DT fMeanGradient; /**< mean gradient over element */
	dMatrixT fF_tmp; /**< F workspace */
	dMatrixT fNEEmat; /**< dimension of stiffness matrix */
	dMatrixT fdiff_b;
	dMatrixT fb_bar;
	dMatrixT fb_sig;
	/*@}*/

	dArrayT   Na_0;
	dArray2DT DNa_0;
	dMatrixT  fGrad_U_0;
	dMatrixT  fGrad_UU_0;

	//dArrayT A;

	/** \name work space - from UpdatedLagrangianT.h */
	/*@{*/
	dMatrixT fCauchyStress;	/**< matrix for Cauchy stress tensor: [nsd] x [nsd] */
	dMatrixT fStressStiff;	/**< "compact" stress stiffness contribution: [nen] x [nen] */
	dMatrixT fGradNa;       /**< shape function gradients matrix: [nsd] x [nen] */

	dMatrixT fG_0; /* discrete spatial gradient operator at centroid */
	dMatrixT fQ;
	/*@}*/

  private:

	LocalArrayT fLocScalarPotential;	// electric potential
    FSDEMatQ1P02DT* fCurrMaterial;
   
    // Stiffness storage
    dMatrixT fAmm_mat;	// mechanical material part of Hessian matrix
    dMatrixT fAmm_geo;	// mechanical geometric part of Hessian matrix
    dMatrixT fAmm_neto; // mechanical part according to Neto et. al formulation
    dMatrixT fAme;	// mechanical-electrical coupling part of Hessian matrix
    dMatrixT fAem;	// electrical-mechanical coupling part of Hessian matrix
    dMatrixT fAee;	// electrical-electrical coupling part of Hessian matrix
    dMatrixT fMassMatrix;	// mass matrix for LHS
    
    /* Electric potential */
    const FieldT* fElectricScalarPotentialField;

  };

} // namespace Tahoe

#endif // _FSDielectricElastomerQ1P02DT_
