/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.h,v 1.8 2006-07-25 16:30:00 a-kopacz Exp $ */
/* created: a-kopacz (07/04/2006) */
#ifndef _FLUID_ELEMENT_H_
#define _FLUID_ELEMENT_H_

/* base class */
#include "ContinuumElementT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class FluidMaterialT;
class FluidMatSupportT;

/* Fluid element; 4 DOF's per node. Fourth degree of freedom
 * is managed by the FieldT
 */
class FluidElementT: public ContinuumElementT
{
public:

  /** list/index of nodal outputs */
  enum NodalOutputCodeT {
    iNodalCrd = 0, /**< (reference) coordinates */
    iNodalVel = 1, /**< velocitiets */
    iNodalAcc = 2, /**< accelerations */
    iNodalPrs = 3 /**< pressures */
  };

  /** list/index of element outputs */
  enum ElementOutputCodeT {
    iNONE = 0 /**< NONE */
  };

  /** list/index of stabilization parameters */
  enum StabParamCodeT {
    iStabParamOne = 0, /**< \tau_m = \tau_c = \tau_PSPG = \tau_SUPG */
    iStabParamNone = 1 /**< \tau_m = \tau_c = \tau_PSPG = \tau_SUPG = 0 */
  };

  /** constructor */
  FluidElementT(const ElementSupportT& support);

  /** destructor */
  virtual ~FluidElementT(void);

  /** \name access to nodal values */
  /*@{*/
  const LocalArrayT& OldVelocities(void) const;
  const LocalArrayT& Velocities(void) const;
  const LocalArrayT& Accelerations(void) const;
  const LocalArrayT& Pressures(void) const;
  /*@}*/

  /** compute nodal force */
  virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

  /** returns the stored energy */
  virtual double InternalEnergy(void);

  /** compute specified output parameter and send for smoothing */
  virtual void SendOutput(int kincode);

protected:

  /** parameters */
  static const int NumNodalOutputCodes;
  static const int NumElementOutputCodes;
  static const int NumStabParamCodes;
  static const int kPressureNDOF;

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** allocate and initialize shape function objects */
	virtual void SetShape(void);

  /** set the \e B matrix at the specified integration point */
  void B(int ip, dMatrixT& B_matrix) const;
  
	/** set the \e B matrix using the given shape function derivatives
	 * Set strain displacement matrix as in Hughes (2.8.20)
	 * \param derivatives of shape function derivatives: [nsd] x [nen]
	 * \param B destination for B */
	void Set_B(const dArray2DT& derivatives, dMatrixT& B) const;

	/** set the \e B matrix for 2D axysymmetric problems using the given shape functions
	 * and derivative using the y-axis as the axis or revolution.
	 * \param shapes shape function values: [nen]
	 * \param derivatives of shape function derivatives: [nsd] x [nen]
	 * \param r distance from the axis of revolution
	 * \param B destination for B */
	void Set_B_axi(const dArrayT& shapes, const dArray2DT& derivatives, double r, dMatrixT& B) const;

	/** increment current element */
	virtual bool NextElement(void);

  /** form shape functions and derivatives */
  virtual void SetGlobalShape(void);

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
    
	/** form the element stiffness matrix
	 * Compute the linearization of the force calculated by SolidElementT::FormKd */
	virtual void FormStiffness(double constK);

	/** internal force */
	virtual void FormKd(double constK);

  /** drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
  virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
  virtual void RHSDriver(void);

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;

	/** return a pointer to a new material list. Recipient is responsible for freeing
	 * the pointer.
	 * \param name list identifier
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

  /** driver for calculating output values */
  virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
    const iArrayT& e_codes, dArray2DT& e_values);

  /**************************************************************/
  /******* implementing interface of ParameterInterfaceT ********/
  /**************************************************************/
  /* virtual void DefineParameters(ParameterListT& list) const;
   * virtual void DefineSubs(SubListT& sub_list) const;
   * virtual ParameterInterfaceT* NewSub(const StringT& name) const;
   * virtual virtual void TakeParameterList(const ParameterListT& list);
  /**************************************************************/
  /**************************************************************/

  /** \name implementation of the ParameterInterfaceT interface */
  /*@{*/
  /** describe the parameters needed by the interface */
  virtual void DefineParameters(ParameterListT& list) const;

  /** information about subordinate parameter lists */
  virtual void DefineSubs(SubListT& sub_list) const;

  /** a pointer to the ParameterInterfaceT of the given subordinate */
  virtual ParameterInterfaceT* NewSub(const StringT& name) const;

  /** accept parameter list */
  virtual void TakeParameterList(const ParameterListT& list);
  /*@}*/

  /** extract the list of material parameters */
  virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

  /** run time */
  FluidMaterialT* fCurrMaterial;

  /*nodal dofs with local ordering.  Includes both velocities and pressures*/
  /*Sets  pressures as the last dof in the array*/
  LocalArrayT fLocDisp;
  LocalArrayT fLocLastDisp;
  LocalArrayT fLocVel;

  /** nodal pressure values with local ordering, shallow copy of fLocDisp */
  LocalArrayT fLocCurPrs;

  /** nodal current/old velocities with local ordering shallow copy of fLocDisp and fLocLastDisp*/
  LocalArrayT fLocCurVel;
  LocalArrayT fLocOldVel;

  /** nodal current accelerations with local ordering shallow copy of fLocVel */
  LocalArrayT fLocCurAcc; // post-processing only

  /** \name work space */
  /*@{*/
  dMatrixT fD; /**< constitutive matrix          */
  dMatrixT fB; /**< "strain-displacement" matrix */
  double tau_m;
  double tau_c;
  /*@}*/

 /** field gradients over the element. The gradients are only computed
	* an integration point at a time and stored.
  * Should we symmetrize? */

  /** pressure */
  dArrayT fPres_list;
  
  /** pressure gradient */
  ArrayT<dArrayT> fGradPres_list;

  /** velocity */
  ArrayT<dArrayT> fVel_list;
  ArrayT<dArrayT> fOldVel_list;

  /** velocity gradient */
  ArrayT<dMatrixT> fGradVel_list;

private:

  /** \name construct output labels array */
  /*@{*/
  virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
    iArrayT& counts) const;
  virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
    iArrayT& counts) const;
  virtual void GenerateOutputLabels(const iArrayT& n_counts,
    ArrayT<StringT>& n_labels, const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;
  /*@}*/

	/** Material Interface/Support */
	FluidMatSupportT* fFluidMatSupport;

	/** pressure index */
	int fPresIndex;

  /** stabilization parameter */
  FluidElementT::StabParamCodeT fStabParam;

  /** FOR DEBUGGING PURPOSES ONLY */
  void WriteCallLocation( char* loc ) const;
};

inline const LocalArrayT& FluidElementT::OldVelocities(void) const { return fLocOldVel; }
inline const LocalArrayT& FluidElementT::Velocities(void) const { return fLocCurVel; }
inline const LocalArrayT& FluidElementT::Accelerations(void) const { return fLocCurAcc; }
inline const LocalArrayT& FluidElementT::Pressures(void) const { return fLocCurPrs; }

} // namespace Tahoe
#endif /* _FLUID_ELEMENT_H_ */