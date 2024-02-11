
#if !defined(_FSThermoMechT_)
#define _FSThermoMechT_

#include <cassert>

#include "UpdatedLagrangianT.h"
#include "FSThermoMechSupportT.h"
#include "dMatrixT.h"

namespace Tahoe {

  /* Forward declarations */
  class FSThermoMechMatT;
  class FSThermoMechT: public UpdatedLagrangianT {

  public:

    /* constructor*/
    FSThermoMechT(const ElementSupportT& support);

    /* destructor*/
    ~FSThermoMechT();

      /* information about subordinate parameter lists */
    virtual void DefineSubs(SubListT& sub_list) const;

      /* a pointer to the ParameterInterfaceT of the given subordinate */
      virtual ParameterInterfaceT* NewSub(const StringT& name) const;

      /* specify parameters needed by the interface */
      virtual void DefineParameters(ParameterListT& list) const;

      /* accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);
      
      /*TDN: initialize step.  Te from last step is extrapolated from ip to nodes*/
      virtual void InitStep(void);

      /** extract the list of material parameters */
     virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

      /* define total # of DOFs/node, i.e. 4 (3 mech, 1 electric) */
	virtual int TotalNumDOF() const;

    virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
        AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;


    // increment current element
    virtual bool NextElement();

	// Temperature Gradient  at current integration point
	//	FiniteStrainT::ComputeGradient() computes the deformation gradient from the element displacements.
	const dArrayT& TemperatureGradient() const;
	  
	// Temperature Gradient  at given integration point
	const dArrayT& TemperatureGradient(int ip) const;
	  
	/** returns nodal temperaturerates.*/
//	const LocalArrayT TemperatureRate(void) const { return *fLocTemperatureRate; };
	  	
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
      
      /** construct the effective mass matrix */
      virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
      
      /** form the residual force vector */
      virtual void RHSDriver(void);

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

    // driver for calculating output values
 /*   virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
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
      /** \name mixed boundary conditions */
      /*@{*/
      void TakeTractionBC(const ParameterListT& list);
      
      /** compute contribution to RHS from mixed BC's */
      void TractionBC_RHS(void);
      
      /** compute contribution to LHS from mixed BC's */
      void TractionBC_LHS(void);
      /*@}*/

      /** TDN: functions to extrapolate element ip effective temperature to nodes*/
      void Extrapolate(void);
      /*utility function to assemble lumped mass values*/
      void AssembleArray(const dArrayT& elem_val, dArrayT& global_val, const iArrayT& elem_nodes);
      /*utility function to assemble ip values*/
      void AssembleArray2D(const dArray2DT& elem_val, dArray2DT& global_val, const iArrayT& elem_nodes);
  protected:

    // temperature gradient
    ArrayT<dArrayT> fTempGrad_List;
	dArrayT         fTempGrad_all;       /**< grouped memory for all temperature gradients */

    // The material support used to construct materials lists. This
    // pointer is only set the first time
    // FSThermoMechT::NewMaterialList is called.
    FSThermoMechSupportT* fFSThermoMechSupport;
	  
  private:
	  
	  /** time integrator */
	  const eIntegratorT* fTempIntegrator;
	 //  eIntegratorT* fTempIntegrator;
	  //temperature field
	  const FieldT* fTemperatureField;	

	  LocalArrayT* fLocTemperatureRate;	// temperature time derivative at current time
	  
      FSThermoMechMatT* fCurrMaterial;
   
    // Stiffness storage
    dMatrixT fAmm_mat;	// mechanical material part of Hessian matrix
    dMatrixT fAmm_geo;	// mechanical geometric part of Hessian matrix
    dMatrixT fAme;	// mechanical-temperature coupling part of Hessian matrix
    dMatrixT fAem;	// temperature-mechanical coupling part of Hessian matrix
    dMatrixT fAee;	// temperature-temperature coupling part of Hessian matrix
    dMatrixT fMassMatrix;	// thermal mass matrix for LHS
	  
	dMatrixT fgrad;  //workspace for calculating temperature gradient.
      dMatrixT fdij;
      dMatrixT fHij;

      /*TDN: Define workspace to extrapolate effective temperature from ips to nodes*/
      /*For a forward euler integration scheme, we only need to extrapolate the previous time step.  Forward Euler is the only way I can think of to solve the diffusion equation for Te locally at the IP*/
      int fNumTe;
      dArrayT fGlobalMass;
      dArray2DT fGlobalVal_last;
      dArrayT felem_mass; 
      dArray2DT felem_val_last;
      //      dArray2DT felem_val;
      //      dArray2DT fGlobalVal;

      dMatrixT fGradEffectiveTemp_last;
      LocalArrayT fLocTe_last;
      ArrayT<dArray2DT> fGradTe_last_List; /**< last deformation gradient */
      //      dMatrixT fGradEffectiveTemp;
      //      LocalArrayT fLocTe;
      //      ArrayT<dArray2DT> fGradTe_List;      /**< deformation gradient */
      
      
/*Boundary conditions*/
      iArray2DT fBCFaces;
      iArray2DT fBCEqnos;
      
      double feps;
      double fT0;
      double falpha;
  };

} // namespace Tahoe

#endif // _FSThermoMechT_
