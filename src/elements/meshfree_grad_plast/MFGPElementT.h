/* $Id: MFGPElementT.h,v 1.2 2005-11-16 23:05:41 kyonten Exp $ */ 
//DEVELOPMENT
#ifndef _MFGP_ELEMENT_T_H_ 
#define _MFGP_ELEMENT_T_H_ 

/* base class */
#include "MFGPAssemblyT.h"

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ElementCardT.h"
#include "D3MeshfreeShapeFunctionT.h"
#include "eIntegratorT.h"
#include "iAutoArrayT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "OutputSetT.h"

/* direct members */
#include "nVariMatrixT.h"

/* exception codes */
#include "ExceptionT.h"

namespace Tahoe 
{
/* forward declarations */
class MFGPMaterialT;
class D3MeshFreeSupportT;
class D3MeshFreeShapeFunctionT;
class MeshFreeFractureSupportT;

/** class for deformation of solids: meshfree small strain gradient plasticity */
class MFGPElementT: public MFGPAssemblyT
{
	
public:
	/** list/index of nodal outputs */
	enum NodalOutputCodeT {
		iNodalCoord = 0, /**< (reference) coordinates */
 	     iNodalDisp = 1, /**< displacements */
 	   iNodalLambda = 2, /**< plastic multipliers */
       iNodalStress = 3, /**< extrapolated stresses */
       iNodalStrain = 4, /**< extrapolated strains */
    iNodalLapStrain = 5, /**< extrapolated laplacian of strains */
    iNodalLapLambda = 6, /**< extrapolated laplacian of plastic multipliers */
      iMaterialData = 7 /**< extrapolated model output */
		};
	
	/** list/index of element outputs */
	enum ElementOutputCodeT {
          iIPLambda = 0, /**< integration point plastic multipliers */
          iIPStress = 1, /**< integration point stresses */
          iIPStrain = 2, /**< integration point strains */
       iIPLapStrain = 3, /**< integration point laplacian of strains */
       iIPLapLambda = 4, /**< integration point laplacian of plastic multipliers */
    iIPMaterialData = 5  /**< integration point material model output */
      	};
      	
	/** constructor */
	MFGPElementT(const ElementSupportT& support);

	/** destructor */
	virtual ~MFGPElementT(void);
	
	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
	
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);
	
	/** contribution to the nodal residual forces. Return the contribution of this element
	 * group to the residual for the given solver group. ParticleT::InternalForce 
	 * returns the internal force calculated with the latest call to ElementBaseT::FormRHS. */
	virtual const dArray2DT& InternalForce(int group);
	
	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	
	/** \name total strain */
	/*@{*/
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;
	/*@}*/
	
	/** total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;
	/*@}*/
	
	/** \name laplacian of total strain */
	/*@{*/
	const dSymMatrixT& LapLinearStrain(void) const;
	const dSymMatrixT& LapLinearStrain(int ip) const;
	/*@}*/
	
	/** laplacian of total strain from the end of the previous time step */
	/*@{*/
	const dSymMatrixT& LapLinearStrain_last(void) const;
	const dSymMatrixT& LapLinearStrain_last(int ip) const;
	/*@}*/
	
	/** \name total lambda */
	/*@{*/
	const dArrayT& Lambda(void) const;
	const dArrayT& Lambda(int ip) const;
	/*@}*/
	
	/** total lambda from the end of the previous time step */
	/*@{*/
	const dArrayT& Lambda_last(void) const;
	const dArrayT& Lambda_last(int ip) const;
	/*@}*/
	
	/** \name laplacian of total lambda */
	/*@{*/
	const dArrayT& LapLambda(void) const;
	const dArrayT& LapLambda(int ip) const;
	/*@}*/
	
	/** laplacian of total lambda from the end of the previous time step */
	/*@{*/
	const dArrayT& LapLambda_last(void) const;
	const dArrayT& LapLambda_last(int ip) const;
	/*@}*/
	
	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations( AutoArrayT<const iArray2DT*>& eq_u,
							AutoArrayT<const RaggedArray2DT<int>*>& eq_lambda);

	/* appends group connectivities to the array */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	
	/** write element output */
	virtual void WriteOutput(void);	
	
	/** returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
	
	/* retrieve nodal unknowns */
	virtual void NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const;
	
	/* weight the computational effort of every node */
	virtual void WeightNodalCost(iArrayT& weight) const;
	
	/* initialize/finalize time increment */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state
	
	/** accessors */
	D3MeshFreeSupportT& D3MeshFreeSupport(void) const;
	
	/* extrapolate the integration point stresses and internal variables, 
   	 * check the yield condition on the nodes of each background grid (element), 
     * and pass the flag whether the nodes are elastically or plastically loaded
	*/
	/*@{*/  
	virtual void CheckNodalYield(void);
	
	/* calculate yield condition at the nodes */ 
	const double& ComputeNodalYield(const dArrayT& qn, double& ff);
	/*@}*/
	
	/*@}*/
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

protected:
	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MFGPMatSupportT to be initialized. If NULL, allocate
	 *        a new MFGPMatSupportT and initialize it. */
	virtual MFGPMatSupportT* NewMFGPMatSupport(MFGPMatSupportT* p = NULL) const;
	
	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MFGPMatListT* NewMFGPMatList(const StringT& name, int size);
	 
	/** set the \e B1 matrix using the given shape function derivatives
	 * \param first derivatives of shape function derivatives: [nsd] x [nen]
	 * \param B1 destination for B1: [nstr]x[nsd*nnd] */
	void Set_B1(const dArray2DT& DNa, dMatrixT& B1);
		
	/** set the \e B3 matrix using the given shape function derivatives
	 * \param third derivatives of shape function derivatives: [nsd*nsd] x [nen]
	 * \param B3 destination for B3: [nstr]x[nsd*nnd] */
	void Set_B3(const dArray2DT& DDDNa, dMatrixT& B3);
                        
    /** set the \e psi_lam matrix using the given shape function
	 * \param shape function: [1] x [nnd]
	 * \param psi_lam destination for psi_lam: [1]x[nnd] */
    /* psi_lam: [1]x[nnd] */
    void Set_PsiLam(const double* Na, dMatrixT& psi_lam);
        
    /** set the \e B4 matrix using the given shape function derivatives
	 * \param second derivatives of shape function derivatives: [nstr] x [nen]
	 * \param B4 destination for B4: [1]x[nnd] */
	void Set_B4(const dArray2DT& DDNa, dMatrixT& B4);
	
	/** calculate the internal force contribution ("-k*d") */
	void FormKd(double constK);

	/** form the element stiffness matrix */
	void FormStiffness(double constK);
	
	/** compute shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/** allocate and initialize shape function objects */
	virtual void SetShape(void);

	/* increment current element */
	virtual bool NextElement(void);

	/* driver for calculating output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

	/** \name construct output labels array */
	/*@{*/
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void GenerateOutputLabels(const iArrayT& n_counts, ArrayT<StringT>& n_labels, 
		const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;
	/*@}*/
	
	/** \name return values */
	/*@{*/
  	ArrayT<dSymMatrixT> fStrain_List;
  	ArrayT<dSymMatrixT> fStrain_last_List;
  	ArrayT<dSymMatrixT> fLapStrain_List;
  	ArrayT<dSymMatrixT> fLapStrain_last_List;
  	ArrayT<dArrayT> fLambda_List;
  	ArrayT<dArrayT> fLambda_last_List;
  	ArrayT<dArrayT> fLapLambda_List;
  	ArrayT<dArrayT> fLapLambda_last_List;
	/*@}*/
	
protected:
	/** \name work space */
	/*@{*/
	dMatrixT fGradU;
	/*@}*/
	
	/** mass type */
	MassTypeT fMassType;
	
	/** \name total force 
	 * Storage for the internal force calculated by this element group. */
	/*@{*/
	bool fStoreInternalForce;
	dArray2DT fForce;
	/*@}*/
	
	/* parameters */
	static const int NumNodalOutputCodes;
	static const int NumElementOutputCodes;
	
	/* flags for stress smoothing */
	bool qNoExtrap;
	
	/** \name work space */
	/*@{*/
	dMatrixT fB1, fB3, fB4, fPsiLam; /**< strain-displacement matrix */
	dMatrixT fCUU1, fCUU2, fCULam1, fCULam2; /**< constitutive matrix */
	dMatrixT fCLamU1, fCLamU2, fCLamLam1, fCLamLam2; /**< constitutive matrix */
	/*@}*/
	
	MFGPMatSupportT* fMFGPMatSupport;
	MFGPMaterialT*  fCurrMaterial;

private:
	/* dynamic wrapper */
	nVariMatrixT<double> fB1_wrap, fB3_wrap, fB4_wrap, fPsiLam_wrap;
	nVariMatrixT<double> fKulambda_wrap, fKulambda_temp_wrap;
	nVariMatrixT<double> fKlambdau_wrap, fKlambdau_temp_wrap;
	
	/** make field at bounding nodes nodally exact */
	bool fAutoBorder;
	
	/* field ready flag */
	bool fFieldSet;

	/** support for meshless calculations */
	MeshFreeFractureSupportT* fMFFractureSupport_displ;
	MeshFreeFractureSupportT* fMFFractureSupport_plast;
	
	/** pointer to list parameters needed to construct meshless shape functions. This
	 * pointer is set during MFGPElementT::TakeParamaterListT and used during
	 * MFGPElementT::SetShape */
	const ParameterListT* fMeshfreeParameters;
	
	/** \name equations */
	/*@{*/
	ArrayT<iArray2DT> fEqnos_displ;
	ArrayT<iArray2DT> fEqnos_plast;
	RaggedArray2DT<int> fMonolithicEquations;
	/*@}*/
};

/* return total strains */
inline const dSymMatrixT& MFGPElementT::LinearStrain(void) const
{
	return fStrain_List[CurrIP()];
}

inline const dSymMatrixT& MFGPElementT::LinearStrain(int ip) const
{
	return fStrain_List[ip];
}

inline const dSymMatrixT& MFGPElementT::LinearStrain_last(void) const
{
	return fStrain_last_List[CurrIP()];
}

inline const dSymMatrixT& MFGPElementT::LinearStrain_last(int ip) const
{
	return fStrain_last_List[ip];
}

/* return laplacian of total strains */
inline const dSymMatrixT& MFGPElementT::LapLinearStrain(void) const
{
	return fLapStrain_List[CurrIP()];
}

inline const dSymMatrixT& MFGPElementT::LapLinearStrain(int ip) const
{
	return fLapStrain_List[ip];
}

inline const dSymMatrixT& MFGPElementT::LapLinearStrain_last(void) const
{
	return fLapStrain_last_List[CurrIP()];
}

inline const dSymMatrixT& MFGPElementT::LapLinearStrain_last(int ip) const
{
	return fLapStrain_last_List[ip];
}

/* return total lambdas */
inline const dArrayT& MFGPElementT::Lambda(void) const
{
	return fLambda_List[CurrIP()];
}

inline const dArrayT& MFGPElementT::Lambda(int ip) const
{
	return fLambda_List[ip];
}

inline const dArrayT& MFGPElementT::Lambda_last(void) const
{
	return fLambda_last_List[CurrIP()];
}

inline const dArrayT& MFGPElementT::Lambda_last(int ip) const
{
	return fLambda_last_List[ip];
}

/* return laplacian of total lambdas */
inline const dArrayT& MFGPElementT::LapLambda(void) const
{
	return fLapLambda_List[CurrIP()];
}

inline const dArrayT& MFGPElementT::LapLambda(int ip) const
{
	return fLapLambda_List[ip];
}

inline const dArrayT& MFGPElementT::LapLambda_last(void) const
{
	return fLapLambda_last_List[CurrIP()];
}

inline const dArrayT& MFGPElementT::LapLambda_last(int ip) const
{
	return fLapLambda_last_List[ip];
}

} // namespace Tahoe 
#endif /* _MFGP_ELEMENT_T_H_ */



