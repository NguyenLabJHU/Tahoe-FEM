/* $Id: GradSmallStrainT.h,v 1.13 2004-08-04 22:02:13 rdorgan Exp $ */ 
#ifndef _GRAD_SMALL_STRAIN_T_H_ 
#define _GRAD_SMALL_STRAIN_T_H_ 

/* base classes */
#include "SmallStrainT.h"

/* direct members */
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSSSolidMatT;
class GradSSMatSupportT;
class ShapeFunctionT;
class ShapeTools; 
 
/** element formulation with gradient plasticity constitutive model */
class GradSmallStrainT: public SmallStrainT
{
public:

	/** constructor */
	GradSmallStrainT(const ElementSupportT& support);

	/** destructor */
	~GradSmallStrainT(void);
	
	/** \name field */
	/*@{*/
	const double& LinearPMultiplier(void) const;
	const double& LinearPMultiplier(int ip) const;
	const double& LinearPMultiplier_last(void) const;
	const double& LinearPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name gradient field */
	/*@{*/
	const double& LinearGradPMultiplier(void) const;
	const double& LinearGradPMultiplier(int ip) const;
	const double& LinearGradPMultiplier_last(void) const;
	const double& LinearGradPMultiplier_last(int ip) const;
	/*@}*/
	
	/** \name Laplacian field */
	/*@{*/
	const double& LinearLapPMultiplier(void) const;
	const double& LinearLapPMultiplier(int ip) const;
	const double& LinearLapPMultiplier_last(void) const;
	const double& LinearLapPMultiplier_last(int ip) const;
	/*@}*/

	/** return the number of degrees of freedom for lambda per node */
	int NumDOF_PMultiplier(void) const { return fPMultiplier->NumDOF();} ;
	
	/** number of element integration points for lambda */
	int NumIP_PMultiplier(void) const { return fNumIP_PMultiplier;} ;
	
	/** number of nodes per element for the lambda.
	 * This value will initially be taken to be the number
	 * of nodes per element for the displacement field. */
	int NumElementNodes_PMultiplier(void) const { return fNumElementNodes_PMultiplier;} ;

	/** reference to element shape functions */
	const ShapeTools& ShapeFunction(void) const;
	
	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
						   AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** return the description of the given inline subordinate parameter list. */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** extract the list of material parameters */
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;

protected:
	
	/* element degree of continuity types */
	enum TypeT {kC0 = 0,	kC1 = 1};

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *	a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;
	
	/** return a pointer to a new material list. Recipient is responsible for freeing 
	 * the pointer. 
	 * \param name list identifier
	 * \param size length of the list */
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);

	/** define the elements blocks for the element group */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);

	/** initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
	
	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);
	
	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);
	
	/** increment current element */
	virtual bool NextElement(void);	
	
	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
private:
	/** set the shape function matrices using the given shape functions */
	virtual void Set_h(dMatrixT& h) const;
	virtual void Set_p(dMatrixT& p) const;
	virtual void Set_q(dMatrixT& q) const;

protected:
	/** \name return values */
	/*@{*/
	dArrayT fPMultiplier_List;
	dArrayT fPMultiplier_last_List;

	dArrayT fGradPMultiplier_List;
	dArrayT fGradPMultiplier_last_List;

	dArrayT fLapPMultiplier_List;
	dArrayT fLapPMultiplier_last_List;

	dArrayT fYield_List;
	/*@}*/
	  
	/** \name element lambda in local ordering for current element */
	/*@{*/
	LocalArrayT fLocPMultiplier;      /**< hardness: for 1d arranged as { r1, r2; r1x, r2x } */
	LocalArrayT fLocLastPMultiplier;  /**< hardness from last time increment */
  	dArrayT fLocPMultiplierTranspose; /**< hardness: for 1d arranged as { r1, r1x; r2, r2x } */
	/*@}*/
	
  	/** the material support used to construct materials lists. This pointer
  	 * is only set the first time SmallStrainT::NewMaterialList is called. */
	GradSSMatSupportT* fGradSSMatSupport;

	//	/* run time */
	GradSSSolidMatT*  fCurrMaterial_Grad;
	
 private:

	/** connectivities for the multiplier */
	ArrayT<iArray2DT> fConnectivities_PMultiplier;
	iArray2DT fConnectivities_All;

	/* \name fields */
	/*@{*/
	const FieldT* fDisplacement; /**< displacement field */
	const FieldT* fPMultiplier;        /**< hardening parameter field */
	ArrayT<KBC_ControllerT*> fFixedPMultiplier; /**< fixed conditions block-by-block */
	/*@}*/
	
	/** \name shape functions for lambda */
	ShapeTools* fShapes_PMultiplier;

	/** \name work space */
	/*@{*/
	/** shape functions for PMultiplier */
	dMatrixT fh, fhT; /**<  shape functions */
	dMatrixT fp;      /**<  gradient of shape functions */
	dMatrixT fq;      /**<  Laplacian of shape functions */
	
	/** stiffnesses */
	ElementMatrixT fK_bb;               /**< elastic stiffness matrix */
	ElementMatrixT fK_bh;               /**< off-diagonal matrices */
	ElementMatrixT fK_hb;               /**< off-diagonal matrices */
	ElementMatrixT fK_hh, fK_hp, fK_hq; /**< Gradient dependent matrices */
	ElementMatrixT fK_ct;               /**< plastic multiplier constraint matrix */

	/** returned matrices obtained from material model */
	dMatrixT fDM_bb;                    /**< elastic stiffness modulus */
	dMatrixT fODM_hb, fODM_bh;            /**< off-diagonal moduli */
	dMatrixT fGM_hh, fGM_hp, fGM_hq;    /**< gradient dependent moduli */
	dMatrixT fI;
	/*@}*/
	
	/** array of nodes for the hardening field */
	iArrayT fNodesPMultiplier;
	/*@}*/
	
	/** \name dimensions */
	/*@{*/
	int fNumSD;                 /**< number of spatial dimensions */
	int fNumIP_Disp;            /**< number of integration points for displacement field*/
	int fNumElementNodes_Disp;  /**< number of nodes per element for displacement field */
	int fNumDOF_Disp;           /**< number of degrees of freedom for displacement field */
	int fNumDOF_PMultiplier;          /**< number of degrees of freedom for lambda */
	int fNumEQ_Total;           /**< number of total equations */
	/*@}*/

	/** \name input data for PMultiplier */
	/*@{*/
	int fNumIP_PMultiplier;                     /**< number of integration points for lambda */
	int fNumElementNodes_PMultiplier;           /**< number of nodes per element for lambda */
	int fDegreeOfContinuity_PMultiplier;        /**< degree of continuity of PMultiplier shape functions */
	double fNodalConstraint;              /**< constraint constants */
	/*@}*/
		
	/** \name print debug information */
	/*@{*/
	bool fprint_GlobalShape;
	bool fprint_Kd;
	bool fprint_KdMatrix;	
	bool fprint_Stiffness;
	bool fprint_StiffnessMatrix;
	bool fprint_All;
	/*@}*/
};

/* inlines */

inline const double& GradSmallStrainT::LinearPMultiplier(void) const
{
	return fPMultiplier_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearPMultiplier(int ip) const
{
	return fPMultiplier_List[ip]; 
}

inline const double& GradSmallStrainT::LinearPMultiplier_last(void) const
{
	return fPMultiplier_last_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearPMultiplier_last(int ip) const
{
	return fPMultiplier_last_List[ip]; 
}

inline const double& GradSmallStrainT::LinearGradPMultiplier(void) const
{
	return fGradPMultiplier_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearGradPMultiplier(int ip) const
{
	return fGradPMultiplier_List[ip]; 
}

inline const double& GradSmallStrainT::LinearGradPMultiplier_last(void) const
{
	return fGradPMultiplier_last_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearGradPMultiplier_last(int ip) const
{
	return fGradPMultiplier_last_List[ip]; 
}

inline const double& GradSmallStrainT::LinearLapPMultiplier(void) const
{
	return fLapPMultiplier_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearLapPMultiplier(int ip) const
{
	return fLapPMultiplier_List[ip]; 
}

inline const double& GradSmallStrainT::LinearLapPMultiplier_last(void) const
{
	return fLapPMultiplier_last_List[CurrIP()]; 
}

inline const double& GradSmallStrainT::LinearLapPMultiplier_last(int ip) const
{
	return fLapPMultiplier_last_List[ip]; 
}

/* accessors */
inline const ShapeTools& GradSmallStrainT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes_PMultiplier)
	{
		cout << "\n GradSmallStrainT::ShapeFunction: no shape functions" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif
	return *fShapes_PMultiplier;
}

} // namespace Tahoe 

#endif /* _GRAD_SMALL_STRAIN_T_H_ */
