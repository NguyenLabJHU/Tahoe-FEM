/* $Id: GradSmallStrainT.h,v 1.4 2004-04-23 18:44:36 rdorgan Exp $ */ 
#ifndef _GRAD_SMALL_STRAIN_T_H_ 
#define _GRAD_SMALL_STRAIN_T_H_ 

/* base classes */
#include "SmallStrainT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSSSolidMatT;
class GradSSMatSupportT;
class ShapeFunctionT;
class C1ShapeTools; 
 
/** element formulation with gradient plasticity constitutive model */
class GradSmallStrainT: public SmallStrainT
{
public:

	/** constructor */
	GradSmallStrainT(const ElementSupportT& support, const FieldT& disp, 
					   const FieldT& field);

	/** destructor */
	~GradSmallStrainT(void);
	
	/** data initialization */
	virtual void Initialize(void); 
	
	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
						   AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/** \name isotropic hardening */
	/*@{*/
	const double& LinearField(void) const { return fField_List[CurrIP()]; };
	const double& LinearField(int ip) const { return fField_List[ip]; };
	/*@}*/
	
	/** \name isotropic hardening from the end of the previous time step */
	/*@{*/
	const double& LinearField_last(void) const { return fField_last_List[CurrIP()]; };
	const double& LinearField_last(int ip) const { return fField_last_List[ip]; };
	/*@}*/
	
	/** return the number of degrees of freedom for field per node */
	int NumDOF_Field(void) const { return fField.NumDOF();} ;
	
	/** number of element integration points for field */
	int NumIP_Field(void) const { return fNumIP_Field;} ;
	
	/** number of nodes per element for the field.
	 * This value will initially be taken to be the number
	 * of nodes per element for the displacement field. */
	int NumElementNodes_Field(void) const { return NumElementNodes();} ;

	/** reference to element shape functions */
	const C1ShapeTools& ShapeFunction(void) const;
	
	/** return the geometry code */
	GeometryT::CodeT GeometryCode_Field(void) const;
	
protected:
	
	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *	a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;
	
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
	
	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;
	
private:
	/** set the shape function matrices using the given shape functions */
	virtual void Set_h(dMatrixT& h) const;
	virtual void Set_p(dMatrixT& p) const;
	virtual void Set_q(dMatrixT& q) const;

protected:
	/** \name return values */
	/*@{*/
	dArrayT fField_List;
	dArrayT fField_last_List;

	dArrayT fYield_List;
	/*@}*/
	  
	/** \name element field in local ordering for current element */
	/*@{*/
	LocalArrayT fLocField;      /**< hardness */
	LocalArrayT fLocLastField;  /**< hardness from last time increment */

	dArrayT fLocFieldTranspose;      /**< hardness */
	dArrayT fLocLastFieldTranspose;  /**< hardness from last time increment */
	/*@}*/
	
private:
	/* \name fields */
	/*@{*/
	const FieldT& fDisplacement; /**< displacement field */
	const FieldT& fField;        /**< hardening parameter fiel */
	/*@}*/
	
	/** \name shape functions for field */
	C1ShapeTools* fShapes_Field;
	
	/** the material support used to construct materials lists. This pointer
	 * is only set the first time GradSmallStrainT::NewMaterialList is called. */
	GradSSMatSupportT* fGradSSMatSupport;
	
	/* run time */
	GradSSSolidMatT*  fCurrMaterial_Grad;
	
	/** \name work space */
	/*@{*/
	/** shape functions for Field */
	dMatrixT fh;      /**<  shape functions */
	dMatrixT fhT;     /**<  shape functions (Transpose) */
	dMatrixT fp;      /**<  gradient shape functions */
	dMatrixT fpT;     /**<  gradient shape functions (Transpose) */
	dMatrixT fq;      /**<  Laplacian shape functions */
	dMatrixT fqT;     /**<  Laplacian shape functions (Transpose) */
	
	/** stiffnesses */
	ElementMatrixT fK_bb;               /**< elastic stiffness matrix */
	ElementMatrixT fK_bh;               /**< off-diagonal matrices */
	ElementMatrixT fK_hb;               /**< off-diagonal matrices */
	ElementMatrixT fK_hh, fK_hp, fK_hq; /**< Gradient matrices */
	ElementMatrixT fK_ct;  /**< constraint matrix */
       
	/** returned matrices */
	dMatrixT fDM_bb;
	dMatrixT fOM_hb, fOM_bh;
	dMatrixT fGM_hh, fGM_hp, fGM_hq;
	dMatrixT fI;
	/*@}*/
	
	/** \name dimensions */
	/*@{*/
	int fNumIP_Disp;            /**< number of integration points */
	int fNumSD;                 /**< number of spatial dimensions */
	int fNumDOF_Disp;           /**< number of degrees of freedom for displacement field */
	int fNumDOF_Field;          /**< number of degrees of freedom for field */
	int fNumEQ_Total;           /**< number of total equations */
	/*@}*/

	/** \name input data for Field */
	/*@{*/
	GeometryT::CodeT fGeometryCode_Field; /**< element parameter */
	int fNumIP_Field;                     /**< number of integration points for field */
	int fNumElementNodes_Field;           /**< number of integration points */
	int fDegreeOfContinuity_Field;        /**< degree of continuity of Field shape functions */
	double fNodalConstraint;              /**< constraint constants */
	/*@}*/
		
	/** \name print debug information */
	/*@{*/
	bool print_GlobalShape;
	bool print_Kd;
	bool print_Stiffness;
	/*@}*/

	/* element degree of continuity types */
	enum TypeT {
		C0 = 0,
		C1 = 1};
};

/* inlines */

/* return the geometry code */
inline GeometryT::CodeT GradSmallStrainT::GeometryCode_Field(void) const
{ return fGeometryCode_Field; }

/* accessors */
inline const C1ShapeTools& GradSmallStrainT::ShapeFunction(void) const
{
#if __option(extended_errorcheck)
	if (!fShapes_Field)
	{
		cout << "\n GradSmallStrainT::ShapeFunction: no shape functions" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif
	return *fShapes_Field;
}

} // namespace Tahoe 

#endif /* _GRAD_SMALL_STRAIN_T_H_ */
