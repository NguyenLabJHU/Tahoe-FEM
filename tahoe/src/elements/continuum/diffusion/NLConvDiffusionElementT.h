/* $Id: NLConvDiffusionElementT.h,v 1.2 2016-11-05 15:46:19 tdnguye Exp $ */
#ifndef _NL_DIFFUSE_CONV_T_H_
#define _NL_DIFFUSE_CONV_T_H_

/* base class */
#include "DiffusionElementT.h"
//#include "FiniteStrainAxiT.h"

namespace Tahoe {

/* forward declarations */
class NLDiffusionMaterialT;
class ScheduleT;

/** diffusion element with nonlinear properties and (mixed) boundary conditions */
class NLConvDiffusionElementT: public DiffusionElementT
{
public:
	
	/** constructor */
	NLConvDiffusionElementT(const ElementSupportT& support);

	/** collecting element group equation numbers. */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** form of tangent matrix, symmetric by default. The tangent for nonlinear
	 * diffusion is generally nonsymmetric */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kNonSymmetric; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:
    // added by RXiao
    virtual bool Axisymmetric(void) const { return true; };
    
    /** allocate and initialize local arrays */
    virtual void SetLocalArrays(void);
    // form shape functions and derivatives
    virtual void SetGlobalShape(void);

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form the residual force vector */
	virtual void RHSDriver(void);

	/** calculate the internal force contribution ("-k*d") */
	virtual void FormKd(double constK);

	/** form the element stiffness matrix */
	virtual void FormStiffness(double constK);

	/** construct a new material support and return a pointer. Recipient is responsible for
	 * for freeing the pointer.
	 * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
	 *        a new MaterialSupportT and initialize it. */
	virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;
    
////    virtual void SetLocalArrays(void);

private:

	/** \name mixed boundary conditions */
	/*@{*/
	void TakeTractionBC(const ParameterListT& list);

	/** compute contribution to RHS from mixed BC's */
	void TractionBC_RHS(void);

	/** compute contribution to LHS from mixed BC's */
	void TractionBC_LHS(void);
	/*@}*/

protected:

	/** field values over the element. The interpolated values are only computed
	 * an integration point at a time and stored. */
	dArrayT fField_list;
    
 ////   LocalArrayT fLocLastDisp;

	/** \name mixed boundary condition parameters 
	 * Nonlinear, mixed boundary conditions of the form
	   \f[
	   		q_i n_i = \epsilon (\Theta - \Theta_0)^{\alpha}
	   \f]
	 */
	/*@{*/
     LocalArrayT* fLocDisplacement;      /**< (optional) nodal temperatures */
 	LocalArrayT* fLocDisplacement_last; /**< (optional) last nodal temperatures */
	/** list of nodes on each face */
	iArray2DT fBCFaces;
	iArray2DT fBCEqnos;
	
	dArrayT feps;
	dArrayT falpha;
    dArrayT fschedulenum;
    
    /** returns true if the material requires the deformation gradient */
 //   bool Needs_F(int material_number) const;
    
    /** returns true if the material requires the deformation gradient
     * from the end of the last time increment */
 //   bool Needs_F_last(int material_number) const;
	/*@}*/
private:
    const ScheduleT* fSchedule;/**< schedule for mixed boundary conditions*/
    
    /** \name radius to integration points computed during FiniteStrainAxiT::SetGlobalShape */
    /*@{*/
    /** current coords with local ordering */
    LocalArrayT fLocCurrCoords;
    /** integration point radius in undeformed configuration */
    dArrayT fRadius_X;
    
    /** integration point radius in current configuration */
    dArrayT fRadius_x;
    /*@}*/



};

} /* namespace Tahoe */

#endif /* _NL_DIFFUSE_CONV_T_H_ */
