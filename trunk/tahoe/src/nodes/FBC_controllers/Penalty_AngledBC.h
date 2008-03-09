/*  $Id:Penalty_AngledBC.h v 1.0 
 * 
 *
 *  Created by vicky on 3/5/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
 #ifndef _ANGLED_KBCT_H_
 #define _ANGLED_KBCT_H_
 
 /*base class*/
 #include "FBC_ControllerT.h"
 #include "Traction_CardT.h"
 #include "ScheduleT.h"
 #include "ElementMatrixT.h"

 #include "dArray2DT.h"
/*Angled kinematic boundary conditions 
  The class implements kinematic boundary conditions normal to an angled surface using the penalty method.
  Given d is the vector of nodal displacements of dimensions neq x 1, define the Lagrangian  as,
  L(d) = 1/2 d^T K d - F^T d + k^2 g(d)^T g(d).
  The vector g(d) has dimensions nkbc x 1 where nkbc is the number of nodes with angled KBC,
  and g(d)_A = sum_i^nsd d_iA ni  + u0, where n is the normal vector defined at a node.
  The resulting matrix equation is (K+Kbar) d = F + Fbar */
  

namespace Tahoe {

class  DomainIntegrationT;

class Penalty_AngledBC: public FBC_ControllerT
{
	public:
	
	Penalty_AngledBC();
  
	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const {return(GlobalT::kSymmetric);};

	/** set to initial conditions */
	virtual void InitialCondition(void);

	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);
	virtual void ApplyRHS(void);

	/* initialize/finalize step */
	virtual void InitStep(void){};
	virtual void CloseStep(void){};

	/* reset to the last known solution */
	virtual void Reset(void){};

	/** \name writing results */
	/*@{*/
	/** register data for output */
	virtual void RegisterOutput(void){};

	/** write results */
	virtual void WriteOutput(ostream& out) const{};

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
	protected:
	virtual void ReadAngledBC(const ParameterListT& list);	
	
	
	private:
	
	int fNumFacets_Tot;
	int fNumFacetNodes;
	int fNumElementNodes;
	
	double fK;   /*penalty parameter*/
	/*time schedule*/
	const ScheduleT* fSchedule;
	double fValue;
	
	/*integration domain*/
	ArrayT<DomainIntegrationT*> fDomain;

	/*Use Traction Card to store angled BC.  The applied normal displacement is analogous to the pressure*/
	ArrayT<StringT> fside_set_IDs;
	ArrayT<iArray2DT> fBC_sides;
	ArrayT<iArray2DT> fBC_global_nodes;
	ArrayT<iArray2DT> fBC_local_nodes;
	ArrayT<iArray2DT> fBC_eqnos;
	
	/*LHS and RHS*/
	ElementMatrixT fLHS; //stiffness matrix;
	dArrayT fRHS; // force vector
	dArrayT fnormal;
	dMatrixT fQ;
	dMatrixT fjacobian;
	iArrayT fnodes;
	iArrayT feqnos;
	dArrayT fip_disp;
	LocalArrayT fInitCoords;
	LocalArrayT fDisp;

};

} //namespace Tahoe

#endif /*ANGLED_KCBT_H_*/