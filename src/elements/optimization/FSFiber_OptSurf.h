/* $Id: FSFiber_OptSurf.h,v 1.1 2009-04-23 03:03:43 thao Exp $ */
/*Class to calculate objective function and gradient for inverse elasticity problems using the adjoint method */
/*A Oberai, NH Gakhale, GR Feijoo (2003) Inverse Problems 19:297-313*/

/*currently parameters are read in and acted on at the element level only.*/
/*assumes a homogeneous distribution*/
#ifndef _FSFiber_OptSurf_
#define _FSFiber_OptSurf_

/* base class */
#include "FSFiber_Optimize_Dual.h"
#include "DomainIntegrationT.h"


namespace Tahoe {

/** Interface for linear strain deformation and field gradients */
class FSFiber_OptSurf: public FSFiber_Optimize_Dual
{
  public:
      
	/** constructor */
	FSFiber_OptSurf(const ElementSupportT& support);


	virtual void InitialCondition(void);
	
	virtual void InitStep(void);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/


  protected:
  
	/** \name form the residual force vector */
	virtual void RHSDriver(void);

	/** calculate the internal force contribution*/
	virtual void FormKd(double constK);

	/** calculate objective function*/
	virtual void Compute_Cost(void);

	/** calculate objective function*/
	virtual void Compute_Gradients(void);

	/** form shape functions and derivatives */
	virtual void SetGlobalShape(void);
		
	void ApplySurfaceForce(void);
		
	private:
	
	int fNumFacets;
	int fNumFacetNodes;
	/*Use traction card scheme to store angled BC.  The applied normal displacement is analogous to the pressure*/
	ArrayT<StringT> fside_set_IDs;
	ArrayT<StringT> fblock_IDs;
	ArrayT<StringT> fField_Params;
	iArrayT fparam_map;
	iArrayT fblock_index;
	
	ArrayT<iArray2DT> fdata_sides;
	ArrayT<iArray2DT> fdata_global_nodes;
	ArrayT<iArray2DT> fdata_local_nodes;
	ArrayT<iArray2DT> fdata_eqnos;
	
	LocalArrayT fLocSurfDisp;
	LocalArrayT fLocSurfCoords;
		
	/*workspaces*/
	iArrayT fsurf_nodes;
	iArrayT feqnos;
	dArrayT fnormal;
	dMatrixT fQ;
	dMatrixT fjacobian;
	
	dArrayT fsurf_rhs;
};

}
#endif /* _FSFiber_OptSurf_Dual_ */

