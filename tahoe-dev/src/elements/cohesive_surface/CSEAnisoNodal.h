/* $Id: CSEAnisoNodal.h,v 1.1 2006-08-30 17:32:01 thao Exp $ */
/* created: paklein (11/19/1997) */
#ifndef _CSE_ANISO_NODAL_H_
#define _CSE_ANISO_NODAL_H_

/* base class */
#include "CSEAnisoT.h"

namespace Tahoe {

/** Cohesive surface elements with vector argument cohesive relations. */
class CSEAnisoNodal: public CSEAnisoT
{
public:

	/** constructors */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	CSEAnisoNodal(const ElementSupportT& support);
#else
	CSEAnisoNodal(ElementSupportT& support, bool rotate);
#endif
	
	/** close current time increment */
	virtual void CloseStep(void);

	/** restore the element group to its state at the beginning of the
	 * current time step. */
	virtual GlobalT::RelaxCodeT ResetStep(void); 

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** write restart data to the output stream. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. */
	virtual void ReadRestart(istream& in);
#endif

	/** state variable array */
	dArray2DT& StateVariables(void) { return fStateVariables; };
	/** state variable array */
	dArray2DT& StateVariables_n(void) { return fStateVariables_n; };

	virtual const int NumIP() const;
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** set the active elements.
	 * \param array of status flags for all elements in the group */
	virtual void SetStatus(const ArrayT<ElementCardT::StatusT>& status);

	void InitializeStateVars(const iArrayT& elem_nodes, const dArray2DT& traction);
	void ResetStateVars(const iArrayT& elem_nodes, dArray2DT& traction);
	int ComputeTraction(const iArrayT& elem_nodes, const dArray2DT& nodal_stress, dArray2DT& traction);
		
protected:

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

private:
	double ComputeRotation(const dMatrixT& jacobian, dMatrixT& Q); 
	double ComputeRotation(const dArray2DT& DNa, const dMatrixT& jacobian, dMatrixT& Q, ArrayT<dMatrixT>& dQ); 

protected:

	int fNumIP;

	/** \name state variable storage arrays. 
	 * arrays have dimensions: [npairs] x [nvar] */
	/*@{*/
	dArray2DT fStateVariables;
	
	/** previous converged solution */
	dArray2DT fStateVariables_n;
	/*@}*/
	int fNumStateVars;

	LocalArrayT fFacetCoords; /*nodal coordinates of mid plane in local ordering, nfn x nsd*/
	dArrayT fdelta_glob;  /*nodal jump vectors in local ordering, nfn x nsd*/		

	dArray2DT fjump;  /*matrix to transform nodal displacements (nen) to nodal jump vectors (nfn)*/
	dArray2DT fmid;   /*matrix to calcualate mid plane nodal coords (nfn)*/
	dArrayT fBa;      /*shape function of jump vector (nen)*/
	dArray2DT fDHa;	  /*shape function derivatives of jump vector*/
	dMatrixT fId;
	/*@}*/
	
	/*work space*/
	dArrayT fperm_dx;
	dMatrixT fperm_dx_dxi;
	dMatrixT fperm_dx_deta;
	
	dArrayT fNodalArea;
	dMatrixT fJacobian; 	/*jacobian calculated at nodes (nsd x nsd-1)*/
	dArrayT fNa;  /*nodal shape functions for a given IP (nfn)*/
	dArray2DT fDNa;  /*nodal shape functions for a given IP (nfn x nsd)*/
	dArrayT fCoords; /*parent domain coords at nodes (nsd-1)*/	
	iArrayT facet1;
	iArrayT facet2;
	iArrayT fNodePairMap;
	

};

inline const int CSEAnisoNodal::NumIP(void) const
{
	return(fNumIP);
}
 
} // namespace Tahoe 
#endif /* _CSE_ANISO_NODAL_H_ */
