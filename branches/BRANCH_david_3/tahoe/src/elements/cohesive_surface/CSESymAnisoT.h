	/* $Id: CSESymAnisoT.h,v 1.6 2004-09-10 21:40:05 paklein Exp $ */
#ifndef _CSE_SYM_ANISO_T_H_
#define _CSE_SYM_ANISO_T_H_

/* base class */
#include "CSEAnisoT.h"

namespace Tahoe {

/* forward declarations */
class SurfacePotentialT;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
class TiedPotentialBaseT;
#endif

/** Cohesive surface elements with only 1 "active" facet. The other
 *  facet's displacements are presumed symmetry-related to the active
 *  facet. Consequently, these elements have half as many nodes as
 *  a CSE with two active facets. Currently, enforced symmetry is 
 *  hardcoded to be mode I, but this can be changed.  
 */
class CSESymAnisoT: public CSEAnisoT
{
public:

	/* constructors */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	CSESymAnisoT(const ElementSupportT& support);
#else
	CSESymAnisoT(ElementSupportT& support, bool rotate);
#endif

	/* writing output */
	virtual void RegisterOutput(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name construction of connectivities */
	/*@{*/
	/** extract element block info from parameter list to be used. Read element 
	 * connectivity data as a sideset since only half the number of element nodes 
	 * are expected from a single face of the continuum elements*/
	virtual void CollectBlockInfo(const ParameterListT& list, ArrayT<StringT>& block_ID,  
		ArrayT<int>& mat_index) const;
	/*@}*/

	/** read element connectivity data as a sideset since only half the
	    number of element nodes are expected from a single face of the
	    continuum elements */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	virtual void ReadConnectivity(ifstreamT& in, ostream& out);
#else
	virtual void ReadConnectivity(void);
#endif 

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);
	
	/* extrapolate the integration point stresses and strains and extrapolate */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
						const iArrayT& e_codes, dArray2DT& e_values);

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;
	
	/* number of facet nodes as a function of number of element nodes */
	virtual int NumFacetNodes(void) { return NumElementNodes(); }

private:

	/* operations with pseudo rank 3 (list in j) matrices */
	void u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
		dMatrixT& Qu);

	void Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
		dMatrixT& Qu);

private:

	/** stiffness matrix */
	dMatrixT fK;
	
	/** information about the side sets for OutputSetT */
	ArrayT<StringT> fSideSet_ID;

};

} // namespace Tahoe 
#endif /* _CSE_SYM_ANISO_T_H_ */
