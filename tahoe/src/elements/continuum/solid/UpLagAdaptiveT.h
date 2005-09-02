/* $Id: UpLagAdaptiveT.h,v 1.3 2005-02-13 22:18:02 paklein Exp $ */
#ifndef _UPDATED_LAGRANGIAN_ADAPTIVE_T_H_
#define _UPDATED_LAGRANGIAN_ADAPTIVE_T_H_

/* configuration */
#include "ElementsConfig.h"
#ifdef COHESIVE_SURFACE_ELEMENT

/* base class */
#include "UpdatedLagrangianT.h"

/* direct members */
#include "InverseMapT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/** forward declarations */
class CSEAnisoT;
class TiedNodesT;

/** update Lagrangian, finite strain solid */
class UpLagAdaptiveT: public UpdatedLagrangianT
{
public:

	/** constructors */
	UpLagAdaptiveT(const ElementSupportT& support, const FieldT& field);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** compute residual */
	virtual void RHSDriver(void);

	/** calculate the internal force contribution ("-k*d"). Accumulates nodal
	 * stresses. */
	virtual void FormKd(double constK);

	/** element level reconfiguration for the current time increment. This
	 * provides an interface for element-level adaptivity. The nature of
	 * the changes are indicated by the GlobalT::RelaxCodeT return value. */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

private:

	/** determine the tied nodes and reset the constraints based on the list of
	 * active elements */
	void SetNetwork(const ArrayT<ElementCardT::StatusT>& active_elements);

	/** generate the list of leaders for all nodes based on the active element list */
	void FindLeaders(const iArray2DT& connects, const ArrayT<ElementCardT::StatusT>& active, iArrayT& same_as) const;

protected:

	/** release value */
	double fReleaseThreshold;

	/** cohesive element group */
	CSEAnisoT* fCSE;

	/** nodes used by cohesive elements */
	iArrayT fCSENodesUsed;
	
	/** inverse map of global to local number for CSE nodes used */
	InverseMapT fGlobalToLocal;
	
	/** connectivity of cohesive elements in local numbering */
	iArray2DT fConnectivitiesCSELocal;
	
	/** active flag */
	ArrayT<ElementCardT::StatusT> fCSEActive;

	/** inverse connectivities of CSE, elements per node */
	RaggedArray2DT<int> fInverseConnectivitiesCSELocal;

	/** leaders for all nodes used, -1 if none */
	iArrayT fSameAs;
	
	/** enforc tied nodes */
	TiedNodesT* fTied;

	/** \name nodal stresses */
	/*@{*/
	iArrayT   fAvgCount;
	dArray2DT fNodalValues;
	dArray2DT fNodalExtrapolation;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* COHESIVE_SURFACE_ELEMENT */

#endif /* _UPDATED_LAGRANGIAN_ADAPTIVE_T_H_ */
