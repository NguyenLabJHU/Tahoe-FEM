/* $Id: TotalLagrangianCBSurfaceT.h,v 1.4 2005-07-14 05:27:28 paklein Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_
#define _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_

/* base class */
#include "TotalLagrangianT.h"

namespace Tahoe {

class FCC3D_Surf;

/** total Lagrangian, finite strain element for working with Cauchy-Born approach
 * for modeling surface effects */
class TotalLagrangianCBSurfaceT: public TotalLagrangianT
{
public:

	/** constructor */
	TotalLagrangianCBSurfaceT(const ElementSupportT& support);

	/** destructor */
	virtual ~TotalLagrangianCBSurfaceT(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/
		
protected:

	/** list of elements on the surface */
	iArrayT fSurfaceElements;

	/** elements neighbors */
	iArray2DT fSurfaceElementNeighbors;

	/** surface model number */
	iArray2DT fSurfaceElementFacesType;

	/** surface normals */
	ArrayT<dArrayT> fNormal;

	/** surface Cauchy-Born models */
	ArrayT<FCC3D_Surf*> fSurfaceCB;

	/** support for the surface models */
	FSMatSupportT* fSurfaceCBSupport;

	/** deformation gradients at the surface integration points */
	ArrayT<dMatrixT> fF_Surf_List;
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_ */
