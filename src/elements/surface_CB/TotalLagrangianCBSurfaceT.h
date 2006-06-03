/* $Id: TotalLagrangianCBSurfaceT.h,v 1.8 2006-06-03 23:04:13 hspark Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_
#define _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_

/* base class */
#include "TotalLagrangianT.h"

namespace Tahoe {

class FCC3D_Surf;
class EAMFCC3DMatT_surf;

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

	/** reduce the coordinates to a surface layer on the given face
	 *\param coords should enter with the coordinates of entire element and
	 *       returns with the coordinates defining the surface layer */
	void SurfaceLayer(LocalArrayT& coords, int face, double thickness) const;

protected:

	/** list of elements on the surface */
	iArrayT fSurfaceElements;

	/** elements neighbors */
	iArray2DT fSurfaceElementNeighbors;

	/** surface model number */
	iArray2DT fSurfaceElementFacesType;

	/** surface normals */
	ArrayT<dArrayT> fNormal;

    /** EAM Surface CB model */
	ArrayT<EAMFCC3DMatT_surf*> fEAMSurfaceCB;

	/** surface Cauchy-Born models */
	ArrayT<FCC3D_Surf*> fSurfaceCB;

	/** support for the surface models */
	FSMatSupportT* fSurfaceCBSupport;

	/** deformation gradients at the surface integration points */
	ArrayT<dMatrixT> fF_Surf_List;

	/** indicator for EAM_CB or FCC_CB */
	StringT fIndicator;

	/** \name split integration */
	/*@{*/
	LocalArrayT fSplitInitCoords;
	ShapeFunctionT* fSplitShapes;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_ */
