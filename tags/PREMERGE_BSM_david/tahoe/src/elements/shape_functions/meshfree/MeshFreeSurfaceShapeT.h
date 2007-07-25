/* $Id: MeshFreeSurfaceShapeT.h,v 1.3 2002-07-05 22:28:37 paklein Exp $ */
/* created: paklein (06/03/2000)                                          */
/* Class to manage CSE integrals, where the dimension of                  */
/* the field variable is 1 greater than the dimension of the parent       */
/* domain. Jump quantities imply jump between any field variable          */
/* across the CSE.                                                        */

#ifndef _MF_SURFACE_SHAPE_T_H_
#define _MF_SURFACE_SHAPE_T_H_

/* direct members */
#include "SurfaceShapeT.h"
#include "GeometryT.h"
#include "nMatrixGroupT.h"
#include "VariArrayT.h"
#include "nArrayGroupT.h"
#include "VariLocalArrayT.h"
#include "nVariArray2DT.h"
#include "nArray2DGroupT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class MeshFreeSurfaceSupportT;
template <class TYPE> class RaggedArray2DT;

class MeshFreeSurfaceShapeT
{
public:

	/* constructors */	
	MeshFreeSurfaceShapeT(GeometryT::CodeT geometry_code, int num_ip,
		MeshFreeSupportT& mf_support, const LocalArrayT& loc_disp,
		const dArray2DT& facet_coords, int num_facet_nodes,
		bool store_shape);

	/* destructors */	
	~MeshFreeSurfaceShapeT(void);

	/* set all local parameters */
	void Initialize(void);

	/* recompute selected facet shape functions */
	void ResetFacets(const ArrayT<int>* reset_facets = NULL);

	/* access to full neighbors database */
	const ArrayT<int>& NeighborCounts(int side) const;
	const RaggedArray2DT<int>& Neighbors(int side) const;
	void NeighborCounts(ArrayT<int>& counts) const;
	void Neighbors(RaggedArray2DT<int>& neighbors) const;
		// neighbors must be configured first using counts
	
	/* initialize data for the specified facet - returns total
	 * number of facet nodes */
	int SetFacet(int facet);
	
	/* local node numbers on each side of the current facet */
	const iArrayT& NodesOnFacet(void) const;
	
	/* integration management */
	void TopIP(void);
	int  NextIP(void);
	
/**** for the current integration point ****/
	const int& CurrIP(void) const;
	double IPWeight(void) const;

	/* (undeformed) integration point coordinates */
	const dArrayT& IPCoords(void);

	/* jump in the nodal values (from side 1 to 2) */
	const dArrayT& InterpolateJumpU(const LocalArrayT& nodal_1,
		const LocalArrayT& nodal_2); //KEEP?? - not with access to
//         nodes on facets separately		
	const dArrayT& InterpolateJumpU(const LocalArrayT& nodalU);
	
	/* jump gradient tables:
	 *
	 *     fgrad_d = d delta_i/d u_j	[i] = FieldDim
	 *                              	[j] = NumNodes*FieldDim
*
	 *     fgrad_dTgrad_d = d delta_k/d u_i	d delta_k/d u_j
	 *                          	[i],[j] = NumNodes*FieldDim
	 */
	const dMatrixT& Grad_d(void) const;
	const dMatrixT& Grad_dTGrad_d(void) const;

	/* jacobian of the area transformation */
	void Jacobian(double& j0, double& j);
	void Jacobian(double& j0, double& j, dMatrixT& Q);
	void Jacobian(double& j0, double& j, dMatrixT& Q, ArrayT<dMatrixT>& dQ);

/*******************************************/

private:

	/* configure work space arrays */
	void Construct(void);

	/* compute the jacobian of the mapping to the current coordinates
	 * at the current integration point*/
	void SetDomainJacobian(void);

	/* compute Q derivative for the given side */
	void Set_dQ(const dMatrixT& Q, double j, ArrayT<dMatrixT>& dQ);

	/* set shape function tables */
	void SetShapeFunctionTables(void);

	/* set Jacobian derivatives */
	void SetJacobianDerivatives(void);

private:

	/* local displacement arrays */
	int fFieldDim;
	const LocalArrayT& fLocDisp; // [ndof] x [nnd_tot]

	/* midline facet coordinates */
	LocalArrayT fFacetCoords;
	dArray2DT fRefFacetCoords;

	/* integration shape functions */
	SurfaceShapeT fRefSurfaceShape;
	const int& fCurrIP;

	/* meshfree surface support */
	MeshFreeSurfaceSupportT *fMFSurfaceSupport;

	iArrayT fneighbors;      // [nnd_tot]
	dArray2DT fjump_phi;     // [nip] x [nnd_tot]
	ArrayT<dArray2DT> fDphi; // [nip] : [nsd] x [nnd_tot]
	ArrayT<dArray2DT> fDphi_tmp; // [nip] : [nsd] x [nnd_1,2]

	VariArrayT<int> fneighbors_man;   // [nnd_tot]
	nVariArray2DT<double> fphi_man;   // [nip] x [nnd_tot]
	nArray2DGroupT<double> fDphi_man; // [nsd] x [nnd_tot]
	

	/* shape function tables */
	ArrayT<dMatrixT> fgrad_d;        // [nip] : [fFieldDim] x [fFieldDim*nnd_tot]
	ArrayT<dMatrixT> fgrad_dTgrad_d; // [nip] : [fFieldDim*nnd_tot] x [fFieldDim*nnd_tot]
	nMatrixGroupT<double> fMatrixManager_1; // [fFieldDim] x [fFieldDim*nnd_tot]
	nMatrixGroupT<double> fMatrixManager_2; // [fFieldDim*nnd_tot] x [fFieldDim*nnd_tot]

	/* shape function jacobian derivative tables (current ip) */
	ArrayT<dMatrixT> fdx_dsdu; // [fFieldDim - 1] : [fFieldDim] x [fFieldDim*nnd_tot]

	/* return values */
	dArrayT fInterp; // [fFieldDim]
	
	/* coordinate transformation */
	dMatrixT fRefJacobian; // [fFieldDim] x [fFieldDim - 1]
	dMatrixT fJacobian;    // [fFieldDim] x [fFieldDim - 1]
	dMatrixT fJ_tmp;       // [fFieldDim] x [fFieldDim - 1]
	dMatrixT fJ_1, fJ_2;   // [fFieldDim] x [fFieldDim]
	
	/* work space */
	dArrayT fu_vec; // [nnd*fFieldDim]
	VariArrayT<double> fVectorManager_1; // [nnd*fFieldDim]
	
	ArrayT<dArrayT> fnnd_vec; // [fFieldDim - 1] : [nnd]
	nArrayGroupT<double> fVectorManager_2; // [nnd]

	dArrayT fx_vec; // shallow
	
	/* work space for 3D jacobian and derivatives */
	dMatrixT fM1, fM2; // [fFieldDim] x [nnd*fFieldDim]
	nMatrixGroupT<double> fMatrixManager_3; // [fFieldDim] x [nnd*fFieldDim]
		//share by both sides during shape function calculations, i.e., not
		//kept constant during an element calculation
};

/* inlines */

/* integration management */
inline void MeshFreeSurfaceShapeT::TopIP(void) { fRefSurfaceShape.TopIP(); }
inline int MeshFreeSurfaceShapeT::NextIP(void) { return fRefSurfaceShape.NextIP(); }
inline const int& MeshFreeSurfaceShapeT::CurrIP(void) const { return fCurrIP; }
inline double MeshFreeSurfaceShapeT::IPWeight(void) const { return fRefSurfaceShape.IPWeight(); }

/* local node numbers on each side */
inline const iArrayT& MeshFreeSurfaceShapeT::NodesOnFacet(void) const
{
	return fneighbors;
}

inline const dMatrixT& MeshFreeSurfaceShapeT::Grad_d(void) const
{
	return fgrad_d[fCurrIP];
}

inline const dMatrixT& MeshFreeSurfaceShapeT::Grad_dTGrad_d(void) const
{
	return fgrad_dTgrad_d[fCurrIP];
}

/* (undeformed) integration point coordinates */
inline const dArrayT& MeshFreeSurfaceShapeT::IPCoords(void)
{
	return fRefSurfaceShape.IPCoords();
}

} // namespace Tahoe 
#endif /* _MF_SURFACE_SHAPE_T_H_ */
