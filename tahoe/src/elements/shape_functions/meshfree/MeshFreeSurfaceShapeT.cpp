/* $Id: MeshFreeSurfaceShapeT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (06/03/2000)                                          */

#include "MeshFreeSurfaceShapeT.h"

#include "Constants.h"
#include "ExceptionCodes.h"

#include "SurfaceShapeT.h"
#include "MeshFreeSurfaceSupportT.h"
#include "MeshFreeSupportT.h"

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

/* constructor */
MeshFreeSurfaceShapeT::MeshFreeSurfaceShapeT(GeometryT::CodeT geometry_code,
	int num_ip, MeshFreeSupportT& mf_support, const LocalArrayT& loc_disp,
	const dArray2DT& facet_coords, int num_facet_nodes,
	bool store_shape):
	/* local displacement arrays */
	fFieldDim(loc_disp.MinorDim()),
	fLocDisp(loc_disp),

	fFacetCoords(LocalArrayT::kInitCoords, num_facet_nodes,
		GeometryT::GeometryToNumSD(geometry_code) + 1),
	fRefSurfaceShape(geometry_code, num_ip, 2*num_facet_nodes, fFieldDim, fFacetCoords),
	fCurrIP(fRefSurfaceShape.CurrIP()),

	/* dynamics work space managers */
	fDphi_man(0),
	fMatrixManager_1(0),
	fMatrixManager_2(0),

	fVectorManager_2(0),
	fMatrixManager_3(0)
{
	/* meshfree surface support */
	fMFSurfaceSupport = new MeshFreeSurfaceSupportT(mf_support,
		fRefSurfaceShape, fFacetCoords, facet_coords, num_facet_nodes,
		store_shape);	
	if (!fMFSurfaceSupport) throw eOutOfMemory;

	/* dimension arrays */
	Construct();
}

/* destructors */	
MeshFreeSurfaceShapeT::~MeshFreeSurfaceShapeT(void)
{
	delete fMFSurfaceSupport;
	fMFSurfaceSupport = NULL;
}

/* set all local parameters */
void MeshFreeSurfaceShapeT::Initialize(void)
{
	/* initialize integration domain */
	fRefSurfaceShape.Initialize();

	/* initialize (all) surface shape functions */
	fMFSurfaceSupport->ResetFacets();
}

/* recompute selected facet shape functions */
void MeshFreeSurfaceShapeT::ResetFacets(const ArrayT<int>* reset_facets)
{
	fMFSurfaceSupport->ResetFacets(reset_facets);
}

/* access to full neighbors database */
const ArrayT<int>& MeshFreeSurfaceShapeT::NeighborCounts(int side) const
{
	return fMFSurfaceSupport->NeighborCounts(side);
}

const RaggedArray2DT<int>& MeshFreeSurfaceShapeT::Neighbors(int side) const
{
	return fMFSurfaceSupport->Neighbors(side);
}

void MeshFreeSurfaceShapeT::NeighborCounts(ArrayT<int>& counts) const
{
	const ArrayT<int>& counts_1 = fMFSurfaceSupport->NeighborCounts(0);
	const ArrayT<int>& counts_2 = fMFSurfaceSupport->NeighborCounts(1);
	counts.Allocate(counts_1.Length());
	for (int i = 0; i < counts.Length(); i++)
		counts[i] = counts_1[i] + counts_2[i];
}

void MeshFreeSurfaceShapeT::Neighbors(RaggedArray2DT<int>& neighbors) const
{
	const RaggedArray2DT<int>& neighbors_1 = fMFSurfaceSupport->Neighbors(0);
	const RaggedArray2DT<int>& neighbors_2 = fMFSurfaceSupport->Neighbors(1);

	/* overall check */
	if (neighbors.Length() != neighbors_1.Length() + neighbors_2.Length())
	{
		cout << "\n MeshFreeSurfaceShapeT::Neighbors: destination array must be configured" << endl;
		throw eSizeMismatch;
	}
	
	int nnd = neighbors_1.MajorDim();
	for (int i = 0; i < nnd; i++)
	{
		int nnd_1 = neighbors_1.MinorDim(i);
		int nnd_2 = neighbors_2.MinorDim(i);
		int nnd   = neighbors.MinorDim(i);

		/* row-by-row check */
		if (nnd != nnd_1 + nnd_2)
		{
			cout << "\n MeshFreeSurfaceShapeT::Neighbors: row " << i
			     << " of destination is length "
			     << nnd << "\n" <<   "     instead of length "
			     << nnd_1 + nnd_2 << endl;
			throw eSizeMismatch;
		}
	
		/* copy data */
		memcpy(neighbors(i), neighbors_1(i), sizeof(int)*nnd_1);
		memcpy(neighbors(i) + nnd_1, neighbors_2(i), sizeof(int)*nnd_2);
	}
}

/* initialize data for the specified facet - returns total
* number of facet nodes */
int MeshFreeSurfaceShapeT::SetFacet(int facet)
{
	/* total number of neighborhood nodes */
	int nnd = fMFSurfaceSupport->NumberOfNeighbors(facet);

	/* dimension joined shape function storage */
	fneighbors_man.SetLength(nnd, false);
	fphi_man.Dimension(fDphi.Length(), nnd);
	fDphi_man.Dimension(fFieldDim, nnd);

	/* (shallow) temp space */
	iArrayT neigh;
	dArray2DT phi;

	/* copy data from side 1 */
	fjump_phi = 0.0;
	fMFSurfaceSupport->LoadData(facet, 0, neigh, phi, fDphi_tmp);
	int nnd_1 = neigh.Length();
	fneighbors.CopyPart(0, neigh, 0, nnd_1);
	fjump_phi.BlockColumnCopyAt(phi, 0);
	fjump_phi *= -1.0;
	for (int i0 = 0; i0 < fDphi_tmp.Length(); i0++)
		fDphi[i0].BlockColumnCopyAt(fDphi_tmp[i0], 0);

	/* copy data from side 2 */
	fMFSurfaceSupport->LoadData(facet, 1, neigh, phi, fDphi_tmp);
	int nnd_2 = neigh.Length();
	fneighbors.CopyPart(nnd_1, neigh, 0, nnd_2);
	fjump_phi.BlockColumnCopyAt(phi, nnd_1);
	for (int i1 = 0; i1 < fDphi_tmp.Length(); i1++)
		fDphi[i1].BlockColumnCopyAt(fDphi_tmp[i1], nnd_1);

	/* dimension arrays */	
	fMatrixManager_1.Dimension(fFieldDim, fFieldDim*nnd);
	fMatrixManager_2.Dimension(fFieldDim*nnd, fFieldDim*nnd);

	/* set facet reference coordinates */
	const dArray2DT& all_facets = fMFSurfaceSupport->FacetCoords();
	fRefFacetCoords.Set(fRefSurfaceShape.NumFacetNodes(),
		fRefSurfaceShape.NumSD() + 1, all_facets(facet));
	fFacetCoords.FromTranspose(fRefFacetCoords);
	
	/* set shape function tables */
	SetShapeFunctionTables();

	/* total number of facet nodes */
	return nnd;
}

/**** for the current integration point ***/

/* jump in the nodal values (from side 1 to 2) */
const dArrayT& MeshFreeSurfaceShapeT::InterpolateJumpU(const LocalArrayT& nodal_1,
	const LocalArrayT& nodal_2)
{
#if __option(extended_errorcheck)
	if (nodal_1.NumberOfNodes() + nodal_2.NumberOfNodes() !=
	    fjump_phi.MinorDim()) throw eSizeMismatch;
#endif

	int nnd_1 = nodal_1.NumberOfNodes();
	int nnd_2 = nodal_2.NumberOfNodes();
	for (int i = 0; i < fFieldDim; i++)
	{
		double* phi = fjump_phi(fCurrIP);
	
		double jump = 0.0;
		double* u = nodal_1(i);
		for (int j1 = 0; j1 < nnd_1; j1++)
			jump += (*phi++)*(*u++);

		u = nodal_2(i);
		for (int j2 = 0; j2 < nnd_2; j2++)
			jump += (*phi++)*(*u++);
		
		fInterp[i] = jump;
	}
	return fInterp;
}

const dArrayT& MeshFreeSurfaceShapeT::InterpolateJumpU(const LocalArrayT& nodalU)
{
	for (int i = 0; i < fFieldDim; i++)
		fInterp[i] = fjump_phi.DotRow(fCurrIP, nodalU(i));
	return fInterp;
}

/* jacobian of the area transformation at the current integration
* using the nodes on the specified facet */
void MeshFreeSurfaceShapeT::Jacobian(double& j0, double& j)
{
	/* Jacobian matrix of the surface transformation */
	SetDomainJacobian();

	/* parent domain */
	const ParentDomainT& parent_domain = fRefSurfaceShape.ParentDomain();

	/* surface map jacobian */
	j0 = parent_domain.SurfaceJacobian(fRefJacobian);
	j  = parent_domain.SurfaceJacobian(fJacobian);
}

void MeshFreeSurfaceShapeT::Jacobian(double& j0, double& j, dMatrixT& Q)
{
	/* Jacobian matrix of the surface transformation */
	SetDomainJacobian();	

	/* parent domain */
	const ParentDomainT& parent_domain = fRefSurfaceShape.ParentDomain();

	/* surface map jacobian */
	j0 = parent_domain.SurfaceJacobian(fRefJacobian);
	j  = parent_domain.SurfaceJacobian(fJacobian, Q);
}

void MeshFreeSurfaceShapeT::Jacobian(double& j0, double& j, dMatrixT& Q, ArrayT<dMatrixT>& dQ)
{
	/* Jacobian matrix of the surface transformation */
	SetDomainJacobian();	

	/* parent domain */
	const ParentDomainT& parent_domain = fRefSurfaceShape.ParentDomain();

	/* get Jacobian matrix of the surface transformation */
	j0 = parent_domain.SurfaceJacobian(fRefJacobian);
	j  = parent_domain.SurfaceJacobian(fJacobian, Q);

	/* set Jacobian derivatives */
	SetJacobianDerivatives();

	/* components for each side */
	Set_dQ(Q, j, dQ);
}

/**** for the current integration point ***/

/***********************************************************************
* Private
***********************************************************************/

/* configure work space arrays */
void MeshFreeSurfaceShapeT::Construct(void)
{
	/* dimensions */
	int num_ip = fRefSurfaceShape.NumIP();

	/* shape function derivatives array */
	fDphi.Allocate(num_ip);
	fDphi_tmp.Allocate(num_ip);

	fneighbors_man.SetWard(0, fneighbors);
	fphi_man.SetWard(0, fjump_phi, 0);
	for (int ii = 0; ii < num_ip; ii++)
		fDphi_man.Register(fDphi[ii]);

	/* shape function and derivatives tables */
	fgrad_d.Allocate(num_ip);
	fgrad_dTgrad_d.Allocate(num_ip);

	/* shape function jacobian derivative tables (current ip) */
	fdx_dsdu.Allocate(fFieldDim-1);

	/* register with dynamic workspace managers */
	for (int i = 0; i < num_ip; i++)
	{		
		fMatrixManager_1.Register(fgrad_d[i]);
		fMatrixManager_2.Register(fgrad_dTgrad_d[i]);
	}

	for (int k = 0; k < fFieldDim-1; k++)
		fMatrixManager_1.Register(fdx_dsdu[k]);

	/* return value */
	fInterp.Allocate(fFieldDim);
	
	/* coordinate transformations */	
	fRefJacobian.Allocate(fFieldDim, fFieldDim-1);
	fJacobian.Allocate(fFieldDim, fFieldDim-1);
	fJ_tmp.Allocate(fFieldDim, fFieldDim-1);
	fJ_1.Allocate(fFieldDim, fFieldDim);
	fJ_2.Allocate(fFieldDim, fFieldDim);
	
	/* work space */
	fVectorManager_1.SetWard(0, fu_vec);

	fnnd_vec.Allocate(fFieldDim-1);
	for (int j = 0; j < fnnd_vec.Length(); j++)
		fVectorManager_2.Register(fnnd_vec[j]);
	
	/* 3D work space */
	if (fFieldDim == 3)
	{
		fMatrixManager_3.Register(fM1);
		fMatrixManager_3.Register(fM2);
	}
}

/* compute the jacobian of the mapping to the current coordinates
* at the current integration point*/
void MeshFreeSurfaceShapeT::SetDomainJacobian(void)
{
	/* parent domain */
	const ParentDomainT& parent_domain = fRefSurfaceShape.ParentDomain();

	/* reference transformation */
	parent_domain.DomainJacobian(fFacetCoords, fCurrIP, fRefJacobian);

	/* reference displacement gradients */
	parent_domain.Jacobian(fLocDisp, fDphi[fCurrIP], fJ_1);
	
	/* chain rule */
	fJ_tmp.MultAB(fJ_1, fRefJacobian);
	
	/* final jacobian */
	fJacobian.SetToCombination(1.0, fRefJacobian, 0.5, fJ_tmp);
}

void MeshFreeSurfaceShapeT::Set_dQ(const dMatrixT& Q, double j, ArrayT<dMatrixT>& dQ)
{
	/* checks */
	if (dQ.Length() != fFieldDim) throw eSizeMismatch;

	/* dimensions */
	int nnd = fneighbors.Length();
	int nu = nnd*fFieldDim;

	/* resize */
	fVectorManager_1.SetLength(nu, false);

	/* 2D */
	if (dQ.Length() == 2)
	{
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];
	
		/* unit tangent */
		fx_vec.Set(2, Q(0));
		dMatrixT& dtan_du = fdx_dsdu[0];
		dtan_du.MultTx(fx_vec, fu_vec);
	
		/* first component */
		dQ1.Outer(fx_vec, fu_vec);
		dQ1 -= dtan_du;
		dQ1 /= -j;
		
		/* second component - multiply dQ1 by Q^(pi/2) */
		int num_cols = dQ1.Cols();
		double* row1 = dQ1.Pointer();
		double* row2 = row1 + 1;
		double* pdQ2 = dQ2.Pointer();
		for (int i = 0; i < num_cols; i++)
		{
			*pdQ2++ =-(*row2);
			*pdQ2++ =  *row1;
			
			row1 += 2;
			row2 += 2;	
		}
	}
	else /* 3D */
	{
		/* resize */
		fMatrixManager_3.Dimension(fFieldDim, nu);
	
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];
		dMatrixT& dQ3 = dQ[2];

		/* tangent vectors */
		double* v_m1 = fJacobian(0);
		double* v_m2 = fJacobian(1);
		double    m1 = sqrt(v_m1[0]*v_m1[0] + v_m1[1]*v_m1[1] + v_m1[2]*v_m1[2]);
		if (m1 < kSmall) throw eBadJacobianDet;

		/* tangent gradients */
		ArrayT<dMatrixT>& grad_dd = fdx_dsdu;
		dMatrixT& dm1_du = grad_dd[0];
		dMatrixT& dm2_du = grad_dd[1];

		/* first component */
		fx_vec.Set(3, Q(0));
		dm1_du.MultTx(fx_vec, fu_vec);
		dQ1.Outer(fx_vec, fu_vec);
		dQ1 -= dm1_du;
		dQ1 /= -m1;
		
		/* compute d_normal/d_u */
		for (int i = 0; i < nu; i++)
		{
			CrossProduct(dm1_du(i), v_m2, fM1(i));
			CrossProduct(v_m1, dm2_du(i), fM2(i));
		}
		fM1 += fM2;
		
		/* third component */
		fx_vec.Set(3, Q(2));
		fM1.MultTx(fx_vec, fu_vec);
		dQ3.Outer(fx_vec, fu_vec);
		dQ3 -= fM1;
		dQ3 /= -j;
		
		/* second component */
		double* t1 = Q(0);
		double*  n = Q(2);
		for (int k = 0; k < nu; k++)
		{
			CrossProduct(dQ3(k), t1, dQ2(k));
			CrossProduct(n , dQ1(k), fM1(k));
		}
		dQ2 += fM1;
	}	
}

/* set shape function tables */
void MeshFreeSurfaceShapeT::SetShapeFunctionTables(void)
{
	/* dimensions */
	int nnd   = fneighbors.Length();

	/* set tables */
	dMatrixT shNaMat;
	int num_ip = fRefSurfaceShape.NumIP();
	for (int i = 0; i < num_ip; i++)
	{
		shNaMat.Set(1, nnd, fjump_phi(i));
		fgrad_d[i] = 0.0;
		fgrad_d[i].Expand(shNaMat, fFieldDim);
		fgrad_dTgrad_d[i].MultATB(fgrad_d[i], fgrad_d[i]);
	}
}

/* set Jacobian derivatives */
void MeshFreeSurfaceShapeT::SetJacobianDerivatives(void)
{
	/* work space */
	int nx = fJacobian.Rows();
	int ns = fJacobian.Cols();
	dArrayT dphi_dx(nx); //TEMP - make workspace variable
	dArrayT dphi_ds(ns); //TEMP - make workspace variable
	dMatrixT shNaMat;

	/* facet parameters */
	int nnd = fneighbors.Length();
	const dArray2DT& Dphi = fDphi[fCurrIP];
		
	/* dimension work space */
	fVectorManager_2.Dimension(nnd, false);
		
	/* loop over facet nodes */
	for (int i = 0; i < nnd; i++)
	{
		/* chain rule derivative */
		Dphi.ColumnCopy(i, dphi_dx);
		fJacobian.MultTx(dphi_dx, dphi_ds);
		
		/* loop over domain dimensions */
		for (int j = 0; j < ns; j++)
			fnnd_vec[j][i] = 0.5*dphi_ds[j];
	}
	
	/* loop over domain dimensions */
	for (int j = 0; j < ns; j++)
	{
		/* shallow matrix */
		shNaMat.Set(1, nnd, fnnd_vec[j].Pointer());
		
		/* expand */
		fdx_dsdu[j] = 0.0;
		fdx_dsdu[j].Expand(shNaMat, fFieldDim);	
	}
}
