/* $Id: SurfaceShapeT.cpp,v 1.7.2.3 2002-10-19 17:46:44 paklein Exp $ */
/* created: paklein (11/21/1997) */
#include "SurfaceShapeT.h"

#include "toolboxConstants.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

/* constructor */
SurfaceShapeT::SurfaceShapeT(GeometryT::CodeT geometry_code, int num_ip,
	int num_nodes, int field_dim, const LocalArrayT& coords):
	DomainIntegrationT(geometry_code, num_ip, num_nodes/2),
	fTotalNodes(num_nodes),
	fNumFacetNodes(coords.NumberOfNodes()),
	fFieldDim(field_dim),
	fCoords(coords),
	fFacetCoords(fCoords.Type())
{
	/* coordinate mode */
	if (fCoords.NumberOfNodes() == fNumFacetNodes)
		fFacetCoords.Alias(fCoords);
	else
		fFacetCoords.Allocate(fNumFacetNodes, fCoords.MinorDim());
	
	/* dimension arrays */
	Construct();

	/* element node numbering */
	SetNodesOnFacets(fFacetNodes);
}

SurfaceShapeT::SurfaceShapeT(const SurfaceShapeT& link, const LocalArrayT& coords):
	DomainIntegrationT(link),
	fTotalNodes(link.fTotalNodes),
	fNumFacetNodes(link.fNumFacetNodes),
	fFieldDim(link.fFieldDim),
	fCoords(coords),
	fFacetCoords(fCoords.Type())
{
	/* coordinate mode */
	if (fCoords.NumberOfNodes() == fNumFacetNodes)
		fFacetCoords.Alias(fCoords);
	else
		fFacetCoords.Allocate(fNumFacetNodes, fCoords.MinorDim());

	/* dimension arrays */
	Construct();

	/* element node numbering */
	SetNodesOnFacets(fFacetNodes);
}	

/* set all local parameters */
void SurfaceShapeT::Initialize(void)
{
	/* inherited */
	DomainIntegrationT::Initialize();
	
	/* get jump vector */
	iArrayT jump(fTotalNodes);
	SetJumpVector(jump);
	
	/* work space */
	dArrayT dshape(fTotalNodes);
		
	/* set tables */
	for (int i = 0; i < fNumIP; i++)
	{
		double* Na = fNa(i);
		double* jumpNa = fjumpNa(i);
		int* sign = jump.Pointer();
		
		/* loop over facets */
		for (int k = 0; k < 2; k++)
		{
	    	const double* shape = fDomain->Shape(i);
	    	int* pfacet = fFacetNodes(k);
			for (int j = 0; j < fNumFacetNodes; j++) //ISO
			{
				int dex = *pfacet++;
				Na[dex] = *shape;
				jumpNa[dex] = sign[dex]*(*shape++);
			}
		}

		/* shape function tables */
		dMatrixT shNaMat(1, fTotalNodes, fjumpNa(i));
		fgrad_d[i].Expand(shNaMat, fFieldDim, dMatrixT::kOverwrite);
		fgrad_dTgrad_d[i].MultATB(fgrad_d[i], fgrad_d[i]);
		
		/* shape function derivative tables */
		for (int j = 0; j < (fFieldDim-1); j++)
		{
			/* expand over facets */
			for (int k = 0; k < 2; k++)
			{
				const double* pdshape = fDomain->DShape(i,j);
	    		int* pfacet = fFacetNodes(k);
				for (int l = 0; l < fNumFacetNodes; l++) //ISO
				{
					int dex = *pfacet++;
					dshape[dex] = 0.5*(*pdshape++);
				}	
			}
			
			shNaMat.Set(1, fTotalNodes, dshape.Pointer());
			fgrad_dd(i,j).Expand(shNaMat, fFieldDim, dMatrixT::kOverwrite);
		}
	}	
}

/**** for the current integration point ***/

void SurfaceShapeT::InterpolateJump(const LocalArrayT& nodal, dArrayT& jump) const
{
	int dim = jump.Length();
	for (int i = 0; i < dim; i++)
		jump[i] = fjumpNa.DotRow(fCurrIP, nodal(i));
}

/* interpolate field values to the current integration point */
void SurfaceShapeT::Interpolate(const LocalArrayT& nodal, dArrayT& u) const
{
#if __option(extended_errorcheck)
	if (u.Length() != nodal.MinorDim()) throw ExceptionT::kSizeMismatch;
	if (nodal.NumberOfNodes() != TotalNodes() &&
	    nodal.NumberOfNodes() != NumFacetNodes()) throw ExceptionT::kSizeMismatch;
#endif

	/* average across both sides if all values given */
	bool both_sides = nodal.NumberOfNodes() == TotalNodes();
	double scale = (both_sides) ? 0.5 : 1.0;

	/* reference to the shape functions for one face */
	const dArray2DT& shapes = Na();

	/* a little tricky here because node numbering across
	 * both faces is inconsistent between 2D and 3D */
	int* face_nodes = fFacetNodes(1); /* nodes on 2nd face */
	for (int i = 0; i < u.Length(); i++)
	{	
		/* first face */
		double* p = nodal(i);
		u[i] = scale*shapes.DotRow(fCurrIP, p);
		
		/* second face */
		if (both_sides)
		{
			for (int j = 0; j < NumFacetNodes(); j++)
				u[i] += scale*p[face_nodes[j]]*shapes(fCurrIP, j);
		}
	}	
}

/* extrapolate integration point values to the nodes
*    IPvalues[numvals] : values from a single integration point
*    nodalvalues[fNumNodes x numvals] : extrapolated values */
void SurfaceShapeT::Extrapolate(const dArrayT& IPvalues,
	dArray2DT& nodalvalues)
{
	/* resize workspace */
	fNodalValues.Allocate(fNumFacetNodes, IPvalues.Length());

	/* initialize */
	fNodalValues = 0.0;

	/* extrapolate for one facet */
	fDomain->NodalValues(IPvalues, fNodalValues, CurrIP());

	/* copy values to nodes on the facets */
	for (int i = 0; i < IPvalues.Length(); i++)
		for (int k = 0; k < 2; k++) /* loop over facets */
		{
	    	int* dex = fFacetNodes(k);
			for (int j = 0; j < fNumFacetNodes; j++)
				nodalvalues(*dex++, i) += fNodalValues(j, i);
		}
}	

double SurfaceShapeT::Jacobian(dMatrixT& Q, ArrayT<dMatrixT>& dQ)
{
#if __option(extended_errorcheck)
	if (dQ.Length() != fFieldDim) throw ExceptionT::kSizeMismatch;
#endif

	/* compute facet coordinates */
	if (fCurrIP == 0 && fCoords.NumberOfNodes() != fNumFacetNodes)
		ComputeFacetCoords();

	/* get Jacobian matrix of the surface transformation */
	fDomain->DomainJacobian(fFacetCoords, fCurrIP, fJacobian);	
	double j = fDomain->SurfaceJacobian(fJacobian, Q);
	if (j <= 0.0) throw ExceptionT::kBadJacobianDet;

//NOTE: everything from here down depends only on Q

	if (dQ.Length() == 2)
	{
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];
	
		/* unit tangent */
		fx_vec.Set(2, Q(0));
		dMatrixT& dtan_du = fgrad_dd(CurrIP(), 0);
		dtan_du.MultTx(fx_vec, fu_vec);
	
		/* first component */
		dQ1.Outer(fx_vec, fu_vec);
		dQ1 -= dtan_du;
		dQ1 /= -j;
		
		/* second component - multiply dQ1 by Q^(-pi/2) */
		int num_cols = dQ1.Cols();
		double* row1 = dQ1.Pointer();
		double* row2 = row1 + 1;
		double* pdQ2 = dQ2.Pointer();
		for (int i = 0; i < num_cols; i++)
		{
			*pdQ2++ = *row2;
			*pdQ2++ =-(*row1);
			
			row1 += 2;
			row2 += 2;	
		}
	}
	else
	{
		/* dimensions */
		int nu = fu_vec.Length();
	
		/* components of the rank 3 tensor */
		dMatrixT& dQ1 = dQ[0];
		dMatrixT& dQ2 = dQ[1];
		dMatrixT& dQ3 = dQ[2];

		/* tangent vectors */
		double* v_m1 = fJacobian(0);
		double* v_m2 = fJacobian(1);
		double    m1 = sqrt(v_m1[0]*v_m1[0] + v_m1[1]*v_m1[1] + v_m1[2]*v_m1[2]);
		if (m1 <= 0.0) throw ExceptionT::kBadJacobianDet;

		/* tangent gradients */
		dMatrixT& dm1_du = fgrad_dd(CurrIP(), 0);
		dMatrixT& dm2_du = fgrad_dd(CurrIP(), 1);

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

	return j;
}

void SurfaceShapeT::Shapes(dArrayT& Na) const
{
	/* nodes from both faces */
	if (Na.Length() == fNa.MinorDim())
		fNa.RowAlias(fCurrIP, Na);
	else if (Na.Length() == fFacetNodes.MinorDim()) /* just nodes from the first face */
	{
		double* pNa = fNa(fCurrIP);
		int* lnd = fFacetNodes(0);
		for (int i = 0; i < Na.Length(); i++)
			Na[i] = pNa[lnd[i]];
	}
	else throw ExceptionT::kBadInputValue;
}

/*******************************************/

/* local node numbers on each facet */
void SurfaceShapeT::SetNodesOnFacets(iArray2DT& facetnodes)
{
	/* check */
	if (facetnodes.MajorDim() != 2 &&
	    facetnodes.MinorDim() != fNumFacetNodes)  throw ExceptionT::kSizeMismatch;

	int num_nodes_error = 0;
	int geometry = fDomain->GeometryCode();
	switch (geometry)
	{
		case GeometryT::kLine:
		{
			if (fTotalNodes == 4)
			{
//				int dat_4[4] = {0, 1,
//				                3, 2};
				int dat_4[4] = {1, 0,
				                2, 3};
				iArray2DT temp(2, 2, dat_4);
				facetnodes = temp;
			}
			else if (fTotalNodes == 6)
			{
//				int dat_6[6] = {0, 1, 4,
//				                3, 2, 5};
				int dat_6[6] = {1, 0, 4,
				                2, 3, 5};
				iArray2DT temp(2, 3, dat_6);
				facetnodes = temp;
			}
			else
				num_nodes_error = 1;
			
			break;
		}
		case GeometryT::kTriangle:
		{
			if (fTotalNodes == 6)
			{
				int dat_6[6] = {0, 1, 2,
				                3, 4, 5};
				iArray2DT temp(2, 3, dat_6);
				facetnodes = temp;
			}
			else if (fTotalNodes == 12)
			{
				int dat_12[12] = {0, 1, 2, 6, 7, 8,
				                  3, 4, 5, 9, 10, 11};
				iArray2DT temp(2, 6, dat_12);
				facetnodes = temp;
			}
			else
				num_nodes_error = 1;
				
			break;
		}
		case GeometryT::kQuadrilateral:
		{
			if (fTotalNodes == 8)
			{
				int dat_8[8] = {0, 1, 2, 3,
				                4, 5, 6, 7};
				iArray2DT temp(2, 4, dat_8);
				facetnodes = temp;
			}
			else if (fTotalNodes == 16)
			{
				int dat_16[16] = {0, 1, 2, 3, 8, 9, 10, 11,
				                  4, 5, 6, 7, 12, 13, 14, 15};
				iArray2DT temp(2, 8, dat_16);
				facetnodes = temp;
			}
			else
				num_nodes_error = 1;
				
			break;
		}
		default:
		
			cout << "\n SurfaceShapeT::NodesOnFacets: unsupported geometry: ";
			cout << geometry << endl;
			throw ExceptionT::kGeneralFail;
	}
	
	if (num_nodes_error)
	{
		cout << "\n SurfaceShapeT::NodesOnFacets: " << fTotalNodes;
		cout << " nodes with geometry " << geometry << " is\n";
		cout <<   "      not supported" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* configure work space arrays */
void SurfaceShapeT::Construct(void)
{  	
	/* check dimensions */
	if (fFacetCoords.NumberOfNodes()*2 != fTotalNodes) throw ExceptionT::kSizeMismatch;

	/* shape functions */
	fNa.Allocate(fNumIP, fTotalNodes);

	/* jump shape functions */
	fjumpNa.Allocate(fNumIP, fTotalNodes);

	/* shape function and derivatives tables */
	fgrad_d.Allocate(fNumIP);
	fgrad_dTgrad_d.Allocate(fNumIP);
	fgrad_dd.Allocate(fNumIP, fFieldDim-1);

	int dim = fTotalNodes*fFieldDim;
	for (int i = 0; i < fNumIP; i++)
	{
		fgrad_d[i].Allocate(fFieldDim, dim);
		fgrad_dTgrad_d[i].Allocate(dim, dim);
		
		for (int j = 0; j < fFieldDim-1; j++)
			fgrad_dd(i,j).Allocate(fFieldDim, dim);
	}

	/* return value */
	fInterp.Allocate(fFieldDim);
	
	/* coordinate transformation */	
	int nsd = fFacetCoords.MinorDim();
	fJacobian.Allocate(nsd, nsd-1);
	
	/* surface node numbering */
	fFacetNodes.Allocate(2, fNumFacetNodes);
	
	/* work space */
	fu_vec.Allocate(fTotalNodes*fFieldDim);	
	
	/* 3D work space */
	if (fFieldDim == 3)
	{
		fM1.Allocate(fFieldDim, fu_vec.Length());
		fM2.Allocate(fFieldDim, fu_vec.Length());
	}
}

/* set jump vector, i.e., assignment of element nodes to
* facets: +1 => facet_2, -1 => facet_1. Displacement jump
* is: u_2 - u_1. */
void SurfaceShapeT::SetJumpVector(iArrayT& jump) const
{
	/* check */
	if (jump.Length() != fTotalNodes) throw ExceptionT::kSizeMismatch;

	int num_nodes_error = 0;
	int geometry = fDomain->GeometryCode();
	switch (geometry)
	{
		case GeometryT::kLine:
		{
			if (fTotalNodes == 4)
			{
				iArray2DT temp(2, 2, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
			}
			else if (fTotalNodes == 6)
			{
				iArray2DT temp(3, 2, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
				temp(2,0) =-1;
				temp(2,1) = 1;
			}
			else
				num_nodes_error = 1;
				
			break;
		}
		case GeometryT::kTriangle:
		{
			if (fTotalNodes == 6)
			{
				iArray2DT temp(2, 3, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
			}
			else if (fTotalNodes == 12)
			{
				iArray2DT temp(4, 3, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
				temp.SetRow(2,-1);
				temp.SetRow(3, 1);
			}
			else
				num_nodes_error = 1;
				
			break;
		}
		case GeometryT::kQuadrilateral:
		{
			if (fTotalNodes == 8)
			{
				iArray2DT temp(2, 4, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
			}
			else if (fTotalNodes == 16)
			{
				iArray2DT temp(4, 4, jump.Pointer());
				temp.SetRow(0,-1);
				temp.SetRow(1, 1);
				temp.SetRow(2,-1);
				temp.SetRow(3, 1);
			}
			else
				num_nodes_error = 1;
			
			break;
		}
		default:
		
			cout << "\n SurfaceShapeT::SetJumpVector: unsupported geometry: ";
			cout << geometry << endl;
			throw ExceptionT::kGeneralFail;
	}
	
	if (num_nodes_error)
	{
		cout << "\n SurfaceShapeT::SetJumpVector: " << fTotalNodes;
		cout << " nodes with geometry " << geometry << " is\n";
		cout <<   "      not supported" << endl;
		throw ExceptionT::kGeneralFail
		;
	}
}

/* compute average of the facet coordinates */
void SurfaceShapeT::ComputeFacetCoords(void)
{
#if __option(extended_errorcheck)
	if (fCoords.NumberOfNodes() != 2*fFacetCoords.NumberOfNodes())
		throw ExceptionT::kSizeMismatch;
#endif

	for (int i = 0; i < fFieldDim; i++)
	{
		int* facet1 = fFacetNodes(0);
		int* facet2 = fFacetNodes(1);
		
		double* px    = fCoords(i);
		double* pxmid = fFacetCoords(i);
		for (int j = 0; j < fNumFacetNodes; j++)
			*pxmid++ = 0.5*(px[*facet1++] + px[*facet2++]);
	}
}
