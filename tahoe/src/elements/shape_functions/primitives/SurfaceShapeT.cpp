/* $Id: SurfaceShapeT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (11/21/1997)                                          */
/* Class to manage CSE integrals, where the dimension of                  */
/* the field variable is 1 greater than the dimension of the parent       */
/* domain. Jump quantities imply jump between any field variable          */
/* across the CSE. Revised for new shape function object model            */
/* PAK (09/04/98)                                                         */
/* NOTE:                                                                  */
/* numnodes = total number of element nodes                               */
/* coords = coordinates of nodes on 1st facet (numnodes/2)                */

#include "SurfaceShapeT.h"

#include "Constants.h"
#include "ExceptionCodes.h"

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
	dArrayT	 shNafacets;
	dMatrixT shNaMat;
	for (int i = 0; i < fNumIP; i++)
	{
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
				jumpNa[dex] = sign[dex]*(*shape++);
			}
		}

		/* shape function tables */
		shNaMat.Set(1, fTotalNodes, fjumpNa(i));
		fgrad_d[i] = 0.0;
		fgrad_d[i].Expand(shNaMat, fFieldDim);
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
			fgrad_dd(i,j) = 0.0;
			fgrad_dd(i,j).Expand(shNaMat, fFieldDim);
		}
	}	
}

/**** for the current integration point ***/

/* jump in the nodal values */
const dArrayT& SurfaceShapeT::InterpolateJumpU(const LocalArrayT& nodalU)
{
	for (int i = 0; i < fFieldDim; i++)
		fInterp[i] = fjumpNa.DotRow(fCurrIP, nodalU(i));
	return fInterp;
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
	if (dQ.Length() != fFieldDim) throw eSizeMismatch;
#endif

	/* compute facet coordinates */
	if (fCurrIP == 0 && fCoords.NumberOfNodes() != fNumFacetNodes)
		ComputeFacetCoords();

	/* get Jacobian matrix of the surface transformation */
	fDomain->DomainJacobian(fFacetCoords, fCurrIP, fJacobian);	
	double j = fDomain->SurfaceJacobian(fJacobian, Q);
	if (j < kSmall) throw eBadJacobianDet;

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
		if (m1 < kSmall) throw eBadJacobianDet;

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

/*******************************************/

/* local node numbers on each facet */
void SurfaceShapeT::SetNodesOnFacets(iArray2DT& facetnodes)
{
	/* check */
	if (facetnodes.MajorDim() != 2 &&
	    facetnodes.MinorDim() != fNumFacetNodes)  throw eSizeMismatch;

	int num_nodes_error = 0;
	int geometry = fDomain->GeometryCode();
	switch (geometry)
	{
		case GeometryT::kLine:
		{
			if (fTotalNodes == 4)
			{
				int dat_4[4] = {0, 1,
				                3, 2};
				iArray2DT temp(2, 2, dat_4);
				facetnodes = temp;
			}
			else if (fTotalNodes == 6)
			{
				int dat_6[6] = {0, 1, 4,
				                3, 2, 5};
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
			throw eGeneralFail;
	}
	
	if (num_nodes_error)
	{
		cout << "\n SurfaceShapeT::NodesOnFacets: " << fTotalNodes;
		cout << " nodes with geometry " << geometry << " is\n";
		cout <<   "      not supported" << endl;
		throw eGeneralFail;
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* configure work space arrays */
void SurfaceShapeT::Construct(void)
{  	
	/* check dimensions */
	if (fFacetCoords.NumberOfNodes()*2 != fTotalNodes) throw eSizeMismatch;

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
	fJacobian.Allocate(fFieldDim, fFieldDim-1);
	
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
	if (jump.Length() != fTotalNodes) throw eSizeMismatch;

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
			throw eGeneralFail;
	}
	
	if (num_nodes_error)
	{
		cout << "\n SurfaceShapeT::SetJumpVector: " << fTotalNodes;
		cout << " nodes with geometry " << geometry << " is\n";
		cout <<   "      not supported" << endl;
		throw eGeneralFail;
	}
}

/* compute average of the facet coordinates */
void SurfaceShapeT::ComputeFacetCoords(void)
{
#if __option(extended_errorcheck)
	if (fCoords.NumberOfNodes() != 2*fFacetCoords.NumberOfNodes())
		throw eSizeMismatch;
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
