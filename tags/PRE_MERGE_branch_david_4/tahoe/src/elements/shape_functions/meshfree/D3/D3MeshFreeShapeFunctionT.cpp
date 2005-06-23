/* $Id: D3MeshFreeShapeFunctionT.cpp,v 1.5 2005-04-06 01:03:55 kyonten Exp $ */
/* created: paklein (10/23/1999) */
#include "D3MeshFreeShapeFunctionT.h"
#include "D3MeshFreeSupport2DT.h"
//#include "MeshFreeSupport3DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
D3MeshFreeShapeFunctionT::D3MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
	const LocalArrayT& coords, const dArray2DT& all_coords,
	const iArray2DT& connects, const iArrayT& nongridnodes,
	const int& currelement, const ParameterListT& mf_support_params):
	D2MeshFreeShapeFunctionT(geometry_code, numIP, coords, all_coords, connects,
		nongridnodes, currelement, mf_support_params),
	fDDDNaU(numIP),
	fDDDNa_tmp(numIP),
	fD3MFSupport(NULL)
{
	const char caller[] = "D3MeshFreeShapeFunctionT::D3MeshFreeShapeFunctionT";

	/* construct MLS support - otherwise must be set by constructor of derived classes */
	if (mf_support_params.Name() == "D3_meshfree_support_2D")
		fD3MFSupport = new D3MeshFreeSupport2DT(fDomain, all_coords, connects, nongridnodes);

	/* initialize */
	if (fD3MFSupport) fD3MFSupport->TakeParameterList(mf_support_params);

	/* delete MLS support for base class */
	delete fMFSupport;
	fMFSupport = fD2MFSupport = fD3MFSupport;
	
	/* set as field shape function */
	//SetUShapeFunctions(fNaU, fDNaU, fDDNaU, fDDDNaU);
}

/* class-dependent initializations */
void D3MeshFreeShapeFunctionT::Initialize(void)
{
	/* inherited */
	//D2MeshFreeShapeFunctionT::Initialize();
	ShapeFunctionT::Initialize();
	
	/* check */
	if (!fD3MFSupport)
		ExceptionT::GeneralFail("D3MeshFreeShapeFunctionT::Initialize",
			"meshfree support not set");
}

/* compute local shape functions and derivatives */ 	
void D3MeshFreeShapeFunctionT::SetDerivatives(void)
{
	/* inherited (set geometry shape functions) */
	//D2MeshFreeShapeFunctionT::SetDerivatives();
	ShapeFunctionT::SetDerivatives();
	//TEMP - need to redesign things here. Skipping base
	//       class and going to base-base class

	/* load MLS field shape functions */
	fD3MFSupport->LoadElementData(fCurrElement, fNeighbors, fNaU,
		fDNaU, fDDNaU, fDDDNaU);
	//TEMP - if this weren't specific to the number of derivatives,
	//       wouldn't need to call the base class functions out-of-order

	/* blend for interpolant nodes */
	if (fExactNodes.Length() > 0) BlendElementData();
}

int D3MeshFreeShapeFunctionT::SetDerivativesAt(const dArrayT& x, AutoArrayT<int>& nodes)
{
	/* compute derivatives */
	if (fD3MFSupport->SetFieldAt(x))
	{
		const dArray2DT& Grad_x = fD3MFSupport->DFieldAt();
	
		/* copy nodal neighor data */
		fNeighbors.Alias(fMFSupport->NeighborsAt());
		nodes.Dimension(fNeighbors.Length());
		nodes = fNeighbors;
		
		/* set next calls to GradU */
		SetGrad_x(Grad_x);
		return 1;
	}
	else
		return 0;
}

/* 3rd order shape function gradients matrix */
void D3MeshFreeShapeFunctionT::D3GradNa(dMatrixT& D3_grad_Na) const
{
	/* current integration point data */
	const dArray2DT& DDDNa = fDDDNaU[fCurrIP];
	int numderiv = DDDNa.MajorDim();
	int numnodes = DDDNa.MinorDim();
	
	for (int i = 0; i < numderiv; i++)	
		for (int a = 0; a < numnodes; a++)	
			D3_grad_Na(i,a) = DDDNa(i,a); // looks return of transpose?
}

/* reconstruct displacement field and all derivatives */
void D3MeshFreeShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
	dArray2DT& Dfield, dArray2DT& DDfield, dArray2DT& DDDfield, iArrayT& nodes)
{
	/* fetch list of nodes to compute */
	nodes.Alias(fD3MFSupport->NodesUsed());
	
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = DOF.MinorDim();
	int nsd = NumSD();
	int nxx = dSymMatrixT::NumValues(nsd);
	
	/* allocate output space */
	field.Dimension(nnd, ndf);
	Dfield.Dimension(nnd, ndf*nsd);
	DDfield.Dimension(nnd, ndf*nxx);
	DDDfield.Dimension(nnd, ndf*(nsd*nsd));  

	/* MLS nodal data */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;
	dArray2DT DDphi;
	dArray2DT DDDphi;

	/* "local" data for each node */
	int maxneighbors = (fD3MFSupport->NodeNeighbors()).MaxMinorDim();
	dArrayT space(maxneighbors*ndf);
	LocalArrayT locdisp(LocalArrayT::kDisp);
	locdisp.Set(0, ndf, NULL); // must have minor dim to set global
	locdisp.SetGlobal(DOF);
	dArrayT dof;

	/* loop over nodes in the set */
	dMatrixT Du, DDu, DDDu;
	dArrayT tmp;
	for (int i = 0; i < nnd; i++)
	{
		int node = nodes[i];

		/* fetch MLS data */
		fD3MFSupport->LoadNodalData(node, neighbors, phi, Dphi, DDphi, DDDphi);
		
		/* fetch neighbor data */
		int len = neighbors.Length();
		locdisp.Set(len, ndf, space.Pointer());
		locdisp.SetLocal(neighbors);

		/* compute nodal values */
		Du.Set(ndf, nsd, Dfield(i));
		DDu.Set(ndf, nxx, DDfield(i));
		DDDu.Set(ndf, nsd*nsd, DDDfield(i)); 
		for (int j = 0; j < ndf; j++)
		{
			dof.Set(len, locdisp(j));

			/* displacement field */
			field(i, j) = dArrayT::Dot(dof, phi);
			
			/* first derivatives */
			for (int k = 0; k < nsd; k++)
			{
				Dphi.RowAlias(k, tmp);
				Du(j, k) = dArrayT::Dot(dof, tmp);
			}

			/* second derivatives */
			for (int l = 0; l < nxx; l++)
			{
				DDphi.RowAlias(l, tmp);
				DDu(j, l) = dArrayT::Dot(dof, tmp);
			}
			
			/* third derivatives */
			for (int l = 0; l < nsd*nsd; l++)  
			{
				DDDphi.RowAlias(l, tmp);
				DDDu(j, l) = dArrayT::Dot(dof, tmp);
			}
		}
	}
}
