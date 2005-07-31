/* $Id: D2MeshFreeShapeFunctionT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#include "D2MeshFreeShapeFunctionT.h"
#include "D2MeshFreeSupport2DT.h"
//#include "MeshFreeSupport3DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

/* constructor */
D2MeshFreeShapeFunctionT::D2MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
	const LocalArrayT& coords, const dArray2DT& all_coords,
	const iArray2DT& connects, const iArrayT& nongridnodes,
	MeshFreeT::FormulationT code, double dextra, int complete, int store_shape,
	const int& currelement):
	MeshFreeShapeFunctionT(geometry_code, numIP, coords, all_coords, connects,
		nongridnodes, code, dextra, complete, store_shape, currelement),
	fDDNaU(numIP),
	fDDNa_tmp(numIP)
{
	if (all_coords.MinorDim() == 2)
		fD2MFSupport = new D2MeshFreeSupport2DT(*fDomain, all_coords, connects,
							nongridnodes, code, dextra, complete, store_shape);
	else
	{
		cout << "\n D2MeshFreeShapeFunctionT::D2MeshFreeShapeFunctionT: no 3D yet" << endl;
		throw eBadInputValue;
	}
	if (!fD2MFSupport) throw eOutOfMemory;

	/* delete MLS support for base class */
	delete fMFSupport;
	fMFSupport = fD2MFSupport;
	//TEMP - need to destruct the MLS support from the base class
}

/* compute local shape functions and derivatives */ 	
void D2MeshFreeShapeFunctionT::SetDerivatives(void)
{
	/* inherited (set geometry shape functions) */
	ShapeFunctionT::SetDerivatives();
	//TEMP - need to redesign things here. Skipping base
	//       class and going to base-base class

	/* load MLS field shape functions */
	fD2MFSupport->LoadElementData(fCurrElement, fNeighbors, fNaU,
		fDNaU, fDDNaU);
	//TEMP - if this weren't specific to the number of derivatives,
	//       wouldn't need to call the base class functions out-of-order

	/* blend for interpolant nodes */
	if (fExactNodes.Length() > 0) BlendElementData();
}

int D2MeshFreeShapeFunctionT::SetDerivativesAt(const dArrayT& x, AutoArrayT<int>& nodes)
{
	/* compute derivatives */
	if (fD2MFSupport->SetFieldAt(x, nodes))
	{
		const dArray2DT& Grad_x = fD2MFSupport->DFieldAt();
	
		/* keep neighbors data */
		fNeighbors.Alias(nodes);
		
		/* set next calls to GradU */
		SetGrad_x(Grad_x);

		return 1;
	}
	else
		return 0;
}

/* 2nd order shape function gradients matrix */
void D2MeshFreeShapeFunctionT::D2GradNa(dMatrixT& D2_grad_Na) const
{
	/* current integration point data */
	const dArray2DT& DDNa = fDDNaU[fCurrIP];
int numderiv = DDNa.MajorDim();
	int numnodes = DDNa.MinorDim();
	
	for (int i = 0; i < numderiv; i++)	
		for (int a = 0; a < numnodes; a++)	
			D2_grad_Na(i,a) = DDNa(i,a); // looks return of transpose?
}

/* reconstruct displacement field and all derivatives */
void D2MeshFreeShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
	dArray2DT& Dfield, dArray2DT& DDfield, iArrayT& nodes)
{
	/* fetch list of nodes to compute */
	nodes.Alias(fD2MFSupport->NodesUsed());
	
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = DOF.MinorDim();
	int nsd = NumSD();
	int nxx = dSymMatrixT::NumValues(nsd);
	
	/* allocate output space */
	field.Allocate(nnd, ndf);
	Dfield.Allocate(nnd, ndf*nsd);
	DDfield.Allocate(nnd, ndf*nxx);

	/* MLS nodal data */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;
	dArray2DT DDphi;

	/* "local" data for each node */
	int maxneighbors = (fD2MFSupport->NodeNeighbors()).MaxMinorDim();
	dArrayT space(maxneighbors*ndf);
	LocalArrayT locdisp(LocalArrayT::kDisp);
	locdisp.Set(0, ndf, NULL); // must have minor dim to set global
	locdisp.SetGlobal(DOF);
	dArrayT dof;

	/* loop over nodes in the set */
	dMatrixT Du, DDu;
	dArrayT tmp;
	for (int i = 0; i < nnd; i++)
	{
		int node = nodes[i];

		/* fetch MLS data */
		fD2MFSupport->LoadNodalData(node, neighbors, phi, Dphi, DDphi);
		
		/* fetch neighbor data */
		int len = neighbors.Length();
		locdisp.Set(len, ndf, space.Pointer());
		locdisp.SetLocal(neighbors);

		/* compute nodal values */
		Du.Set(ndf, nsd, Dfield(i));
		DDu.Set(ndf, nxx, DDfield(i));
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
		}
	}
}
