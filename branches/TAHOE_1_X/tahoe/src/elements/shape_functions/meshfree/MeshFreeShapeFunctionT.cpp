/* $Id: MeshFreeShapeFunctionT.cpp,v 1.15 2004-06-26 18:39:57 paklein Exp $ */
/* created: paklein (09/10/1998) */
#include "MeshFreeShapeFunctionT.h"

#include "toolboxConstants.h"
#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* constructor */
MeshFreeShapeFunctionT::MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
	const LocalArrayT& coords, const dArray2DT& all_coords,
	const iArray2DT& connects, const iArrayT& nongridnodes,
	const int& currelement, ifstreamT& in):
	ShapeFunctionT(geometry_code, numIP, coords),
	fMFSupport(NULL),
	fCurrElement(currelement),
	fDNaU(numIP),
	fXConnects(connects)
{
	/* construct MLS support */
	if (all_coords.MinorDim() == 2)
		fMFSupport = new MeshFreeSupport2DT(fDomain, all_coords, connects,
							nongridnodes, in);
	else
		fMFSupport = new MeshFreeSupport3DT(fDomain, all_coords, connects,
							nongridnodes, in);

	if (!fMFSupport) throw ExceptionT::kOutOfMemory;

	/* set as field shape function */
	SetUShapeFunctions(fNaU, fDNaU);
}

/* destructor */
MeshFreeShapeFunctionT::~MeshFreeShapeFunctionT(void) { delete fMFSupport; }

/* initialization - modifications to the support size must
* occur before setting the neighbor data. Coordinates and
* connecitivies must be set */
void MeshFreeShapeFunctionT::SetSupportSize(void)
{
	/* initialize MLS data */
	fMFSupport->InitSupportParameters();
}

void MeshFreeShapeFunctionT::SetNeighborData(void)
{
	/* initialize MLS data */
	fMFSupport->InitNeighborData();
}

void MeshFreeShapeFunctionT::SetExactNodes(const iArrayT& exact_nodes)
{
	/* keep copy of the node list */
	fExactNodes = exact_nodes;
	
	/* initialize shape function blending data */
	if (fExactNodes.Length() > 0) InitBlend();
}

/* specify nodes/cells to skip when doing MLS calculations */
void MeshFreeShapeFunctionT::SetSkipNodes(const iArrayT& skip_nodes)
{
	fMFSupport->SetSkipNodes(skip_nodes);
}

void MeshFreeShapeFunctionT::SetSkipElements(const iArrayT& skip_elements)
{
	fMFSupport->SetSkipElements(skip_elements);
}

/* compute local shape functions and derivatives */ 	
void MeshFreeShapeFunctionT::SetDerivatives(void)
{
	/* inherited (set geometry shape functions) */
	ShapeFunctionT::SetDerivatives();

	/* load MLS field shape functions */
	fMFSupport->LoadElementData(fCurrElement, fNeighbors, fNaU, fDNaU);

//TEMP - see if we're using only meshfree nodes
#if 0
ofstream MLS_out("MLS.out", ios::app);
MLS_out << "element: " << fCurrElement + 1 << '\n';
MLS_out << "  count: " << fNeighbors.Length() << '\n';
for (int i = 0; i < fNeighbors.Length(); i++)
{
	MLS_out << setw(kIntWidth) << i + 1;
	MLS_out << setw(kIntWidth) << fNeighbors[i] + 1;
	for (int j = 0; j < NumIP(); j++)
		MLS_out << setw(kDoubleWidth) << fNaU(j,i);
	MLS_out << '\n';
}
MLS_out << endl;
#endif
//END

	/* blend for interpolant nodes */
	if (fExactNodes.Length() > 0) BlendElementData();
}

int MeshFreeShapeFunctionT::SetDerivativesAt(const dArrayT& x)
{
	/* compute derivatives */
	if (fMFSupport->SetFieldAt(x))
	{
		const dArray2DT& Grad_x = fMFSupport->DFieldAt();
	
		/* keep neighbors data */
		fNeighbors.Alias(fMFSupport->NeighborsAt());
		
		/* set next calls to GradU */
		SetGrad_x(Grad_x);
		
		return 1;
	}
	else
		return 0;
}

void MeshFreeShapeFunctionT::UseDerivatives(const iArrayT& neighbors,
	const dArray2DT& Dfield) // load external values
{
	/*set neighbors */
	fNeighbors.Alias(neighbors);

	/* set next calls to GradU */
	SetGrad_x(Dfield);
}

/* cutting facet functions */
void MeshFreeShapeFunctionT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	fMFSupport->SetCuttingFacets(facet_coords, num_facet_nodes);
}

void MeshFreeShapeFunctionT::ResetFacets(const ArrayT<int>& facets)
{
	fMFSupport->ResetFacets(facets);
}

const ArrayT<int>& MeshFreeShapeFunctionT::ResetNodes(void) const
{
	return fMFSupport->ResetNodes();
}

const ArrayT<int>& MeshFreeShapeFunctionT::ResetCells(void) const
{
	return fMFSupport->ResetCells();
}

/* access to MLS field neighbor data */
const iArrayT& MeshFreeShapeFunctionT::ElementNeighborsCounts(void) const
{
	return fMFSupport->ElementNeighborsCounts();
}

const RaggedArray2DT<int>& MeshFreeShapeFunctionT::ElementNeighbors(void) const
{
	return fMFSupport->ElementNeighbors();
}

const RaggedArray2DT<int>& MeshFreeShapeFunctionT::NodeNeighbors(void) const
{
	return fMFSupport->NodeNeighbors();
}

/* reconstruct displacement field */
void MeshFreeShapeFunctionT::SelectedNodalField(const dArray2DT& all_DOF,
	const iArrayT& nodes, dArray2DT& field)
{
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = all_DOF.MinorDim();
	
	/* allocate output space */
	field.Dimension(nnd, ndf);

	/* MLS nodal data */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;

	/* "local" data for each node */
	int maxneighbors = (fMFSupport->NodeNeighbors()).MaxMinorDim();
	dArrayT space(maxneighbors*ndf);
	LocalArrayT locdisp(LocalArrayT::kDisp);
	locdisp.Set(0, ndf, NULL); // must have minor dim to set global
	locdisp.SetGlobal(all_DOF);
	dArrayT dof;

	/* loop over nodes in the set */
	for (int i = 0; i < nnd; i++)
	{
		int node = nodes[i];

		/* fetch MLS data */
		fMFSupport->LoadNodalData(node, neighbors, phi, Dphi);
		
		/* blend for interpolant nodes */
		if (fExactNodes.Length() > 0) BlendNodalData(node, neighbors, phi);

		/* fetch neighbor data */
		int len = neighbors.Length();
		locdisp.Set(len, ndf, space.Pointer());
		locdisp.SetLocal(neighbors);

		/* compute displacements */
		for (int j = 0; j < ndf; j++)
		{
			dof.Set(len,locdisp(j));
			field(i,j) = dArrayT::Dot(dof, phi);	
		}
	}
}

void MeshFreeShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
	iArrayT& nodes)
{
	/* fetch list of nodes to compute */
	nodes.Alias(fMFSupport->NodesUsed());
	
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = DOF.MinorDim();
	
	/* allocate output space */
	field.Dimension(nnd, ndf);

	/* MLS nodal data */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;

	/* "local" data for each node */
	int maxneighbors = (fMFSupport->NodeNeighbors()).MaxMinorDim();
	dArrayT space(maxneighbors*ndf);
	LocalArrayT locdisp(LocalArrayT::kDisp);
	locdisp.Set(0, ndf, NULL); // must have minor dim to set global
	locdisp.SetGlobal(DOF);
	dArrayT dof;

	/* loop over nodes in the set */
	for (int i = 0; i < nnd; i++)
	{
		int node = nodes[i];

		/* fetch MLS data */
		fMFSupport->LoadNodalData(node, neighbors, phi, Dphi);
		
		/* blend for interpolant nodes */
		if (fExactNodes.Length() > 0) BlendNodalData(node, neighbors, phi);

		/* fetch neighbor data */
		int len = neighbors.Length();
		locdisp.Set(len, ndf, space.Pointer());
		locdisp.SetLocal(neighbors);

		/* compute displacements */
		for (int j = 0; j < ndf; j++)
		{
			dof.Set(len,locdisp(j));
			field(i,j) = dArrayT::Dot(dof, phi);	
		}
	}
}

void MeshFreeShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
	dArray2DT& Dfield, iArrayT& nodes)
{
	/* interpolant nodes produce discontinuous derivatives */
	if (fExactNodes.Length() > 0)
	{
		cout << "\n MeshFreeShapeFunctionT::NodalField: derivatives not continuous with\n"
		     <<   "     interpolant field nodes" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* fetch list of nodes to compute */
	nodes.Alias(fMFSupport->NodesUsed());
	
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = DOF.MinorDim();
	int nsd = NumSD();
	
	/* allocate output space */
	field.Dimension(nnd, ndf);
	Dfield.Dimension(nnd, ndf*nsd);

	/* MLS nodal data */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;

	/* "local" data for each node */
	int maxneighbors = (fMFSupport->NodeNeighbors()).MaxMinorDim();
	dArrayT space(maxneighbors*ndf);
	LocalArrayT locdisp(LocalArrayT::kDisp);
	locdisp.Set(0, ndf, NULL); // must have minor dim to set global
	locdisp.SetGlobal(DOF);
	dArrayT dof;

	/* loop over nodes in the set */
	dMatrixT Du;
	dArrayT tmp;
	for (int i = 0; i < nnd; i++)
	{
		int node = nodes[i];

		/* fetch MLS data */
		fMFSupport->LoadNodalData(node, neighbors, phi, Dphi);
		
		/* fetch neighbor data */
		int len = neighbors.Length();
		locdisp.Set(len, ndf, space.Pointer());
		locdisp.SetLocal(neighbors);

		/* compute nodal values */
		Du.Set(ndf, nsd, Dfield(i));
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
		}
	}
}

/* print the shape function values to the output stream */
void MeshFreeShapeFunctionT::Print(ostream& out) const
{
/* inherited */
	ShapeFunctionT::Print(out);

	out << "\n Field shape functions:\n";
	for (int i = 0; i < fNeighbors.Length(); i++)
	{
		out << setw(kIntWidth) << fNeighbors[i];
		
		for (int j = 0; j < NumIP(); j++)
			out << setw(kDoubleWidth) << fNaU(j,i);
		out << '\n';
	}

	out << "\n Field shape function derivatives:\n";
//	ArrayT<dArray2DT> fDNaU;

// print out associated node numbers with the values and
// derivatives
//	if (fDNaX.Pointer() == fDNaU.Pointer())
//	    out << " isoparametric \n";
//	else	
//	    for (int i = 0; i < fDNaU.Length(); i++)
//		    out << fDNaU[i] << '\n';
}

void MeshFreeShapeFunctionT::PrintAt(ostream& out) const
{
	/* shape functions array */
	const dArrayT& phi = fMFSupport->FieldAt();

	out << "\n Field shape functions:\n";
	for (int i = 0; i < fNeighbors.Length(); i++)
	{
		out << setw(kIntWidth)    << fNeighbors[i];
		out << setw(kDoubleWidth) << phi[i];
		out << '\n';
	}
}

/* write MLS statistics */
void MeshFreeShapeFunctionT::WriteParameters(ostream& out) const
{
	fMFSupport->WriteParameters(out);
}

void MeshFreeShapeFunctionT::WriteStatistics(ostream& out) const
{
	fMFSupport->WriteStatistics(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* initialize blending database */
void MeshFreeShapeFunctionT::InitBlend(void)
{
	/* initialize element flags */
	fElemHasExactNode.Dimension(fXConnects.MajorDim());
	fElemHasExactNode = -1;

	/* map to exact nodes */
	int shift, max;
	fExactNodes.MinMax(shift, max);
	int range = max - shift + 1;
	iArrayT exact_map(range);
	exact_map = 0;
	for (int jj = 0; jj < fExactNodes.Length(); jj++)
		exact_map[fExactNodes[jj] - shift] = 1;
		
	/* mark elements and check */
	int nel = fXConnects.MajorDim();
	int nen = fXConnects.MinorDim();
	iArrayT hit_map(range);
	hit_map = 0;
	for (int j = 0; j < nel; j++)
	{
		const int* pelem = fXConnects(j);
		for (int k = 0; k < nen; k++)
		{
			int shifted_node = *pelem++ - shift;
			if (shifted_node > -1 &&
			    shifted_node < range &&
			    exact_map[shifted_node] == 1)
			{
				fElemHasExactNode[j] = 1;
				hit_map[shifted_node] = 1;
			}
		}
	}

	/* off-grid nodes */
	bool first_hit = true;
	for (int k = 0; k < hit_map.Length(); k++)
		if (exact_map[k] == 1 && hit_map[k] == 0)
		{
			if (first_hit)
			{
				cout << "\n MeshFreeShapeFunctionT::InitBlend: interpolant nodes should lie\n"
				     << "     on the integration grid. Exceptions may exist for parallel\n"
				     << "     execution. Suspect nodes:\n"
				     << setw(2*kIntWidth) << "node (loc)"
				     << setw(kIntWidth) << "is skip" << '\n';
					first_hit = false;
			}
		
			int node = k + shift;
			const iArrayT& skip_nodes = fMFSupport->SkipNodes();
			cout << setw(2*kIntWidth) << node + 1
			     << setw(kIntWidth) << skip_nodes.HasValue(node) << '\n';		
		}
			
	/* free work space */
	exact_map.Free();
	hit_map.Free();

	/* initialize element data */
	int num_flagged = fElemHasExactNode.Count(1);
	fElemFlags.Dimension(num_flagged, nen);
	fElemFlags = 0;
	
	/* create element flagged data and make map */
	int flagged_count = 0;
	for (int ii = 0; ii < fElemHasExactNode.Length(); ii++)
		if (fElemHasExactNode[ii] == 1)
		{
			/* convert to map of element data */
			fElemHasExactNode[ii] = flagged_count;

			const int* pelem = fXConnects(ii);
			int* pflags = fElemFlags(flagged_count);
			for (int j = 0; j < nen; j++)
			{
				if (fExactNodes.HasValue(*pelem)) *pflags = 1;
				pelem++; pflags++;
			}

			flagged_count++;
		}

//TEMP - write flag data
#if 0
for (int i = 0; i < fElemHasExactNode.Length(); i++)
{
	cout << i+1 << '\t' << fElemHasExactNode[i] << '\n';
	if (fElemHasExactNode[i] > -1)
		fElemFlags.PrintRow(fElemHasExactNode[i], cout);
}	
throw ExceptionT::kStop;
#endif
//END
	
	/* allocate ramp function */
	fR.Dimension(NumIP());
	fDR.Dimension(NumIP(), NumSD());
		
	/* allocate work space for blended shape functions */
	fDNa_tmp.Dimension(NumIP());
	int el_maxsize = (fMFSupport->ElementNeighbors()).MaxMinorDim();
	int nd_maxsize = (fMFSupport->NodeNeighbors()   ).MaxMinorDim();
	felSpace.Dimension(NumIP()*el_maxsize*(1 + NumSD()));
	fndSpace.Dimension(nd_maxsize*(1 + NumSD()));
}

/* blend FE/MLS shape functions for interpolant nodes */
void MeshFreeShapeFunctionT::BlendElementData(void)
{
	/* need to remove any contribution from interpolant nodes */
	if (fElemHasExactNode[fCurrElement] == -1) return;
	
	/* construct ramp function */
	const dArray2DT& rNa = Na();
	const ArrayT<dArray2DT>& rDNa = DNaX();
	
	int nip = NumIP();
	int nsd = NumSD();
	int nen = rNa.MinorDim();
	int nnd = fNeighbors.Length();
	
	fR  = 0.0;
	fDR = 0;			

	int* pelem_flags = fElemFlags(fElemHasExactNode[fCurrElement]);
	for (int ii = 0; ii < nip; ii++)
	{
		/* ramp function */
		int* pflag = pelem_flags;
		const double* pNa = rNa(ii);
		for (int j = 0; j < nen; j++)
			fR[ii] += (1 - *pflag++)*(*pNa++);
		
		/* ramp function derivatives */
		for (int k = 0; k < nsd; k++)
		{
			int*   pflag = pelem_flags;
			double* pDR  = fDR(ii);
			const double* pDNa = (rDNa[ii])(k);
			for (int j = 0; j < nen; j++)
				pDR[k] += (1 - *pflag++)*(*pDNa++);
		}
	}

//TEMP
#if 0
	MLS_out << " element: " << fCurrElement+1 << '\n';
	MLS_out << fR  << '\n';
	MLS_out << fDR << '\n';
#endif
//END

	/* initialize/configure temp space */
	felSpace = 0.0;
	double* pelspace = felSpace.Pointer();
	fNa_tmp.Set(nip, nnd, pelspace);
	pelspace += fNa_tmp.Length();
	for (int j = 0; j < nip; j++)
	{
		fDNa_tmp[j].Set(nsd, nnd, pelspace);
		pelspace += fDNa_tmp[j].Length();
	}	
	
	/* identify the cell nodes in the neighbors list */
	fNeighExactFlags.Dimension(nnd);
	for (int k = 0; k < nnd; k++)
	{
		int loc = -1;
		const int* pelem = fXConnects(fCurrElement);
		int listnode = fNeighbors[k];
		for (int i = 0; i < nen && loc < 0; i++)
			if (*pelem++ == listnode) loc = i;
			
		fNeighExactFlags[k] = loc; // -1 if not grid node, else local number
	}

//TEMP
#if 0
	iArrayT n_tmp(fNeighExactFlags.Length(), fNeighExactFlags.Pointer());
	MLS_out << " element: " << fCurrElement+1 << '\n';
	fNeighbors++;
	MLS_out << fNeighbors.wrap(10) << '\n';
	fNeighbors--;
	MLS_out << n_tmp.wrap(10);
#endif
//END
						
	/* compute blended shape functions */
	for (int i = 0; i < nip; i++)
	{
		/* shape function */
		int*   pflag = fNeighExactFlags.Pointer();
		double*  pNa = fNa_tmp(i);
		double* pPhi = fNaU(i);
		double     R = fR[i];
		for (int j = 0; j < nnd; j++)
		{
			if (*pflag > -1)
				*pNa++ = (1.0 - R)*rNa(i,*pflag) + R*(*pPhi++);
			else
				*pNa++ = R*(*pPhi++);
		
			pflag++;
		}
				
		/* shape function derivatives */
		for (int s = 0; s < nsd; s++)
		{
			/* local number of interpolant node in current grid */
			int*    pflag = fNeighExactFlags.Pointer();
		
			double*  pPhi = fNaU(i);
			double*  pDNa = (fDNa_tmp[i])(s);
			double* pDPhi = (fDNaU[i])(s);
			double*   pDR = fDR(i);
			for (int j = 0; j < nnd; j++)
			{
				if (*pflag > -1)
					*pDNa++ = pDR[s]*((*pPhi++) - rNa(i,*pflag)) +
					          (1.0 - R)*(rDNa[i])(s,*pflag) +
					          R*(*pDPhi++);
				else
					*pDNa++ = pDR[s]*(*pPhi++) + R*(*pDPhi++);
			
				pflag++;
			}
		}
	}
	
	/* reset pointers */
	fNaU.Alias(fNa_tmp);
	for (int l = 0; l < nip; l++)
		fDNaU[l].Alias(fDNa_tmp[l]);
}

void MeshFreeShapeFunctionT::BlendNodalData(int node, const iArrayT& nodes, dArrayT& phi)
{
// NOTES:
// (1) only the ramp function is continuous across element boundaries
// (2) ramp function derivatives are discontinuous
// (3) implies cannot evaluate ramp function derivatives, or the shape
//     function derivatives at the element boundaries.	
// (4) function also assumes that no "extra" points have been included
//     in the integration cells containing exact nodes

	/* check if at an interpolant node */
	int dex;
	if (fExactNodes.HasValue(node, dex))
	{
		/* temp space */
		dArrayT phi_temp(phi.Length(), fndSpace.Pointer());

		/* no contribution from other nodes */
		phi_temp = 0.0;
		
		/* interpolant at the current node */
		int loc_dex;
		nodes.HasValue(node, loc_dex);
		phi_temp[loc_dex] = 1.0;

		/* reset pointer */
		phi.Alias(phi_temp);
	}
}

/* read/write nodal meshfree parameters */
void MeshFreeShapeFunctionT::SetNodalParameters(const iArrayT& node, const dArray2DT& nodal_params)
{
	fMFSupport->SetSupportParameters(node, nodal_params);
}

void MeshFreeShapeFunctionT::GetNodalParameters(const iArrayT& node, dArray2DT& nodal_params) const
{
	fMFSupport->GetSupportParameters(node, nodal_params);
}

dArray2DT& MeshFreeShapeFunctionT::NodalParameters(void) 
{ 
	return fMFSupport->NodalParameters(); 
}
