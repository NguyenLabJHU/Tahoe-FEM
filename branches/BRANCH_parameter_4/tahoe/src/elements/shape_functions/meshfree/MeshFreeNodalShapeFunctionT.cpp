/* $Id: MeshFreeNodalShapeFunctionT.cpp,v 1.3 2004-06-24 21:03:50 cjkimme Exp $ */
#include "MeshFreeNodalShapeFunctionT.h"
#include "toolboxConstants.h"
#include "MeshFreeSupport2DT.h"
#include "MeshFreeSupport3DT.h"
#include "LocalArrayT.h"
#include "dMatrixT.h"

using namespace Tahoe;

/* constructor */
MeshFreeNodalShapeFunctionT::MeshFreeNodalShapeFunctionT(int numSD, 
	const dArray2DT& all_coords, const iArray2DT& connects,
	const dArray2DT& nonNodes, ifstreamT& in):
	fSD(numSD),
	fMFSupport(NULL),
	fDNaU(1),
	fNonGridNodes()
{
#pragma unused(nonNodes)
	/* store auxiliary points for shape function evaluation */

	/* construct MLS support */
	if (all_coords.MinorDim() == 2)
		fMFSupport = new MeshFreeSupport2DT(NULL, all_coords, connects,
							fNonGridNodes, in);
	else
		fMFSupport = new MeshFreeSupport3DT(NULL, all_coords, connects,
							fNonGridNodes, in);

	if (!fMFSupport) 
		ExceptionT::OutOfMemory("MeshFreeNodalShapeFunctionT::MeshFreeNodalShapeFunctionT","Cannot create support\n");

}

/* destructor */
MeshFreeNodalShapeFunctionT::~MeshFreeNodalShapeFunctionT(void) { delete fMFSupport; }

/* initialization - modifications to the support size must
* occur before setting the neighbor data. Coordinates and
* connecitivies must be set */
void MeshFreeNodalShapeFunctionT::SetSupportSize(void)
{
	/* initialize MLS data */
	fMFSupport->InitSupportParameters();
}

void MeshFreeNodalShapeFunctionT::SetNeighborData(void)
{
	/* initialize MLS data */
	fMFSupport->InitNeighborData();
}

/* specify nodes/cells to skip when doing MLS calculations */
void MeshFreeNodalShapeFunctionT::SetSkipNodes(const iArrayT& skip_nodes)
{
	fMFSupport->SetSkipNodes(skip_nodes);
}

/* compute shape function at arbitrary point */
int MeshFreeNodalShapeFunctionT::SetFieldAt(const dArrayT& x, const dArrayT* shift)
{
	/* compute derivatives */
	if (fMFSupport->SetFieldAt(x, shift))
	{
		/* keep neighbors data */
		fNeighbors.Alias(fMFSupport->NeighborsAt());
		
		return 1;
	}
	else
		return 0;
}

const dArrayT& MeshFreeNodalShapeFunctionT::FieldAt()
{
	return fMFSupport->FieldAt();
}

/* compute local shape functions and derivatives */ 	
void MeshFreeNodalShapeFunctionT::SetDerivatives(void)
{

	/* load MLS field shape functions */
//	fMFSupport->LoadElementData(fCurrElement, fNeighbors, fNaU, fDNaU);

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

}

int MeshFreeNodalShapeFunctionT::SetDerivativesAt(const dArrayT& x)
{
	/* compute derivatives */
	if (fMFSupport->SetFieldAt(x))
	{
		const dArray2DT& Grad_x = fMFSupport->DFieldAt();
	
		/* keep neighbors data */
		fNeighbors.Alias(fMFSupport->NeighborsAt());
		
		/* set next calls to GradU */
//		SetGrad_x(Grad_x);
		
		return 1;
	}
	else
		return 0;
}

const dArray2DT& MeshFreeNodalShapeFunctionT::DFieldAt(void)
{
	return fMFSupport->DFieldAt();
}

void MeshFreeNodalShapeFunctionT::UseDerivatives(const iArrayT& neighbors,
	const dArray2DT& Dfield) // load external values
{
#pragma unused(Dfield)
	/*set neighbors */
	fNeighbors.Alias(neighbors);

	/* set next calls to GradU */
//	SetGrad_x(Dfield);
}

/* cutting facet functions */
void MeshFreeNodalShapeFunctionT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	fMFSupport->SetCuttingFacets(facet_coords, num_facet_nodes);
}

void MeshFreeNodalShapeFunctionT::ResetFacets(const ArrayT<int>& facets)
{
	fMFSupport->ResetFacets(facets);
}

const ArrayT<int>& MeshFreeNodalShapeFunctionT::ResetNodes(void) const
{
	return fMFSupport->ResetNodes();
}

const ArrayT<int>& MeshFreeNodalShapeFunctionT::ResetCells(void) const
{
	return fMFSupport->ResetCells();
}

/* access to MLS field neighbor data */
/*const iArrayT& MeshFreeNodalShapeFunctionT::ElementNeighborsCounts(void) const
{
	return fMFSupport->ElementNeighborsCounts();
}

const RaggedArray2DT<int>& MeshFreeNodalShapeFunctionT::ElementNeighbors(void) const
{
	return fMFSupport->ElementNeighbors();
}*/

const RaggedArray2DT<int>& MeshFreeNodalShapeFunctionT::NodeNeighbors(void) const
{
	return fMFSupport->NodeNeighbors();
}

/* reconstruct displacement field */
void MeshFreeNodalShapeFunctionT::SelectedNodalField(const dArray2DT& all_DOF,
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

void MeshFreeNodalShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
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

void MeshFreeNodalShapeFunctionT::NodalField(const dArray2DT& DOF, dArray2DT& field,
	dArray2DT& Dfield, iArrayT& nodes)
{

	/* fetch list of nodes to compute */
	nodes.Alias(fMFSupport->NodesUsed());
	
	/* dimensions */
	int nnd = nodes.Length();
	int ndf = DOF.MinorDim();
	
	/* allocate output space */
	field.Dimension(nnd, ndf);
	Dfield.Dimension(nnd, ndf*fSD);

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
		Du.Set(ndf, fSD, Dfield(i));
		for (int j = 0; j < ndf; j++)
		{
			dof.Set(len, locdisp(j));

			/* displacement field */
			field(i, j) = dArrayT::Dot(dof, phi);
			
			/* first derivatives */
			for (int k = 0; k < fSD; k++)
			{
				Dphi.RowAlias(k, tmp);
				Du(j, k) = dArrayT::Dot(dof, tmp);
			}
		}
	}
}

/* print the shape function values to the output stream */
void MeshFreeNodalShapeFunctionT::Print(ostream& out) const
{

	out << "\n Field shape functions:\n";
	for (int i = 0; i < fNeighbors.Length(); i++)
	{
		out << setw(kIntWidth) << fNeighbors[i];
		
		for (int j = 0; j < 1; j++)
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

void MeshFreeNodalShapeFunctionT::PrintAt(ostream& out) const
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
void MeshFreeNodalShapeFunctionT::WriteParameters(ostream& out) const
{
	fMFSupport->WriteParameters(out);
}

void MeshFreeNodalShapeFunctionT::WriteStatistics(ostream& out) const
{
	fMFSupport->WriteStatistics(out);
}

dArray2DT& MeshFreeNodalShapeFunctionT::NodalParameters(void) 
{ 
	return fMFSupport->NodalParameters(); 
}

dArrayT& MeshFreeNodalShapeFunctionT::NodalVolumes(void)
{
	return fMFSupport->NodalVolumes();
}

/***********************************************************************
* Private
***********************************************************************/

/* read/write nodal meshfree parameters */
void MeshFreeNodalShapeFunctionT::SetNodalParameters(const iArrayT& node, const dArray2DT& nodal_params)
{
	fMFSupport->SetSupportParameters(node, nodal_params);
}

void MeshFreeNodalShapeFunctionT::GetNodalParameters(const iArrayT& node, dArray2DT& nodal_params) const
{
	fMFSupport->GetSupportParameters(node, nodal_params);
}

