/* $Id: D2MeshFreeSupportT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#include "D2MeshFreeSupportT.h"

#include <math.h>
#include <string.h>

#include "ExceptionCodes.h"
#include "Constants.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "ParentDomainT.h"
#include "iGridManagerT.h"

/* variable length memory managers */
#include "nVariArray2DT.h"
#include "VariArrayT.h"
#include "nArrayGroupT.h"

/* D2 MLS solvers */
#include "D2OrthoMLS2DT.h"
#include "MLSSolverT.h"

/* constructor */
D2MeshFreeSupportT::D2MeshFreeSupportT(const ParentDomainT& domain, const dArray2DT& coords,
	const iArray2DT& connects, const iArrayT& nongridnodes, FormulationT code,
	double dextra, int complete, bool store_shape):
	MeshFreeSupportT(domain, coords, connects, nongridnodes, code, dextra,
		complete, store_shape),
	fD2EFG(NULL)
{
	/* only EFG solver is different for D2 */
	if (code == kEFG)
	{
		/* construct D2 MLS solver */
		if (fCoords.MinorDim() == 2)
			fD2EFG = new D2OrthoMLS2DT(fComplete);
		else
		{
			cout << "\n D2MeshFreeSupportT::D2MeshFreeSupportT: no 3D yet" << endl;
			throw eBadInputValue;
		}
		if (!fD2EFG) throw eOutOfMemory;	
		fD2EFG->Initialize();
	
	//TEMP - this will be better later
	
		/* set inherited */
		delete fEFG;
		fEFG = fD2EFG;
	}
}

/* steps to initialization - modifications to the support size must
* occur before setting the neighbor data */
void D2MeshFreeSupportT::SetNeighborData(void)
{
	/* inherited */
	MeshFreeSupportT::SetNeighborData();

//TEMP - reset nodal work space for higher order derivatives
//       this process could be redesigned

	int nip        = fDomain.NumIP();
	int nsd        = fCoords.MinorDim();
	int stress_dim = dSymMatrixT::NumValues(nsd);

	/* space for nodal calculations */
	int max_n_size = fnNeighborData.MaxMinorDim();
	fndShapespace.Allocate(max_n_size*(1 + nsd + stress_dim));

	/* data and element integration point shape functions */
	int max_e_size = feNeighborData.MaxMinorDim();
	felShapespace.Allocate(nip*max_e_size*(1 + nsd + stress_dim));
}

/* "load" data for the specified node (global numbering) */
void D2MeshFreeSupportT::LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
	dArray2DT& Dphi, dArray2DT& DDphi)
{
	int tag = node ; //TEMP - leftover from OFFSET
	
	/* fetch neighbor data */
	fnNeighborData.RowAlias(tag, neighbors);

	/* dimensions */
	int nsd = fCoords.MinorDim();
	int nst = dSymMatrixT::NumValues(nsd);
	int nnd = neighbors.Length();
	
	if (fStoreShape)
	{
		/* recompute */
		if (fReformNode != kNoReform)
		{
			fReformNode = kNoReform;
			SetNodalShapeFunctions();
		}

#if __option (extended_errorcheck)
		if (fnPhiData.MinorDim(node) < 1)
		{
			cout << "\n D2MeshFreeSupportT::LoadNodalData: requesting empty data for node: ";
			cout << node + 1 << endl;
			throw eGeneralFail;
		}
#endif

		/* set shallow copies */
		fnPhiData.RowAlias(tag, phi);
		Dphi.Set(nsd, nnd, fnDPhiData(tag));
		DDphi.Set(nst, nnd, fnDDPhiData(tag));
	}
	else
	{
		/* set shallow data */
		double* pdata = fndShapespace.Pointer();
		phi.Set(nnd, pdata);
		pdata += phi.Length();
		
		Dphi.Set(nsd, nnd, pdata);
		pdata += Dphi.Length();
		
		DDphi.Set(nst, nnd, pdata);
	
		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi, DDphi);
	}
}

/* "load" data for the specified element (0...)
* for all integration points in the element */
void D2MeshFreeSupportT::LoadElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi)
{
#if __option(extended_errorcheck)
	if (Dphi.Length() != fDomain.NumIP()) throw eSizeMismatch;
#endif

	/* element neighbors */
	feNeighborData.RowAlias(element, neighbors);

	/* dimensions */
	int nip = fDomain.NumIP();
	int nsd = fCoords.MinorDim();
	int nst = dSymMatrixT::NumValues(nsd);
	int nnd = neighbors.Length();

	if (fStoreShape)
	{
		/* recompute */
		if (fReformElem != kNoReform)
		{
			fReformElem = kNoReform;
			SetElementShapeFunctions();
		}
	
		/* load functions */
		phi.Set(nip, nnd, fePhiData(element));
			
		/* load derivatives */
		double*  Dphi_ptr = feDPhiData(element);
		double* DDphi_ptr = feDDPhiData(element);
		for (int i = 0; i < nip; i++)
		{
			/* 1st derivatives */
			Dphi[i].Set(nsd, nnd, Dphi_ptr);
			Dphi_ptr += Dphi[i].Length();

			/* 2nd derivatives */
			DDphi[i].Set(nst, nnd, DDphi_ptr);
			DDphi_ptr += DDphi[i].Length();
		}
	}
	else
	{
		/* set shallow space */
		double* pelspace = felShapespace.Pointer();
		phi.Set(nip, nnd, pelspace);
		pelspace += phi.Length();
		
		/* loop over integration points */
		for (int i = 0; i < nip; i++)
		{
			/* pointers for 1st derivatives */
			Dphi[i].Set(nsd, nnd, pelspace);
			pelspace += Dphi[i].Length();
			
			/* pointers for 2nd derivatives */
			DDphi[i].Set(nst, nnd, pelspace);
			pelspace += DDphi[i].Length();
		}
		
		/* compute */
		ComputeElementData(element, neighbors, phi, Dphi, DDphi);
	}
}

/* return the field derivatives at the specified point */
int D2MeshFreeSupportT::SetFieldAt(const dArrayT& x, AutoArrayT<int>& nodes)
//const dArray2DT& D2MeshFreeSupportT::FieldDerivativesAt(const dArrayT& x, AutoArrayT<int>& tags)
{
// collect nodes within a cellspan and then collect list of
// the union of the nodes and neighbors and filter those to
// find "all" the nodes_i within dmax_i of x. There is now
// mechanism by which to change the local nodes dmax if the
// search for a neighborhood fails.

	/* check nodes and their neighbors */
	BuildNeighborhood(x, nodes);
			
	/* check */
	int dim = (fD2EFG) ? fD2EFG->NumberOfMonomials() :
	                      fRKPM->BasisDimension();
	if (nodes.Length() < dim)
	{
		cout << "\n D2MeshFreeSupportT::SetFieldAt: could not build neighborhood at:\n";
		cout << x << endl;
		return 0;
	}
	else
	{
		/* dimension */
		fcoords_man.SetMajorDimension(nodes.Length(), false);	
		fdmax_man.Dimension(nodes.Length(), false);
	
		/* collect local lists */
		fcoords.RowCollect(nodes, fCoords);
		fdmax.Collect(nodes, fDmax);
	
		/* compute MLS field */
		if (fD2EFG)
			fD2EFG->SetField(fcoords, fdmax, x);
		else
		{
			fvolume.Collect(nodes, fVolume);
			fRKPM->SetField(fcoords, fdmax, fvolume, x, 2);
		}
		
		return 1;
	}
}

/* return values */
const dArray2DT& D2MeshFreeSupportT::DDFieldAt(void) const
{	
	return (fD2EFG) ? fD2EFG->DDphi() : fRKPM->DDphi();
}

/*************************************************************************
* Protected
*************************************************************************/

/* (re-)compute some/all nodal shape functions and derivatives */
void D2MeshFreeSupportT::SetNodalShapeFunctions(void)
{
	/* initialize database */
	InitNodalShapeData();

	/* shallow space */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;
	dArray2DT DDphi;

	/* selectively or all */
	iArrayT nodes;
	if (fResetNodes.Length() > 0)
		nodes.Alias(fResetNodes);
	else
		nodes.Alias(fNodesUsed);

	/* fill data tables */
	int len = nodes.Length();
	for (int j = 0; j < len; j++)
	{
		/* current node (global numbering) */
		int node = nodes[j];

		/* load space */
		LoadNodalData(node, neighbors, phi, Dphi, DDphi);

		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi, DDphi);
	}
	
	/* clear */
	fResetNodes.Allocate(0);
}

/* compute all integration point shape functions and derivatives */
void D2MeshFreeSupportT::SetElementShapeFunctions(void)
{
	/* initialize database */
	InitElementShapeData();

	/* dimensions */
	int nip = fDomain.NumIP();
	int nel = fConnects.MajorDim();

	/* work space */
	iArrayT    neighbors;
	dArray2DT phi;
	ArrayT<dArray2DT> Dphi(nip);
	ArrayT<dArray2DT> DDphi(nip);

	/* selectively or all */
	ArrayT<int>* elems  = (fResetElems.Length() > 0) ? &fResetElems : NULL;
	int lim = (elems != NULL) ? fResetElems.Length() : nel;
	
	for (int j = 0; j < lim; j++)
	{
		int elem = (elems != NULL) ? fResetElems[j] : j; // 0,...
	
		/* fetch element data and space */
		LoadElementData(elem, neighbors, phi, Dphi, DDphi);

		/* compute */
		ComputeElementData(elem, neighbors, phi, Dphi, DDphi);
	}

	/* clear */
	fResetElems.Allocate(0);
}

/*************************************************************************
* Private
*************************************************************************/

void D2MeshFreeSupportT::ComputeNodalData(int node, const iArrayT& neighbors,
	dArrayT& phi, dArray2DT& Dphi, dArray2DT& DDphi)
{
	/* set dimensions */
	int count = neighbors.Length();
	fdmax_man.Dimension(count, false);
	fcoords_man.SetMajorDimension(count, false);
	
	/* collect local lists */
	fcoords.RowCollect(neighbors, fCoords);
	fdmax.Collect(neighbors, fDmax);
		
	/* coords of current node */
	dArrayT x_node;
	fCoords.RowAlias(node, x_node);
	
	/* process boundaries */
	fdmax_ip = fdmax;
	ProcessBoundaries(fcoords, x_node, fdmax);
	// set dmax = -1 for nodes that are inactive at x_node
		
	/* compute MLS field */
	if (fD2EFG)
	{
		fD2EFG->SetField(fcoords, fdmax_ip, x_node);
			
		/* copy field data */
		phi   = fD2EFG->phi();
		Dphi  = fD2EFG->Dphi();
		DDphi = fD2EFG->DDphi();
	}
	else
	{
		fvolume.Collect(neighbors, fVolume);
		fRKPM->SetField(fcoords, fdmax_ip, fvolume, x_node, 2);
			
		/* copy field data */
		phi   = fRKPM->phi();
		Dphi  = fRKPM->Dphi();
		DDphi = fRKPM->DDphi();
	}
}

void D2MeshFreeSupportT::ComputeElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi)
{
	/* dimensions */
	int nsd = fCoords.MinorDim();
	int nip = fDomain.NumIP();
	int nnd = neighbors.Length();
	int nen = fConnects.MinorDim();

	/* set dimensions */
	fdmax_man.Dimension(nnd, false);
	fcoords_man.SetMajorDimension(nnd, false);

	/* collect neighbor data */
	fcoords.RowCollect(neighbors, fCoords);
	fdmax.Collect(neighbors, fDmax);

	/* workspace */
	iArrayT     elementnodes;
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, nsd);
	loccoords.SetGlobal(fCoords);
		
	/* integration point coordinates */
	fConnects.RowAlias(element, elementnodes);
	loccoords.SetLocal(elementnodes);
	fDomain.Interpolate(loccoords, fx_ip_table);
		
	/* loop over integration points */
	dArrayT x_ip;
	for (int i = 0; i < nip; i++)
	{
		/* fetch ip coordinates */
		fx_ip_table.RowAlias(i, x_ip);

		/* process boundaries */
		fdmax_ip = fdmax;
		ProcessBoundaries(fcoords, x_ip, fdmax_ip);
		// set dmax = -1 for nodes that are inactive at x_node
	
		/* compute MLS field */
		if (fD2EFG)
		{
			fD2EFG->SetField(fcoords, fdmax_ip, x_ip);
		
			/* store field data */
			phi.SetRow(i, fD2EFG->phi());
			Dphi[i]  = fD2EFG->Dphi();
			DDphi[i] = fD2EFG->DDphi();
		}
		else
		{
			fvolume.Collect(neighbors, fVolume);
			fRKPM->SetField(fcoords, fdmax_ip, fvolume, x_ip, 2);
		
			/* store field data */
			phi.SetRow(i, fRKPM->phi());
			Dphi[i]  = fRKPM->Dphi();
			DDphi[i] = fRKPM->DDphi();
		}
	}
}

/* allocate and set pointers for shape function databases */
void D2MeshFreeSupportT::InitNodalShapeData(void)
{
	/* inherited */
	MeshFreeSupportT::InitNodalShapeData();

	/* dimensions */
	int nst = dSymMatrixT::NumValues(fCoords.MinorDim());

	/* configure nodal storage */
	fnDDPhiData.Configure(fnNeighborCount, nst);	
}

void D2MeshFreeSupportT::InitElementShapeData(void)
{
	/* inherited */
	MeshFreeSupportT::InitElementShapeData();

	/* dimensions */
	int nst = dSymMatrixT::NumValues(fCoords.MinorDim());
	int nip = fDomain.NumIP();

	/* configure element storage */
	feDDPhiData.Configure(feNeighborCount, nip*nst);
}
