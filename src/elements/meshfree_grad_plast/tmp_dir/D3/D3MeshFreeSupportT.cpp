/* $Id: D3MeshFreeSupportT.cpp,v 1.3 2004-10-05 15:10:33 kyonten Exp $ */
/* created: paklein (10/23/1999)                                          */

#include "D3MeshFreeSupportT.h"

#include <math.h>
#include <string.h>

#include "ExceptionT.h"
#include "toolboxConstants.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "ParentDomainT.h"
#include "iGridManagerT.h"

/* variable length memory managers */
#include "nVariArray2DT.h"
#include "VariArrayT.h"
#include "nArrayGroupT.h"

/* D3 MLS solvers */
#include "MLSSolverT.h"

using namespace Tahoe;

/* constructor */
D3MeshFreeSupportT::D3MeshFreeSupportT(const ParentDomainT* domain, const dArray2DT& coords,
	const iArray2DT& connects, const iArrayT& nongridnodes):
	D2MeshFreeSupportT(domain, coords, connects, nongridnodes)
{
	/* only EFG solver is different for D3 */
	if (fMeshfreeType == kEFG)
	{
		cout << "\n D3MeshFreeSupportT::D3MeshFreeSupportT: no EFG implemented" << endl;
		throw ExceptionT::kBadInputValue;
	}
}

/* steps to initialization - modifications to the support size must
* occur before setting the neighbor data */
void D3MeshFreeSupportT::InitNeighborData(void)
{
	/* inherited */
	D2MeshFreeSupportT::InitNeighborData();

//TEMP - reset nodal work space for higher order derivatives
//       this process could be redesigned

	int nip        = fDomain->NumIP();
	int nsd        = fCoords->MinorDim();
	int stress_dim = dSymMatrixT::NumValues(nsd);

	/* space for nodal calculations */
	int max_n_size = fnNeighborData.MaxMinorDim();
	fndShapespace.Dimension(max_n_size*(1 + nsd + stress_dim));

	/* data and element integration point shape functions */
	int max_e_size = feNeighborData.MaxMinorDim();
	felShapespace.Dimension(nip*max_e_size*(1 + nsd + stress_dim));
}

/* "load" data for the specified node (global numbering) */
void D3MeshFreeSupportT::LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
	dArray2DT& Dphi, dArray2DT& DDphi, dArray2DT& DDDphi)
{
	int tag = node ; //TEMP - leftover from OFFSET
	
	/* fetch neighbor data */
	fnNeighborData.RowAlias(tag, neighbors);

	/* dimensions */
	int nsd = fCoords->MinorDim();
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
			cout << "\n D3MeshFreeSupportT::LoadNodalData: requesting empty data for node: ";
			cout << node + 1 << endl;
			throw ExceptionT::kGeneralFail;
		}
#endif

		/* set shallow copies */
		fnPhiData.RowAlias(tag, phi);
		Dphi.Set(nsd, nnd, fnDPhiData(tag));
		DDphi.Set(nst, nnd, fnDDPhiData(tag));
		DDDphi.Set(nsd*nsd, nnd, fnDDDPhiData(tag)); //kyonten
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
		pdata += DDphi.Length();
		
		DDDphi.Set(nsd*nsd, nnd, pdata); //kyonten
	
		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi, DDphi, DDDphi);
	}
}

/* "load" data for the specified element (0...)
* for all integration points in the element */
void D3MeshFreeSupportT::LoadElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi, ArrayT<dArray2DT>& DDDphi)
{
#if __option(extended_errorcheck)
	if (Dphi.Length() != fDomain->NumIP()) throw ExceptionT::kSizeMismatch;
#endif

	/* element neighbors */
	feNeighborData.RowAlias(element, neighbors);

	/* dimensions */
	int nip = fDomain->NumIP();
	int nsd = fCoords->MinorDim();
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
		double* DDDphi_ptr = feDDDPhiData(element);
		for (int i = 0; i < nip; i++)
		{
			/* 1st derivatives */
			Dphi[i].Set(nsd, nnd, Dphi_ptr);
			Dphi_ptr += Dphi[i].Length();

			/* 2nd derivatives */
			DDphi[i].Set(nst, nnd, DDphi_ptr);
			DDphi_ptr += DDphi[i].Length();
			
			/* 3rd derivatives */
			DDDphi[i].Set(nsd*nsd, nnd, DDDphi_ptr); //kyonten
			DDDphi_ptr += DDDphi[i].Length();
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
			
			/* pointers for 3rd derivatives */
			DDDphi[i].Set(nsd*nsd, nnd, pelspace); //kyonten
			pelspace += DDDphi[i].Length();
		}
		
		/* compute */
		ComputeElementData(element, neighbors, phi, Dphi, DDphi, DDDphi);
	}
}


/* return values */
const dArray2DT& D3MeshFreeSupportT::DDDFieldAt(void) const
{	
	return fRKPM->DDDphi();  //kyonten (uncommented)
}

/*************************************************************************
* Protected
*************************************************************************/

/* (re-)compute some/all nodal shape functions and derivatives */
void D3MeshFreeSupportT::SetNodalShapeFunctions(void)
{
	/* initialize database */
	InitNodalShapeData();

	/* shallow space */
	iArrayT   neighbors;
	dArrayT   phi;
	dArray2DT Dphi;
	dArray2DT DDphi;
	dArray2DT DDDphi;

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
		LoadNodalData(node, neighbors, phi, Dphi, DDphi, DDDphi);

		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi, DDphi, DDDphi);
	}
	
	/* clear */
	fResetNodes.Dimension(0);
}

/* compute all integration point shape functions and derivatives */
void D3MeshFreeSupportT::SetElementShapeFunctions(void)
{
	/* initialize database */
	InitElementShapeData();

	/* dimensions */
	int nip = fDomain->NumIP();
	int nel = fConnects->MajorDim();

	/* work space */
	iArrayT    neighbors;
	dArray2DT phi;
	ArrayT<dArray2DT> Dphi(nip);
	ArrayT<dArray2DT> DDphi(nip);
	ArrayT<dArray2DT> DDDphi(nip);

	/* selectively or all */
	ArrayT<int>* elems  = (fResetElems.Length() > 0) ? &fResetElems : NULL;
	int lim = (elems != NULL) ? fResetElems.Length() : nel;
	
	for (int j = 0; j < lim; j++)
	{
		int elem = (elems != NULL) ? fResetElems[j] : j; // 0,...
	
		/* fetch element data and space */
		LoadElementData(elem, neighbors, phi, Dphi, DDphi, DDDphi);

		/* compute */
		ComputeElementData(elem, neighbors, phi, Dphi, DDphi, DDDphi);
	}

	/* clear */
	fResetElems.Dimension(0);
}

/*************************************************************************
* Private
*************************************************************************/

void D3MeshFreeSupportT::ComputeNodalData(int node, const iArrayT& neighbors,
	dArrayT& phi, dArray2DT& Dphi, dArray2DT& DDphi, dArray2DT& DDDphi)
{
	/* set dimensions */
	int count = neighbors.Length();
	fnodal_param_man.SetMajorDimension(count, false);
	fcoords_man.SetMajorDimension(count, false);
	
	/* collect local lists */
	fcoords.RowCollect(neighbors, *fCoords);
	fnodal_param.Collect(neighbors, fNodalParameters);
		
	/* coords of current node */
	dArrayT x_node;
	fCoords->RowAlias(node, x_node);
	
	/* process boundaries */
	fnodal_param_ip = fnodal_param;
	ProcessBoundaries(fcoords, x_node, fnodal_param_ip);
	// set dmax = -1 for nodes that are inactive at x_node
		
	/* compute MLS field */
	if (fD2EFG)
	{
		cout << "\n D3MeshFreeSupportT::D3MeshFreeSupportT: no EFG implemented" << endl;
		throw ExceptionT::kBadInputValue;
	}
	else
	{
		fvolume_man.SetLength(count, false);
		fvolume.Collect(neighbors, fVolume);
		fRKPM->SetField(fcoords, fnodal_param_ip, fvolume, x_node, 2);
			
		/* copy field data */
		phi   = fRKPM->phi();
		Dphi  = fRKPM->Dphi();
		DDphi = fRKPM->DDphi();
		DDDphi = fRKPM->DDDphi(); //kyonten (uncommented)
	}
}

void D3MeshFreeSupportT::ComputeElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi, ArrayT<dArray2DT>& DDDphi)
{
	/* dimensions */
	int nsd = fCoords->MinorDim();
	int nip = fDomain->NumIP();
	int nnd = neighbors.Length();
	int nen = fConnects->MinorDim();

	/* set dimensions */
	fnodal_param_man.SetMajorDimension(nnd, false);
	fcoords_man.SetMajorDimension(nnd, false);
	fvolume_man.SetLength(nnd, false);

	/* collect neighbor data */
	fcoords.RowCollect(neighbors, *fCoords);
	fnodal_param.RowCollect(neighbors, fNodalParameters);

	/* workspace */
	iArrayT     elementnodes;
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, nsd);
	loccoords.SetGlobal(*fCoords);
		
	/* integration point coordinates */
	fConnects->RowAlias(element, elementnodes);
	loccoords.SetLocal(elementnodes);
	fDomain->Interpolate(loccoords, fx_ip_table);
		
	/* loop over integration points */
	dArrayT x_ip;
	for (int i = 0; i < nip; i++)
	{
		/* fetch ip coordinates */
		fx_ip_table.RowAlias(i, x_ip);

		/* process boundaries */
		fnodal_param_ip = fnodal_param;
		ProcessBoundaries(fcoords, x_ip, fnodal_param_ip);
		// set dmax = -1 for nodes that are inactive at x_node
	
		/* compute MLS field */
		if (fD2EFG)
		{
			cout << "\n D3MeshFreeSupportT::D3MeshFreeSupportT: no EFG implemented" << endl;
			throw ExceptionT::kBadInputValue;
		}
		else
		{
			fvolume.Collect(neighbors, fVolume);
			fRKPM->SetField(fcoords, fnodal_param_ip, fvolume, x_ip, 2);
		
			/* store field data */
			phi.SetRow(i, fRKPM->phi());
			Dphi[i]  = fRKPM->Dphi();
			DDphi[i] = fRKPM->DDphi();
			DDDphi[i] = fRKPM->DDDphi();  //kyonten (uncommented) 
		}
	}
}

/* allocate and set pointers for shape function databases */
void D3MeshFreeSupportT::InitNodalShapeData(void)
{
	/* inherited */
	D2MeshFreeSupportT::InitNodalShapeData();

	/* dimensions */
	int nst = dSymMatrixT::NumValues(fCoords->MinorDim());
	int nsd = fCoords->MinorDim(); //kyonten

	/* configure nodal storage */
	fnDDDPhiData.Configure(fnNeighborCount, nsd*nsd); //kyonten	
}

void D3MeshFreeSupportT::InitElementShapeData(void)
{
	/* inherited */
	D2MeshFreeSupportT::InitElementShapeData();

	/* dimensions */
	int nst = dSymMatrixT::NumValues(fCoords->MinorDim());
	int nsd = fCoords->MinorDim(); //kyonten
	int nip = fDomain->NumIP();

	/* configure element storage */
	feDDDPhiData.Configure(feNeighborCount, nip*(nsd*nsd)); //kyonten
}
