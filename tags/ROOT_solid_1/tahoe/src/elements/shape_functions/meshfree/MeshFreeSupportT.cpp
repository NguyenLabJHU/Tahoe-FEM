/* $Id: MeshFreeSupportT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (09/07/1998)                                          */
/* supporting functions for EFG elements                                  */

#include "MeshFreeSupportT.h"

#include <math.h>
#include <string.h>
#include "ExceptionCodes.h"
#include "Constants.h"
#include "dArray2DT.h"

/* variable length memory managers */
#include "nVariArray2DT.h"
#include "VariArrayT.h"
#include "nArrayGroupT.h"

/* search grid */
#include "iGridManagerT.h"
#include "iNodeT.h"

/* MLS solvers */
#include "OrthoMLS2DT.h" // 2D EFG
#include "OrthoMLS3DT.h" // 3D EFG
#include "MLSSolverT.h"  // RKPM 2D/3D

/* element integration domain operations */
#include "ParentDomainT.h"
#include "LocalArrayT.h"

/* parameters */
const int kMaxNumGrid   = 100;
const int kListInitSize =  10; // initial neighbor list size
const int kListSizeinc  =   5; // list size increment on overflow

static    int Max(int a, int b) { return (a > b) ? a : b; };
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
MeshFreeSupportT::MeshFreeSupportT(const ParentDomainT& domain, const dArray2DT& coords,
	const iArray2DT& connects, const iArrayT& nongridnodes, FormulationT code,
	double dextra, int complete, bool store_shape):
	fDomain(domain),
	fEFG(NULL),
	fRKPM(NULL),
	fGrid(NULL),
	fCoords(coords),
	fConnects(connects),
	fNonGridNodes(nongridnodes),
	fDextra(dextra),
	fComplete(complete),
	fStoreShape(store_shape),
	fNumFacetNodes(0),
	fCutCoords(NULL),
	fdmax_man(25),
	fcoords_man(25, fcoords, fCoords.MinorDim()),
	fx_ip_table(fDomain.NumIP(), fCoords.MinorDim()),
	fReformNode(kNotInit),
	fReformElem(kNotInit)	
{
	/* checks */
	if (fDomain.NumSD()    !=   fCoords.MinorDim() ||
	    fDomain.NumNodes() != fConnects.MinorDim()) throw eBadInputValue;
	if (fDextra < 0.0) throw eBadInputValue;

	/* construct MLS solver */
	if (code == kEFG)
	{
		if (fDomain.NumSD() == 2)
			fEFG = new OrthoMLS2DT(fComplete);
		else
			fEFG = new OrthoMLS3DT(fComplete);
		if (!fEFG) throw eOutOfMemory;	
		fEFG->Initialize();
	}
	else if (code == kRKPM)
	{
		fRKPM = new MLSSolverT(fDomain.NumSD(), fComplete);
		if (!fRKPM) throw eOutOfMemory;	
		fRKPM->Initialize();
	}
	else
		throw eBadInputValue;

	/* variable size managers */
	fdmax_man.Register(fvolume );
	fdmax_man.Register(fdmax   );
	fdmax_man.Register(fdmax_ip);
}

/* destructor */
MeshFreeSupportT::~MeshFreeSupportT(void)
{
	delete fEFG;
	delete fRKPM;
	delete fGrid;
}

/* steps to initialization - modifications to the support size must
* occur before setting the neighbor data */
void MeshFreeSupportT::SetSupportSize(void)
{
	/* collect numbers for used nodes */
	SetNodesUsed();

	/* set search grid */
	SetSearchGrid();
//NOTE: do this every time SetSupportSize is called???
	
	/* initialize Dmax for nodes in connectivity set */
	cout << "\n MeshFreeSupportT::SetSupportSize: setting nodal support" << endl;
	SetDmax();
}

void MeshFreeSupportT::SetNeighborData(void)
{
	/* data and nodal shape functions */
	cout << " MeshFreeSupportT::SetNeighborData: setting nodal data" << endl;
	SetNodeNeighborData(fCoords);

	/* data and element integration point shape functions */
	cout << " MeshFreeSupportT::SetNeighborData: setting integration point data" << endl;
	SetElementNeighborData(fConnects);
}

/* specify nodes/cells to skip when doing MLS calculations */
void MeshFreeSupportT::SetSkipNodes(const iArrayT& skip_nodes)
{
	if (skip_nodes.Length() == 0)
		fSkipNode.Free();
	else
	{
		fSkipNode.Allocate(fCoords.MajorDim());
		fSkipNode = 0;
		for (int i = 0; i < skip_nodes.Length(); i++)
			fSkipNode[skip_nodes[i]] = 1;
	}
}

void MeshFreeSupportT::SetSkipElements(const iArrayT& skip_elements)
{
	if (skip_elements.Length() == 0)
		fSkipElement.Free();
	else
	{
		fSkipElement.Allocate(fConnects.MajorDim());
		fSkipElement = 0;
		for (int i = 0; i < skip_elements.Length(); i++)
			fSkipElement[skip_elements[i]] = 1;
	}
}

/* synchronize Dmax with another set (of active EFG nodes) */
void MeshFreeSupportT::SynchronizeDmax(dArrayT& synchDmax)
{
	/* should be over the same global node set (marked by length) */
	if (fDmax.Length() != synchDmax.Length())
	{
		cout << "\n MeshFreeSupportT::SynchronizeDmax: must be over same global node set" << endl;
		throw eSizeMismatch;
	}

	double* pthis = fDmax.Pointer();
	double* pthat = synchDmax.Pointer();
	int length = fDmax.Length();
	for (int i = 0; i < length; i++)
	{
		*pthis = *pthat = Max(*pthis,*pthat);
		pthis++; pthat++;
	}
}

void MeshFreeSupportT::SetDmax(const iArrayT& node, const dArrayT& Dmax)
{
	if (node.Length() != Dmax.Length()) throw eSizeMismatch;

	/* Dmax initialized externally */
	if (fDmax.Length() == 0)
	{
		if (node.Length() != fCoords.MajorDim())
		{
			cout << "\n MeshFreeSupportT::SetDmax: must initialize d_max for ALL nodes" << endl;
			throw eGeneralFail;
		}
	
		/* initialize */
		fDmax.Allocate(fCoords.MajorDim());
		fDmax = -1;

		/* assume not following "normal" route to initialization */
		if (!fGrid)
		{
			fNodesUsed.Allocate(fCoords.MajorDim());
			fNodesUsed.SetValueToPosition();
			SetSearchGrid();
		}

		//TEMP - set nodal volume all to 1.0
		if (fRKPM)
		{
			fVolume.Allocate(fCoords.MajorDim());
			fVolume = 1.0;
		}
	}

	/* map in */
	for (int i = 0; i < node.Length(); i++)
		fDmax[node[i]] = Dmax[i];
}

void MeshFreeSupportT::GetDmax(const iArrayT& node, dArrayT& Dmax) const
{
	Dmax.Collect(node, fDmax);
}

/* cutting facet functions */
void MeshFreeSupportT::SetCuttingFacets(const dArray2DT& facet_coords,
	int num_facet_nodes)
{
	/* set cutting facet data */
	fNumFacetNodes = num_facet_nodes;
	fCutCoords = &facet_coords;
}

void MeshFreeSupportT::ResetFacets(const ArrayT<int>& facets)
{
	/* shape functions not yet initialized */
	if (fReformNode == -1 && fReformElem == -1) return;

	/* check */
	if (facets.Length() == 0)
		return;
	else if (fCutCoords == NULL)
	{
		cout << "\n MeshFreeSupportT::ResetFacets: facet coordinates not set" << endl;
		throw eGeneralFail;
	}
	
	/* collect nodes to reset */
	int nsd = fDomain.NumSD();
	dArray2DT facet_coords;
	AutoArrayT<int> resetnodes;
	for (int ii = 0; ii < facets.Length(); ii++)
	{
		/* alias to facet data */
		facet_coords.Set(fNumFacetNodes, nsd, (*fCutCoords)(facets[ii]));

		/* compute characteristic facet size */
		double sqr_max = 0.0;
		for	(int k = 1; k < fNumFacetNodes; k++)
		{
			double sqr_d = 0.0;
			double* A = facet_coords(k-1);
			double* B = facet_coords(k);
			for (int i = 0; i < fDomain.NumSD(); i++)
			{
				double dx = B[i] - A[i];
				sqr_d += dx*dx;
			}
			
			sqr_max = (sqr_d > sqr_max) ? sqr_d : sqr_max;
		}
		double facet_size = sqrt(sqr_max);
	
		/* collect nodes */
		for (int j = 0; j < fNumFacetNodes; j++)
		{
			const AutoArrayT<iNodeT>& inodes = fGrid->HitsInRegion(facet_coords(j), facet_size);
			for (int i = 0; i < inodes.Length(); i++)
			{
				int nd = inodes[i].Tag();
	
				/* must be meshfree node */
				if (fDmax[nd] > 0.0)			
				{
					/* the node */
					resetnodes.AppendUnique(nd);
					
					/* its meshfree neighbors */
					int  nnd = fnNeighborData.MinorDim(nd);
					int* pnd = fnNeighborData(nd);
					for (int k = 0; k < nnd; k++)
					{
						if (fDmax[*pnd] > 0.0) resetnodes.AppendUnique(*pnd);
						pnd++;
					}
				}
			}
		}
	}

	/* found nodes to reset */
	if (resetnodes.Length() > 0)
	{
		/* element shape functions not yet initialized */
		if (fReformElem != kNotInit)
		{
			/* collect elements containing any reset nodes */
			iArrayT elem;
			int nel = fConnects.MajorDim();
			for (int k = 0; k < resetnodes.Length(); k++)
				/* check elements for presence of reset node */
				for (int i = 0; i < nel; i++)
					if (!fResetElems.HasValue(i))
					{
							/* connectivity */
						feNeighborData.RowAlias(i, elem);
					
						/* check element */
						if (elem.HasValue(resetnodes[k])) fResetElems.Append(i);
					}
	
			/* set flag */
			fReformElem = kReform;
		}
		// worthwhile to keep nodes per element neighborhood map, i.e., inverse of
		// feNeighborData, to save time here?

		/* node flag */
		if (fReformNode != kNotInit)
		{
			fResetNodes.AppendUnique(resetnodes); // add to any existing nodes
			fReformNode = kReform;
		}
	}
}

/* "load" data for the specified node (global numbering) */
void MeshFreeSupportT::LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
	dArray2DT& Dphi)
{
	int tag = node; //TEMP - from before OFFSET fix
	
	/* fetch neighbor data */
	fnNeighborData.RowAlias(tag, neighbors);

	/* dimensions */
	int nsd = fCoords.MinorDim();
	int nnd = neighbors.Length();
	
	if (fStoreShape)
	{
		/* recompute */
		if (fReformNode != kNoReform)
		{
			fReformNode = kNoReform;
			cout << " MeshFreeSupportT::LoadNodalData: computing nodal shape functions" << endl;
			SetNodalShapeFunctions();
		}

#if __option (extended_errorcheck)
		if (fnPhiData.MinorDim(node) < 1)
		{
			cout << "\n MeshFreeSupportT::LoadNodalData: requesting empty data for node: ";
			cout << node << endl;
			throw eGeneralFail;
		}
#endif

		/* set shallow copies */
		fnPhiData.RowAlias(tag, phi);
		Dphi.Set(nsd, nnd, fnDPhiData(tag));
	}
	else
	{
		/* set shallow data */
		double* pdata = fndShapespace.Pointer();
		phi.Set(nnd, pdata);
		Dphi.Set(nsd, nnd, pdata + nnd);
	
		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi);
	}
}

/* "load" data for the specified element (0...)
* for all integration points in the element */
void MeshFreeSupportT::LoadElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi)
{
#if __option(extended_errorcheck)
	if (Dphi.Length() != fDomain.NumIP()) throw eSizeMismatch;
#endif

	/* element neighbors */
	feNeighborData.RowAlias(element, neighbors);

	/* dimensions */
	int nip = fDomain.NumIP();
	int nsd = fCoords.MinorDim();
	int nnd = neighbors.Length();

	if (fStoreShape)
	{
		/* recompute */
		if (fReformElem != kNoReform)
		{
			fReformElem = kNoReform;
			cout << " MeshFreeSupportT::LoadElementData: computing int. pt. shape functions" << endl;
			SetElementShapeFunctions();

			/* determine storage efficiency */
			if (fePhiData.Length() > 0)
			{
				int zero_count = 0;
				double* phi = fePhiData.Pointer();
				for (int i = 0; i < fePhiData.Length(); i++)
					if (fabs(*phi++) < kSmall)
						zero_count++;
				cout << " MeshFreeSupportT::LoadElementData: int. pt. storage efficiency: "
				     << 1.0 - double(zero_count)/fePhiData.Length() << endl;
			}
		}
	
		/* load functions */
		phi.Set(nip, nnd, fePhiData(element));
			
		/* load derivatives */
		double* Dphi_ptr = feDPhiData(element);
		for (int i = 0; i < nip; i++)
		{
			Dphi[i].Set(nsd, nnd, Dphi_ptr);
			Dphi_ptr += nsd*nnd;
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
			Dphi[i].Set(nsd, nnd, pelspace);
			pelspace += Dphi[i].Length();
		}
		
		/* compute */
		ComputeElementData(element, neighbors, phi, Dphi);
	}
}

/* return the field derivatives at the specified point */
int MeshFreeSupportT::SetFieldUsing(const dArrayT& x, const ArrayT<int>& nodes)
{
	/* check */
	int dim = (fEFG) ? fEFG->NumberOfMonomials() :
	                  fRKPM->BasisDimension();
	if (nodes.Length() < dim)
	{
		cout << "\n MeshFreeSupportT::SetFieldUsing: could not build neighborhood at:\n";
		cout << x << '\n';
		cout << " insufficient number of nodes: " << nodes.Length() << "/" << dim << '\n';
		iArrayT tmp;
		tmp.Alias(nodes);
		tmp++;
		cout << tmp.wrap(5) << endl;
		tmp--;
		return 0;
	}
	else
	{
		/* copy neighor nodes */
		fneighbors = nodes;
	
		/* dimension */
		fcoords_man.SetMajorDimension(fneighbors.Length(), false);	
		fdmax_man.Dimension(fneighbors.Length(), false);
	
		/* collect local lists */
		fcoords.RowCollect(fneighbors, fCoords);
		fdmax.Collect(fneighbors, fDmax);
	
		/* compute MLS field */
		int OK;
		if (fEFG)
			OK = fEFG->SetField(fcoords, fdmax, x);
		else
		{
			fvolume.Collect(fneighbors, fVolume);
			OK = fRKPM->SetField(fcoords, fdmax, fvolume, x, 1);
		}

		/* error */
		if (!OK)
		{
			int d_width = cout.precision() + kDoubleExtra;
			cout << "\n MeshFreeSupportT::SetFieldUsing: could not compute:\n";
			cout << " coordinates :" << x.no_wrap() << '\n';
			cout << " neighborhood: " << fneighbors.Length() << '\n';
			cout << setw(kIntWidth) << "node"
			     << setw(  d_width) << "dmax"
			     << setw(  d_width) << "dist"
			     << setw(  d_width) << "x" << '\n';
			dArrayT dist(x.Length());			
			for (int i = 0; i < fneighbors.Length(); i++)
			{
				cout << setw(kIntWidth) << fneighbors[i] + 1;
				cout << setw(  d_width) << fdmax[i];
				fcoords.RowCopy(i, dist);
				dist -= x;		
				cout << setw(  d_width) << dist.Magnitude();
				fcoords.PrintRow(i, cout);
			}
			cout.flush();
			return 0;
		}
		else
			return 1;
	}
}

int MeshFreeSupportT::SetFieldAt(const dArrayT& x, const dArrayT* shift)
{
	/* collect all nodes covering x */
	int result;
	if (shift != NULL)
	{
		dArrayT x_shift(x.Length());
		x_shift.SumOf(x, *shift);
		result = BuildNeighborhood(x_shift, fneighbors);
	}
	else
		result = BuildNeighborhood(x, fneighbors);

	/* insufficient neighborhood */
	if (!result)
		cout << "\n MeshFreeSupportT::SetFieldAt: BuildNeighborhood: failed" << endl;
	else
	{
		/* set field */
		result = SetFieldUsing(x, fneighbors);
		if (!result)
			cout << "\n MeshFreeSupportT::SetFieldAt: SetFieldUsing: failed" << endl;
	}
	return result;
}

/* return values */
const dArrayT& MeshFreeSupportT::FieldAt(void) const
{	
	return (fEFG) ? fEFG->phi() : fRKPM->phi();
}

const dArray2DT& MeshFreeSupportT::DFieldAt(void) const
{	
	return (fEFG) ? fEFG->Dphi() : fRKPM->Dphi();
}

/* write MLS statistics */
void MeshFreeSupportT::WriteStatistics(ostream& out) const
{
	/* neighbor count data */
	int min = fnNeighborData.MinMinorDim(0);
	int max = fnNeighborData.MaxMinorDim();
	int used = fNodesUsed.Length();
	out << "\n MLS shape function data:\n";
	out << " Minimum number of nodal neighbors . . . . . . . = " << min << '\n';
	out << " Maximum number of nodal neighbors . . . . . . . = " << max << '\n';
	out << " Average number of nodal neighbors . . . . . . . = ";
	if (used != 0)
		out << fnNeighborData.Length()/used << '\n';
	else
		out << "-\n";

	/* collect neighbor distribution */
	iArrayT counts(max + 1);
	counts = 0;

	/* skip nodes */
	if (fSkipNode.Length() == fnNeighborData.MajorDim())
	{
		for (int j = 0; j < fnNeighborData.MajorDim(); j++)
			if (!fSkipNode[j])
				counts[fnNeighborData.MinorDim(j)]++;
	}	
	else
		for (int j = 0; j < fnNeighborData.MajorDim(); j++)
			counts[fnNeighborData.MinorDim(j)]++;

	out << " Nodal neighbor number distribution:\n";
	out << setw(kIntWidth) << "number"
	    << setw(kIntWidth) << "count" << '\n';
	for (int k = 0; k < counts.Length(); k++)
		out << setw(kIntWidth) << k
	    	<< setw(kIntWidth) << counts[k] << '\n';	

	/* neighbor distance data */
	double dmin, dmax;
	fDmax.MinMax(dmin, dmax, 1);

	double*  pd = fDmax.Pointer();
	double dsum = 0.0;
	for (int i = 0; i < fDmax.Length(); i++)
	{
		if (*pd > 0.0) dsum += *pd;
		pd++;
	}

	out << " Largest support radius (dmax) . . . . . . . . . = " << dmin << '\n';
	out << " Smallest support radius (dmax). . . . . . . . . = " << dmax << '\n';
	out << " Average support radius (dmax) . . . . . . . . . = ";
	if (used != 0)
		out << dsum/fNodesUsed.Length() << '\n';
	else
		out << "-\n";

	/* memory requirements */
	int nsd = fCoords.MinorDim();
	int nip = fDomain.NumIP();
	out << "\n MLS storage requirements:\n";
	int n_count = fnNeighborCount.Sum();
	out << " Total number of nodal neighbors . . . . . . . . = " << n_count << '\n';
	out << " Nodal shape function storage. . . . . . . . . . = " << n_count*sizeof(double) << " bytes\n";
	out << " Nodal shape function derivatives storage. . . . = " << nsd*n_count*sizeof(double) << " bytes\n";
	int e_count = feNeighborCount.Sum();
	out << " Total number of integration point neighbors . . = " << e_count << '\n';
	out << " i.p. shape function storage . . . . . . . . . . = " << nip*e_count*sizeof(double) << " bytes\n";
	out << " i.p. shape function derivatives storage . . . . = " << nip*nsd*e_count*sizeof(double) << " bytes\n";

	/* search grid statistics */
	fGrid->WriteStatistics(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* initialize search grid */
void MeshFreeSupportT::SetSearchGrid(void)
{
	/* must set nodes used first */
	if (fNodesUsed.Length() == 0)
	{
		cout << "\n MeshFreeSupportT::SetSearchGrid: must set used nodes first" << endl;
		throw eGeneralFail;
	}

	/* try to get roughly least 10 per grid */
	int ngrid = int(pow(fNodesUsed.Length()/10.0, 1.0/fDomain.NumSD())) + 1;
	ngrid = (ngrid < 2) ? 2 : ngrid;
	ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

	/* free any existing */
	delete fGrid;

	/* construct a search grid */
	iArrayT n_grid(fDomain.NumSD());
	n_grid = ngrid;
	fGrid = new iGridManagerT(n_grid, fCoords, &fNodesUsed);	
	if (!fGrid) throw eOutOfMemory;

	/* reset search grid bounds */
	fGrid->Reset();
}

/* generate lists of all nodes that fall within Dmax of the
* nodal coords (self included) */
void MeshFreeSupportT::SetNodeNeighborData(const dArray2DT& coords)
{
#if __option(extended_errorcheck)
	if (coords.MajorDim() != fDmax.Length()) throw eSizeMismatch;
#endif

	int numnodes = fDmax.Length();
	int nsd      = coords.MinorDim();
	
	/* initialize neighbor counts */
	fnNeighborCount.Allocate(numnodes);
	fnNeighborCount = 0;

	/* work space */
	iArray2DT a1(numnodes, kListInitSize), a2;
	iArray2DT* pthis = &a1;
	iArray2DT* pnext = &a2;

	/* loop over "active" nodes */
	AutoArrayT<int> neighbors;
	for (int ndex = 0; ndex < fNodesUsed.Length(); ndex++)
	{
		int i = fNodesUsed[ndex];
		double* target = coords(i);
		double  tol    = fDmax[i];
		
		/* find nodes in neighborhood */
		fGrid->Neighbors(i, tol, neighbors);
		neighbors.Append(i); // add self to the list
			
		/* need more space */
		int count = neighbors.Length();
		if (count > (*pthis).MinorDim())
		{
			(*pnext).Allocate(numnodes, count + kListSizeinc);
					
			/* copy and swap */
			SwapData(fnNeighborCount, &pthis, &pnext);
		}
					
		/* set data */
		fnNeighborCount[i] = count;
		memcpy((*pthis)(i), neighbors.Pointer(), sizeof(int)*count);
	}					

	/* configure node neighbor data */
	fnNeighborData.Configure(fnNeighborCount);
	
	/* copy data in */
	for (int j = 0; j < numnodes; j++)
		fnNeighborData.SetRow(j, (*pthis)(j));

	/* space for nodal calculations */
	int maxsize = fnNeighborData.MaxMinorDim();
	fndShapespace.Allocate(maxsize*(1 + nsd));
}

/* generate lists of all nodes that fall within Dmax of the
* element integration points */
void MeshFreeSupportT::SetElementNeighborData(const iArray2DT& connects)
{
	/* dimensions */
	int numelems = connects.MajorDim();
	int nen      = connects.MinorDim();
	int nsd      = fCoords.MinorDim();
	int nip      = fDomain.NumIP();
	
	/* initialize neighbor counts */
	feNeighborCount.Allocate(numelems);
	feNeighborCount = 0;

	/* work space */
	iArray2DT a1(numelems, kListInitSize), a2;
	iArray2DT* pthis = &a1;
	iArray2DT* pnext = &a2;
	
	/* domain coords */
	iArrayT     elementnodes;
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, nsd);
	loccoords.SetGlobal(fCoords);

	/* integration point coordinate tables */	
	dArray2DT x_ip_table(nip, nsd);
	dArrayT   x_ip, x_node;
	dArrayT   x_ab(nsd);

	/* "big" system */
	const int big = 25000;
	bool big_system = numelems >= big;
	if (big_system)
		cout << setw(2*kIntWidth) << "count"
		     << setw(kIntWidth) << "min"
		     << setw(kIntWidth) << "max"
		     << setw(kDoubleWidth) << "avg" << endl;
	int max_count = -1, min_count = -1, sum_count = 0;

	/* loop over elements */
	AutoArrayT<int> nodeset;
	iArrayT neighbors;
	for (int i = 0; i < numelems; i++)
	{
		/* "big" system */
		if (big_system && fmod(i + 1, big) == 0)
		{
			cout << setw(2*kIntWidth) << i + 1
			     << setw(kIntWidth) << min_count
			     << setw(kIntWidth) << max_count
			     << setw(kDoubleWidth) << double(sum_count)/big <<endl;

			/* reset */
			max_count = -1, min_count = -1, sum_count = 0;
		}

		/* connectivity */
		int* pelem = connects(i);
		
		/* integration point coordinates */
		fConnects.RowAlias(i, elementnodes);
		loccoords.SetLocal(elementnodes);
		fDomain.Interpolate(loccoords, x_ip_table);
	
		/* collect unique list of neighboring nodes */
		nodeset.Allocate(0);
		for (int j = 0; j < nen; j++)
		{
			int node = pelem[j];
		
			/* fetch nodal neighors */
			fnNeighborData.RowAlias(node, neighbors);
			
			/* add to element list */
			for (int n = 0; n < neighbors.Length(); n++)
			{
				int neigh = neighbors[n];
			
				if (!nodeset.HasValue(neigh))
				{
					/* fetch nodal coords */
					fCoords.RowAlias(neigh, x_node);
				
					/* loop over integration points */
					int hit = 0;
					for (int k = 0; k < nip && !hit; k++)
					{
						/* fetch ip coordinates */
						x_ip_table.RowAlias(k, x_ip);
				
						/* distance */
						x_ab.DiffOf(x_node, x_ip);
				
						/* will contribute */
						if (x_ab.Magnitude() < fDmax[neigh])
							hit = 1;
					}
				
					/* within the neighborhood */	
					if (hit) nodeset.Append(neigh);
				}
			}
		}
	
		/* need more space */
		int setlength = nodeset.Length();
		if (setlength > (*pthis).MinorDim())
		{
			int newsize = Max(setlength,(*pthis).MinorDim() + kListSizeinc);

			/* allocate */		
			(*pnext).Allocate(numelems, newsize);

			/* copy and swap */
			SwapData(feNeighborCount, &pthis, &pnext);		
		}

		/* statistics for "big" systems */
		if (big_system)
		{
			if (max_count == -1 || setlength > max_count)
				max_count = setlength;			
			if (min_count == -1 || setlength < min_count)
				min_count = setlength;
			sum_count += setlength;
		}
		
		/* add to database */
		feNeighborCount[i] = setlength;
		memcpy((*pthis)(i), nodeset.Pointer(), sizeof(int)*setlength);	
	}

	/* configure element neighbor data */
	feNeighborData.Configure(feNeighborCount);

	/* copy data in */
	for (int j = 0; j < numelems; j++)
		feNeighborData.SetRow(j, (*pthis)(j));
		
	/* space element calculation */
	int maxsize = feNeighborData.MaxMinorDim();
	felShapespace.Allocate(nip*maxsize*(1 + nsd));
}

/* (re-)compute some/all nodal shape functions and derivatives */
void MeshFreeSupportT::SetNodalShapeFunctions(void)
{
	/* initialize database */
	InitNodalShapeData();

	/* shallow space */
	iArrayT    neighbors;
	dArrayT    phi;
	dArray2DT  Dphi;

	/* selectively or all */
	iArrayT nodes;
	if (fResetNodes.Length() > 0)
	{
		nodes.Alias(fResetNodes);
//DEBUG
//nodes.Alias(fNodesUsed);
//cout << "\n MeshFreeSupportT::SetNodalShapeFunctions: **** no selective computation ****" << endl;
//DEBUG
	}
	else
		nodes.Alias(fNodesUsed);
		
	/* message */
	cout << " MeshFreeSupportT::SetNodalShapeFunctions: node count: "
	     << nodes.Length() << endl;		

	/* fill data tables */
	int len = nodes.Length();
	for (int j = 0; j < len; j++)
	{
		/* current node (global numbering) */
		int node = nodes[j];

		/* load space */
		LoadNodalData(node, neighbors, phi, Dphi);

		/* compute */
		ComputeNodalData(node, neighbors, phi, Dphi);
	}
	
	/* clear */
	fResetNodes.Allocate(0);
}

/* compute all integration point shape functions and derivatives */
void MeshFreeSupportT::SetElementShapeFunctions(void)
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

	/* selectively or all */
	ArrayT<int>* elems  = (fResetElems.Length() > 0) ? &fResetElems : NULL;

//DEBUG
//ArrayT<int>* elems = NULL;
//if (fResetElems.Length() > 0)
//cout << "\n MeshFreeSupportT::SetElementShapeFunctions: **** no selective computation ****" << endl;
//DEBUG

	int lim = (elems != NULL) ? fResetElems.Length() : nel;
	for (int j = 0; j < lim; j++)
	{
		int elem = (elems != NULL) ? fResetElems[j] : j; // 0,...
	
		/* fetch element data and space */
		LoadElementData(elem, neighbors, phi, Dphi);

		/* compute */
		ComputeElementData(elem, neighbors, phi, Dphi);
	}

	/* message */
	cout << " MeshFreeSupportT::SetElementShapeFunctions: cell count: "
	     << lim << endl;		

	/* clear */
	fResetElems.Allocate(0);	
}

/*************************************************************************
* Private
*************************************************************************/

void MeshFreeSupportT::ComputeElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi)
{
	/* skip nodes */
	if (fSkipElement.Length() > 0 && fSkipElement[element] == 1)
	{
		/* zero values */
		phi = 0.0;
		for (int i = 0; i < Dphi.Length(); i++)
			Dphi[i] = 0.0;
		return ;
	}

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
		int OK;
		if (fEFG)
		{
			OK = fEFG->SetField(fcoords, fdmax_ip, x_ip);
	
			/* store field data */
			phi.SetRow(i, fEFG->phi());
			Dphi[i] = fEFG->Dphi();
		}
		else
		{
			fvolume.Collect(neighbors, fVolume);
			OK = fRKPM->SetField(fcoords, fdmax_ip, fvolume, x_ip, 1);

			/* store field data */
			phi.SetRow(i, fRKPM->phi());
			Dphi[i] = fRKPM->Dphi();
		}

		/* error */
		ostream& err = cout;
		if (!OK)
		{
			int d_width = err.precision() + kDoubleExtra;
			err << "\n MeshFreeSupportT::ComputeElementData: could not compute field at \n"
			     <<   "     integration point " << i+1 << " of element " << element+1<< '\n';
			err << " coordinates: " << x_ip.no_wrap() << '\n';
			err << " neighborhood: " << neighbors.Length() << '\n';
			err << setw(kIntWidth) << "node"
			     << setw(  d_width) << "dmax"
			     << setw(  d_width) << "dist"
			     << setw(  d_width) << "x" << '\n';
			dArrayT dist(x_ip.Length());
			for (int j = 0; j < neighbors.Length(); j++)
			{
				err << setw(kIntWidth) << neighbors[j] + 1;
				err << setw(  d_width) << fdmax_ip[j];
				fcoords.RowCopy(j, dist);
				dist -= x_ip;		
				err << setw(  d_width) << dist.Magnitude();
				fcoords.PrintRow(j, err);				
			}
			err.flush();
			throw eGeneralFail;	
		}	
	}
}

void MeshFreeSupportT::ComputeNodalData(int node, const iArrayT& neighbors,
	dArrayT& phi, dArray2DT& Dphi)
{
	/* skip nodes */
	if (fSkipNode.Length() > 0 && fSkipNode[node] == 1)
	{
		/* zero values */
		phi = 0.0;
		Dphi = 0.0;
		return;
	}
	
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
	ProcessBoundaries(fcoords, x_node, fdmax);
	// set dmax = -1 for nodes that are inactive at x_node

	/* compute MLS field */
	int OK;
	if (fEFG)
	{
		OK = fEFG->SetField(fcoords, fdmax, x_node);
		
		/* copy field data */
		phi  = fEFG->phi();
		Dphi = fEFG->Dphi();
	}
	else
	{
		fvolume.Collect(neighbors, fVolume);
		OK = fRKPM->SetField(fcoords, fdmax, fvolume, x_node, 1);
		
		/* copy field data */
		phi  = fRKPM->phi();
		Dphi = fRKPM->Dphi();
	}
	
	/* error */
	if (!OK)
	{
		int d_width = cout.precision() + kDoubleExtra;
		cout << "\n MeshFreeSupportT::ComputeNodalData: could not compute field at node: "
		     << node + 1 << '\n';
		cout << " coordinates: " << x_node.no_wrap() << '\n';
		cout << " neighborhood: " << neighbors.Length() << '\n';
		cout << setw(kIntWidth) << "node"
		     << setw(  d_width) << "dmax"
		     << setw(  d_width) << "dist"
		     << setw(  d_width) << "x" << '\n';
		dArrayT dist(x_node.Length());
		for (int i = 0; i < neighbors.Length(); i++)
		{
			cout << setw(kIntWidth) << neighbors[i] + 1;
			cout << setw(  d_width) << fdmax[i];
			fcoords.RowCopy(i, dist);
			dist -= x_node;		
			cout << setw(  d_width) << dist.Magnitude();
			fcoords.PrintRow(i, cout);
		}
		cout.flush();		
		throw eGeneralFail;	
	}
}

/* set dmax for each node in the connectivities */
void MeshFreeSupportT::SetDmax(void)
{
	/* dimensions */
	int nnd = fCoords.MajorDim();
	int nsd = fCoords.MinorDim();

	//TEMP - set nodal volume all to 1.0
	if (fRKPM)
	{
		fVolume.Allocate(nnd);
		fVolume = 1.0;
	}

	/* initialize max distance to "inactive" */
	fDmax.Allocate(nnd);
	fDmax = -1.0;

	int min_neighbors = (fEFG) ? fEFG->NumberOfMonomials() :
	                            fRKPM->BasisDimension();
	//TEMP
	min_neighbors += 1; // don't count self

	/* quick check */
	if (fNodesUsed.Length() < min_neighbors)
	{
		cout << " MeshFreeSupportT::SetDmax: not enough meshfree nodes" << endl;
		throw eGeneralFail;
	}

	/* "big" system */
	const int big = 25000;
	bool big_system = fNodesUsed.Length() >= big;
	if (big_system)
		cout << " MeshFreeSupportT::SetDmax: setting nodal neighborhoods:\n"
		     << setw(2*kIntWidth) << "count"
		     << setw(kIntWidth) << "min"
		     << setw(kIntWidth) << "max"
		     << setw(kDoubleWidth) << "avg" << endl;
	int max_count = -1, min_count = -1, sum_count = 0;

	/* loop over all nodes */
	AutoArrayT<double> dists;
	AutoArrayT<int> nodes;
//	dArrayT dists_sh;
	iArrayT nodes_sh;
	double dmax = 0.0;
	for (int ndex = 0; ndex < fNodesUsed.Length(); ndex++)
	{
		/* "big" system */
		if (big_system && fmod(ndex + 1, big) == 0)
		{
			cout << setw(2*kIntWidth) << ndex + 1
			     << setw(kIntWidth) << min_count
			     << setw(kIntWidth) << max_count
			     << setw(kDoubleWidth) << double(sum_count)/big <<endl;

			/* reset */
			max_count = -1, min_count = -1, sum_count = 0;
		}

		/* resolve tag */
		int i = fNodesUsed[ndex];
		double* target = fCoords(i);

		/* search parameters */
		bool cell_search = true;
		int cell_span = 0;
		int num_neighbors = 0;
		double tol = 0.0;

		int iteration = 0;
		int max_iterations = 3;
		while (num_neighbors < min_neighbors && ++iteration <= max_iterations)
		{
			/* get neighborhood nodes */
			const AutoArrayT<iNodeT>& inodes = (cell_search) ?
				fGrid->HitsInRegion(target, ++cell_span) :
				fGrid->HitsInRegion(target, tol);

			/* set valid neighborhood */
			dists.Allocate(0);
			nodes.Allocate(0);
		
			/* collect visible nodes */
			num_neighbors = inodes.Length();
			for (int j = 0; num_neighbors >= min_neighbors && j < inodes.Length(); j++)
			{
				const double* coords = inodes[j].Coords();
				
				/* check visibility */
				if (Visible(target, coords))
				{
					/* compute distances^2 */
					double dist = 0.0;
					for (int k = 0; k < nsd; k++)
					{
						double dx_k = coords[k] - target[k];
						dist += dx_k*dx_k;
					}
					dists.Append(dist);
					nodes.Append(inodes[j].Tag());
				}
				else
					num_neighbors--;
			}
			
			/* sort values */
			if (num_neighbors >= min_neighbors)
			{
				//dists_sh.Alias(dists);
				//dists_sh.SortAscending(nodes);
				nodes_sh.Alias(nodes);
				nodes_sh.SortAscending(dists);

				/* support size */						
				//dmax = 1.01*sqrt(dists_sh[min_neighbors - 1]);
				dmax = 1.01*sqrt(dists[min_neighbors - 1]);

				/* check validity of search */
				if (cell_search && fGrid->CellSpan(cell_span) < dmax)
				{
					/* wipe results */
					num_neighbors = 0;
					
					/* do distance search */
					cell_search = false;
					tol = 1.01*dmax;
				}
			}
		}

		/* neighborhood growing too big */
		if (num_neighbors < min_neighbors && iteration > max_iterations)
		{
			cout << "\n MeshFreeSupportT::SetDmax: failed to find valid neighborhood around point\n";
			cout << "     " << i << " after " << max_iterations << " iterations:\n";
			cout << "      x: ";
			fCoords.PrintRow(i, cout);
			throw eGeneralFail;
		}
		else
		{
			/* include nodes at same distance */
			bool repeat = true;
			int num_support = min_neighbors;
			for (int j = min_neighbors - 1; j < num_neighbors - 1 && repeat; j++)
				//if (fabs(dists_sh[j] - dists_sh[j+1]) < kSmall)
				if (fabs(dists[j] - dists[j+1]) < kSmall)
					num_support++;
				else
					repeat = false;
		
			/* insurance for constructing the nodal fields */
			for (int i = 0; i < num_support; i++)
			{
				double& d_i = fDmax[nodes[i]];
				d_i = (d_i < dmax) ? dmax : d_i;
			}
			
			/* statistics for "big" systems */
			if (big_system)
			{
				if (max_count == -1 || num_neighbors > max_count)
					max_count = num_neighbors;			
				if (min_count == -1 || num_neighbors < min_count)
					min_count = num_neighbors;
				sum_count += num_neighbors;
			}
		}
	}
	
	/* scale neighborhoods */
	fDmax *= fDextra;

//TEMP
#if 0
cout << " MeshFreeSupportT::SetDmax:\n";
for (int i = 0; i < fDmax.Length(); i++)
	cout << setw(kIntWidth) << i+1 << setw(kDoubleWidth) << fDmax[i] << '\n';
#endif
//END		
						
}

/* set dmax for each node in the connectivities */
void MeshFreeSupportT::SetDmaxFromConnects(void)
{
// NOTE: this version of finding Dmax only works if all the nodes
//       are members of the connectivities defining the integration
//       grid cells. The Dmax computed here is also strictly not what
//       its definition says: the distance to the farthest node in the
//       smallest neighborhood that's required for the MLS fit.
	
	/* dimensions */
	int nnd = fCoords.MajorDim();
	int nsd = fCoords.MinorDim();
	int nel = fConnects.MajorDim();
	int nen = fConnects.MinorDim();

	/* loop over elements */
	dArrayT x_0, x_i;
	dArrayT r(nsd);
	for (int i = 0; i < nel; i++)
	{
		int* pelem = fConnects(i);
		for (int j = 0; j < nen; j++)
		{
			/* current node */
			int node = pelem[j];

			/* fetch origin coords */
			fCoords.RowAlias(node, x_0);
			
			/* find max neighbor distance */
			double& maxdist = fDmax[node];
			for (int k = 0; k < nen; k++)
				if (k != j)
				{
					/* fetch origin coords */
					fCoords.RowAlias(pelem[k], x_i);
			
					/* distance vector */
					r.DiffOf(x_i, x_0);
					double dist = r.Magnitude();
					
					/* keep max */
					maxdist = Max(dist, maxdist);
				}
		}
	}
}

void MeshFreeSupportT::SetNodesUsed(void)
{
	/* dimensions */
	int nnd = fCoords.MajorDim();

	/* markers */
	iArrayT used(nnd);
	used = 0;

	/* mark nodes used on the integration grid */
	int  tot = fConnects.Length();
	int* pnd = fConnects.Pointer();
	for (int i = 0; i < tot; i++)
		used[*pnd++] = 1;
	
	/* mark other EFG nodes */
	tot = fNonGridNodes.Length();
	pnd = fNonGridNodes.Pointer();
	for (int j = 0; j < tot; j++)
		used[*pnd++] = 1;

	/* count nodes used */
	int numused = used.Count(1);
	fNodesUsed.Allocate(numused);
		
	/* collect node numbers (ascending order) */	
	int* pnused = fNodesUsed.Pointer();
	int*  pused = used.Pointer();
	for (int k = 0; k < nnd; k++)
		if (*pused++) *pnused++ = k;
}

/* collect all nodes covering the point x */
int MeshFreeSupportT::BuildNeighborhood(const dArrayT& x, AutoArrayT<int>& nodes)
{
/* NOTE: relies on the fact that the support of any node is big
*       enough to ensure coverage of all neighboring nodes */

	/* collect nodes in neighborhood of point x */
	int cell_span = 0;
	double* target = x.Pointer();
	const AutoArrayT<iNodeT>* inodes = &fGrid->HitsInRegion(target, ++cell_span);
	while (inodes->Length() < 1 && cell_span <= 3)
		inodes = &fGrid->HitsInRegion(target, ++cell_span);
	if (inodes->Length() < 1)
	{
		cout << "\n MeshFreeSupportT::BuildNeighborhood: failed to find any nodes around:\n";
		cout << x << "\n     after 4 iterations" << endl;
		throw eGeneralFail;
	}
	
	/* find biggest support */
	double support_max = 0.0;
	for (int ii = 0; ii < inodes->Length(); ii++)
	{
		int tag = ((*inodes)[ii]).Tag();
		double dmax = fDmax[tag];
		support_max = (dmax > support_max) ? dmax : support_max;
	}
	
	/* need to re-collect nodes */
	if (fGrid->CellSpan(cell_span) < support_max)
		inodes = &fGrid->HitsInRegion(target, support_max);

	/* work space */	
	int nsd = fCoords.MinorDim();
	AutoArrayT<int> out_nodes; // make as class workspace?

	/* initialize */
	nodes.Allocate(0);
	out_nodes.Allocate(0);
				
	/* run through nodes */
	for (int i = 0; i < inodes->Length(); i++)
	{
		int tag = ((*inodes)[i]).Tag();
		double* x_node = fCoords(tag);

//NOTE - reverse this. finding the r^2 is probably cheaper than
//       determining whether the node is visible
	
		/* check node */
		if (!nodes.HasValue(tag) &&     // redundant?
		    !out_nodes.HasValue(tag) && // redundant?
		     Visible(target, x_node))
		{
			double dist = 0.0;
			for (int j = 0; j < nsd; j++)
			{
				double dx = x_node[j] - target[j];
				dist += dx*dx;
			}
				
			/* add to list */
			double dmax_i = fDmax[tag];
			if (dist < dmax_i*dmax_i)
				nodes.Append(tag);
			else
				out_nodes.Append(tag);
		}		
	}

	/* verify */
	int dim = (fEFG) ? fEFG->NumberOfMonomials() :
	                  fRKPM->BasisDimension();
	if (nodes.Length() < dim)
	{
		int d_width = OutputWidth(cout, x.Pointer());
		cout << "\n MeshFreeSupportT::SetFieldUsing: insufficient number of nodes: "
		     << nodes.Length() << "/" << dim << '\n';
		cout <<   "         x: " << x.no_wrap() << '\n'
		     <<   "     nodes: " << inodes->Length() << '\n'
		     <<   " cell span: " << cell_span << '\n'
		     <<   "    radius: " << fGrid->CellSpan(cell_span) << "\n\n";
		cout << setw(kIntWidth) << "node"
		     << setw(5) << "OK"
		     << setw(5) << "vis"
		     << setw(5) << "cov"
		     << setw(d_width) << "dist"
		     << setw(d_width) << "dmax"
		     << setw(d_width) << "x" << '\n';

		/* work space */
		dArrayT xn, r(nsd);
		
		/* write nodes */
		for (int i = 0; i < inodes->Length(); i++)
		{
			int tag = ((*inodes)[i]).Tag();
			fCoords.RowAlias(tag, xn);
			r.DiffOf(x, xn);
			double dist = r.Magnitude();
			cout << setw(kIntWidth) << tag+1
		         << setw(5) << nodes.HasValue(tag)
		         << setw(5) << Visible(x.Pointer(), xn.Pointer())
		         << setw(5) << (dist < fDmax[tag])
		         << setw(  d_width) << dist
		         << setw(  d_width) << fDmax[tag]
		         << setw(  d_width) << xn.no_wrap() << '\n';		
		}
		return 0;
	}
	else
		return 1;
}

/* allocate and set pointers for shape function databases */
void MeshFreeSupportT::InitNodalShapeData(void)
{
	/* dimensions */
	int nsd = fCoords.MinorDim();

	/* configure nodal storage */
	int step;
	try {
		step = 0;
		fnPhiData.Configure(fnNeighborCount);
		step = 1;
		fnDPhiData.Configure(fnNeighborCount, nsd);
	}
	
	/* write memory information */
	catch (int error)
	{
		if (error == eOutOfMemory)
		{
			cout << "\n MeshFreeSupportT::InitNodalShapeData: out of memory for nodal shape data: "
			     << ((step == 0) ? "phi" : "Dphi") << '\n';
			int total_count = fnNeighborCount.Sum();
			cout << "     phi: " << total_count*sizeof(double)     << " bytes\n";
			cout << "    Dphi: " << nsd*total_count*sizeof(double) << " bytes\n";
		}
		throw error;
	}
}

void MeshFreeSupportT::InitElementShapeData(void)
{
	/* dimensions */
	int nsd = fCoords.MinorDim();
	int nip = fDomain.NumIP();

	/* configure element storage */
	int step;
	try {
		step = 0;
		fePhiData.Configure(feNeighborCount, nip);
		step = 1;
		feDPhiData.Configure(feNeighborCount, nip*nsd);
	}

	/* write memory information */
	catch (int error)
	{
		if (error == eOutOfMemory)
		{
			cout << "\n MeshFreeSupportT::InitNodalShapeData: out of memory for nodal shape data: "
			     << ((step == 0) ? "phi" : "Dphi") << '\n';
			int total_count = feNeighborCount.Sum();
			cout << "     phi: " << nip*total_count*sizeof(double)     << " bytes\n";
			cout << "    Dphi: " << nip*nsd*total_count*sizeof(double) << " bytes\n";
		}
		throw error;
	}
}

/* swap data */
void MeshFreeSupportT::SwapData(const iArrayT& counts, iArray2DT** pfrom,
	iArray2DT** pto)
{
#if __option(extended_errorcheck)
	if ((**pfrom).MajorDim() != (**pto).MajorDim()) throw eSizeMismatch;
	if ((**pfrom).MinorDim() >= (**pto).MinorDim()) throw eGeneralFail;
#endif

	iArray2DT& from = **pfrom;
	iArray2DT& to   = **pto;

	/* copy non-empty rows */
	int numrows = counts.Length();
	for (int i = 0; i < numrows; i++)
		if (counts[i] > 0)
			memcpy(to(i), from(i), sizeof(int)*counts[i]);		

	/* swap pointers */
	iArray2DT* temp = *pfrom;
	*pfrom = *pto;
	*pto   = temp;
}
