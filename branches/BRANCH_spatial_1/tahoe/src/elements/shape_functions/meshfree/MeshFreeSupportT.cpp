/* $Id: MeshFreeSupportT.cpp,v 1.32 2004-10-14 20:24:56 paklein Exp $ */
/* created: paklein (09/07/1998) */
#include "MeshFreeSupportT.h"

#include <math.h>
#include <string.h>
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArray2DT.h"
#include "ParameterContainerT.h"

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

/* use disk to reduce memory usage */
#undef _USE_DISK_

using namespace Tahoe;

/* parameters */
const int kMaxNumGrid   = 100;
const int kListInitSize =  10; // initial neighbor list size
const int kListSizeinc  =   5; // list size increment on overflow

static    int Max(int a, int b) { return (a > b) ? a : b; };
static double Max(double a, double b) { return (a > b) ? a : b; };
static double AbsMax(double a, double b) { return (fabs(a) > fabs(b)) ? fabs(a) : fabs(b); };

/* constructor */
MeshFreeSupportT::MeshFreeSupportT(const ParentDomainT* domain,  
	const dArray2DT& coords, const iArray2DT& connects, const iArrayT& nongridnodes):
	ParameterInterfaceT("meshfree_support"),
	fDomain(domain),
	fDextra(0.0),
	fEFG(NULL),
	fRKPM(NULL),
	fGrid(NULL),
	fCoords(&coords),
	fConnects(&connects),
	fNonGridNodes(&nongridnodes),
	fNumFacetNodes(0),
	fCutCoords(NULL),
	fvolume_man(25, fvolume),
	fnodal_param_man(25, true),
	fcoords_man(25, fcoords, fCoords->MinorDim()),
	fReformNode(kNotInit),
	fReformElem(kNotInit)	
{
	if (fDomain)
	{
		fIP = fDomain->NumIP();
		fx_ip_table.Dimension(fIP, fCoords->MinorDim());
		fSD = fDomain->NumSD();
		
		/* checks */
		if (fSD != fCoords->MinorDim() || fDomain->NumNodes() != fConnects->MinorDim()) 
	    		ExceptionT::BadInputValue("MeshFreeSupportT::MeshFreeSupportT","Dimension mismatch\n");
	}
	else // Presumably nodally integrated elements are calling this
	{
		fx_ip_table.Dimension(1, fCoords->MinorDim());
		fSD = fCoords->MinorDim();
		fIP = 1;
	}
}

MeshFreeSupportT::MeshFreeSupportT(void):
	ParameterInterfaceT("meshfree_support"),
	fDomain(NULL),
	fDextra(0.0),
	fEFG(NULL),
	fRKPM(NULL),
	fGrid(NULL),
	fCoords(NULL),
	fConnects(NULL),
	fNonGridNodes(NULL),
	fNumFacetNodes(0),
	fCutCoords(NULL),
	fnodal_param_man(25, true),
	fReformNode(kNotInit),
	fReformElem(kNotInit)	
{

}

/* destructor */
MeshFreeSupportT::~MeshFreeSupportT(void)
{
	delete fEFG;
	delete fRKPM;
	delete fGrid;
}

/* write parameters */
void MeshFreeSupportT::WriteParameters(ostream& out) const
{
	/* common parameters */
	out << "\n Meshfree support parameters:\n";
	out << " Store shape functions . . . . . . . . . . . . . = " << ((fStoreShape) ? "TRUE" : "FALSE") << '\n';
	out << " Meshfree formulation. . . . . . . . . . . . . . = " << fMeshfreeType << '\n';
	out << "    [" << kEFG << "]: Element-free Galerkin (EFG)\n";
	out << "    [" << kRKPM << "]: Reproducing Kernel Particle Method (RPKM)\n";	
	if (fMeshfreeType == kEFG)
	{
		out << " Meshfree formulation. . . . . . . . . . . . . . = " << fDextra << '\n';
		out << " Basis function completeness . . . . . . . . . . = " << fEFG->Completeness() << '\n';
		out << " Number of basis function monomials. . . . . . . = " << fEFG->NumberOfMonomials() << '\n';
	}
	else if (fMeshfreeType == kRKPM)
	{
		fRKPM->WriteParameters(out);
	}
	else throw ExceptionT::kGeneralFail;
	out << '\n';
}

/* steps to initialization - modifications to the support size must
* occur before setting the neighbor data */
void MeshFreeSupportT::InitSupportParameters(void)
{
	/* collect numbers for used nodes */
	SetNodesUsed();

	/* set search grid */
	SetSearchGrid();
	//NOTE: do this every time SetSupportParameters is called???
	
	/* initialize support size for nodes in connectivity set */
	cout << "\n MeshFreeSupportT::InitSupportParameters: setting nodal support" << endl;
	if (fMeshfreeType == kEFG || fRKPM->SearchType() == WindowT::kSpherical)
	{
		cout << "\n MeshFreeSupportT::InitSupportParameters: spherical search" << endl;
		SetSupport_Spherical_Search();
	}
	else if (fRKPM->SearchType() == WindowT::kConnectivity)
	{
		cout << "\n MeshFreeSupportT::InitSupportParameters: connectivity search" << endl;
		SetSupport_Cartesian_Connectivities();
	} 
	else throw ExceptionT::kGeneralFail;
}

void MeshFreeSupportT::InitNeighborData(void)
{
	/* data and nodal shape functions */
	cout << " MeshFreeSupportT::InitNeighborData: setting nodal data" << endl;
//	SetNodeNeighborData(fCoords);
	SetNodeNeighborData_2(*fCoords);

	/* data and element integration point shape functions */
	cout << " MeshFreeSupportT::InitNeighborData: setting integration point data" << endl;
//	SetElementNeighborData(fConnects);
	if (fDomain)
		SetElementNeighborData_2(*fConnects);
}

/* specify nodes/cells to skip when doing MLS calculations */
void MeshFreeSupportT::SetSkipNodes(const iArrayT& skip_nodes)
{
	if (skip_nodes.Length() == 0)
		fSkipNode.Free();
	else
	{
		fSkipNode.Dimension(fCoords->MajorDim());
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
		fSkipElement.Dimension(fConnects->MajorDim());
		fSkipElement = 0;
		for (int i = 0; i < skip_elements.Length(); i++)
			fSkipElement[skip_elements[i]] = 1;
	}
}

/* synchronize Dmax with another set (of active EFG nodes) */
void MeshFreeSupportT::SynchronizeSupportParameters(dArray2DT& nodal_params)
{
	if (fMeshfreeType == kEFG)
	{
		/* should be over the same global node set (marked by length) */
		if (fNodalParameters.Length() != nodal_params.Length())
		{
			cout << "\n MeshFreeSupportT::SynchronizeSupportParameters: list of nodal\n"
			     << " parameters must the same length over the global node set" << endl;
			throw ExceptionT::kSizeMismatch;
		}
		
		/* "synchronize" means take max of dmax */
		double* pthis = fNodalParameters.Pointer();
		double* pthat = nodal_params.Pointer();
		int length = fNodalParameters.Length();
		for (int i = 0; i < length; i++)
		{
			*pthis = *pthat = Max(*pthis,*pthat);
			pthis++; pthat++;
		}
	}
	else if (fMeshfreeType == kRKPM)
		/* handled by the MLS solver */
		fRKPM->SynchronizeSupportParameters(fNodalParameters, nodal_params);
	else throw ExceptionT::kGeneralFail;
}

void MeshFreeSupportT::SetSupportParameters(const iArrayT& node, const dArray2DT& nodal_params)
{
	/* make sure grid is set - not a good place for this */
	if (!fGrid)
	{
		fNodesUsed.Dimension(fCoords->MajorDim());
		fNodesUsed.SetValueToPosition();
		SetSearchGrid();
	}

	/* map in */
	fNodalParameters.Assemble(node, nodal_params);
}

void MeshFreeSupportT::GetSupportParameters(const iArrayT& node, dArray2DT& nodal_params) const
{
	/* copy selected rows */
	nodal_params.RowCollect(node, fNodalParameters);
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
		throw ExceptionT::kGeneralFail;
	}
	
	/* collect nodes to reset.
	 * The set of nodes to reset is determined by first calculating
	 * a characteristic facet size. Then, all nodes within this
	 * distance from the facet nodes are collected. Any node with
	 * a non-zero support size is appended to the list for recomputation.
	 * This approach does NOT strictly guarantee that affected nodes
	 * are found. More correctly, we would need to find all nodes whose
	 * support is cut by the facet. */
	dArray2DT facet_coords;
	AutoArrayT<int> resetnodes;
	for (int ii = 0; ii < facets.Length(); ii++)
	{
		/* alias to facet data */
		facet_coords.Alias(fNumFacetNodes, fSD, (*fCutCoords)(facets[ii]));

		/* compute characteristic facet size */
		double sqr_max = 0.0;
		for	(int k = 1; k < fNumFacetNodes; k++)
		{
			double sqr_d = 0.0;
			double* A = facet_coords(k-1);
			double* B = facet_coords(k);
			for (int i = 0; i < fSD; i++)
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
	
				/* must an "active" meshfree node - determined by non-zero
				 * nodal parameter, i.e., the support size */
				if (fNodalParameters(nd,0) > 0.0)			
				{
					/* the node */
					resetnodes.AppendUnique(nd);
					
					/* its meshfree neighbors */
					int  nnd = fnNeighborData.MinorDim(nd);
					int* pnd = fnNeighborData(nd);
					for (int k = 0; k < nnd; k++)
					{
						if (fNodalParameters(*pnd,0) > 0.0) resetnodes.AppendUnique(*pnd);
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
			int nel = fConnects->MajorDim();
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
	int nsd = fCoords->MinorDim();
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
			throw ExceptionT::kGeneralFail;
		}
#endif

		/* set shallow copies */
		fnPhiData.RowAlias(tag, phi);
		Dphi.Set(nsd, nnd, fnDPhiData(tag));
	}
	else
	{
		/* check work space */
		if (fndShapespace.Length() < nnd*(nsd + 1))
		{
			cout << " MeshFreeSupportT::LoadNodalData: work space is not allocated" << endl;
			throw ExceptionT::kGeneralFail;
		}
	
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
	if (fDomain && Dphi.Length() != fDomain->NumIP()) throw ExceptionT::kSizeMismatch;
#endif

	/* element neighbors */
	feNeighborData.RowAlias(element, neighbors);

	/* dimensions */
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
		phi.Set(fIP, nnd, fePhiData(element));
			
		/* load derivatives */
		double* Dphi_ptr = feDPhiData(element);
		for (int i = 0; i < fIP; i++)
		{
			Dphi[i].Set(fSD, nnd, Dphi_ptr);
			Dphi_ptr += fSD*nnd;
		}
	}
	else
	{
		/* set shallow space */
		double* pelspace = felShapespace.Pointer();
		phi.Set(fIP, nnd, pelspace);
		pelspace += phi.Length();
		
		/* loop over integration points */
		for (int i = 0; i < fIP; i++)
		{
			Dphi[i].Set(fSD, nnd, pelspace);
			pelspace += Dphi[i].Length();
		}

		/* check */
		if (pelspace - felShapespace.Pointer() > felShapespace.Length())
		{
			cout << " MeshFreeSupportT::LoadElementData: element work space is not allocated" << endl;
			throw ExceptionT::kGeneralFail;
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
		fnodal_param_man.SetMajorDimension(fneighbors.Length(), false);
	
		/* collect local lists */
		fcoords.RowCollect(fneighbors, *fCoords);
		fnodal_param.RowCollect(fneighbors, fNodalParameters);
	
		/* compute MLS field */
		int OK;
		if (fEFG)
			OK = fEFG->SetField(fcoords, fnodal_param, x);
		else
		{
			/* nodal volumes */
			fvolume_man.SetLength(fneighbors.Length(), false);
			fvolume.Collect(fneighbors, fVolume);
			
			/* compute field */
			OK = fRKPM->SetField(fcoords, fnodal_param, fvolume, x, 1);
		}

		/* error */
		if (!OK)
		{
			int d_width = cout.precision() + kDoubleExtra;
			cout << "\n MeshFreeSupportT::SetFieldUsing: could not compute:\n";
			cout << " coordinates :" << x.no_wrap() << '\n';
			cout << " neighborhood: " << fneighbors.Length() << '\n';
			cout << setw(kIntWidth) << "node"
			     << setw(  d_width) << "dist"
			     << setw(fnodal_param.MinorDim()*d_width) << "nodal parameters"
			     << setw(fcoords.MinorDim()*d_width) << "x" << '\n';
			dArrayT dist(x.Length());			
			for (int i = 0; i < fneighbors.Length(); i++)
			{
				cout << setw(kIntWidth) << fneighbors[i] + 1;
				fcoords.RowCopy(i, dist);
				dist -= x;		
				cout << setw(  d_width) << dist.Magnitude();
				fnodal_param.PrintRow(i, cout);
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
	out << "\n MLS shape function data:\n";

	/* nodal neighbor count data */
	if (true) {
	int min = fnNeighborData.MinMinorDim(0);
	int max = fnNeighborData.MaxMinorDim();
	int used = fNodesUsed.Length();
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
	}

	/* nodal support data */
	int d_width = OutputWidth(out, fNodalParameters.Pointer());
	out << "\n Support size distribution (unscaled):\n";
	out << setw(d_width) << "min"
	    << setw(d_width) << "max"
	    << setw(d_width) << "avg" << '\n';
	for (int i = 0; i < fNodalParameters.MinorDim(); i++)
	{
		double min, max, sum;
		min = -1.0; 
		max = sum = 0.0;
		int count = 0;
		for (int j = 0; j < fNodalParameters.MajorDim(); j++)
		{
			double value = fNodalParameters(j,i);
			if (value > 0.0)
			{
				if (min == -1 || value < min)
					min = value;
				else if (value > max)
					max = value;
				sum += value;
				count++;
			}
		}
		
		/* results */
		out << setw(d_width) << min
	        << setw(d_width) << max;
		if (fNodalParameters.MajorDim() > 0)
			out << setw(d_width) << ((count > 0) ? sum/count : 0.0);
		else
			out << setw(d_width) << 0.0;
		out << '\n';
	}

	/* element neighbor count data */
	if (true) {
	int min = feNeighborData.MinMinorDim(0);
	int max = feNeighborData.MaxMinorDim();

	/* collect neighbor distribution */
	iArrayT counts(max + 1);
	counts = 0;

	/* skip nodes */
	int used = 0;
	if (fSkipElement.Length() == feNeighborData.MajorDim())
	{
		for (int j = 0; j < feNeighborData.MajorDim(); j++)
			if (!fSkipElement[j]) {
				counts[feNeighborData.MinorDim(j)]++;
				used++;
			}
	}	
	else
	{
		used = feNeighborData.MajorDim();
		for (int j = 0; j < feNeighborData.MajorDim(); j++)
			counts[feNeighborData.MinorDim(j)]++;		
	}

	out << " Minimum number of element neighbors . . . . . . . = " << min << '\n';
	out << " Maximum number of element neighbors . . . . . . . = " << max << '\n';
	out << " Average number of element neighbors . . . . . . . = ";
	if (used != 0)
		out << feNeighborData.Length()/used << '\n';
	else
		out << "-\n";

	out << " Element neighbor number distribution:\n";
	out << setw(kIntWidth) << "number"
	    << setw(kIntWidth) << "count" << '\n';
	for (int k = 0; k < counts.Length(); k++)
		out << setw(kIntWidth) << k
	    	<< setw(kIntWidth) << counts[k] << '\n';
	}

	/* memory requirements */
	int nsd = fCoords->MinorDim();
	
	out << "\n MLS storage requirements:\n";
	int n_count = fnNeighborCount.Sum();
	out << " Total number of nodal neighbors . . . . . . . . = " << n_count << '\n';
	out << " Nodal shape function storage. . . . . . . . . . = " << n_count*sizeof(double) << " bytes\n";
	out << " Nodal shape function derivatives storage. . . . = " << nsd*n_count*sizeof(double) << " bytes\n";
	int e_count = feNeighborCount.Sum();
	out << " Total number of integration point neighbors . . = " << e_count << '\n';
	out << " i.p. shape function storage . . . . . . . . . . = " << fIP*e_count*sizeof(double) << " bytes\n";
	out << " i.p. shape function derivatives storage . . . . = " << fIP*nsd*e_count*sizeof(double) << " bytes\n";

	/* search grid statistics */
	fGrid->WriteStatistics(out);
}

/* describe the parameters needed by the interface */
void MeshFreeSupportT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	ParameterT store_shape(ParameterT::Boolean, "store_shapefunctions");
	store_shape.SetDefault(true);
	list.AddParameter(store_shape);
}

/* information about subordinate parameter lists */
void MeshFreeSupportT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* window function */
	sub_list.AddSub("meshfree_formulation", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* MeshFreeSupportT::NewSub(const StringT& name) const
{
	if (name == "meshfree_formulation")
	{
		ParameterContainerT* window_choice = new ParameterContainerT(name);
		window_choice->SetSubSource(this);
		window_choice->SetListOrder(ParameterListT::Choice);
		window_choice->AddSub("EFG");
		window_choice->AddSub("RKPM");
		
		return window_choice;
	}
	else if (name == "EFG") 
	{
		ParameterContainerT* efg = new ParameterContainerT(name);
	
		ParameterT d_extra(ParameterT::Double, "support_scaling");
		d_extra.AddLimit(0.0, LimitT::Lower);
		d_extra.SetDefault(1.0);
		efg->AddParameter(d_extra);
		
		ParameterT completeness(ParameterT::Integer, "completeness");
		completeness.AddLimit(1, LimitT::LowerInclusive);
		completeness.SetDefault(1);
		efg->AddParameter(completeness);
	
		return efg;
	}
	else if (name == "RKPM")
	{
		ParameterContainerT* rkpm = new ParameterContainerT(name);
		rkpm->SetSubSource(this);
	
		ParameterT completeness(ParameterT::Integer, "completeness");
		completeness.AddLimit(1, LimitT::LowerInclusive);
		completeness.SetDefault(1);
		rkpm->AddParameter(completeness);
	
		/* window function choice */
		rkpm->AddSub("window_function_choice", ParameterListT::Once, true);
	
		return rkpm;
	}
	else if (name == "window_function_choice")
	{
		ParameterContainerT* window_choice = new ParameterContainerT("window_function_choice");
		window_choice->SetListOrder(ParameterListT::Choice);
		window_choice->SetSubSource(this);
		window_choice->AddSub("gaussian_window");
		window_choice->AddSub("rect_gaussian_window");
		window_choice->AddSub("cubic_spline_window");
		window_choice->AddSub("rect_cubic_spline_window");
		return window_choice;
	}
	else if (name == "gaussian_window")
	{
		ParameterContainerT* window = new ParameterContainerT(name);
		
		ParameterT support_scaling(ParameterT::Double, "support_scaling");
		support_scaling.AddLimit(0.0, LimitT::Lower);
		support_scaling.SetDefault(1.0);
		window->AddParameter(support_scaling);

		ParameterT sharpening_factor(ParameterT::Double, "sharpening_factor");
		sharpening_factor.AddLimit(0.0, LimitT::Lower);
		sharpening_factor.SetDefault(0.4);
		window->AddParameter(sharpening_factor);

		ParameterT cutoff_factor(ParameterT::Double, "cutoff_factor");
		cutoff_factor.AddLimit(0.0, LimitT::Lower);
		cutoff_factor.SetDefault(2.0);
		window->AddParameter(cutoff_factor);
		
		return window;
	}
	else if (name == "rect_gaussian_window")
	{
		ParameterContainerT* window = new ParameterContainerT(name);
		window->SetSubSource(this);
		window->SetDescription("1 or nsd sets of Gaussian window parameters");
		
		/* separate parameters for each coordinate direction */
		window->AddSub("gaussian_window", ParameterListT::OnePlus);

		return window;
	}
	else if (name == "cubic_spline_window")
	{
		ParameterContainerT* window = new ParameterContainerT(name);
		window->AddParameter(ParameterT::Double, "support_scaling");
		return window;
	}
	else if (name == "rect_cubic_spline_window")
	{
		ParameterContainerT* window = new ParameterContainerT(name);
		window->SetSubSource(this);
		window->SetDescription("1 or nsd sets of spline parameters");
		
		/* separate parameters for each coordinate direction */
		window->AddSub("cubic_spline_window", ParameterListT::OnePlus);

		return window;	
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void MeshFreeSupportT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MeshFreeSupportT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fStoreShape = list.Parameter("store_shapefunctions");
	
	/* coordinates */
	const dArray2DT& coordinates = NodalCoordinates();
	int nsd = coordinates.MinorDim();

	/* meshfree formulation */
	const ParameterListT& formulation = list.GetListChoice(*this, "meshfree_formulation");
	if (formulation.Name() == "EFG") {

		/* meshfree formulation flag */
		fMeshfreeType = kEFG;
	
		/* EFG parameters */
		fDextra = formulation.GetParameter("support_scaling");
		int completeness = formulation.GetParameter("completeness");

		/* construct MLS solver */	
		if (nsd == 2)
			fEFG = new OrthoMLS2DT(completeness);
		else if (nsd == 3)
			fEFG = new OrthoMLS3DT(completeness);
		else
			ExceptionT::GeneralFail(caller, "%dD not supported", nsd); 
		
		/* initialize */
		fEFG->Initialize();
		
		/* just one nodal field parameter */
		fnodal_param_man.SetMinorDimension(1);
		fNodalParameters.Dimension(fCoords->MajorDim(), 1);
		fNodalParameters = -1;	
	}
	else if (formulation.Name() == "RKPM") {

		/* meshfree formulation flag */
		fMeshfreeType = kRKPM;

		/* construct MLS solver */
		int completeness = formulation.GetParameter("completeness");
		const ParameterListT& window = formulation.GetListChoice(*this, "window_function_choice");		
		fRKPM = New_MLSSolverT(nsd, completeness, window);

		/* initialize */
		fRKPM->Initialize();

		/* dimension nodal field parameter */
		fnodal_param_man.SetMinorDimension(fRKPM->NumberOfSupportParameters());
		
		/* allocate work space */
		fVolume.Dimension(fCoords->MajorDim());
		fVolume = 1.0;
		fNodalParameters.Dimension(fCoords->MajorDim(), fRKPM->NumberOfSupportParameters());
		fNodalParameters = -1;
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized meshfree formulation \"%s\"",
			formulation.Name().Pointer());

	/* memory manager */
	fnodal_param_man.Register(fnodal_param   );
	fnodal_param_man.Register(fnodal_param_ip);
}

/* construct a new MLSSolverT with the given parameters */
MLSSolverT* MeshFreeSupportT::New_MLSSolverT(int nsd, int completeness, const ParameterListT& window)
{
	const char caller[] = "MeshFreeSupportT::New_MLSSolverT";

	WindowTypeT window_type;
	dArrayT window_params;
	if (window.Name() == "gaussian_window") 
	{
		window_type = kGaussian;
		window_params.Dimension(3);
		window_params[0] = window.GetParameter("support_scaling");
		window_params[1] = window.GetParameter("sharpening_factor");
		window_params[2] = window.GetParameter("cutoff_factor");
	}
	else if (window.Name() == "rect_gaussian_window")
	{
		window_type = kBrick;
		int num_params = window.NumLists("gaussian_window");
		if (num_params != 1 && num_params != nsd)
			ExceptionT::GeneralFail(caller, "expecting 1 or %d \"gaussian_window\" in \"rect_gaussian_window\"", nsd);

		window_params.Dimension(3*nsd);
		int index = 0;
		for (int i = 0; i < nsd; i++) {
			const ParameterListT& param = window.GetList("gaussian_window", (num_params == 1) ? 0 : i);
			window_params[index++] = param.GetParameter("support_scaling");
			window_params[index++] = param.GetParameter("sharpening_factor");
			window_params[index++] = param.GetParameter("cutoff_factor");
		}		
	}
	else if (window.Name() == "cubic_spline_window")
	{
		window_type = kCubicSpline;
		window_params.Dimension(1);
		window_params[0] = window.GetParameter("support_scaling");
	}		
	else if (window.Name() == "rect_cubic_spline_window")
	{
		window_type = kRectCubicSpline;
		int num_params = window.NumLists("cubic_spline_window");
		if (num_params != 1 && num_params != nsd)
			ExceptionT::GeneralFail(caller, "expecting 1 or %d \"cubic_spline_window\" in \"rect_cubic_spline_window\"", nsd);

		window_params.Dimension(nsd);
		for (int i = 0; i < nsd; i++) {
			const ParameterListT& param = window.GetList("cubic_spline_window", (num_params == 1) ? 0 : i);
			window_params[i] = param.GetParameter("support_scaling");
		}
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized window function choice \"%s\"",
			window.Name().Pointer());

	/* construct MLS solver */
	return new MLSSolverT(nsd, completeness, window_type, window_params);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* initialize search grid */
void MeshFreeSupportT::SetSearchGrid(void)
{
	const char caller[] = "MeshFreeSupportT::SetSearchGrid";

	/* must set nodes used first */
	if (fNodesUsed.Length() == 0)
		ExceptionT::GeneralFail(caller, "must set used nodes first");

	/* try to get roughly least 10 per grid */
	int ngrid = int(pow(fNodesUsed.Length()/10.0, 1.0/fSD)) + 1;
	ngrid = (ngrid < 2) ? 2 : ngrid;
	ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

	/* free any existing */
	delete fGrid;

	/* construct a search grid */
	iArrayT n_grid(fSD);
	n_grid = ngrid;
	fGrid = new iGridManagerT(n_grid, *fCoords, &fNodesUsed);	
	if (!fGrid) ExceptionT::OutOfMemory(caller);

	/* reset search grid bounds */
	fGrid->Reset();
}

/* generate lists of all nodes that fall within support of the
* nodal coords (self included) */
void MeshFreeSupportT::SetNodeNeighborData(const dArray2DT& coords)
{
	const char caller[] = "MeshFreeSupportT::SetNodeNeighborData";
#if __option(extended_errorcheck)
	if (coords.MajorDim() != fNodalParameters.MajorDim())
		ExceptionT::SizeMismatch(caller);
#endif

	int numnodes = fNodalParameters.MajorDim();
	int nsd      = coords.MinorDim();
	
	/* initialize neighbor counts */
	fnNeighborCount.Dimension(numnodes);
	fnNeighborCount = 0;

	/* work space */
	iArray2DT a1(numnodes, kListInitSize), a2;
	iArray2DT* pthis = &a1;
	iArray2DT* pnext = &a2;

	/* determine search type and size */
	dArrayT support_size;
	dArray2DT rect_support_size;
	WindowT::SearchTypeT search_type = WindowT::kNone;
	if (fMeshfreeType == kEFG) {
		search_type = WindowT::kSpherical;
		support_size.Dimension(fNodalParameters.MajorDim());
		support_size.SetToScaled(fDextra, fNodalParameters);
	}
	else if (fMeshfreeType == kRKPM) {
		search_type = fRKPM->SearchType();
		if (search_type == WindowT::kSpherical) {
			support_size.Dimension(fNodalParameters.MajorDim());
			fRKPM->SphericalSupportSize(fNodalParameters, support_size);
		}
		else {
			rect_support_size.Dimension(fNodalParameters.MajorDim(), nsd);
			fRKPM->RectangularSupportSize(fNodalParameters, rect_support_size);		
		}
	}

	/* loop over "active" nodes */
	dArrayT nodal_params;
	AutoArrayT<int> neighbors;
	for (int ndex = 0; ndex < fNodesUsed.Length(); ndex++)
	{
		int i = fNodesUsed[ndex];
		const double* target = coords(i);
		
		/* find nodes in neighborhood */
		if (search_type == WindowT::kSpherical)
			fGrid->Neighbors(i, support_size[i], neighbors);
		else {
			rect_support_size.RowAlias(i, nodal_params);
			fGrid->Neighbors(i, nodal_params, neighbors);
		}		

		/* add self to the list */
		neighbors.Append(i);
			
		/* need more space */
		int count = neighbors.Length();
		if (count > (*pthis).MinorDim())
		{
			(*pnext).Dimension(numnodes, count + kListSizeinc);
					
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
	fndShapespace.Dimension(maxsize*(1 + nsd));
}

void MeshFreeSupportT::SetNodeNeighborData_2(const dArray2DT& coords)
{
	const char caller[] = "MeshFreeSupportT::SetNodeNeighborData_2";
#if __option(extended_errorcheck)
	if (coords.MajorDim() != fNodalParameters.MajorDim()) 
		ExceptionT::SizeMismatch(caller);
#endif

	int numnodes = fNodalParameters.MajorDim();
	int nsd      = coords.MinorDim();
	
	/* initialize neighbor counts */
	fnNeighborCount.Dimension(numnodes);
	fnNeighborCount = 0;

	/* work space */
	iArrayT allocator;
	ArrayT<int*> pointers(numnodes);
	pointers = NULL;

	/* determine search type and size */
	dArrayT support_size;
	dArray2DT rect_support_size;
	WindowT::SearchTypeT search_type = WindowT::kNone;
	if (fMeshfreeType == kEFG) {
		search_type = WindowT::kSpherical;
		support_size.Dimension(fNodalParameters.MajorDim());
		support_size.SetToScaled(fDextra, fNodalParameters);
	}
	else if (fMeshfreeType == kRKPM) {
		search_type = fRKPM->SearchType();
		if (search_type == WindowT::kSpherical) {
			support_size.Dimension(fNodalParameters.MajorDim());
			fRKPM->SphericalSupportSize(fNodalParameters, support_size);
		}
		else {
			rect_support_size.Dimension(fNodalParameters.MajorDim(), nsd);
			fRKPM->RectangularSupportSize(fNodalParameters, rect_support_size);		
		}
	}

	/* loop over "active" nodes */
	dArrayT nodal_params;
	AutoArrayT<int> neighbors;
	for (int ndex = 0; ndex < fNodesUsed.Length(); ndex++)
	{
		int i = fNodesUsed[ndex];
		const double* target = coords(i);
		
		/* find nodes in neighborhood */
		if (search_type == WindowT::kSpherical)
			fGrid->Neighbors(i, support_size[i], neighbors);
		else
		{
			rect_support_size.RowAlias(i, nodal_params);
			fGrid->Neighbors(i, nodal_params, neighbors);
		}		

		/* add self to the list */
		neighbors.Append(i);
		int count = neighbors.Length();

		/* store count */
		fnNeighborCount[i] = count;

		/* store neighbors */
		allocator.Dimension(count);
		neighbors.CopyInto(allocator);
		allocator.ReleasePointer(pointers.Pointer(i));
		allocator.Free();
	}					

	/* use disk for bigger problems */
#if defined(_USE_DISK_) && !defined(__CPLANT__)
	bool use_disk = (numnodes > 250000);
#else
	bool use_disk = false;
#endif

	/* write data to disk before allocating new memory */
	if (use_disk)
	{
		StringT file = "nodesupport.tmp";
		cout << "\n MeshFreeSupportT::SetNodeNeighborData_2: system is big. flushing\n"
		     <<   "     node neighbor data to disk file: \"" << file << '\"' << endl;

		/* write data to disk */
		ArrayT<int> dump;
		ofstreamT tmp_out(file);
		if (!tmp_out.is_open()) throw ExceptionT::kGeneralFail;
		for (int j = 0; j < numnodes; j++)
		{
			dump.Set(fnNeighborCount[j], pointers[j]);
			dump.WriteBinary(tmp_out);
		}
		tmp_out.close();

		/* free all temp space */
		for (int k = 0; k < numnodes; k++)
			delete[] pointers[k];
		pointers.Free();

		/* configure node neighbor data */
		fnNeighborData.Configure(fnNeighborCount);

		/* read contents of */
		ifstreamT tmp_in(file);
		if (!tmp_in.is_open()) throw ExceptionT::kGeneralFail;
		fnNeighborData.ReadDataBinary(tmp_in);
		tmp_in.close();
	}
	else
	{
		/* configure node neighbor data */
		fnNeighborData.Configure(fnNeighborCount);
	
		/* copy data in */
		for (int j = 0; j < numnodes; j++)
			fnNeighborData.SetRow(j, pointers[j]);

		/* free all temp space */
		for (int k = 0; k < numnodes; k++)
			delete[] pointers[k];
	}
		
	/* space for nodal calculations */
	int maxsize = fnNeighborData.MaxMinorDim();
	fndShapespace.Dimension(maxsize*(1 + nsd));
}

/* generate lists of all nodes that fall within Dmax of the
* element integration points */
void MeshFreeSupportT::SetElementNeighborData(const iArray2DT& connects)
{
	/* verify that this routine should have been called */
	if (!fDomain)
		ExceptionT::GeneralFail("MeshFreeSupportT::SetElementNeighborData","No Domain available\n");


	/* dimensions */
	int numelems = connects.MajorDim();
	int nen      = connects.MinorDim();
		
	/* initialize neighbor counts */
	feNeighborCount.Dimension(numelems);
	feNeighborCount = 0;

	/* work space */
	iArray2DT a1(numelems, kListInitSize), a2;
	iArray2DT* pthis = &a1;
	iArray2DT* pnext = &a2;
	
	/* domain coords */
	iArrayT     elementnodes;
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, fSD);
	loccoords.SetGlobal(*fCoords);

	/* integration point coordinate tables */	
	dArray2DT x_ip_table(fIP, fSD);
	dArrayT   x_ip, x_node, nodal_params;

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
		if (big_system && ((i + 1) % big) == 0)
		{
			cout << setw(2*kIntWidth) << i + 1
			     << setw(kIntWidth) << min_count
			     << setw(kIntWidth) << max_count
			     << setw(kDoubleWidth) << double(sum_count)/big <<endl;

			/* reset */
			max_count = -1, min_count = -1, sum_count = 0;
		}

		/* connectivity */
		const int* pelem = connects(i);
		
		/* integration point coordinates */
		fConnects->RowAlias(i, elementnodes);
		loccoords.SetLocal(elementnodes);
		fDomain->Interpolate(loccoords, x_ip_table);
		
		/* collect unique list of neighboring nodes */
		nodeset.Dimension(0);
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
					fCoords->RowAlias(neigh, x_node);
				
					/* loop over integration points */
					int hit = 0;
					for (int k = 0; k < fIP && !hit; k++)
					{
						/* fetch ip coordinates */
						x_ip_table.RowAlias(k, x_ip);
						if (Covers(x_ip, x_node, neigh)) hit = 1;
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
			(*pnext).Dimension(numelems, newsize);

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
	felShapespace.Dimension(fIP*maxsize*(1 + fSD));
}

/* generate lists of all nodes that fall within Dmax of the
* element integration points */
void MeshFreeSupportT::SetElementNeighborData_2(const iArray2DT& connects)
{
	/* verify that this routine should have been called */
	if (!fDomain)
		ExceptionT::GeneralFail("MeshFreeSupportT::SetElementNeighborData_2","No Domain available\n");


	/* dimensions */
	int numelems = connects.MajorDim();
	int nen      = connects.MinorDim();
	int nsd      = fCoords->MinorDim();
		
	/* initialize neighbor counts */
	feNeighborCount.Dimension(numelems);
	feNeighborCount = 0;

	/* work space */
	iArrayT allocator;
	ArrayT<int*> pointers(numelems);
	pointers = NULL;

	/* domain coords */
	iArrayT     elementnodes;
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, fSD);
	loccoords.SetGlobal(*fCoords);

	/* integration point coordinate tables */	
	dArray2DT x_ip_table(fIP, fSD);
	dArrayT   x_ip, x_node, nodal_params;

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
		if (big_system && ((i + 1) % big) == 0)
		{
			cout << setw(2*kIntWidth) << i + 1
			     << setw(kIntWidth) << min_count
			     << setw(kIntWidth) << max_count
			     << setw(kDoubleWidth) << double(sum_count)/big <<endl;

			/* reset */
			max_count = -1, min_count = -1, sum_count = 0;
		}

		/* connectivity */
		const int* pelem = connects(i);
		
		/* integration point coordinates */
		fConnects->RowAlias(i, elementnodes);
		loccoords.SetLocal(elementnodes);
		fDomain->Interpolate(loccoords, x_ip_table);
		
		/* collect unique list of neighboring nodes */
		nodeset.Dimension(0);
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
					fCoords->RowAlias(neigh, x_node);
				
					/* loop over integration points */
					int hit = 0;
					for (int k = 0; k < fIP && !hit; k++)
					{
						/* fetch ip coordinates */
						x_ip_table.RowAlias(k, x_ip);
						if (Covers(x_ip, x_node, neigh)) hit = 1;
					}
				
					/* within the neighborhood */	
					if (hit) nodeset.Append(neigh);
				}
			}
		}
		int setlength = nodeset.Length();

		/* statistics for "big" systems */
		if (big_system)
		{
			if (max_count == -1 || setlength > max_count)
				max_count = setlength;			
			if (min_count == -1 || setlength < min_count)
				min_count = setlength;
			sum_count += setlength;
		}
		
		/* store count */
		feNeighborCount[i] = setlength;

		/* store neighbors */
		allocator.Dimension(setlength);
		nodeset.CopyInto(allocator);
		allocator.ReleasePointer(pointers.Pointer(i));
		allocator.Free();
	}

	/* use disk for bigger problems */
#if defined(_USE_DISK_) && !defined(__CPLANT__)
	bool use_disk = (numelems > 250000);
#else
	bool use_disk = false;
#endif

	/* write data to disk before allocating new memory */
	if (use_disk)
	{
		StringT file = "elemsupport.tmp";
		cout << "\n MeshFreeSupportT::SetElementNeighborData_2: system is big. flushing\n"
		     <<   "     element neighbor data to disk file: \"" << file << '\"' << endl;

		/* write data to disk */
		ArrayT<int> dump;
		ofstreamT tmp_out(file);
		if (!tmp_out.is_open()) throw ExceptionT::kGeneralFail;
		for (int j = 0; j < numelems; j++)
		{
			dump.Set(feNeighborCount[j], pointers[j]);
			dump.WriteBinary(tmp_out);
		}
		tmp_out.close();

		/* free all temp space */
		for (int k = 0; k < numelems; k++)
			delete[] pointers[k];
		pointers.Free();

		/* configure element neighbor data */
		feNeighborData.Configure(feNeighborCount);

		/* read contents of */
		ifstreamT tmp_in(file);
		if (!tmp_in.is_open()) throw ExceptionT::kGeneralFail;
		feNeighborData.ReadDataBinary(tmp_in);
		tmp_in.close();
	}
	else
	{
		/* configure element neighbor data */
		feNeighborData.Configure(feNeighborCount);

		/* copy data in */
		for (int j = 0; j < numelems; j++)
			feNeighborData.SetRow(j, pointers[j]);
		
		/* free all temp space */
		for (int k = 0; k < numelems; k++)
			delete[] pointers[k];
	}
		
	/* space element calculation */
	int maxsize = feNeighborData.MaxMinorDim();
	felShapespace.Dimension(fIP*maxsize*(1 + fSD));
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
	fResetNodes.Dimension(0);
}

/* compute all integration point shape functions and derivatives */
void MeshFreeSupportT::SetElementShapeFunctions(void)
{
	/* initialize database */
	InitElementShapeData();

	/* dimensions */
	int nel = fConnects->MajorDim();

	/* work space */
	iArrayT    neighbors;
	dArray2DT phi;
	ArrayT<dArray2DT> Dphi(fIP);

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
	fResetElems.Dimension(0);	
}

/*************************************************************************
* Private
*************************************************************************/

/* test coverage - temporary function until OrthoMLSSolverT class
 * is brought up to date */
bool MeshFreeSupportT::Covers(const dArrayT& field_x, const dArrayT& node_x, 
	int node) const
{
	bool covers = false;
	if (fMeshfreeType == kEFG)
	{
		double dist = 0.0;
		for (int j = 0; j < field_x.Length(); j++)
		{
			double dx = node_x[j] - field_x[j];
			dist += dx*dx;
		}
				
		/* covers */
		double dmax_i = fNodalParameters(node, 0);
		if (dist < dmax_i*dmax_i) covers = true;
	}
	else if (fMeshfreeType == kRKPM)
	{
		dArrayT nodal_params;
		fNodalParameters.RowAlias(node, nodal_params);
		covers = fRKPM->Covers(node_x, field_x, nodal_params);
	}
	else
	{
		cout << "\n MeshFreeSupportT::Covers: unexpected meshfree type: " 
		     << fMeshfreeType << endl;
		throw ExceptionT::kGeneralFail;
	}
	return covers;
}

void MeshFreeSupportT::ComputeElementData(int element, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi)
{
	/* Verify that there is an element */
	if (!fDomain)
		ExceptionT::GeneralFail("MeshFreeSupportT::ComputeElementData","No Domain available\n");

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
	LocalArrayT loccoords(LocalArrayT::kUnspecified, nen, fSD);
	loccoords.SetGlobal(*fCoords);
		
	/* integration point coordinates */
	fConnects->RowAlias(element, elementnodes);
	loccoords.SetLocal(elementnodes);
	fDomain->Interpolate(loccoords, fx_ip_table);

	/* loop over integration points */
	dArrayT x_ip;
	for (int i = 0; i < fIP; i++)
	{
		/* fetch ip coordinates */
		fx_ip_table.RowAlias(i, x_ip);

		/* process boundaries */
		fnodal_param_ip = fnodal_param;
		ProcessBoundaries(fcoords, x_ip, fnodal_param_ip);
		// set dmax = -1 for nodes that are inactive at x_node
	
		/* compute MLS field */
		int OK;
		if (fMeshfreeType == kEFG)
		{
			OK = fEFG->SetField(fcoords, fnodal_param_ip, x_ip);
	
			/* store field data */
			phi.SetRow(i, fEFG->phi());
			Dphi[i] = fEFG->Dphi();
		}
		else if (fMeshfreeType == kRKPM)
		{
			fvolume.Collect(neighbors, fVolume);
			OK = fRKPM->SetField(fcoords, fnodal_param_ip, fvolume, x_ip, 1);

			/* store field data */
			phi.SetRow(i, fRKPM->phi());
			Dphi[i] = fRKPM->Dphi();
		}
		else throw ExceptionT::kGeneralFail;

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
			     << setw(  d_width) << "dist"
			     << setw(fnodal_param_ip.MinorDim()*d_width) << "nodal parameters"
			     << setw(fcoords.MinorDim()*d_width) << "x" << '\n';
			dArrayT dist(x_ip.Length());
			for (int j = 0; j < neighbors.Length(); j++)
			{
				err << setw(kIntWidth) << neighbors[j] + 1;
				fcoords.RowCopy(j, dist);
				dist -= x_ip;		
				err << setw(  d_width) << dist.Magnitude();
				fnodal_param_ip.PrintRow(j, err);
				fcoords.PrintRow(j, err);				
			}
			err.flush();
			throw ExceptionT::kGeneralFail;	
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
	fnodal_param_man.SetMajorDimension(count, false);
	fcoords_man.SetMajorDimension(count, false);
	
	/* collect local lists */
	fcoords.RowCollect(neighbors, *fCoords);
	fnodal_param.RowCollect(neighbors, fNodalParameters);
		
	/* coords of current node */
	dArrayT x_node;
	fCoords->RowAlias(node, x_node);
	
	/* process boundaries */
	ProcessBoundaries(fcoords, x_node, fnodal_param);
	// set dmax = -1 for nodes that are inactive at x_node

	/* compute MLS field */
	int OK;
	if (fEFG)
	{
		OK = fEFG->SetField(fcoords, fnodal_param, x_node);
		
		/* copy field data */
		phi  = fEFG->phi();
		Dphi = fEFG->Dphi();
	}
	else
	{
		/* dimension */
		fvolume_man.SetLength(count, false);
		fvolume.Collect(neighbors, fVolume);
		
		/* set field */
		OK = fRKPM->SetField(fcoords, fnodal_param, fvolume, x_node, 1);
		
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
		     << setw(  d_width) << "dist"
		     << setw(fnodal_param.MinorDim()*d_width) << "dmax"
		     << setw(fcoords.MinorDim()*d_width) << "x" << '\n';
		dArrayT dist(x_node.Length());
		for (int i = 0; i < neighbors.Length(); i++)
		{
			cout << setw(kIntWidth) << neighbors[i] + 1;
			fcoords.RowCopy(i, dist);
			dist -= x_node;		
			cout << setw(  d_width) << dist.Magnitude();
			fnodal_param.PrintRow(i, cout);
			fcoords.PrintRow(i, cout);
		}
		cout.flush();		
		throw ExceptionT::kGeneralFail;	
	}
}

/* set support for each node in the connectivities */
void MeshFreeSupportT::SetSupport_Spherical_Search(void)
{
	const char caller[] = "MeshFreeSupportT::SetSupport_Spherical_Search";

	/* dimensions */
	int nnd = fCoords->MajorDim();
	int nsd = fCoords->MinorDim();

	/* MLS solver specific check */
	if (fRKPM && fRKPM->NumberOfSupportParameters() != 1)
		ExceptionT::GeneralFail(caller, "expecting 1 support parameter not %d", 
			fRKPM->NumberOfSupportParameters());

	int min_neighbors = (fEFG) ? fEFG->NumberOfMonomials() :
	                            fRKPM->BasisDimension();
	//TEMP
	min_neighbors += 1; // don't count self

	/* quick check */
	if (fNodesUsed.Length() < min_neighbors)
		ExceptionT::GeneralFail(caller, "not enough meshfree nodes");

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
	iArrayT nodes_sh;
	double dmax = 0.0;
	for (int ndex = 0; ndex < fNodesUsed.Length(); ndex++)
	{
		/* "big" system */
		if (big_system && ((ndex + 1) % big) == 0)
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
		const double* target = (*fCoords)(i);

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
			dists.Dimension(0);
			nodes.Dimension(0);
		
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
			fCoords->PrintRow(i, cout);
			throw ExceptionT::kGeneralFail;
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
				double& d_i = fNodalParameters(nodes[i],0);
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
}

/* set the support for each node in the connectivities */
void MeshFreeSupportT::SetSupport_Spherical_Connectivities(void)
{
// NOTE: this version of finding Dmax only works if all the nodes
//       are members of the connectivities defining the integration
//       grid cells. The Dmax computed here is also strictly not what
//       its definition says: the distance to the farthest node in the
//       smallest neighborhood that's required for the MLS fit.
	
	/* dimensions */
	int nnd = fCoords->MajorDim();
	int nsd = fCoords->MinorDim();
	int nel = fConnects->MajorDim();
	int nen = fConnects->MinorDim();

	/* loop over elements */
	dArrayT x_0, x_i;
	for (int i = 0; i < nel; i++)
	{
		const int* pelem = (*fConnects)(i);
		for (int j = 0; j < nen; j++)
		{
			/* current node */
			int node = pelem[j];

			/* fetch origin coords */
			fCoords->RowAlias(node, x_0);
			
			/* find max neighbor distance */
			double& maxdist = fNodalParameters(node,0);
			for (int k = 0; k < nen; k++)
				if (k != j)
				{
					/* fetch origin coords */
					fCoords->RowAlias(pelem[k], x_i);
			
					/* separating distance */
					double dist = dArrayT::Distance(x_i, x_0);
					
					/* keep max */
					maxdist = Max(dist, maxdist);
				}
		}
	}
}

/* set the support for each node in the connectivities */
void MeshFreeSupportT::SetSupport_Cartesian_Connectivities(void)
{
	/* dimensions */
	int nnd = fCoords->MajorDim();
	int nsd = fCoords->MinorDim();
	int nel = fConnects->MajorDim();
	int nen = fConnects->MinorDim();

	/* loop over elements */
	dArrayT x_0, x_i, r(nsd), support;
	for (int i = 0; i < nel; i++)
	{
		const int* pelem = (*fConnects)(i);
		for (int j = 0; j < nen; j++)
		{
			/* current node */
			int node = pelem[j];

			/* fetch origin coords */
			fCoords->RowAlias(node, x_0);
			
			/* find max neighbor distance */
			for (int k = 0; k < nen; k++)
				if (k != j) /* skip self */
				{
					/* fetch node coords */
					fCoords->RowAlias(pelem[k], x_i);

					/* check visibility */
					if (Visible(x_0.Pointer(), x_i.Pointer()))
					{
						/* fetch support size */
						fNodalParameters.RowAlias(node, support);

						/* separation */
						r.DiffOf(x_i, x_0);
					
						/* keep max in each direction */
						for (int l = 0; l < nsd; l++)
							if (support[l] < 0.0)
								support[l] = fabs(r[l]);
							else
								support[l] = AbsMax(r[l], support[l]);
					}
				}
		}
	}
}

void MeshFreeSupportT::SetNodesUsed(void)
{
	/* dimensions */
	int nnd = fCoords->MajorDim();

	/* markers */
	iArrayT used(nnd);
	used = 0;

	/* mark nodes used on the integration grid */
	int  tot = fConnects->Length();
	const int* pnd = fConnects->Pointer();
	for (int i = 0; i < tot; i++)
		used[*pnd++] = 1;
	
	/* mark other EFG nodes */
	tot = fNonGridNodes->Length();
	pnd = fNonGridNodes->Pointer();
	for (int j = 0; j < tot; j++)
		used[*pnd++] = 1;

	/* count nodes used */
	int numused = used.Count(1);
	fNodesUsed.Dimension(numused);
		
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
	const double* target = x.Pointer();
	const AutoArrayT<iNodeT>* inodes = &fGrid->HitsInRegion(target, ++cell_span);
	while (inodes->Length() < 1 && cell_span <= 3)
		inodes = &fGrid->HitsInRegion(target, ++cell_span);
	if (inodes->Length() < 1)
	{
		cout << "\n MeshFreeSupportT::BuildNeighborhood: failed to find any nodes around:\n";
		cout << x << "\n     after 4 iterations" << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* collect support nodes */
	if (fMeshfreeType == kEFG || fRKPM->SearchType() == WindowT::kSpherical)
	{
		/* find biggest support */
		double support_max = 0.0;
		for (int ii = 0; ii < inodes->Length(); ii++)
		{
			int tag = ((*inodes)[ii]).Tag();
			double dmax = fNodalParameters(tag,0);
			support_max = (dmax > support_max) ? dmax : support_max;
		}
	
		/* need to re-collect nodes */
		if (fGrid->CellSpan(cell_span) < support_max)
			inodes = &fGrid->HitsInRegion(target, support_max);
	}
	else if (fRKPM->SearchType() == WindowT::kConnectivity)
	{
		/* find biggest support */
		dArrayT support_max(x.Length()), support;
		support_max = 0.0;
		for (int ii = 0; ii < inodes->Length(); ii++)
		{
			int tag = ((*inodes)[ii]).Tag();
			fNodalParameters.RowAlias(tag, support);
			for (int i = 0; i < support.Length(); i++)
				support_max[i] = Max(support_max[i], support[i]);
		}
	
		/* re-collect using max support */
		inodes = &fGrid->HitsInRegion(target, support_max);
	}
	else throw ExceptionT::kGeneralFail;

	/* work space */	
	int nsd = fCoords->MinorDim();
	AutoArrayT<int> out_nodes; // make as class workspace?

	/* initialize */
	nodes.Dimension(0);
	out_nodes.Dimension(0);
				
	/* run through nodes */
	dArrayT x_node, nodal_params;
	for (int i = 0; i < inodes->Length(); i++)
	{
		int tag = ((*inodes)[i]).Tag();
		fCoords->RowAlias(tag, x_node);

//NOTE - reverse this. finding the r^2 is probably cheaper than
//       determining whether the node is visible
	
		/* check node */
		if (!nodes.HasValue(tag) &&     // redundant?
		    !out_nodes.HasValue(tag) && // redundant?
		     Visible(target, x_node.Pointer()))
		{
			/* add to list */
			if (Covers(x, x_node, tag))
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
		     << setw(fNodalParameters.MinorDim()*d_width) << "dmax"
		     << setw(x.Length()*d_width) << "x" << '\n';

		/* work space */
		dArrayT xn, nodal_params;
		
		/* write nodes */
		for (int i = 0; i < inodes->Length(); i++)
		{
			int tag = ((*inodes)[i]).Tag();
			fCoords->RowAlias(tag, xn);
			double dist = dArrayT::Distance(x, xn);

			cout << setw(kIntWidth) << tag+1
		         << setw(5) << nodes.HasValue(tag)
		         << setw(5) << Visible(x.Pointer(), xn.Pointer())
		         << setw(5) << Covers(x, xn, tag)
		         << setw(  d_width) << dist;
			fNodalParameters.PrintRow(tag, cout);
			cout << xn.no_wrap() << '\n';		
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
	int nsd = fCoords->MinorDim();

	/* configure nodal storage */
	int step;
	try {
		step = 0;
		fnPhiData.Configure(fnNeighborCount);
		step = 1;
		fnDPhiData.Configure(fnNeighborCount, nsd);
	}
	
	/* write memory information */
	catch (ExceptionT::CodeT error)
	{
		if (error == ExceptionT::kOutOfMemory)
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
	int nsd = fCoords->MinorDim();
	int nip = fIP;

	/* configure element storage */
	int step;
	try {
		step = 0;
		fePhiData.Configure(feNeighborCount, nip);
		step = 1;
		feDPhiData.Configure(feNeighborCount, nip*nsd);
	}

	/* write memory information */
	catch (ExceptionT::CodeT error)
	{
		if (error == ExceptionT::kOutOfMemory)
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
	if ((**pfrom).MajorDim() != (**pto).MajorDim()) throw ExceptionT::kSizeMismatch;
	if ((**pfrom).MinorDim() >= (**pto).MinorDim()) throw ExceptionT::kGeneralFail;
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
