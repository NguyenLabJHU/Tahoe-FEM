/* $Id: SamplingSurfaceT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (10/19/2000)                                          */

#include "SamplingSurfaceT.h"
#include "SurfaceShapeT.h"
#include "MeshFreeSupportT.h"

/* constructor */
SamplingSurfaceT::SamplingSurfaceT(GeometryT::CodeT code, int num_facet_nodes,
	int num_samples, MeshFreeSupportT& mf_support):
	fCode(code),
	fNumFacetNodes(num_facet_nodes),
	fNumSamples(num_samples),
	fMFSupport(mf_support),
	fLocFacetCoords(LocalArrayT::kInitCoords, fNumFacetNodes,
		GeometryT::GeometryToNumSD(fCode)+1)
{
	/* set surface geometry */
	fSurfaceShape = new SurfaceShapeT(fCode, fNumSamples, 2*fNumFacetNodes,
		fLocFacetCoords.MinorDim(), fLocFacetCoords);
	if (!fSurfaceShape) throw eOutOfMemory;
	fSurfaceShape->Initialize();
	
	/* work space */
	fQ.Allocate(fLocFacetCoords.MinorDim());
}

/* destructor */
SamplingSurfaceT::~SamplingSurfaceT(void)
{
	delete fSurfaceShape;
	fSurfaceShape = NULL;
}

/* compute/store shape functions at sampling points on surface */
void SamplingSurfaceT::SetSamplingPoints(const dArray2DT& facet_coords,
	const iArray2DT& facet_connects)
{
	/* copy data */
	fFacetCoords = facet_coords;
	fFacetConnects = facet_connects;

//NOTE: could also
// (1) fFacetCoords.Swap(facet_coords)
// (2) facet_coords.ShallowCopy(fFacetCoords);
// -> less memory movement

	/* dimensions */
	fNumFacets = fFacetConnects.MajorDim();
	if (fNumFacets == 0) return;

	/* checks */
	if (fFacetCoords.MinorDim() != fSurfaceShape->NumSD() + 1)
	{
		cout << "\n SamplingSurfaceT::SetSamplingPoints: dimension of coordinates\n"
		     <<   "     " << fFacetCoords.MinorDim()
		     << " does not match dimension of surface " << fSurfaceShape->NumSD() + 1
		     << endl;
		throw eGeneralFail;
	}
	if (fFacetConnects.MinorDim() != fNumFacetNodes)
	{
		cout << "\n SamplingSurfaceT::SetSamplingPoints: dimension of facet connectivites\n"
		     <<   "     " << fFacetConnects.MinorDim()
		     << " does not match dimension of surface " << fNumFacetNodes << endl;
		throw eGeneralFail;
	}	

	/* set source of facet coordinates */
	fLocFacetCoords.SetGlobal(fFacetCoords);
	
	/* set neighbor data */
	AutoArrayT<int> all_neighbors(20, true);
	AutoArrayT<int> neighbors;
	iArrayT facet_nodes;
	iArrayT neighbor_count(fNumFacets*fNumSamples);
	int dex = 0;
	for (int i = 0; i < fNumFacets; i++)
	{
		/* set facet coordinates */
		fFacetConnects.RowAlias(i, facet_nodes);
		fLocFacetCoords.SetLocal(facet_nodes);
	
		/* loop over integration points */
		fSurfaceShape->TopIP();
		while (fSurfaceShape->NextIP())
		{
			/* field point coordinates */
			const dArrayT& coords = fSurfaceShape->IPCoords();
		
			/* collect neighborhood nodes */
			if (!fMFSupport.BuildNeighborhood(coords, neighbors))
			{
				cout << "\n SamplingSurfaceT::SetSamplingPoints: BuildNeighborhood: failed" << endl;
				throw eGeneralFail;
			}		
			neighbor_count[dex] = neighbors.Length();
			all_neighbors.Append(neighbors);
			
			/* next point */
			dex++;
		}
	}
	
	/* dimension database */
	fNeighborData.Configure(neighbor_count);
	const int* p_neighbors = all_neighbors.Pointer();
	for (int j = 0; j < neighbor_count.Length(); j++)
	{
		int dim = neighbor_count[j];
		fNeighborData.SetRow(j, p_neighbors);
		p_neighbors += dim;
	}
	fPhiData.Configure(neighbor_count);
	fDPhiData.Configure(neighbor_count, facet_coords.MinorDim());
	
	/* compute/store ALL shape function data */
	SetFieldData(NULL);
	
	/* allocate flags */
	fFlag.Allocate(fNumFacets);
	fFlag = 0;
}

/* reset database for facets affected by given nodes */
void SamplingSurfaceT::ResetSamplingPoints(const iArrayT& changed_nodes)
{
	/* node number range */
	int min, max;
	changed_nodes.MinMax(min, max);

	/* check facets */
	AutoArrayT<int> reset_facets;
	iArrayT neighbors;
	for (int i = 0; i < fNumFacets; i++)
	{
		int dex = i*fNumSamples;
		bool hit = false;
		for (int j = 0; !hit && j < fNumSamples; j++)
		{
			fNeighborData.RowAlias(dex, neighbors);
			for (int k = 0; !hit && k < neighbors.Length(); k++)
			{
				int node = neighbors[k];
				if (node >= min &&
				    node <= max &&
				    changed_nodes.HasValue(node))
				{
					reset_facets.Append(i);
					hit = true;
				}
				
				/* next point */
				dex++;
			}
		}
	}
	
	/* selective recomputation */
	SetFieldData(&reset_facets);
}

/* retrieve stored values */
void SamplingSurfaceT::LoadSamplingPoint(int facet, int point, iArrayT& nodes,
	dArrayT& phi, dArray2DT& Dphi)
{
#if __option(extended_errorcheck)
	/* check */
	if (facet < 0 ||
	    facet >= fNumFacets ||
	    point < 0 ||
	    point >= fNumSamples) throw eOutOfRange;
#endif

	int dex = fNumSamples*facet + point;
	fNeighborData.RowAlias(dex, nodes);
	fPhiData.RowAlias(dex, phi);
	Dphi.Set(fLocFacetCoords.MinorDim(), nodes.Length(), fDPhiData(dex));
}

/* reference surface normal */
void SamplingSurfaceT::ReferenceConfiguration(int facet, int point, dMatrixT& Q,
	dArrayT& normal)
{
	/* set facet coordinates */
	iArrayT nodes;
	fFacetConnects.RowAlias(facet, nodes);
	fLocFacetCoords.SetLocal(nodes);
	
	/* set sampling point */
	fSurfaceShape->SetIP(point);
	
	/* get surface jacobian */
	fSurfaceShape->Jacobian(Q);
	
	/* normal is in last column */
	Q.CopyColumn(Q.Cols()-1, normal);
}

/* facet coordinates */
void SamplingSurfaceT::FacetCoordinates(int facet, dArray2DT& coordinates) const
{
	coordinates.RowCollect(fFacetConnects(facet), fFacetCoords);
}

/* write storage statistics */
void SamplingSurfaceT::WriteStatistics(ostream& out) const
{
	/* neighbor count data */
	int  min = fNeighborData.MinMinorDim(0);
	int  max = fNeighborData.MaxMinorDim();
	int used = fNumFacets*fNumSamples;
	out << "\n MLS shape function data:\n";
	out << " Minimum number of nodal neighbors . . . . . . . = " << min << '\n';
	out << " Maximum number of nodal neighbors . . . . . . . = " << max << '\n';
	out << " Average number of nodal neighbors . . . . . . . = ";
	if (used != 0)
		out << fNeighborData.Length()/used << '\n';
	else
		out << "-\n";

	/* collect neighbor distribution */
	iArrayT counts(max + 1);
	counts = 0;

	/* skip points */
	if (fSkipPoints.Length() == fNeighborData.MajorDim())
	{
		for (int j = 0; j < fNeighborData.MajorDim(); j++)
			if (!fSkipPoints[j])
				counts[fNeighborData.MinorDim(j)]++;
	}	
	else
		for (int j = 0; j < fNeighborData.MajorDim(); j++)
			counts[fNeighborData.MinorDim(j)]++;

	out << " Sampling point neighbor number distribution:\n";
	out << setw(kIntWidth) << "number"
	    << setw(kIntWidth) << "count" << '\n';
	for (int k = 0; k < counts.Length(); k++)
		out << setw(kIntWidth) << k
	    	<< setw(kIntWidth) << counts[k] << '\n';	

	/* memory requirements */
	int nsd = fFacetCoords.MinorDim();
	int nip = fNumSamples;
	out << "\n Shape function storage:\n";
	int count = counts.Sum();
	out << " Total number of nodal neighbors . . . . . . . . = " << count << '\n';
	out << " Nodal shape function storage. . . . . . . . . . = " << count*sizeof(double) << " bytes\n";
	out << " Nodal shape function derivatives storage. . . . = " << nsd*count*sizeof(double) << " bytes\n";
}

/*************************************************************************
* Private
*************************************************************************/

/* compute/store shape function data - for set neighbors,
* (facets != NULL) => selective recomputation */
void SamplingSurfaceT::SetFieldData(const ArrayT<int>* facets)
{
	int nsd = fLocFacetCoords.MinorDim();
	iArrayT neighbors;
	iArrayT facet_nodes;
	dArrayT phi;
	dArray2DT Dphi;
	int length = (facets != NULL) ? facets->Length() : fNumFacets;
	for (int k = 0; k < length; k++)
	{
		/* target */
		int facet = (facets != NULL) ? (*facets)[k] : k;
		int dex = facet*fNumSamples;
	
		/* set facet coordinates */
		fFacetConnects.RowAlias(facet, facet_nodes);
		fLocFacetCoords.SetLocal(facet_nodes);
	
		/* loop over integration points */
		fSurfaceShape->TopIP();
		while (fSurfaceShape->NextIP())
		{
			/* field point coordinates */
			const dArrayT& coords = fSurfaceShape->IPCoords();
		
			/* compute field */
			fNeighborData.RowAlias(dex, neighbors);
			if (!fMFSupport.SetFieldUsing(coords, neighbors))
			{
				cout << "\n SamplingSurfaceT::SetSamplingPoints: failed to build field at\n"
				     <<   "     point " << fSurfaceShape->CurrIP()+1
				     << " of facet " << k+1 << endl;
				throw eGeneralFail;
			}
			
			/* store */
			fPhiData.RowAlias(dex, phi);
			int nnd = neighbors.Length();
			Dphi.Set(nsd, nnd, fDPhiData(dex));
			phi = fMFSupport.FieldAt();	
			Dphi = fMFSupport.DFieldAt();
			
			/* next point */
			dex++;
		}
	}
}
