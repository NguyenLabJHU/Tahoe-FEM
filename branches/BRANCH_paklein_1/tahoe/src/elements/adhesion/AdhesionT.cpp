/* $Id: AdhesionT.cpp,v 1.1.2.1 2002-10-17 04:24:24 paklein Exp $ */
#include "AdhesionT.h"

#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "SurfaceShapeT.h"
#include "iArrayT.h"

using namespace Tahoe;

const int kAvgCellNodes = 10;

/* constructor */
AdhesionT::AdhesionT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fGrid(kAvgCellNodes, -1, fFaceCentroids, NULL)
{

}

/* destructor */
AdhesionT::~AdhesionT(void)
{
	for (int i = 0; i < fShapes.Length(); i++)
		delete fShapes[i];
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT AdhesionT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* generate interaction element data */
	bool changed = SetConfiguration();

	/* minimal test of new-ness */
	if (!changed)
		return relax;
	else
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
}

/* initialization after constructor */
void AdhesionT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* set up work space */
	SetWorkSpace();
	
	/* set initial contact configuration */
	SetConfiguration();	
}

/* solution calls */
void AdhesionT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
//not implemented
}

/* Returns the energy as defined by the derived class types */
double AdhesionT::InternalEnergy(void)
{
//not implemented
	return 0.0;
}

/* writing output - nothing to write */
void AdhesionT::RegisterOutput(void) { }
void AdhesionT::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)

}

/* compute specified output parameter and send for smoothing */
void AdhesionT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: contact tractions/forces
}

/* appends group connectivities to the array */
void AdhesionT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	ElementBaseT::ConnectsU(connects_1, connects_2);
}

/* returns no (NULL) geometry connectivies */
void AdhesionT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void AdhesionT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

}

/* echo contact bodies and striker nodes */
void AdhesionT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces;
	in >> num_surfaces;
	out << " Number of surfaces. . . . . . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) throw ExceptionT::kBadInputValue;

	/* read surfaces */
	AutoArrayT<iArray2DT*> surfaces;
	AutoArrayT<GeometryT::CodeT> geom;
	for (int i = 0; i < num_surfaces; i++)
	{
		ArrayT<GeometryT::CodeT> new_geom;
		ArrayT<iArray2DT> new_surfaces;
	
		int spec_mode;
		in >> spec_mode;
		switch (spec_mode)
		{
			case kNodesOnFacet:
			{
				geom.Dimension(1);
				new_surfaces.Dimension(1);
				InputNodesOnFacet(in, new_geom[0], new_surfaces[0]);
				break;
			}
			case kSideSets:
			{
				geom.Dimension(1);
				new_surfaces.Dimension(1);
				InputSideSets(in, new_geom[0], new_surfaces[0]);
				break;
			}
			case kBodyBoundary:
			{
				/* may resize the surfaces array */
				InputBodyBoundary(in, new_geom, new_surfaces);
				num_surfaces = fSurfaces.Length();
				break;
			}
			default:
				cout << "\n AdhesionT::EchoConnectivityData: unknown surface specification\n";
				cout <<   "     mode " << spec_mode << " for surface " << i+1 << '\n';
				throw ExceptionT::kBadInputValue;
		}
		
		/* collect */
		geom.Append(new_geom);
		for (int j = 0; j < new_surfaces.Length(); j++)
		{
			iArray2DT* surf = new iArray2DT;
			surf->Swap(new_surfaces[j]);
			surfaces.Append(surf);
		}
		
		/* next */
		num_surfaces += new_surfaces.Length() - 1;
		i += new_surfaces.Length();
	}
	
	/* set per surface data */
	fSurfaces.Dimension(num_surfaces);
	fLocCurrCoords.Dimension(fSurfaces.Length());
	fShapes.Dimension(fSurfaces.Length());
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/* facets data */
		fSurfaces[i].Swap(*surfaces[i]);

		/* local array for current coordinates */
		fLocCurrCoords[i].SetType(LocalArrayT::kCurrCoords);
		fLocCurrCoords[i].Dimension(fSurfaces[i].MinorDim(), NumSD());
		ElementSupport().RegisterCoordinates(fLocCurrCoords[i]);
	
		/* surface shape functions */
		SurfaceShapeT* shape = new SurfaceShapeT(geom[i], NumIP(geom[i]), 
			fLocCurrCoords[i].NumberOfNodes(), NumSD(), fLocCurrCoords[i]);
	}
	
	/* delete temp space */
	for (int i = 0; surfaces.Length(); i++)
		delete surfaces[i];

	/* echo data and correct numbering offset */
	out << " Contact surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	int surface_count = 0;
	for (int j = 0; j < fSurfaces.Length(); j++)
	{		
	  	iArray2DT& surface = fSurfaces[j];

	  	out << setw(kIntWidth) << j+1
	  	    << setw(kIntWidth) << surface.MajorDim()
	  	    << setw(kIntWidth) << surface.MinorDim() << "\n\n";
	  	
	  	/* set offset for output */
	  	if (ElementSupport().PrintInput())
	  	{
	  		surface++;
	  		surface.WriteNumbered(out);
	  		surface--;
	  		out << '\n';
	  	}
	  	
	  	/* count non-empty */
	  	if (surface.MajorDim() > 0) surface_count++;
	}	
	
	/* remove empty surfaces */
	if (surface_count != fSurfaces.Length())
	{
		out << " Found empty surfaces:\n\n";
		ArrayT<iArray2DT> tmp_surfaces(surface_count);
		surface_count = 0;
		for (int i = 0; i < fSurfaces.Length(); i++)
		{
	  		iArray2DT& surface = fSurfaces[i];
			if (surface.MajorDim() == 0)
				out << " removing surface " << i+1 << '\n';
			else
				tmp_surfaces[surface_count++].Swap(surface);
		}
		
		/* exchange */
		fSurfaces.Swap(tmp_surfaces);
	}
}

void AdhesionT::SetWorkSpace(void)
{
	/* count total number of faces */
	int num_faces = 0;
	for (int i = 0; i < fSurfaces.Length(); i++)
		num_faces += fSurfaces[i].MajorDim();
		
	/* set face index array */
	fFaceIndex.Dimension(num_faces, kFaceIndexDim);
	int* p = fFaceIndex.Pointer();
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		const iArray2DT& surface = fSurfaces[i];
		for (int j = 0; j < surface.MajorDim(); j++)
		{
			*p++ = i; /* surface number */		
			*p++ = j; /* local index on surface */
		}
	}
}

/* generate contact element data - return true if configuration has
* changed since the last call */
bool AdhesionT::SetConfiguration(void)
{
	/* compute the face "centroids" */
	iArrayT face_nodes;
	dArrayT centroid;
	const ElementSupportT& support = ElementSupport();
	for (int i = 0; i < fFaceIndex.MajorDim(); i++)
	{
		/* current facet information */
		int surface_index = fFaceIndex(i,kSurface);
		const iArray2DT& surface = fSurfaces[surface_index];
		LocalArrayT& coords = fLocCurrCoords[surface_index];
		int local_index = fFaceIndex(i,kLocalIndex);
		surface.RowAlias(local_index, face_nodes);
		fFaceCentroids.RowAlias(i, centroid);
		
		/* collect current coordinates */
		coords.SetLocal(face_nodes);
	
		/* compute average coordinates */
		coords.Average(centroid);
	}
	
	/* reset the search grids */
	fGrid.Reset();

	/* search for interacting faces */

#if 0
	bool contact_changed = SetActiveInteractions();
	if (contact_changed)
	{
		/* resize */
		int nel = fActiveStrikers.Length();
		fConnectivities_man.SetMajorDimension(nel, false);
		fEqnos_man.SetMajorDimension(nel, false);

		/* generate connectivities */
		SetConnectivities();	

		/* update dimensions */
		ElementBlockDataT& block = fBlockData[0];
		block.Set(block.ID(), block.StartNumber(), fConnectivities[0]->MinorDim(), block.MaterialID());
	}
	
	return contact_changed;
#endif
	return false;
}

/***********************************************************************
* Private
***********************************************************************/

/* surface input functions */
void AdhesionT::InputNodesOnFacet(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets)
{
	int num_facets = -1;
	int num_facet_nodes = -1;
	in >> geom >> num_facets >> num_facet_nodes;
	if (num_facets < 0 || num_facet_nodes < 0) throw ExceptionT::kBadInputValue;
	
	/* dimension */
	facets.Dimension(num_facets, num_facet_nodes);
	
	/* read */
	in >> facets;
	
	/* correct numbering */
	facets--;
}

void AdhesionT::InputSideSets(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets)
{
	/* read data from parameter file */
	ArrayT<StringT> ss_ID;
	bool multidatabasesets = false; /* change to positive and the parameter file format changes */
	ModelManagerT& model = ElementSupport().Model();
	model.SideSetList(in, ss_ID, multidatabasesets);

	if (ss_ID.Length () != 1) {
		cout << "\n ContactT::InputSideSets: Model Manager read more than one side set, not programmed for this." << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* read side set faces */
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	model.SideSet(ss_ID[0], facet_geom, facet_nodes, facets);
	geom = facet_geom[0];
}

void AdhesionT::InputBodyBoundary(ifstreamT& in, ArrayT<GeometryT::CodeT>& geom,
	ArrayT<iArray2DT>& surfaces)
{
	/* gather element group info */
	int num_blocks = -1;
	in >> num_blocks;
	if (num_blocks < 1) throw ExceptionT::kBadInputValue;
	ArrayT<StringT> IDs(num_blocks);
	for (int i = 0; i < IDs.Length(); i++)
		in >> IDs[i];

	/* get sets of facet */
	GeometryT::CodeT geometry;
	iArrayT surface_nodes;
	ElementSupport().Model().SurfaceFacets(IDs, geometry, surfaces, surface_nodes);
	
	/* face geometries */
	geom.Dimension(surfaces.Length());
	geom = geometry;
}

/* return the number of integration points to use for the given face geometry */
int AdhesionT::NumIP(GeometryT::CodeT code) const
{
	switch (code) 
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
		case GeometryT::kTriangle:
			return 4;
		default:
		{
			cout << "\n AdhesionT::NumIP: unrecognized geometry: " << code << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
	return 0;
}
