/* $Id: ContactT.cpp,v 1.19 2004-06-26 18:39:53 paklein Exp $ */
/* created: paklein (12/11/1997) */
#include "ContactT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ParentDomainT.h"
#include "InverseMapT.h"

using namespace Tahoe;

/* constructor */
ContactT::ContactT(const ElementSupportT& support, const FieldT& field, int numfacetnodes):
	ElementBaseT(support, field),
	fNumFacetNodes(numfacetnodes),
	fnum_contact(-1), fh_max(1)
{

}

/* destructor */
ContactT::~ContactT(void) {	}

/* form of tangent matrix */
GlobalT::SystemTypeT ContactT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT ContactT::RelaxSystem(void)
{
	/* write before reconfiguration since information will be reset */
	if (ElementSupport().WriteOutput())
		WriteContactInfo(ElementSupport().Output());

	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* generate contact element data */
	bool contact_changed = SetContactConfiguration();

	/* minimal test of new-ness */
	if (!contact_changed)
		return relax;
	else
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
}

/* initialization after constructor */
void ContactT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* set up work space */
	SetWorkSpace();
	
	/* set initial contact configuration */
	SetContactConfiguration();	
}

/* initialize current time increment. Reset the contact tracking data. */
void ContactT::InitStep(void)
{
	/* reset tracking data */
	fnum_contact = -1;
	fh_max = 1;
}

/* solution calls */
void ContactT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
//not implemented
}

/* Returns the energy as defined by the derived class types */
double ContactT::InternalEnergy(void)
{
//not implemented
	return 0.0;
}

/* writing output - nothing to write */
void ContactT::RegisterOutput(void) {}
void ContactT::WriteOutput(void) {}

/* compute specified output parameter and send for smoothing */
void ContactT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: contact tractions/forces
}

/* appends group connectivities to the array */
void ContactT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	ElementBaseT::ConnectsU(connects_1, connects_2);

	/* add surface links */
	connects_1.AppendUnique(&fSurfaceLinks);
}

/* returns no (NULL) geometry connectivies */
void ContactT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)
//	connects.Append(NULL);
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void ContactT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Number of facet nodes . . . . . . . . . . . . . = "
	    << fNumFacetNodes << '\n';
}

/* echo contact bodies and striker nodes */
void ContactT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	const char caller[] = "ContactT::EchoConnectivityData";

	int num_surfaces;
	in >> num_surfaces;
	out << " Number of contact surfaces. . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) ExceptionT::BadInputValue(caller);

	/* read contact bodies */
	fSurfaces.Dimension(num_surfaces);
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		int spec_mode;
		in >> spec_mode;
		switch (spec_mode)
		{
			case kNodesOnFacet:	
				InputNodesOnFacet(in, fSurfaces[i]);
				break;
			
			case kSideSets:
				InputSideSets(in, fSurfaces[i]);
				break;
			
			case kBodyBoundary:
				/* may resize the surfaces array */
				InputBodyBoundary(in, fSurfaces, i);
				num_surfaces = fSurfaces.Length();
				break;
		
			default:
				ExceptionT::BadInputValue(caller, "unknown surface specification mode %d for surface %d",
					spec_mode, i+1);
		}
	}
	
	// after the read section, should have valid nodes/facet connectivities
	// for the local database

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
		out << " Found empty contact surfaces:\n\n";
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
	
	/* striker nodes */
	int striker_spec_mode;
	in >> striker_spec_mode;
	switch (striker_spec_mode)
	{
		case kNodeSetList: /* read strikers from node sets */
			ReadStrikers(in, out);
			break;
		
		case kSurfaceNodes: /* collect nodes from contact surfaces */
			StrikersFromSurfaces();
			break;
	
		case kAllStrikers:  /* shallow striker coords */
			fStrikerCoords.Alias(ElementSupport().CurrentCoordinates());
			out << "\n Striker nodes: ALL\n";	

			//TEMP			
			ExceptionT::GeneralFail(caller, "all nodes as strikers not tested");
			
			break;

		case kSideSetList: /* collect strikers from side sets */
			StrikersFromSideSets(in, out);
			break;
	
		default:		
			ExceptionT::BadInputValue(caller, "unknown striker specification mode %d", striker_spec_mode);
	}

	/* echo */
	if (ElementSupport().PrintInput())
	{
		out << "\n Striker nodes:\n";
		fStrikerTags++;
		out << fStrikerTags.wrap(8) << '\n';
		fStrikerTags--;	
	}
	
	/* allocate striker coords */
	fStrikerCoords.Dimension(fStrikerTags.Length(), NumSD());
	
	/* set connectivity name */
	ModelManagerT& model = ElementSupport().Model();
	StringT name ("Contact");
	name.Append (ElementSupport().ElementGroupNumber(this) + 1);

	/* register with the model manager and let it set the ward */
	int nen = fNumFacetNodes + 1; /* facet nodes + 1 striker */
	if (!model.RegisterElementGroup(name, GeometryT::kLine, nen)) 
		ExceptionT::GeneralFail(caller, "could not register contact facets");

	/* set up fConnectivities */
	fConnectivities.Dimension(1);
	fConnectivities[0] = model.ElementGroupPointer(name);

	/* set up fBlockData to store block ID */
	fBlockData.Dimension(1);
	fBlockData[0].Set(name, 0, fConnectivities[0]->MajorDim(), -1);

	/* set managed equation numbers array */
	fEqnos.Dimension(1);
	fEqnos_man.SetWard(0, fEqnos[0], nen*NumDOF());
}

void ContactT::SetWorkSpace(void)
{
	/* allocate map to active strikers data */
	fActiveMap.Dimension(fStrikerTags.Length());
	fActiveMap = -1;

	/* make pseudo-element list to link surfaces in case
	 * bodies are not otherwise interacting (for the bandwidth
	 * reduction) */
	int num_surfaces = fSurfaces.Length();
	if (num_surfaces > 1)
	{
		fSurfaceLinks.Dimension(num_surfaces - 1, 2);
		for (int i = 0; i < num_surfaces - 1; i++)
		{
			fSurfaceLinks(i,0) = (fSurfaces[i  ])[0];
			fSurfaceLinks(i,1) = (fSurfaces[i+1])[0];
		}
	}
}

void ContactT::WriteContactInfo(ostream& out) const
{
	/* contact statistics */
	int d_width = OutputWidth(out, fActiveStrikersForce.Pointer());
	out << "\n Contact tracking: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	out << " Active strikers                = " << fActiveStrikers.Length() << '\n';
	if (fActiveStrikers.Length() > 0)
	{
		out << setw(kIntWidth) << "striker" 
		    << setw(kIntWidth) << "surface"
		    << setw(kIntWidth) << "facet"
		    << setw(fNumFacetNodes*kIntWidth) << "facet nodes" 
		    << setw(d_width) << "force" << '\n';
		for (int i = 0; i < fActiveStrikers.Length(); i++)
		{
			out << setw(kIntWidth) << fActiveStrikers[i] + 1;
			out << setw(kIntWidth) << fHitSurface[i] + 1; int hit_facet = fHitFacets[i];
			out << setw(kIntWidth) << hit_facet + 1;
			for (int j = 0; j < fNumFacetNodes; j++)
				out << setw(kIntWidth) << fSurfaces[fHitSurface[i]](hit_facet, j) + 1;
			out << setw(d_width) << fActiveStrikersForce[i];				
			out << '\n';
		}
		out << endl;
	}

	/* write tracking data */
	if (fnum_contact != -1) {
		out << " Number of contact interactions = " << fnum_contact << '\n';
		out << " Maximum penetration depth      = " << fh_max << '\n';
	} else { /* not set */
		out << " Number of contact interactions = --\n";
		out << " Maximum penetration depth      = --\n";	
	}
}

/* generate contact element data - return true if configuration has
* changed since the last call */
bool ContactT::SetContactConfiguration(void)
{
	int last_num_active = fActiveStrikers.Length();
	bool contact_changed = SetActiveInteractions();
	if (contact_changed)
	{
		/* resize */
		int nel = fActiveStrikers.Length();
		fActiveStrikersForce.Dimension(nel);
		fActiveStrikersForce = 0.0;
		fEqnos_man.SetMajorDimension(nel, false);

		/* update dimensions */
		ElementBlockDataT& block = fBlockData[0];
		block.Set(block.ID(), block.StartNumber(), fConnectivities[0]->MinorDim(), block.MaterialID());
		
		/* reset the model manager */
		ModelManagerT& model = ElementSupport().Model();
		model.ResizeElementGroup(block.ID(), nel);

		/* generate connectivities */
		SetConnectivities();
	}

	/* write list of active strikers */
	iArrayT tmp;
	tmp.Alias(fActiveStrikers);	
	ostream& out = ElementSupport().Output();
	out << "\n            time: " << ElementSupport().Time() << '\n';
	out <<   " previous active: " << last_num_active << '\n';
	out <<   "  current active: " << fActiveStrikers.Length() << '\n';
	if (fActiveStrikers.Length() > 0 && ElementSupport().PrintInput()) {
		out << setw(kIntWidth) << "node" 
		    << setw(kIntWidth) << "surface" 
		    << setw(kIntWidth) << "face" << '\n';
		for (int i = 0; i < fActiveStrikers.Length(); i++)
			out << setw(kIntWidth) << fActiveStrikers[i]+1
			    << setw(kIntWidth) << fHitSurface[i]+1
			    << setw(kIntWidth) << fHitFacets[i]+1 << '\n';
	}

	return contact_changed;
}

/***********************************************************************
* Private
***********************************************************************/

/* surface input functions */
void ContactT::InputNodesOnFacet(ifstreamT& in, iArray2DT& facets)
{
	int num_facets;
	in >> num_facets;
	if (num_facets < 0) throw ExceptionT::kBadInputValue;
	
	/* dimension */
	facets.Dimension(num_facets, fNumFacetNodes);
	
	/* read */
	in >> facets;
	
	/* correct numbering */
	facets--;
}

void ContactT::InputSideSets(ifstreamT& in, iArray2DT& facets)
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
}

void ContactT::InputBodyBoundary(ifstreamT& in, ArrayT<iArray2DT>& surfaces,
		int& surface)
{
	/* gather element group info */
	int elem_group;
	in >> elem_group;
	elem_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(elem_group);
	ArrayT<StringT> IDs;
	element.ElementBlockIDs(IDs);

	/* get sets of facet */
	GeometryT::CodeT geometry;
	ArrayT<iArray2DT> surface_facet_sets;
	iArrayT surface_nodes;
	ElementSupport().Model().SurfaceFacets(IDs, geometry, surface_facet_sets, surface_nodes);

	/* just one surface */
	if (surface_facet_sets.Length() == 1)
		surfaces[surface] = surface_facet_sets[0];
	else if (surface_facet_sets.Length() > 1)
	{
		int num_sets = surface_facet_sets.Length();
		
		/* resize surfaces array */
		surfaces.Resize(surfaces.Length() + (num_sets - 1)); //NOTE: not byte copy!
		for (int i = 0; i < num_sets; i++)
		{
			/* copy */
			surfaces[surface] = surface_facet_sets[i];

			/* next */
			surface++;
		}
	}
}

/* generate striker list from surfaces */
void ContactT::StrikersFromSurfaces(void)
{
	//TEMP just make big for now
	int num_nodes = ElementSupport().NumNodes();
	iArrayT counts(num_nodes);
	counts = 0;

	/* tally occurrences */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		iArray2DT& surface = fSurfaces[i];
		int* psurf  = surface.Pointer();
		int  length = surface.Length();
		for (int j = 0; j < length; j++)
		{
			counts[*psurf]++;
			psurf++;
		}
	}

	/* count surface nodes */
	int  node_count = 0;
	int* pcount = counts.Pointer();
	for (int j = 0; j < num_nodes; j++)
		if (*pcount++ > 0)
			node_count++;

	/* collect */
	fStrikerTags.Dimension(node_count);
	pcount = counts.Pointer();
	int* pstrike = fStrikerTags.Pointer();
	for (int k = 0; k < num_nodes; k++)
		if (*pcount++ > 0)
			*pstrike++ = k;
}

void ContactT::ReadStrikers(ifstreamT& in, ostream& out)
{
#pragma unused(out)

  ModelManagerT& model = ElementSupport().Model();

  /* read list of node set id indexes */
  ArrayT<StringT> ns_ID;
  model.NodeSetList(in, ns_ID);

  /* collect nodes from those indexes */
  model.ManyNodeSets(ns_ID, fStrikerTags);
}

void ContactT::StrikersFromSideSets(ifstreamT& in, ostream& out)
{
#pragma unused(out)

	/* read data from parameter file */
	ArrayT<StringT> ss_ID;
	bool multidatabasesets = true;
	ModelManagerT& model = ElementSupport().Model();
	model.SideSetList(in, ss_ID, multidatabasesets);

	/* list node nodes used */
	iArrayT nodes_used(model.NumNodes());
	nodes_used = 0;

	/* mark nodes used in side sets */
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	iArray2DT facets;
	for (int i = 0; i < ss_ID.Length(); i++)
	{
		/* read side set */
		model.SideSet(ss_ID[i], facet_geom, facet_nodes, facets);
	
		/* mark nodes used */
		for (int j = 0; j < facets.Length(); j++)
			nodes_used[facets[j]] = 1;
	}
	
	/* collect nodes */
	fStrikerTags.Dimension(nodes_used.Count(1));
	int dex = 0;
	for (int i = 0; i < nodes_used.Length(); i++)
		if (nodes_used[i] == 1)
			fStrikerTags[dex++] = i;
}

/* compute the nodal area associated with each striker node */
void ContactT::ComputeNodalArea(const ArrayT<StringT>& striker_blocks, 
	dArrayT& nodal_area, InverseMapT& inverse_map)
{
	/* initialize nodal area */
	nodal_area.Dimension(fStrikerTags.Length());
	nodal_area = 0.0;

	/* get surface faces */
	GeometryT::CodeT geometry;
	ArrayT<iArray2DT> surfaces;
	iArrayT surface_nodes;
	ElementSupport().Model().SurfaceFacets(striker_blocks, geometry, surfaces, surface_nodes);

	/* no surfaces */
	if (surfaces.Length() == 0) return;

	/* map to local id of striker nodes */
	inverse_map.SetOutOfRange(InverseMapT::MinusOne);
	inverse_map.SetMap(fStrikerTags);

	/* shape functions over the faces */
	int nip = 1;
	int nfn = surfaces[0].MinorDim();
	ParentDomainT surf_shape(geometry, nip, nfn);
	surf_shape.Initialize();

	/* coordinates over the face */
	int nsd = NumSD();
	LocalArrayT ref_coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(ref_coords);
	dMatrixT jacobian(nsd, nsd-1);

	/* loop over surfaces */
	const double* Na = surf_shape.Shape(0);
	const double* w  = surf_shape.Weight();
	iArrayT facet_nodes;
	for (int i = 0; i < surfaces.Length(); i++)
	{
		const iArray2DT& surface = surfaces[i];

		/* loop over faces */
		for (int j = 0; j < surface.MajorDim(); j++)
		{
			/* face nodes */
			surface.RowAlias(j, facet_nodes);
		
			/* gather coordinates */
			ref_coords.SetLocal(facet_nodes);
		
			/* coordinate mapping */
			surf_shape.DomainJacobian(ref_coords, 0, jacobian);
			double detj = surf_shape.SurfaceJacobian(jacobian);	
		
			/* loop over face nodes */
			for (int k = 0; k < facet_nodes.Length(); k++)
			{
				/* striker node index */
				int index = inverse_map.Map(facet_nodes[k]);
				if (index != -1)
					nodal_area[index] += w[0]*detj*Na[k];
			}
		}
	}
}
