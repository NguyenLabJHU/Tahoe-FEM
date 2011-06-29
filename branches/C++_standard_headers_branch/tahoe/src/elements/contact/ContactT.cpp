/* $Id: ContactT.cpp,v 1.25 2008-07-14 17:41:12 lxmota Exp $ */
/* created: paklein (12/11/1997) */
#include "ContactT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ModelManagerT.h"
#include "ofstreamT.h"
#include "ParentDomainT.h"
#include "InverseMapT.h"
#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

using namespace Tahoe;

/* constructor */
ContactT::ContactT(const ElementSupportT& support, int numfacetnodes):
	ElementBaseT(support),
	fNumFacetNodes(numfacetnodes),
	fnum_contact(-1),
	fh_max(1)
{
	SetName("contact");
}

/* destructor */
ContactT::~ContactT(void) {	}

/* form of tangent matrix */
GlobalT::SystemTypeT ContactT::TangentType(void) const { return GlobalT::kSymmetric; }

/* prepare for a sequence of time steps */
void ContactT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();

	/* look for axisymmetric groups */
	bool axisymmetric = false;
	int num_groups = ElementSupport().NumElementGroups();
	for (int i = 0; i < num_groups; i++)
		axisymmetric = (axisymmetric || ElementSupport().ElementGroup(i).Axisymmetric());

	/* compute nodal tributary areas */
	ElementSupport().ModelManager().ComputeNodalArea(fStrikerTags,
		fStrikerArea, fStrikerTags_map, axisymmetric);
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

/* register output */
void ContactT::RegisterOutput(void)
{
	/* check */
	if (NumDOF() > 3)
		ExceptionT::GeneralFail("ContactT::RegisterOutput", "ndof %d > 3", NumDOF());

	/* output labels */
	ArrayT<StringT> n_labels(2*NumDOF()); /* displacements + force */

	/* displacements */
	int count = 0;
	const ArrayT<StringT>& labels = Field().Labels();
	for (int i = 0; i < labels.Length(); i++)
		n_labels[count++] = labels[i];

	/* forces */
	const char *force[] = {"F_X", "F_Y", "F_Z"};
	for (int i = 0; i < NumDOF(); i++)
		n_labels[count++] = force[i];

	/* output record */
	OutputSetT output_set(fStrikerTags, n_labels);

	/* register */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

/* write output */
void ContactT::WriteOutput(void)
{
	/* work space */
	dArray2DT n_values(fStrikerTags.Length(), 2*NumDOF());

	/* displacements */
	const dArray2DT& disp = (Field())[0];
	dArray2DT disp_tmp(fStrikerTags.Length(), NumDOF());
	disp_tmp.RowCollect(fStrikerTags, disp);
	n_values.BlockColumnCopyAt(disp_tmp, 0);

	/* forces */
	for (int i = 0; i < NumDOF(); i++)
		n_values.ColumnCopy(i + NumDOF(), fStrikerForce2D, i);

	/* send it */
	ElementSupport().WriteOutput(fOutputID, n_values);
}

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

/* information about subordinate parameter lists */
void ContactT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* surfaces */
	sub_list.AddSub("contact_surface", ParameterListT::OnePlus);

	/* striker nodes */
	sub_list.AddSub("contact_nodes");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ContactT::NewSub(const StringT& name) const
{
	if (name == "contact_surface") {

		ParameterContainerT* contact_surface = new ParameterContainerT(name);
		contact_surface->SetListOrder(ParameterListT::Choice);

		/* surface from side set  */
		ParameterContainerT surface_side_set("surface_side_set");
		surface_side_set.AddParameter(ParameterT::Word, "side_set_ID");
		contact_surface->AddSub(surface_side_set);

		/* surfaces from body boundary */
		ParameterContainerT body_boundary("body_boundary");
		body_boundary.AddParameter(ParameterT::Integer, "body_element_group");
		contact_surface->AddSub(body_boundary);

		return contact_surface;
	}
	else if (name == "contact_nodes") {

		ParameterContainerT* contact_nodes = new ParameterContainerT(name);
		contact_nodes->SetListOrder(ParameterListT::Choice);

		/* strikers from node sets */
		contact_nodes->AddSub("node_ID_list");

		/* strikers from side sets */
		contact_nodes->AddSub("side_set_ID_list");

		/* strikers from surfaces */
		contact_nodes->AddSub(ParameterContainerT("all_surface_nodes"));

		/* all nodes as strikers */
		contact_nodes->AddSub(ParameterContainerT("all_nodes_as_strikers"));

		return contact_nodes;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void ContactT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* dimension */
	int neq = (fNumFacetNodes + 1)*NumDOF();
	fLHS.Dimension(neq);
	fRHS.Dimension(neq);

	/* extract contact geometry from the list */
	ExtractContactGeometry(list);

	/* set up work space */
	SetWorkSpace();

	/* set initial contact configuration */
	SetContactConfiguration();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* echo contact bodies and striker nodes. After the read section, should have valid
 * nodes/facet connectivities for the local database. */
void ContactT::ExtractContactGeometry(const ParameterListT& list)
{
	const char caller[] = "ContactT::ExtractContactGeometry";

	/* output stream */
	ofstreamT& out = ElementSupport().Output();
	bool print_input = ElementSupport().PrintInput();

	/* get surfaces */
	int num_surfaces = list.NumLists("contact_surface");
	fSurfaces.Dimension(num_surfaces);

	/* read contact bodies */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		const ParameterListT& surface_spec = list.GetListChoice(*this, "contact_surface", i);

		if (surface_spec.Name() == "surface_side_set")
			InputSideSets(surface_spec, fSurfaces[i]);
		else if (surface_spec.Name() == "body_boundary")
		{
			/* may resize the surfaces array */
			InputBodyBoundary(surface_spec, fSurfaces, i);
			num_surfaces = fSurfaces.Length();
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized contact surface \"%s\"",
				surface_spec.Name().Pointer());
	}

	/* echo data  */
	out << " Contact surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	for (int j = 0; j < fSurfaces.Length(); j++)
	{
	  	iArray2DT& surface = fSurfaces[j];

	  	out << setw(kIntWidth) << j+1
	  	    << setw(kIntWidth) << surface.MajorDim()
	  	    << setw(kIntWidth) << surface.MinorDim() << "\n\n";

		/* verbose */
		if (print_input) {
			surface++;
			surface.WriteNumbered(out);
			surface--;
			out << '\n';
	  	}
	}

	/* look for empty surfaces */
	int surface_count = 0;
	for (int j = 0; j < fSurfaces.Length(); j++)
		if (fSurfaces[j].MajorDim() > 0)
			surface_count++;

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

	/* get strikers */
	const ParameterListT& striker_spec = list.GetListChoice(*this, "contact_nodes");
	if (striker_spec.Name() == "node_ID_list")
		StrikersFromNodeSets(striker_spec);
	else if (striker_spec.Name() == "side_set_ID_list")
		StrikersFromSideSets(striker_spec);
	else if (striker_spec.Name() == "all_surface_nodes")
		StrikersFromSurfaces();
	else if (striker_spec.Name() == "all_nodes_as_strikers")
	{
		fStrikerCoords.Alias(ElementSupport().CurrentCoordinates());

		//TEMP - not tested
		ExceptionT::GeneralFail("ContactT::EchoConnectivityData", "all nodes as strikers not tested");
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized contact node specification \"\"",
			striker_spec.Name().Pointer());

	/* echo */
	if (print_input) {
		out << "\n Striker nodes:\n";
		fStrikerTags++;
		out << fStrikerTags.wrap(8) << '\n';
		fStrikerTags--;
	}

	/* allocate striker coords */
	fStrikerCoords.Dimension(fStrikerTags.Length(), NumSD());

	/* space for output */
	fStrikerForce2D.Dimension(fStrikerTags.Length(), NumDOF());
	fStrikerForce2D = 0.0;

	/* set connectivity name */
	ModelManagerT& model = ElementSupport().ModelManager();
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
	int d_width = OutputWidth(out, fStrikerForce2D.Pointer());
	out << "\n Contact tracking: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	out << " Active strikers                = " << fActiveStrikers.Length() << '\n';
	dArrayT f_c;
	if (fActiveStrikers.Length() > 0 && ElementSupport().Logging() == GlobalT::kVerbose)
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
			fStrikerForce2D.RowAlias(fActiveStrikers[i], f_c);
			out << setw(d_width) << f_c.Magnitude();
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
		fEqnos_man.SetMajorDimension(nel, false);

		/* update dimensions */
		ElementBlockDataT& block = fBlockData[0];
		block.Set(block.ID(), block.StartNumber(), fConnectivities[0]->MinorDim(), block.MaterialID());

		/* reset the model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
		model.ResizeElementGroup(block.ID(), nel);

		/* generate connectivities */
		SetConnectivities();
	}

	/* write list of active strikers */
	if (ElementSupport().Logging() != GlobalT::kSilent) {
		iArrayT tmp;
		tmp.Alias(fActiveStrikers);
		ostream& out = ElementSupport().Output();
		out << "\n            time: " << ElementSupport().Time() << '\n';
		out <<   " previous active: " << last_num_active << '\n';
		out <<   "  current active: " << fActiveStrikers.Length() << '\n';
		if (fActiveStrikers.Length() > 0 && ElementSupport().Logging() == GlobalT::kVerbose) {
			out << setw(kIntWidth) << "node"
			    << setw(kIntWidth) << "surface"
			    << setw(kIntWidth) << "face" << '\n';
			for (int i = 0; i < fActiveStrikers.Length(); i++)
				out << setw(kIntWidth) << fActiveStrikers[i]+1
				    << setw(kIntWidth) << fHitSurface[i]+1
				    << setw(kIntWidth) << fHitFacets[i]+1 << '\n';
		}
	}

	return contact_changed;
}

/***********************************************************************
 * Private
 ***********************************************************************/

void ContactT::InputSideSets(const ParameterListT& list, iArray2DT& facets)
{
	const char caller[] = "ContactT::InputSideSets";

	/* extract side set ID */
	StringT ss_ID;
	ss_ID = list.GetParameter("side_set_ID");

	/* read side set faces */
	ModelManagerT& model = ElementSupport().ModelManager();
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	model.SideSet(ss_ID, facet_geom, facet_nodes, facets);
}

void ContactT::InputBodyBoundary(const ParameterListT& list, ArrayT<iArray2DT>& surfaces,
	int& surface)
{
	/* gather element group info */
	int elem_group = list.GetParameter("body_element_group");
	elem_group--;
	ElementBaseT& element = ElementSupport().ElementGroup(elem_group);
	ArrayT<StringT> IDs;
	element.ElementBlockIDs(IDs);

	/* get sets of facet */
	GeometryT::CodeT geometry;
	ArrayT<iArray2DT> surface_facet_sets;
	iArrayT surface_nodes;
	ElementSupport().ModelManager().SurfaceFacets(IDs, geometry, surface_facet_sets, surface_nodes);

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

void ContactT::StrikersFromNodeSets(const ParameterListT& list)
{
	/* collect node set id indexes */
	ArrayT<StringT> ns_ID;
	StringListT::Extract(list, ns_ID);

	/* collect nodes from those indexes */
	ModelManagerT& model = ElementSupport().ModelManager();
	iArrayT tmp;
	model.ManyNodeSets(ns_ID, fStrikerTags);
}

void ContactT::StrikersFromSideSets(const ParameterListT& list)
{
	/* read data from parameter file */
	ArrayT<StringT> ss_ID;
	StringListT::Extract(list, ss_ID);

	/* geometry database */
	ModelManagerT& model = ElementSupport().ModelManager();

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

#if 0
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
	ElementSupport().ModelManager().SurfaceFacets(striker_blocks, geometry, surfaces, surface_nodes);

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
#endif
