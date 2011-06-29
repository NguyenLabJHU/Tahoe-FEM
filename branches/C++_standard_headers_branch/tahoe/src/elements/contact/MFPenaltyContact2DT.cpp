/* $Id: MFPenaltyContact2DT.cpp,v 1.17 2005-04-06 01:28:00 paklein Exp $ */
#include "MFPenaltyContact2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "iGridManager2DT.h"
#include "ModelManagerT.h"
#include "OutputSetT.h"

/* meshfree element group types */
#include "MeshFreeSSSolidT.h"
#include "MeshFreeFSSolidT.h"
#include "MeshFreeFSSolidAxiT.h"
#include "MeshFreeSupportT.h"

/* meshless methods usinf SCNI */
#include "SCNIMFT.h"

/* parameters (duplicated from Contact2DT) */
const int kNumFacetNodes = 2;
const int kMaxNumGrid    = 75;

using namespace Tahoe;

/* constructor */
MFPenaltyContact2DT::MFPenaltyContact2DT(const ElementSupportT& support):
	PenaltyContact2DT(support),
	fElementGroup(NULL),
	fMeshFreeSupport(NULL),
	fSCNI(NULL),
	fdvT_man(0, true),
	fRHS_man(0, fRHS),
	fOutputID(-1),
	fOutputForce(false)
{
	SetName("meshfree_contact_2D_penalty");	
}

/* register element for output */
void MFPenaltyContact2DT::RegisterOutput(void)
{
	/* inherited */
//	PenaltyContact2DT::RegisterOutput();
	
	/* write contact forces */
	if (fOutputForce) {
	
		/* output labels */
		ArrayT<StringT> n_labels(2*NumSD());
		n_labels[0] = "D_X";
		n_labels[1] = "D_Y";
		n_labels[2] = "F_X";
		n_labels[3] = "F_Y";

		/* all meshless nodes */
		const iArrayT& all_nf_nodes = (fSCNI) ? fSCNI->NodesUsed() : fMeshFreeSupport->NodesUsed();
	
		/* register output */
		OutputSetT output_set(all_nf_nodes, n_labels);
		fOutputID = ElementSupport().RegisterOutput(output_set);	
	}
}

/* write output */
void MFPenaltyContact2DT::WriteOutput(void)
{
	/* inherited */
//	PenaltyContact2DT::WriteOutput();

	/* write contact forces */
	if (fOutputForce) {
	
		/* work space */
		dArray2DT n_values(fForce.MajorDim(), 2*fForce.MinorDim()), disp(fForce.MajorDim(), fForce.MinorDim());

		/* reconstruct displacement field */
		if (fSCNI) {
			iArrayT all(fForce.MajorDim());
			all.SetValueToPosition();
			fSCNI->InterpolatedFieldAtNodes(all, disp);
		}
		else {
			fElementGroup->NodalDOFs(fMeshFreeSupport->NodesUsed(), disp);
		}

		/* write in forces */
		n_values.BlockColumnCopyAt(disp, 0);				
		n_values.BlockColumnCopyAt(fForce, NumSD());
		
		/* write output */
		ElementSupport().WriteOutput(fOutputID, n_values);
	}
}

/* describe the parameters needed by the interface */
void MFPenaltyContact2DT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PenaltyContact2DT::DefineParameters(list);

	/* the meshless element group */
	list.AddParameter(ParameterT::Integer, "meshless_group");
	
	/* write contact forces to output */
	ParameterT output_force(ParameterT::Boolean, "output_force");
	output_force.SetDefault(fOutputForce);
	list.AddParameter(output_force);
}

/* accept parameter list */
void MFPenaltyContact2DT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "MFPenaltyContact2DT::TakeParameterList";

	/* NOTE: fMeshFreeSupport must be resolved before calling PenaltyContact2DT::TakeParameterList
	 *       because it is needed during MFPenaltyContact2DT:: ExtractContactGeometry */

	/* resolve meshless element group */
	int group = list.GetParameter("meshless_group");
	group--;
	ElementBaseT& element = ElementSupport().ElementGroup(group);
	fElementGroup = &element;

#ifndef __NO_RTTI__
	/* cast to meshfree element types */
	const MeshFreeSSSolidT* mf_ss_solid = dynamic_cast<const MeshFreeSSSolidT*>(fElementGroup);
	const MeshFreeFSSolidT* mf_fs_solid = dynamic_cast<const MeshFreeFSSolidT*>(fElementGroup);
	const MeshFreeFSSolidAxiT* mf_fs_axi_solid = dynamic_cast<const MeshFreeFSSolidAxiT*>(fElementGroup);
	fSCNI = dynamic_cast<const SCNIMFT*>(fElementGroup);
	if (mf_ss_solid)
		fMeshFreeSupport = &(mf_ss_solid->MeshFreeSupport());
	else if (mf_fs_solid)
		fMeshFreeSupport = &(mf_fs_solid->MeshFreeSupport());
	else if (mf_fs_axi_solid)
		fMeshFreeSupport = &(mf_fs_axi_solid->MeshFreeSupport());
	else if (!fSCNI)
		ExceptionT::GeneralFail(caller, "element group %d is not meshfree", group+1);
#else
	/* use name to resolve meshfree type */
	const StringT& element_name = fElementGroup->Name();
	if (element_name == "small_strain_meshfree") {
		const MeshFreeSSSolidT* mf_ss_solid = (const MeshFreeSSSolidT*) fElementGroup;
		fMeshFreeSupport = &(mf_ss_solid->MeshFreeSupport());
	}
	else if (element_name == "large_strain_meshfree") {
		const MeshFreeFSSolidT* mf_fs_solid = (const MeshFreeFSSolidT*) fElementGroup;
		fMeshFreeSupport = &(mf_fs_solid->MeshFreeSupport());
	}
	else if (element_name == "large_strain_meshfree_axi") {
		const MeshFreeFSSolidAxiT* mf_fs_axi_solid = (const MeshFreeFSSolidAxiT*) fElementGroup;
		fMeshFreeSupport = &(mf_fs_axi_solid->MeshFreeSupport());
	}
	else if (element_name.StringMatch("mfparticle"))
		fSCNI = (const SCNIMFT*) fElementGroup;
	else
		ExceptionT::GeneralFail(caller, "could not resolve meshfree group %d", group+1);
#endif

	/* register arrays with memory manager */
	fStrikerCoords_man.SetWard(0, fStrikerCoords, NumSD());
	fdvT_man.Register(fdv1T);
	fdvT_man.Register(fdv2T);

	/* set map of node ID to meshfree point index */
	if (fMeshFreeSupport)
		fNodeToMeshFreePoint.SetMap(fMeshFreeSupport->NodesUsed());

	/* inherited */
	PenaltyContact2DT::TakeParameterList(list);

	/* write contact forces */
	fOutputForce = list.GetParameter("output_force");
	if (fOutputForce) {
		
		/* all meshless nodes */
		const iArrayT& all_nf_nodes = (fSCNI) ? fSCNI->NodesUsed() : fMeshFreeSupport->NodesUsed();

		/* set work space */
		fForce.Dimension(all_nf_nodes.Length(), NumSD());
		fNodesUsed_inv.SetMap(all_nf_nodes);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called by FormRHS and FormLHS */
void MFPenaltyContact2DT::LHSDriver(GlobalT::SystemTypeT)
{
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (formK)
		ExceptionT::GeneralFail("MFPenaltyContact2DT::LHSDriver", "not implemented");
}

void MFPenaltyContact2DT::RHSDriver(void)
{
	const char caller[] = "MFPenaltyContact2DT::RHSDriver";

	/* time integration parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* references to global nodal data */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& disp = Field()[0]; /* displacements */

	/* compute current configuration of active strikers */
	ComputeStrikerCoordinates(fActiveStrikers);

	/* reset tracking data */
	int num_contact = 0;
	double h_max = 0.0;

	/* check */
	if (!fSCNI && !fMeshFreeSupport) ExceptionT::GeneralFail(caller, "no meshless support");

	/* initialize output array */
	if (fOutputForce) fForce = 0.0;

	/* loop over active elements */
	int nsd = NumSD();
	dArrayT tangent(nsd);
	AutoArrayT<int> nodes;
	AutoArrayT<int> eqnos;
	iArray2DT eqnos2D;
	iArrayT neighbors;
	dArrayT Na;
	dArrayT rhs_tmp;
	dArray2DT DNa;
	const int* pelem = fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		/* collect element configuration */
		fElCoord.RowCollect(pelem, init_coords);
		fElDisp.RowCollect(pelem, disp);

		/* current configuration */
		fElCoord.AddScaled(constKd, fElDisp);
	
		/* get facet coords */
		fElCoord.RowAlias(0, fx1);
		fElCoord.RowAlias(1, fx2);

		/* collect information for current meshless striker */
		int striker_node = pelem[2];
		int striker_meshfree_index = fNodeToMeshFreePoint.Map(striker_node);
		int striker_active_index = fNodeToActiveStriker.Map(striker_node);
		fStrikerCoords.RowAlias(striker_active_index, fStriker);

		/* striker neighborhood */
		if (fMeshFreeSupport)
			fMeshFreeSupport->LoadNodalData(striker_node, neighbors, Na, DNa);
		else /* SCNI */ {
			fSCNI_Support.RowAlias(striker_meshfree_index, neighbors);
			fSCNI_Phi.RowAlias(striker_meshfree_index, Na);
		}	
		SetDerivativeArrays(Na);

		/* penetration vectors */
		fv1.DiffOf(fStriker, fx1);
		fv2.DiffOf(fStriker, fx2);

		/* tangent vector */
		tangent.DiffOf(fx2, fx1);

		/* distance to facet (could store some of this) */
		double magtan = tangent.Magnitude();				
		double      h = (fv2[0]*fv1[1] - fv1[0]*fv2[1])/magtan;
//		double  max_d =-magtan/10; //max penetration

		/* contact */
		if (h < 0.0)
		{
			/* tracking data */
			num_contact++;
			h_max = (h < h_max) ? h : h_max;

			/* penetration force */
			double dphi =-fK*h;
			
			/* initialize */
			fRHS_man.SetLength((neighbors.Length() + kNumFacetNodes)*nsd, false);
			fRHS = 0.0;
					
			/* d_tan contribution */
			fdtanT.Multx(tangent, fNEEvec);
			rhs_tmp.Alias(fNEEvec.Length(), fRHS.Pointer());
			rhs_tmp.AddScaled(-dphi*h/(magtan*magtan), fNEEvec);
						
			/* d_area */
			fColtemp1.Set(fdv1T.Rows(), fdv1T(0));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(1));
			fRHS.AddCombination(-dphi*fv2[1]/magtan, fColtemp1,
				                -dphi*fv1[0]/magtan, fColtemp2);
			
			fColtemp1.Set(fdv1T.Rows(), fdv1T(1));
			fColtemp2.Set(fdv2T.Rows(), fdv2T(0));
			fRHS.AddCombination(dphi*fv2[0]/magtan, fColtemp1,
				                dphi*fv1[1]/magtan, fColtemp2);
					
			/* get equation numbers */
			nodes.Dimension(kNumFacetNodes + neighbors.Length());
			nodes[0] = pelem[0];
			nodes[1] = pelem[1];
			nodes.CopyPart(kNumFacetNodes, neighbors, 0, neighbors.Length());
			eqnos.Dimension(nodes.Length()*NumDOF());
			eqnos2D.Alias(nodes.Length(), NumDOF(), eqnos.Pointer());
			Field().SetLocalEqnos(nodes, eqnos2D);

			/* assemble */
			ElementSupport().AssembleRHS(Group(), fRHS, eqnos2D);
			
			/* accumulate forces for output */
			if (fOutputForce) {
				const double* force = fRHS.Pointer(2*2); /* skip face nodes */
				for (int i = 0; i < neighbors.Length(); i++) {
					int lnd = fNodesUsed_inv.Map(neighbors[i]);
					for (int j = 0; j < 2; j++)
						fForce(lnd,j) += *force++;
				}
			}
		}
	}

	/* set tracking */
	SetTrackingData(num_contact, h_max);
}

/* echo contact bodies and striker nodes. After the read section, should have valid 
 * nodes/facet connectivities for the local database. */
void MFPenaltyContact2DT::ExtractContactGeometry(const ParameterListT& list)
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
		//TEMP
		if (fSCNI) ExceptionT::GeneralFail(caller, "\"all_nodes_as_strikers\" not supported with SCNI");
		fStrikerTags.Alias(fMeshFreeSupport->NodesUsed());
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized contact node specification \"\"",
			striker_spec.Name().Pointer());

	/* check to see that all strikers are meshfree - only apply check with meshless methods
	 * using fMeshFreeSupport. The SCNI classes do this check later. */
	if (striker_spec.Name() != "all_nodes_as_strikers" && fMeshFreeSupport) {
		InverseMapT map;
		map.SetOutOfRange(InverseMapT::MinusOne);
		map.SetMap(fMeshFreeSupport->NodesUsed());
		for (int i = 0; i < fStrikerTags.Length(); i++)
			if (map.Map(fStrikerTags[i]) == -1)
				ExceptionT::GeneralFail(caller, "striker %d is not meshfree", fStrikerTags[i]+1);
	}

	/* echo */
	if (print_input) {
		out << "\n Striker nodes:\n";
		fStrikerTags++;
		out << fStrikerTags.wrap(8) << '\n';
		fStrikerTags--;	
	}

	/* collect SCNI data for all striker tags */
	if (fSCNI)
	{
		/* set map of node ID to meshfree point index */
		fNodeToMeshFreePoint.SetMap(fStrikerTags);

		/* map global ID to local numbering */
		fSCNI_LocalID = fStrikerTags;
		if (!fSCNI->GlobalToLocalNumbering(fSCNI_LocalID))
			ExceptionT::GeneralFail(caller, "SCNI global->local failed");
	
		/* collect nodal neighbors and shape functions */
		fSCNI->NodalSupportAndPhi(fSCNI_LocalID, fSCNI_Support, fSCNI_Phi);
	}

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

/* generate contact element data */
bool MFPenaltyContact2DT::SetActiveInteractions(void)
{
	int last_num_active = fActiveStrikers.Length();

	/* current coords of all strikers */
	ComputeStrikerCoordinates(fStrikerTags);

	/* construct search grid if needed */
	if (!fGrid2D)
	{
		/* try to get roughly least 10 per grid */
		int ngrid = int(pow(fStrikerCoords.MajorDim()/10.0,
		                    1.0/fStrikerCoords.MinorDim())) + 1;

		ngrid = (ngrid < 2) ? 2 : ngrid;
		ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

		fGrid2D = new iGridManager2DT(ngrid, ngrid, fStrikerCoords, 0);
		if (!fGrid2D) ExceptionT::OutOfMemory("MFPenaltyContact2DT::SetActiveInteractions");

		/* search grid statistics */
		ostream& out = ElementSupport().Output();
		out << "\n Search grid: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
		fGrid2D->WriteStatistics(out);
	}
	
	/* (re-)set grid boundaries */
	fGrid2D->Reset();

	/* update by-body stored data */
	SetSurfacesData();
		
	/* set striker/facet data */
	SetActiveStrikers();
	
	/* set map from global ID to active striker index */
	fNodeToActiveStriker.SetMap(fActiveStrikers);
	
	/* assume changed unless last and current step have no active */
	if (last_num_active == 0 && fActiveStrikers.Length() == 0)
		return false;
	else
		return true;
}	

void MFPenaltyContact2DT::ComputeStrikerCoordinates(const ArrayT<int>& strikers)
{
	/* dimension */
	fStrikerCoords_man.SetMajorDimension(strikers.Length(), false);

	/* current striker coords */
	if (strikers.Length() > 0) {
	
		/* reconstruct displacement field */
		if (fSCNI) {
			fSCNI_tmp = strikers;
			if (!fSCNI->GlobalToLocalNumbering(fSCNI_tmp))
				ExceptionT::GeneralFail("MFPenaltyContact2DT::ComputeStrikerCoordinates", 
					"SCNI global->local failed");
			fSCNI->InterpolatedFieldAtNodes(fSCNI_tmp, fStrikerCoords);
		}
		else {
			iArrayT tmp;
			tmp.Alias(strikers);
			fElementGroup->NodalDOFs(tmp, fStrikerCoords);
		}

		/* compute current coordinates */
		const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
		for (int i = 0; i < strikers.Length(); i++)
			fStrikerCoords.AddToRowScaled(i, 1.0, init_coords(strikers[i]));
	}
}

/* set derivative arrays given the array of shape functions for the
 * nodes in the neighborhood of the meshfree striker. */
void MFPenaltyContact2DT::SetDerivativeArrays(const dArrayT& mf_shape)
{
	/* dimenions */
	int nsd = NumSD();
	int ndof = NumDOF();
	fdvT_man.Dimension((mf_shape.Length() + kNumFacetNodes)*ndof, nsd);

	/* facet nodes */	
	fdv1T = 0.0;
	fdv1T(0,0) =-1.0;
	fdv1T(1,1) =-1.0;

	fdv2T = 0.0;
	fdv2T(2,0) =-1.0;
	fdv2T(3,1) =-1.0;

	/* striker node */
	int dex = kNumFacetNodes*ndof;
	for (int i = 0; i < mf_shape.Length(); i++) {
		double shape = mf_shape[i];
		fdv1T(dex  , 0) = shape;
		fdv1T(dex+1, 1) = shape;
		fdv2T(dex  , 0) = shape;
		fdv2T(dex+1, 1) = shape;
		dex += ndof;
	}
}
