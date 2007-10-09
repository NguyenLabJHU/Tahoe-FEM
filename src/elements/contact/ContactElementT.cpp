/* $Id: ContactElementT.cpp,v 1.51 2007-10-09 23:24:46 rjones Exp $ */
#include "ContactElementT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ofstreamT.h"
#include "ifstreamT.h"
#include "IOBaseT.h"
#include "iGridManager2DT.h"
#include "XDOF_ManagerT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

#include "SurfaceT.h"
#include "ContactSearchT.h"
#include "ContactNodeT.h"
#include "PenaltyContactElement2DT.h"
#include "PenaltyContactElement3DT.h"
#include "MultiplierContactElement2DT.h"
#include "MultiplierContactElement3DT.h"

#include "ParabolaT.h"
#include "ModSmithFerrante.h"
#include "GreenwoodWilliamson.h"
#include "MajumdarBhushan.h"
#include "GWPlastic.h"

#include "ofstreamT.h"

/* vector functions */
#include "vector2D.h"

/* constants */
const double PI = 2.0*acos(0.0);


#undef TEXT_OUTPUT
#define TEXT_OUTPUT 0

using namespace Tahoe;

/* parameters */ // unfortunately these are also in the derived classes
static const int kMaxNumFaceNodes = 4; // 4node quads

/* constructor */
ContactElementT::ContactElementT(const ElementSupportT& support):
    ElementBaseT(support),
    LHS(ElementMatrixT::kNonSymmetric),
    tmp_LHS(ElementMatrixT::kNonSymmetric),
    fContactSearch(NULL),
    fXDOF_Nodes(NULL),
	fFirstPass(1)
{
	SetName("Jones_contact");

	fNumEnfParameters = 0;
	fNumMultipliers = 0;

	fNumMaterialModelParameters[kDefault] = 0;
	fNumMaterialModelParameters[kModSmithFerrante] = knSF;
	fNumMaterialModelParameters[kGreenwoodWilliamson] = knGW;
	fNumMaterialModelParameters[kMajumdarBhushan] = knMB;
	fNumMaterialModelParameters[kGWPlastic] = knGP;

//    ReadControlData();
}

#if 0
ContactElementT::ContactElementT
(const ElementSupportT& support, XDOF_ManagerT* xdof_nodes):
    ElementBaseT(support),
    fXDOF_Nodes(xdof_nodes),
    LHS(ElementMatrixT::kNonSymmetric),
    tmp_LHS(ElementMatrixT::kNonSymmetric),
    fContactSearch(NULL)
{
	SetName("Jones_contact");

    fNumEnfParameters = 0;
    if (!fXDOF_Nodes) throw ExceptionT::kGeneralFail;

//    ReadControlData();
}
#endif

/* destructor */
ContactElementT::~ContactElementT(void) 
{ 
	delete fContactSearch;
}

/* form of tangent matrix */
GlobalT::SystemTypeT ContactElementT::TangentType(void) const
{
	return GlobalT::kNonSymmetric; 
}

void ContactElementT::SetWorkspace(void)
{	/* workspace matrices */ 
	int nsd = NumSD();
	n1.Dimension(nsd);
	int size_of_eqnum = kMaxNumFaceNodes*nsd;
   	RHS_man.SetWard    (size_of_eqnum,RHS);
   	tmp_RHS_man.SetWard(size_of_eqnum,tmp_RHS);
   	N1_man.SetWard     (size_of_eqnum,N1);
   	N2_man.SetWard     (size_of_eqnum*nsd,N2);
   	LHS_man.SetWard    (size_of_eqnum*size_of_eqnum,LHS);
   	tmp_LHS_man.SetWard(size_of_eqnum*size_of_eqnum,tmp_LHS);
   	N1n_man.SetWard    (size_of_eqnum,N1n);
   	N2n_man.SetWard    (size_of_eqnum,N2n);
   	eqnums1_man.SetWard(size_of_eqnum,eqnums1,NumSD());
   	eqnums2_man.SetWard(size_of_eqnum,eqnums2,NumSD());
   	weights_man.SetWard(size_of_eqnum,weights);
	
	if (fXDOF_Nodes) {
	int size_of_xeqnum = kMaxNumFaceNodes*fNumMultipliers;
	P1_man.SetWard     (size_of_xeqnum*fNumMultipliers,P1);
	P2_man.SetWard     (size_of_xeqnum*fNumMultipliers,P2);
	P1values_man.SetWard(0,P1values,size_of_xeqnum);
	P2values_man.SetWard(0,P2values,size_of_xeqnum);
	xRHS_man.SetWard   (size_of_xeqnum*fNumMultipliers,xRHS);
	tmp_xRHS_man.SetWard   (size_of_xeqnum*fNumMultipliers,tmp_xRHS);
   	xeqnums1_man.SetWard(size_of_xeqnum,xeqnums1,fNumMultipliers);
   	xeqnums2_man.SetWard(size_of_xeqnum,xeqnums2,fNumMultipliers);
   	xconn1_man.SetWard(kMaxNumFaceNodes,xconn1);
   	xconn2_man.SetWard(kMaxNumFaceNodes,xconn2);
	}
}

/* done once per time-step */
GlobalT::RelaxCodeT ContactElementT::RelaxSystem(void)
{
   if (fXDOF_Nodes) {
	/* override - handled by DOFElement::Reconfigure */
	return GlobalT::kNoRelax;
   }
   else {
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
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int ContactElementT::Reconfigure(void)
{ // this overrides Relax
	return 1; // always reconfigure, since SetCont.Conf. is in Gen.El.Data
#if 0
	/* inherited */
        GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

        /* generate contact element data */
        bool contact_changed = SetContactConfiguration();

        /* minimal test of new-ness */
        if (contact_changed)
                relax = GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);

        if (relax != GlobalT::kNoRelax)
                return 1;
        else
                return 0;
#endif
}

/* this sequence supplants Relax */
/* (1) sizes the DOF tags array needed for the current timestep */
void ContactElementT::SetDOFTags(void)
{ 
	bool changed = fContactSearch->SetInteractions();
	
	// last step before losing the mutliplier values

	/* Initialize */
	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[i].InitializeMultiplierMap();

		/* form potential connectivity for step */
 		fSurfaces[i].SetPotentialConnectivity();
	}

	/* Tag potentially active nodes */
	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[i].DetermineMultiplierExtent();
	}
	
	/* Number active nodes and total */
	/* Store last dof tag and value */
	/* Resize DOF tags array for number of potential contacts */
	for (int i = 0; i < fSurfaces.Length(); i++) {
	    fSurfaces[i].AllocateMultiplierTags();
	}

}

/* (2) this function allows the external manager to set the Tags */
iArrayT& ContactElementT::DOFTags(int tag_set)
{
        return fSurfaces[tag_set].MultiplierTags(); 
}

/* (3) generates connectivity based on current tags */
void ContactElementT::GenerateElementData(void)
{ 
	for (int i = 0; i < fSurfaces.Length(); i++) {
		/* hand off location of multipliers */
		const dArray2DT& multipliers 
			= ElementSupport().XDOF_Manager().XDOF(this, i);

		fSurfaces[i].AliasMultipliers(multipliers);

		/* form potential connectivity for step */
 		fSurfaces[i].SetMultiplierConnectivity();
 	}
}

/* set DOF values to the last converged solution, this is called after SetDOF */
void ContactElementT::ResetDOF(dArray2DT& XDOF, int tag_set) const
{
	fSurfaces[tag_set].ResetMultipliers(XDOF);
}

/* return the displacement-ghost node pairs to avoid pivoting*/
const iArray2DT& ContactElementT::DOFConnects(int tag_set) const
{
	ContactSurfaceT& contact_surface = const_cast<ContactSurfaceT&>(fSurfaces[tag_set]);
	return contact_surface.DisplacementMultiplierNodePairs();
}

/* append element equations numbers to the list */
void ContactElementT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
                AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)
	
  /* send potential connectivity */
  for (int i = 0; i < fSurfaces.Length(); i++) {
	const RaggedArray2DT<int>& connectivities   
		= fSurfaces[i].Connectivities(); 
	RaggedArray2DT<int>& equation_numbers 
		= fSurfaces[i].EqNums();
        /* get local equations numbers for u nodes from NodeManager */
	/* Connectivities generated in SetConfiguration */
	if (!fXDOF_Nodes ) {
		Field().SetLocalEqnos(connectivities, equation_numbers);
	}
	else {
		ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), connectivities, equation_numbers);
	}

        /* add to list */
        eq_2.Append(&equation_numbers);
  }

}


/* appends group connectivities to the array for graph-based algorithms */
void ContactElementT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	/* base class uses fConnectivities to create profile */
//ElementBaseT::ConnectsU(connects_1, connects_2);

	/* link surfaces with fictious node-to-node pairs*/
	/* only necessary for bodies out-of-contact */
	connects_1.AppendUnique(&fSurfaceLinks);
	
	/* add node-face interactions */
	for (int i = 0; i < fSurfaces.Length(); i++) {
	  const RaggedArray2DT<int>& connectivities   
		= fSurfaces[i].Connectivities(); 
	  connects_2.Append(&connectivities);
	}
}

/* returns no (NULL) geometry connectivies */
void ContactElementT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)
	connects.Append(NULL);
}

void ContactElementT::RegisterOutput(void) 
{
	int i = 0;
	ArrayT<StringT> labels(2);
    if (fOutputFlags[kGaps]) {
		labels[i++] = "GAP";
	}
    if (fOutputFlags[kMultipliers]) {
		labels[i++] = "PRE";
	}
	fNumOutputVariables = i;
	ArrayT<StringT> n_labels(fNumOutputVariables); 
	for (i = 0; i < fNumOutputVariables; i++) {n_labels[i] = labels[i];}

	fOutputID.Dimension(fSurfaces.Length());
	fOutputID = -1;

    /* register each surface */
	if (fNumOutputVariables) {
    	for (i = 0; i < fOutputID.Length(); i++)
    	{
        	/* set output specifier */
        	OutputSetT output_set
				(GeometryT::kPoint, fSurfaces[i].GlobalNodeNumbers(), n_labels);

        	/* register and get output ID */
        	fOutputID[i] = ElementSupport().RegisterOutput(output_set);
    	}
    }
}


void ContactElementT::WriteOutput(void)
{
ExceptionT::GeneralFail("ContactElementT::WriteOutput", "out of date");
#if 0
// look at EXODUS output in continuumelementT
	/* contact statistics */
	ostream& out = ElementSupport().Output();
	out << "\n Contact tracking: group "
                << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = "
                << ElementSupport().Time() << '\n';
	if (fNumOutputVariables) {
		for(int s = 0; s < fSurfaces.Length(); s++) {
			const ContactSurfaceT& surface = fSurfaces[s];
			dArray2DT n_values(surface.GlobalNodeNumbers().Length(), 
							fNumOutputVariables);
			n_values = 0.0;
			surface.CollectOutput(fOutputFlags,n_values);
			dArray2DT e_values;
			/* send to output */
			ElementSupport().WriteOutput(fOutputID[s], n_values, e_values);
		}
	}

	/* output files */
	StringT filename;
	filename.Root(ElementSupport().Input().filename());
	filename.Append(".", ElementSupport().StepNumber());
	filename.Append("of", ElementSupport().NumberOfSteps());

	for(int s = 0; s < fSurfaces.Length(); s++) {
		const ContactSurfaceT& surface = fSurfaces[s];

#if TEXT_OUTPUT
		if (fOutputFlags[kGaps]) {
                StringT gap_out;
                gap_out = gap_out.Append(filename,".gap");
                gap_out = gap_out.Append(s);
                ofstream gap_file (gap_out);
                surface.PrintGaps(gap_file);


		}
		if (fOutputFlags[kMultipliers]) {
//		surface.PrintMultipliers(cout);
                StringT pressure_out;
                pressure_out = pressure_out.Append(filename,".pre");
                pressure_out = pressure_out.Append(s);
                ofstream pressure_file (pressure_out);
                surface.PrintMultipliers(pressure_file);
		}

		if (fOutputFlags[kNormals]) {
                StringT normal_out;
                normal_out = normal_out.Append(filename,".normal");
                normal_out = normal_out.Append(s);
                ofstream normal_file (normal_out);
                surface.PrintNormals(normal_file);
		}
#endif

		if (fOutputFlags[kMultipliers]) { surface.PrintMultipliers(cout);}
		if (fOutputFlags[kStatus]) { surface.PrintStatus(cout); }
		if (fOutputFlags[kArea]) { surface.PrintContactArea(cout); }
	}
#endif
}

/* compute specified output parameter and send for smoothing */
void ContactElementT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: contact tractions/forces
}

/* solution calls */
void ContactElementT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
//not implemented
}

/* Returns the energy as defined by the derived class types */
double ContactElementT::InternalEnergy(void)
{
//not implemented
        return 0.0;
}

/* information about subordinate parameter lists */
void ContactElementT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* output flags */
	sub_list.AddSub("Jones_contact_output", ParameterListT::ZeroOrOnce);

	/* surfaces */
	sub_list.AddSub("Jones_contact_surfaces");

	/* surface interactions */
	sub_list.AddSub("Jones_contact_surface_pairs", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ContactElementT::NewSub(const StringT& name) const
{
	if (name == "Jones_contact_output")
	{
		ParameterContainerT* output = new ParameterContainerT(name);
		const char* labels[kNumOutputFlags] = {"gaps", "normals", "status", "multipliers", "contact_area"};
		ParameterT value(ParameterT::Integer, "value");
		value.SetDefault(0);
		for (int i = 0; i < kNumOutputFlags; i++) {
			value.SetName(labels[i]);
			output->AddParameter(value);
		}
		return output;
	}
	else if (name == "Jones_contact_surfaces")
	{
		ParameterContainerT* surface = new ParameterContainerT(name);
		surface->AddSub("side_set_ID_list");
		return surface;
	}
	else if (name == "Jones_contact_surface_pairs")
	{
		ParameterContainerT* pair = new ParameterContainerT(name);
		pair->SetSubSource(this);
		
		/* surface numbers */
		pair->AddParameter(ParameterT::Integer, "surface_1");
		pair->AddParameter(ParameterT::Integer, "surface_2");

		/* general parameters for search */
		pair->AddSub("Jones_search");

		/* parameters specific to enforcement */
		pair->AddSub("Jones_enforcement");

		/* constitutive parameters */
		pair->AddSub("Jones_material");

		return pair;	
	}
	else if (name == "Jones_search") {
		ParameterContainerT* search = new ParameterContainerT(name);
    LimitT lower(0.0, LimitT::Lower);

		ParameterT gap_tol(ParameterT::Double, "gap_tol");
    gap_tol.AddLimit(lower);
    gap_tol.SetDefault(1.0);
		search->AddParameter(gap_tol);

		ParameterT xi_tol(ParameterT::Double, "xi_tol");
    xi_tol.AddLimit(lower);
    xi_tol.SetDefault(0.1);
		search->AddParameter(xi_tol);

		ParameterT pass_type(ParameterT::Enumeration, "pass_type");
    pass_type.AddEnumeration("symmetric", ContactElementT::kSymmetric);
    pass_type.AddEnumeration("primary", ContactElementT::kPrimary);
    pass_type.AddEnumeration("secondary", ContactElementT::kSecondary);
		pass_type.AddEnumeration("deformable", ContactElementT::kDeformable);
    pass_type.AddEnumeration("rigid", ContactElementT::kRigid);
    pass_type.SetDefault(ContactElementT::kSymmetric);
    search->AddParameter(pass_type);

		return search;
	}
	else if (name == "Jones_enforcement") {
		ParameterContainerT* enf = new ParameterContainerT(name);
    LimitT lower(0.0, LimitT::Lower);

		ParameterT tangent(ParameterT::Enumeration, "consistent_tangent");
    tangent.AddEnumeration("true", 0);
    tangent.AddEnumeration("false", 1);
    tangent.SetDefault(1);
		enf->AddParameter(tangent);

		ParameterT penalty(ParameterT::Double, "penalty");
    penalty.AddLimit(lower);
    penalty.SetDefault(1.0);
		enf->AddParameter(penalty);

		ParameterT gscale(ParameterT::Double, "Gscale");
    gscale.AddLimit(lower);
    gscale.SetDefault(1.0);
		enf->AddParameter(gscale);

		ParameterT pscale(ParameterT::Double, "Pscale");
    pscale.AddLimit(lower);
    pscale.SetDefault(1.0);
		enf->AddParameter(pscale);

		ParameterT p_tol(ParameterT::Double, "p_tol");
    p_tol.AddLimit(lower);
    p_tol.SetDefault(1.0);
		enf->AddParameter(p_tol);

		ParameterT mat_type(ParameterT::Enumeration, "material_type");
    mat_type.AddEnumeration("none", ContactElementT::kDefault);
    mat_type.AddEnumeration("mod_SmithFerrante", ContactElementT::kModSmithFerrante);
    mat_type.AddEnumeration("GreenwoodWilliamson", ContactElementT::kGreenwoodWilliamson);
		mat_type.AddEnumeration("MajumdarBhushan", ContactElementT::kMajumdarBhushan);
    mat_type.SetDefault(ContactElementT::kDefault);
    enf->AddParameter(mat_type);

		return enf;
	}
	else if (name == "Jones_material")
		// make default zero
		// or move into contact pair
		return new DoubleListT(name);
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void ContactElementT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ContactElementT::TakeParameterList";

	/* inherited */
	ElementBaseT::TakeParameterList(list);
    
	/* output flags */
	fOutputFlags.Dimension(kNumOutputFlags);
	fOutputFlags = 0;
	const ParameterListT* output = list.List("Jones_contact_output");
	if (output) {
		const char* labels[kNumOutputFlags] = {"gaps", "normals", "status", "multipliers", "contact_area"};
		for (int i = 0; i < kNumOutputFlags; i++)
			fOutputFlags[i] = output->GetParameter(labels[i]);
	}

	/* dimension check */
	const ParameterListT& contact_surface = list.GetList("Jones_contact_surfaces");
	int num_surfaces = contact_surface.GetList("side_set_ID_list").NumLists("String");
	fSurfaces.Dimension(num_surfaces);
	int num_pairs = list.NumLists("Jones_contact_surface_pairs");
	if (num_pairs < 1 || num_pairs > num_surfaces*(num_surfaces-1))
		ExceptionT::BadInputValue(caller);

	/* extract pair data */
	fSearchParameters.Dimension(num_surfaces);
	fEnforcementParameters.Dimension(num_surfaces);
	fMaterialParameters.Dimension(num_surfaces);
	for (int i = 0; i < num_pairs ; i++)
	{
		/* pair parameters */
		const ParameterListT& pair = list.GetList("Jones_contact_surface_pairs", i);
		int s1 = pair.GetParameter("surface_1");
		int s2 = pair.GetParameter("surface_2");
		s1--; s2--;

		/* general parameters for search */
		dArrayT& search_parameters = fSearchParameters(s1,s2);
		const ParameterListT& search = pair.GetList("Jones_search");
		search_parameters.Dimension(kNumSearchParameters);
		search_parameters[kGapTol] = search.GetParameter("gap_tol");
		search_parameters[kXiTol] = search.GetParameter("xi_tol");
		//search_parameters[kPass] = (double) search.GetParameter("pass_type");
		switch ((int) search.GetParameter("pass_type")) {
		case kSymmetric : { search_parameters[kPass] = kSymmetric ;}
		case kPrimary :   { search_parameters[kPass] = kPrimary ;}
		case kSecondary : { search_parameters[kPass] = kSecondary ;}
		}

		/* parameters specific to enforcement */
		dArrayT& enf_parameters = fEnforcementParameters(s1,s2);
		const ParameterListT& enforcement = pair.GetList("Jones_enforcement");
		enf_parameters.Dimension(kNumEnfParameters);
		switch ((int) enforcement.GetParameter("consistent_tangent")) {
		case 0 : { enf_parameters[kConsistentTangent] = 0.0 ;}
		case 1 : { enf_parameters[kConsistentTangent] = 1.0 ;}
		}
		enf_parameters[kPenalty] = enforcement.GetParameter("penalty");
		enf_parameters[kGScale] = enforcement.GetParameter("Gscale");
		enf_parameters[kPScale] = enforcement.GetParameter("Pscale");
		enf_parameters[kTolP] = enforcement.GetParameter("p_tol");
		enf_parameters[kMaterialType] = (int) enforcement.GetParameter("material_type");

		/* material parameters */
		const ParameterListT& material = pair.GetList("Jones_material");
		int material_code = (int) enf_parameters[kMaterialType];
		int NumMatParameters = Num_of_Parameters(material_code);
		if (NumMatParameters != material.NumLists())
			ExceptionT::BadInputValue(caller, "expecting %d values in \"Jones_material\" not %d", NumMatParameters, material.NumLists());
		if (NumMatParameters > 0 ) {
			dArrayT& mat_parameters = fMaterialParameters(s1,s2);
			mat_parameters.Dimension(NumMatParameters);
			for (int ii = 0; ii < mat_parameters.Length(); ii++)
				mat_parameters[ii] = material.GetList(ii).GetParameter("value");
		}
	}
	fSearchParameters.CopySymmetric();
	fEnforcementParameters.CopySymmetric();
	fMaterialParameters.CopySymmetric();

	/* get side set ID's */
	const ParameterListT& surface_params = list.GetList("Jones_contact_surfaces");
	ArrayT<StringT> ss_ID;
	StringListT::Extract(surface_params.GetList("side_set_ID_list"), ss_ID);

	/* initialize surfaces, connect nodes to coordinates */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		ContactSurfaceT& surface = fSurfaces[i];
		surface.SetTag(i);

		/* translate ID list */
		surface.InputSideSets(ElementSupport(), ss_ID[i], ElementSupport().Output());
		surface.PrintConnectivityData(ElementSupport().Output());
	
		/* initialize */
		surface.Initialize(ElementSupport(), fNumMultipliers);
	}

	/* create search object */
	fContactSearch = new ContactSearchT(fSurfaces, fSearchParameters);

	/* workspace matrices */
	SetWorkspace();

	/* for bandwidth reduction in the case of no contact 
	 * make node-to-node pseudo-connectivities to link all bodies */
	if (num_surfaces > 1)
	{
		fSurfaceLinks.Dimension(num_surfaces - 1, 2);
		for (int i = 0; i < num_surfaces - 1; i++)
		{
			fSurfaceLinks(i,0) = fSurfaces[i  ].GlobalNodes()[0];
			fSurfaceLinks(i,1) = fSurfaces[i+1].GlobalNodes()[0];
		}
	}

	if (fXDOF_Nodes) {
		iArrayT numDOF(fSurfaces.Length());// the number of tag-sets
		numDOF = fNumMultipliers;
		/* this calls GenerateElementData */
		/* register with node manager */
		ElementSupport().XDOF_Manager().XDOF_Register(this, numDOF);
	}
	else {
		/* set initial contact configuration */
		bool changed = SetContactConfiguration();	
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/
void ContactElementT::TakePairData(const ParameterListT& list)
{
  /* write out search parameter matrix */
	ofstreamT& out = ElementSupport().Output();
	out << " Interaction parameters ............................\n";
	int num_surfaces = fSearchParameters.Rows();
  for (int i = 0; i < num_surfaces ; i++)
	{
		for (int j = i ; j < num_surfaces ; j++)
		{
			const dArrayT& search_parameters = fSearchParameters(i,j);
			const dArrayT& enf_parameters = fEnforcementParameters(i,j);
			const dArrayT& mat_parameters = fMaterialParameters(i,j);
			if (enf_parameters.Length() != kNumEnfParameters)
				ExceptionT::GeneralFail("PenaltyContactElement2DT::TakeParameterList",
					"expecting %d enforcement parameters not %d",
					kNumEnfParameters, enf_parameters.Length());

			/* only print allocated parameter arrays */
			if (search_parameters.Length() == kNumSearchParameters) {
		  	  out << "  surface pair: ("  << i << "," << j << ")\n" ;
			  out << "  gap tolerance:    "
					<< search_parameters[kGapTol] << '\n';
			  out << "  xi tolerance :    "
					<< search_parameters[kXiTol] << '\n';
			  out << "  pass flag    :    "
					<< (int) search_parameters[kPass] << '\n';
			  out << "  consis. tangent:  "
                    << (int) enf_parameters[kConsistentTangent] << '\n';
			  out << "  penalty :         "
					<< enf_parameters[kPenalty] << '\n';
			  out << "  penalty types:\n"
				  << "     Linear              " 
 				  << PenaltyContactElement2DT::kDefault << "\n"
				  << "     ModSmithFerrante    " 
				  << PenaltyContactElement2DT::kModSmithFerrante << "\n"
				  << "     GreenwoodWilliamson " 
			      << PenaltyContactElement2DT::kGreenwoodWilliamson << "\n"
			      << "     MajumdarBhushan     " 
			      << PenaltyContactElement2DT::kMajumdarBhushan << "\n"
				  << "     GWPlastic           " 
			      << PenaltyContactElement2DT::kGWPlastic           << "\n";
			  out << "  penalty Type :         "
					<< (int) enf_parameters[kMaterialType] << '\n';
			  switch ((int) enf_parameters[kMaterialType]) 
			  {
			  case kDefault: // no other parameters
			    out << "  <no parameters> \n";
				break;	
			  case kModSmithFerrante:
				out << "  Smith-Ferrante A : "
					<< mat_parameters[kSmithFerranteA] << '\n';
				out << "  Smith-Ferrante B : "
					<< mat_parameters[kSmithFerranteB] << '\n';
				break;	
			  case kGreenwoodWilliamson:
				out << "  Average asperity height            : "
					<< mat_parameters[kAsperityHeightMean] << '\n';
				out << "  Asperity height standard deviation : "
					<< mat_parameters[kAsperityHeightStandardDeviation] << '\n';
				out << "  Asperity density                   : "
					<< mat_parameters[kAsperityDensity] << '\n';
				out << "  Asperity Radius                    : "
					<< mat_parameters[kAsperityTipRadius] << '\n';
				out << "  Hertzian Modulus                   : "
					<< mat_parameters[kHertzianModulus] << '\n';
				break;	
			  case kMajumdarBhushan:
			  	out << " Asperity height standard deviation : "
			  		<< mat_parameters[kSigma] << '\n';
			  	out << "  Asperity roughness scale : "
			  		<< mat_parameters[kRoughnessScale] << '\n';
			  	out << "  Fractal dimension : "
			  		<< mat_parameters[kFractalDimension] << '\n';
			  	out << "  Hertzian Modulus                   : "
					<< mat_parameters[kEPrime] << '\n';
				out << "  Area Fraction                   : "
					<< mat_parameters[kAreaFraction] << '\n';
				break;	
			  case kGWPlastic:
				out << "  Mean height                        : "
					<< mat_parameters[kMean] << '\n';
				out << "  Standard deviation                 : "
					<< mat_parameters[kStandardDeviation] << '\n';
				out << "  Asperity density                   : "
					<< mat_parameters[kDensity] << '\n';
				out << "  Elastic modulus                    : "
					<< mat_parameters[kModulus] << '\n';
				out << "  Yield value                        : "
					<< mat_parameters[kYield] << '\n';
				out << "  Length-scale                       : "
					<< mat_parameters[kLength] << '\n';
				out << "  AsperityArea                       : "
					<< mat_parameters[kAsperityArea] << '\n';
				break;	
			  default:
				throw ExceptionT::kBadInputValue;
		  	  }
			}
		}
	}
	out <<'\n';

	/* set up Penalty functions */
	fPenaltyFunctions.Dimension(num_surfaces*(num_surfaces-1));
  for (int i = 0; i < num_surfaces ; i++)
  {
        for (int j = 0 ; j < num_surfaces ; j++)
        {
          dArrayT& enf_parameters = fEnforcementParameters(i,j);
          dArrayT& mat_parameters = fMaterialParameters(i,j);
		  if (enf_parameters.Length()) {
			switch ((int) enf_parameters[kMaterialType]) 
			{
			case PenaltyContactElement2DT::kDefault:
				// Macauley bracket:  <-x> ???  linear force
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] 
						= new ParabolaT(1.0);
				break;
			case PenaltyContactElement2DT::kModSmithFerrante:
				{
                double A = mat_parameters[kSmithFerranteA];
                double B = mat_parameters[kSmithFerranteB];
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] 
						= new ModSmithFerrante(A,B);
				}
				break;
			case PenaltyContactElement2DT::kGreenwoodWilliamson:
				{
                /* parameters for Greenwood-Williamson load formulation */
                double gw_m = mat_parameters[kAsperityHeightMean];
                double gw_s = mat_parameters[kAsperityHeightStandardDeviation];
                double gw_dens = mat_parameters[kAsperityDensity];
                double gw_mod = mat_parameters[kHertzianModulus];
                double gw_rad = mat_parameters[kAsperityTipRadius];
                double material_coeff=(4.0/3.0)*gw_dens*gw_mod*sqrt(gw_rad);
          		double area_coeff = PI*gw_dens*gw_rad;
				enf_parameters[kPenalty] *= material_coeff; // overwrite
				fPenaltyFunctions[LookUp(i,j,num_surfaces)]
                         = new GreenwoodWilliamson(1.5,gw_m,gw_s);
				}
				break;
			case PenaltyContactElement2DT::kMajumdarBhushan:
				{
                /* parameters for Majumdar-Bhushan load formulation */
                double mb_s = mat_parameters[kSigma];
                double mb_g = mat_parameters[kRoughnessScale];
                double mb_mod = mat_parameters[kEPrime];
                double mb_f = mat_parameters[kFractalDimension];
                double mb_c = mat_parameters[kAreaFraction];
                double material_coeff;
                if (mb_f==1.5)
                	material_coeff=sqrt(PI*mb_g)*mb_mod;
                else
                	material_coeff=(4.0/3.0)*sqrt(PI)*mb_mod
							*pow(mb_g,mb_f-1)*mb_f/(3.0-2.0*mb_f);
                double area_coeff = 0.5/mb_f;
				enf_parameters[kPenalty] *= material_coeff;// overwrite
				fPenaltyFunctions[LookUp(i,j,num_surfaces)]
                                        = new MajumdarBhushan(mb_f,mb_s,mb_c);
				}
				break;
			case PenaltyContactElement2DT::kGWPlastic:
				{
                /* parameters for cyclic formulation */
                double gp_mu = mat_parameters[kMean];
                double gp_sigma = mat_parameters[kStandardDeviation];
                double gp_dens = mat_parameters[kDensity];
                double gp_mod = mat_parameters[kModulus];
                double gp_len = mat_parameters[kLength];
                double gp_yld = mat_parameters[kYield];
                double gp_aa = mat_parameters[kAsperityArea];
                double gp_ade = mat_parameters[kAdhesionEnergy];
                double gp_adm = mat_parameters[kAdhesionModulus];
                double material_coeff=gp_dens;
          		double area_coeff = gp_dens;
				enf_parameters[kPenalty] *= material_coeff; // overwrite
				fPenaltyFunctions[LookUp(i,j,num_surfaces)]
                     = new GWPlastic(gp_mu,gp_sigma,gp_mod,gp_yld,gp_len,gp_aa,
									 gp_ade,gp_adm);
				}
				break;
			default:
				throw ExceptionT::kBadInputValue;
			}
		  }
		}
	}
}

#if 0
/* echo contact surfaces */
void ContactElementT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces = fSearchParameters.Rows();
	/* surfaces */
	out << " Surface connectivity data .........................\n";
	fSurfaces.Dimension(num_surfaces); 
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		int spec_mode;
		in >> spec_mode;
		ContactSurfaceT& surface = fSurfaces[i];
		surface.SetTag(i);
		/* read connectivity data */
		switch (spec_mode)
		{
			case kSideSets:
				surface.InputSideSets(ElementSupport(), in, out);
				break;
			
			default:
				cout << "\n ContactElementT::EchoSurfaceData:"
                                     << " unknown surface specification\n";
				cout <<   "     mode " << spec_mode 
                                     << " for surface " << i+1 << '\n';
				throw ExceptionT::kBadInputValue;
		}
		surface.PrintConnectivityData(out);
	}
}
#endif

/* generate contact element data - return true if configuration has
 * changed since the last call */
/* generate connectivity data based on current node-face pairs */
bool ContactElementT::SetContactConfiguration(void)
{
	bool changed = fContactSearch->SetInteractions();
	
	if (changed) { 
		/* form potential connectivity for step */
  		for (int i = 0; i < fSurfaces.Length(); i++) {
			fSurfaces[i].SetPotentialConnectivity();
  		}
	}

	return changed;
}

bool ContactElementT::UpdateContactConfiguration(void)
{
	bool changed = fContactSearch->UpdateInteractions();
	return changed;
}

int ContactElementT::PassType (int s1, int s2) const
{
    const dArrayT& parameters = fSearchParameters(s1,s2);
    int pass_code = (int) parameters[kPass];
    if (s1 == s2 || pass_code == 0) {
        return kSymmetric;
    }
    else if (s1 < s2) {
        switch (pass_code) {
            case  1: return kPrimary; break;
            case  2: return kSecondary; break;
            case -1: return kDeformable; break;
            case -2: return kRigid; break;
        }
    }
    else {//(s1 > s2)
        switch (pass_code) {
            case  2: return kPrimary; break;
            case  1: return kSecondary; break;
            case -2: return kDeformable; break;
            case -1: return kRigid; break;
        }
    }

    /* dummy return value */
    return kPrimary;
}

