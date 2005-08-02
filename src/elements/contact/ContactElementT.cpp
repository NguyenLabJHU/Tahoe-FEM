/* $Id: ContactElementT.cpp,v 1.50 2005-08-02 19:37:24 paklein Exp $ */
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
#include "SurfaceT.h"
#include "ContactSearchT.h"
#include "ContactNodeT.h"
#include "OutputSetT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

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
{ 	/* workspace matrices */ 
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
        pair->AddSub("Jones_properties");

		return pair;	
	}
	else if (name == "Jones_search")
		return new DoubleListT(name);
	else if (name == "Jones_enforcement")
		return new DoubleListT(name);
	else if (name == "Jones_properties")
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
        search_parameters.Dimension(search.NumLists());
        for (int i = 0; i < search_parameters.Length(); i++)
        	search_parameters[i] = search.GetList(i).GetParameter("value");

        /* parameters specific to enforcement */
        dArrayT& enf_parameters = fEnforcementParameters(s1,s2);
        const ParameterListT& enforcement = pair.GetList("Jones_enforcement");
        enf_parameters.Dimension(enforcement.NumLists());
        for (int i = 0; i < enf_parameters.Length(); i++)
        	enf_parameters[i] = enforcement.GetList(i).GetParameter("value");

		/* material parameters */
		const ParameterListT& properties = pair.GetList("Jones_properties");
		int material_code = (int) enf_parameters[enf_parameters.Length()-1];
		int NumMatParameters = Num_of_Parameters(material_code);
		if (NumMatParameters != properties.NumLists())
			ExceptionT::BadInputValue(caller, "expecting %d values in \"Jones_properties\" not %d",
				NumMatParameters, properties.NumLists());
		dArrayT& mat_parameters = fMaterialParameters(s1,s2);
		mat_parameters.Dimension(NumMatParameters);
		for (int i = 0; i < mat_parameters.Length(); i++)
			mat_parameters[i] = properties.GetList(i).GetParameter("value");
    }
	fSearchParameters.CopySymmetric();
	fEnforcementParameters.CopySymmetric();
	fMaterialParameters.CopySymmetric();

	/* initialize surfaces, connect nodes to coordinates */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		ContactSurfaceT& surface = fSurfaces[i];
		surface.SetTag(i);

		/* get side set ID's */
		const ParameterListT& surface_params = list.GetList("Jones_contact_surfaces", i);
		ArrayT<StringT> ss_ID;
		StringListT::Extract(surface_params.GetList("side_set_ID_list"), ss_ID);

		/* translate ID list */
		surface.InputSideSets(ElementSupport(), ss_ID, ElementSupport().Output());
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

