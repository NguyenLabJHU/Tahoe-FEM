/* $Id: ContactElementT.cpp,v 1.23 2001-09-24 20:37:24 rjones Exp $ */

#include "ContactElementT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ofstreamT.h"
#include "fstreamT.h"
#include "IOBaseT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "iGridManager2DT.h"
#include "XDOF_ManagerT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "SurfaceT.h"
#include "ContactSearchT.h"

/* parameters */ // unfortunately these are also in the derived classes
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
ContactElementT::ContactElementT
(FEManagerT& fe_manager, int num_enf_params):
	ElementBaseT(fe_manager),
	LHS(ElementMatrixT::kNonSymmetric),
	tmp_LHS(ElementMatrixT::kNonSymmetric),
	opp_LHS(ElementMatrixT::kNonSymmetric)
{
	fNumEnfParameters = num_enf_params;
	fXDOF_Nodes = NULL;
	fNumMultipliers = 0;
	ReadControlData();
}

ContactElementT::ContactElementT
(FEManagerT& fe_manager, int num_enf_params, XDOF_ManagerT* xdof_nodes):
        ElementBaseT(fe_manager),
	fXDOF_Nodes(xdof_nodes),
	LHS(ElementMatrixT::kNonSymmetric),
	tmp_LHS(ElementMatrixT::kNonSymmetric),
	opp_LHS(ElementMatrixT::kNonSymmetric)
{
	fNumEnfParameters = num_enf_params;
	if (!fXDOF_Nodes) throw eGeneralFail;
	ReadControlData();
}


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

/* initialization after construction */
void ContactElementT::Initialize(void)
{
	/* inherited, calls EchoConnectivityData */
	ElementBaseT::Initialize();

	/* initialize surfaces, connect nodes to coordinates */
	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[i].Initialize(ElementBaseT::fNodes,fNumMultipliers);
	}
#if 0
        /* set console access */
        iAddVariable("penalty_parameter", fpenalty);
#endif

	/* create search object */
	fContactSearch = 
	  new ContactSearchT(fSurfaces, fSearchParameters);

	/* workspace matrices */
	SetWorkspace();

	/* for bandwidth reduction in the case of no contact 
	 * make node-to-node pseudo-connectivities to link all bodies */
	int num_surfaces = fSurfaces.Length();
	if (num_surfaces > 1)
	{
		fSurfaceLinks.Allocate(num_surfaces - 1, 2);
		for (int i = 0; i < num_surfaces - 1; i++)
		{
			fSurfaceLinks(i,0) = fSurfaces[i  ].GlobalNodes()[0];
			fSurfaceLinks(i,1) = fSurfaces[i+1].GlobalNodes()[0];
		}
	}

	if (fXDOF_Nodes) {
		iArrayT numDOF(fSurfaces.Length());
		numDOF = fNumMultipliers;
		/* this calls GenerateElementData */
		/* register with node manager */
		fNodes->XDOF_Register(this, numDOF);
	}
	else {
		/* set initial contact configuration */
		bool changed = SetContactConfiguration();	
	}
}

void ContactElementT::SetWorkspace(void)
{
	// ARE THESE RIGHT?
	/* workspace matrices */
	n1.Allocate(fNumSD);
	l1.Allocate(fNumSD);
   	RHS_man.SetWard(kMaxNumFaceDOF,RHS);
   	tmp_RHS_man.SetWard(kMaxNumFaceDOF,tmp_RHS);
   	LHS_man.SetWard(kMaxNumFaceDOF,LHS);
   	opp_LHS_man.SetWard(kMaxNumFaceDOF,opp_LHS);
   	tmp_LHS_man.SetWard(kMaxNumFaceDOF,tmp_LHS);
   	N1_man.SetWard(kMaxNumFaceDOF,N1);
   	N2_man.SetWard(kMaxNumFaceDOF,N2);
   	N1n_man.SetWard(kMaxNumFaceDOF,N1n);
   	N2n_man.SetWard(kMaxNumFaceDOF,N2n);
   	weights_man.SetWard(kMaxNumFaceNodes,weights);
   	eqnums_man.SetWard(kMaxNumFaceDOF,eqnums,fNumSD);
   	opp_eqnums_man.SetWard(kMaxNumFaceDOF,opp_eqnums,fNumSD);
#if 0
  	/* dynamic work space managers for element arrays */
	fXDOFConnectivities_man.SetWard(0, fXDOFConnectivities, fNumElemNodes + 1);
	fXDOFEqnos_man.SetWard(0, fXDOFEqnos, fNumElemEqnos);

#endif
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

	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[i].InitializeMultiplierMap();
	}

	for (int i = 0; i < fSurfaces.Length(); i++) {
		/* tag potentially active nodes */
		fSurfaces[i].DetermineMultiplierExtent();
	}
	
	for (int i = 0; i < fSurfaces.Length(); i++) {
		/* number active nodes and total */
        /* store last dof tag and value */
        /* resize DOF tags array for number of potential contacts */
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
		const dArray2DT& multipliers = fNodes->XDOF(this, i);
		fSurfaces[i].AliasMultipliers(multipliers);

		/* set multiplier connectivity on faces */
		fSurfaces[i].SetMultiplierConnectivity();

		/* form potential connectivity for step */
 		fSurfaces[i].SetPotentialConnectivity();
 	}
}

/* set DOF values to the last converged solution, this is called after SetDOF */
void ContactElementT::ResetDOF(dArray2DT& XDOF, int tag_set) const
{
	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[tag_set].ResetMultipliers(XDOF);
	}
}

/* return the displacement-ghost node pairs to avoid pivoting*/
const iArray2DT& ContactElementT::DOFConnects(int tag_set) const
{ 
	return fSurfaces[tag_set].RealGhostNodePairs();
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
          ElementBaseT::fNodes->
		SetLocalEqnos(connectivities, equation_numbers);
	}
	else {
          ElementBaseT::fNodes->
	   XDOF_SetLocalEqnos(connectivities, equation_numbers);
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

void ContactElementT::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)

// look at EXODUS output in continuumelementT
        /* contact statistics */
        ostream& out = fFEManager.Output();
        out << "\n Contact tracking: group "
                << fFEManager.ElementGroupNumber(this) + 1 << '\n';
        out << " Time                           = "
                << fFEManager.Time() << '\n';

        /* output files */
        StringT filename;
        filename.Root(fFEManager.Input().filename());
        filename.Append(".", fFEManager.StepNumber());
        filename.Append("of", fFEManager.NumberOfSteps());

        for(int s = 0; s < fSurfaces.Length(); s++) {
            const ContactSurfaceT& surface = fSurfaces[s];

           if (fOutputFlags[kGaps]) {
                StringT gap_out;
                gap_out = gap_out.Append(filename,".gap");
                gap_out = gap_out.Append(s);
                ofstream gap_file (gap_out);
                surface.PrintGap(gap_file);
           }

           if (fOutputFlags[kNormals]) {
                StringT normal_out;
                normal_out = normal_out.Append(filename,".normal");
                normal_out = normal_out.Append(s);
                ofstream normal_file (normal_out);
                surface.PrintNormals(normal_file);
           }

           if (fOutputFlags[kStatus]) {
				surface.PrintStatus(cout);
           }

           if (fOutputFlags[kMultipliers]) {
				surface.PrintMultipliers(cout);
           }

//              surface.PrintContactArea(cout);
  }

}

/* compute specified output parameter and send for smoothing */
void ContactElementT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: contact tractions/forces
}

/* writing output - nothing to write */
void ContactElementT::RegisterOutput(void) {}

/* solution calls */
void ContactElementT::AddNodalForce(int node, dArrayT& force)
{
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


/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void ContactElementT::ReadControlData(void)
{
    /* streams */
    ifstreamT& in = fFEManager.Input(); 
    ostream&  out = fFEManager.Output(); 

    /* print flags */
    fOutputFlags.Allocate(kNumOutputFlags);
    for (int i = 0; i < fOutputFlags.Length(); i++) {
        in >> fOutputFlags[i];
    }
	out << " Print gaps                 = " << fOutputFlags[kGaps] << '\n';
	out << " Print normals              = " << fOutputFlags[kNormals] << '\n';

    int num_surfaces;
    in >> num_surfaces;
    if (num_surfaces < 1) throw eBadInputValue;
	out << " Number of contact surfaces. . . . . . . . . . . = "
	    << num_surfaces << '\n';

    int num_pairs;
    in >> num_pairs;
    if (num_pairs < 1 || num_pairs > num_surfaces*(num_surfaces-1))
        throw eBadInputValue;
	out << " Number of surface pairs with data . . . . . . . = "
	    << num_pairs << '\n';

	/* parameters */
	out << " Number of search parameters . . . . . . . . . . = "
	    << kSearchNumParameters << '\n';
	out << " Number of enforcement parameters. . . . . . . . = "
	    << fNumEnfParameters << '\n';
    fSearchParameters.Allocate(num_surfaces);
    fEnforcementParameters.Allocate(num_surfaces);
    int s1, s2;
    for (int i = 0; i < num_pairs ; i++)
    {
        in >> s1 >> s2;
        s1--; s2--;
        dArrayT& search_parameters = fSearchParameters(s1,s2);
        /* general parameters for search */
        search_parameters.Allocate (kSearchNumParameters);
        dArrayT& enf_parameters    = fEnforcementParameters(s1,s2);
        /* parameters specific to enforcement */
        enf_parameters.Allocate (fNumEnfParameters);
        for (int j = 0 ; j < search_parameters.Length() ; j++)
        {
            in >> search_parameters[j];
        }
        for (int j = 0 ; j < enf_parameters.Length() ; j++)
        {
            in >> enf_parameters[j];
        }
    }
    fSearchParameters.CopySymmetric();
    fEnforcementParameters.CopySymmetric();

}

/* print element group data */
void ContactElementT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

}

/* echo contact surfaces */
void ContactElementT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces = fSearchParameters.Rows();
	/* surfaces */
	out << " Surface connectivity data .........................\n";
	fSurfaces.Allocate(num_surfaces); 
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
				surface.InputSideSets(fFEManager, in, out);
				break;
			
			default:
				cout << "\n ContactElementT::EchoSurfaceData:"
                                     << " unknown surface specification\n";
				cout <<   "     mode " << spec_mode 
                                     << " for surface " << i+1 << '\n';
				throw eBadInputValue;
		}
		surface.PrintConnectivityData(out);
	}
}

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
