/* $Id: ContactElementT.cpp,v 1.17 2001-08-09 15:12:12 rjones Exp $ */

#include "ContactElementT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "IOBaseT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "iGridManager2DT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "SurfaceT.h"
#include "ContactSearchT.h"


/* constructor */
ContactElementT::ContactElementT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager)
{
	fXDOF_Nodes = NULL;
}

ContactElementT::ContactElementT
(FEManagerT& fe_manager, XDOF_ManagerT* xdof_nodes):
        ElementBaseT(fe_manager),
	fXDOF_Nodes(xdof_nodes)
{
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

/* initialization after constructor */
void ContactElementT::Initialize(void)
{
	/* inherited, calls EchoConnectivityData */
	ElementBaseT::Initialize();

	/* initialize surfaces, connect nodes to coordinates */
	for (int i = 0; i < fSurfaces.Length(); i++) {
		SurfaceT& surface = fSurfaces[i];
		surface.Initialize(ElementBaseT::fNodes);
	}
	
	/* create search object */
	fContactSearch = 
	  new ContactSearchT(fSurfaces, fSearchParameters);

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


	/* set initial contact configuration */
	bool changed = SetContactConfiguration();	
}

GlobalT::RelaxCodeT ContactElementT::RelaxSystem(void)
{
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

/* returns the array for the DOF tags needed for the current config */
iArrayT& ContactElementT::SetDOFTags(void)
{
#if 0
        /* store history */
        int old_length = fActiveStrikers.Length();
        fLastActiveMap = fActiveMap;
        dArrayT constraints;
        constraints.Alias(fXDOF_Nodes->XDOF(this));
        fLastDOF = constraints;

        /* resize DOF tags array */
        fContactDOFtags.Allocate(fActiveStrikers.Length());

        /* write list of active strikers */
        iArrayT tmp;
        tmp.Alias(fActiveStrikers);
        ostream& out = fFEManager.Output();
        out << "\nold: " << old_length << '\n';
        out << "new: " << fActiveStrikers.Length() << endl;
        out << "\n            time: " << fFEManager.Time() << '\n';
        out <<   " active strikers: " << tmp.Length()   << '\n';
        tmp++;
        out << tmp.wrap(8) << '\n';
        tmp--;

        return fContactDOFtags;
#endif
}

const iArrayT& ContactElementT::DOFTags(void) const
{
#if 0
        return fContactDOFtags;
#endif
}

/* generate element data (based on current striker/body data) */
void ContactElementT::GenerateElementData(void)
{
#if 0
        /* inherited - set nodal connectivities */
        Contact2DT::SetConnectivities();

        /* dimension */
        int num_active = fConnectivities.MajorDim();

        /* resize work space */
        fXDOFConnectivities_man.SetMajorDimension(num_active, false);
        fXDOFEqnos_man.SetMajorDimension(num_active, false);
        for (int i = 0; i < num_active; i++)
        {
                int*  pelem = fConnectivities(i);
                int* pxelem = fXDOFConnectivities(i);

                /* XDOF element tags */
                pxelem[0] = pelem[0]; // 1st facet node
                pxelem[1] = pelem[1]; // 2nd facet node
                pxelem[2] = pelem[2]; // striker node
                pxelem[3] = fContactDOFtags[i]; // contact DOF tag
        }
#endif
}

/* return the contact elements */
const iArray2DT& ContactElementT::DOFConnects(void) const
{
#if 0
        return fXDOFConnectivities;
#endif
}


/* restore the DOF values to the last converged solution */
void ContactElementT::ResetDOF(dArray2DT& DOF) const
{
#if 0
        /* alias */
        dArrayT constraints;
        constraints.Alias(DOF);
        constraints = 0.0;
        for (int i = 0; i < fLastActiveMap.Length(); i++)
        {
                int old_map = fLastActiveMap[i];
                int new_map = fActiveMap[i];
                if (old_map > -1 && new_map > -1)
                        constraints[new_map] = fLastDOF[old_map];
        }
#endif
}

/* returns 1 if group needs to reconfigure DOF's, else 0 */
int ContactElementT::Reconfigure(void)
{
#if 0
        /* inherited */
        GlobalT::RelaxCodeT relax = Contact2DT::RelaxSystem();
        if (relax != GlobalT::kNoRelax)
                return 1;
        else
                return 0;
#endif
}


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

/* writing output - nothing to write */
void ContactElementT::RegisterOutput(void) {}

void ContactElementT::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)
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


//    		surface.PrintContactArea(cout);
  }

}

/* compute specified output parameter and send for smoothing */
void ContactElementT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: contact tractions/forces
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
        /* get local equations numbers from NodeManager */
	/* Connectivities generated in SetConfiguration */
        ElementBaseT::fNodes->
		SetLocalEqnos(connectivities, equation_numbers);
#if 0
cout << "\n connectivities \n" ;//<< connectivities <<"\n";
for (int k = 0; k < connectivities.MajorDim(); k++) {
cout << "(" << k << ")";
for (int j = 0; j < connectivities.MinorDim(k); j++) {
cout << connectivities(k)[j] <<" ";
if ((j+1) == 1) cout << ",";
}
cout << "\n";
}

cout << "\n eq numbers \n" ; // << equation_numbers << "\n";
for (int k = 0; k < equation_numbers.MajorDim(); k++) {
cout << "(" << k << ")";
for (int j = 0; j < equation_numbers.MinorDim(k); j++) {
cout << equation_numbers(k)[j] <<" ";
if ((j+1) == fNumSD) cout << ",";
}
cout << "\n";
}
#endif
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
	}
}

/* returns no (NULL) geometry connectivies */
void ContactElementT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)
	connects.Append(NULL);
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void ContactElementT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

}

/* echo contact surfaces */
void ContactElementT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces;
	in >> num_surfaces;
	out << " Number of contact surfaces. . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) throw eBadInputValue;

	/* surfaces */
	fSurfaces.Allocate(num_surfaces); 
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		int spec_mode;
		in >> spec_mode;
		ContactSurfaceT& surface = fSurfaces[i];
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
		surface.SetTag(i);
		surface.PrintConnectivityData(out);
		surface.AllocateContactNodes();
	}

	/* print flags */
	fOutputFlags.Allocate(kNumOutputFlags);
	for (int i = 0; i < fOutputFlags.Length(); i++) {
		in >> fOutputFlags[i];
	}

	/* parameters */
	fSearchParameters.Allocate(num_surfaces);
	fEnforcementParameters.Allocate(num_surfaces);

	int num_pairs;
	in >> num_pairs;
	out << " Number of surface pairs with data . . . . . . . . = "
	    << num_pairs << '\n';
	int s1, s2;
	for (int i = 0; i < num_pairs ; i++) 
	{
		in >> s1 >> s2;
		s1--; s2--;
		dArrayT& search_parameters = fSearchParameters(s1,s2);
		search_parameters.Allocate (kSearchNumParameters);
		dArrayT& enf_parameters    = fEnforcementParameters(s1,s2);
		// add parameters specific to enforcement
		enf_parameters.Allocate (kEnfNumParameters);
		for (int j = 0 ; j < kSearchNumParameters ; j++)
		{
			in >> search_parameters[j]; 
		}
		for (int j = 0 ; j < kEnfNumParameters ; j++)
		{
			in >> enf_parameters[j]; 
		}
	}
	
	fSearchParameters.CopySymmetric();
	fEnforcementParameters.CopySymmetric();

	/* write out search parameter matrix */
	for (int i = 0; i < num_surfaces ; i++) 
        {
                for (int j = 0 ; j < num_surfaces ; j++)
                {
			dArrayT& search_parameters = fSearchParameters(i,j);
			dArrayT& enf_parameters = fEnforcementParameters(i,j);
			out << "(" << i << "," << j << ")" ;
			if (search_parameters.Length() 
				== kSearchNumParameters) {
			  for (int k = 0 ; k < kSearchNumParameters ; k++)
			  {
				out << search_parameters[k];
				out << '\n';
				out << enf_parameters[k];
			  }
			}
			out << '\n';
                }
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
			fSurfaces[i].SetContactStatus(fEnforcementParameters);
  		}
	}

	return changed;
}

bool ContactElementT::UpdateContactConfiguration(void)
{
        bool changed = fContactSearch->UpdateInteractions();
  	for (int i = 0; i < fSurfaces.Length(); i++) {
		fSurfaces[i].UpdateContactStatus(fEnforcementParameters);
  	}

        return changed;
}

