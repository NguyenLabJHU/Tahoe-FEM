/* $Id: ContactElementT.cpp,v 1.9 2001-04-27 00:55:25 rjones Exp $ */

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

/* element level reconfiguration for the current solution */
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
	cout << "\nTHROWING EXCEPTION in ContactElementT::Initialize"
	     << " to stop execution before getting into more trouble\n";
	throw; //HACK
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
	out << "\n Contact tracking: group " << fFEManager.ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << fFEManager.Time() << '\n';
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
	ElementBaseT::ConnectsU(connects_1, connects_2);

	/* link surfaces with fictious node-to-node pairs*/
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

	/* read contact surfaces */
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

	fSearchParameters.Allocate(num_surfaces);

	int num_pairs;
	in >> num_pairs;
	out << " Number of surface pairs with data . . . . . . . . = "
	    << num_pairs << '\n';
	int s1, s2;
	int num_param = kNumParameters;
	for (int i = 0; i < num_pairs ; i++) 
	{
		in >> s1 >> s2;
		s1--; s2--;
		dArrayT& search_parameters = fSearchParameters(s1,s2);
		search_parameters.Allocate(num_param);
		for (int j = 0 ; j < num_param ; j++)
		{
			in >> search_parameters[j]; 
		}
	}
	
	fSearchParameters.CopySymmetric();

	/* write out search parameter matrix */
	for (int i = 0; i < num_surfaces ; i++) 
        {
                for (int j = 0 ; j < num_surfaces ; j++)
                {
			dArrayT& search_parameters = fSearchParameters(i,j);
			out << "(" << i << "," << j << ")" ;
			if (search_parameters.Length() == num_param) {
			  for (int k = 0 ; k < num_param ; k++)
			  {
				out << search_parameters[k];
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
  		}
	}

	return changed;
}

bool ContactElementT::UpdateContactConfiguration(void)
{
        bool changed = fContactSearch->UpdateInteractions();

        return changed;
}

