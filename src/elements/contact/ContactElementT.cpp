/* $Id: ContactElementT.cpp,v 1.3 2001-04-11 14:48:57 rjones Exp $ */

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


/* constructor */
ContactElementT::ContactElementT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager)
{
}

/* destructor */
ContactElementT::~ContactElementT(void) {	}

/* form of tangent matrix */
GlobalT::SystemTypeT ContactElementT::TangentType(void) const
{
	return GlobalT::kSymmetric; //HACK
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

	/* set up work space */
//SetWorkSpace();

	/* initialize surfaces */
	for (int i = 0; i < fSurfaces.Length(); i++) {
		SurfaceT& surface = fSurfaces[i];
		surface.Initialize(ElementBaseT::fNodes);
	}
	
	/* set initial contact configuration */
	SetContactConfiguration();	
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

/* appends group connectivities to the array */
void ContactElementT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	ElementBaseT::ConnectsU(connects_1, connects_2);

	/* add surface links */
	connects_1.AppendUnique(&fSurfaceLinks);
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
		SurfaceT& surface = fSurfaces[i];
		switch (spec_mode)
		{
			case kSideSets:
				surface.InputSideSets
				    (ElementBaseT::FEManager(), in, out);
				break;
			
			default:
				cout << "\n ContactElementT::EchoSurfaceData:"
                                     << " unknown surface specification\n";
				cout <<   "     mode " << spec_mode 
                                     << " for surface " << i+1 << '\n';
				throw eBadInputValue;
		}
		surface.PrintData(out);
	}
}

void ContactElementT::SetWorkSpace(void)
{
//HACK what is going on here????
	/* set the managed array - can only be set once */
//fConnectivities_man.SetWard(0, fConnectivities, fNumElemNodes);
//fEqnos_man.SetWard(0, fEqnos, fNumElemEqnos);

	/* make pseudo-element list to link surfaces in case
	 * bodies are not otherwise interacting (for the bandwidth
	 * reduction) */
	int num_surfaces = fSurfaces.Length();
	if (num_surfaces > 1)
	{
		fSurfaceLinks.Allocate(num_surfaces - 1, 2);
		for (int i = 0; i < num_surfaces - 1; i++)
		{
//fSurfaceLinks(i,0) = (fSurfaces[i  ])[0];
//fSurfaceLinks(i,1) = (fSurfaces[i+1])[0];
		}
	}
}

/* generate contact element data - return true if configuration has
* changed since the last call */
bool ContactElementT::SetContactConfiguration(void)
{
	bool contact_changed = 1;
//bool contact_changed = SetActiveInteractions();
	if (contact_changed)
	{
	}
	
	return contact_changed;
}

/***********************************************************************
* Private
***********************************************************************/
