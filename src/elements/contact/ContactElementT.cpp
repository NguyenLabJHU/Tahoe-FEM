/* $Id: ContactElementT.cpp,v 1.1 2001-03-22 17:57:42 rjones Exp $ */

#include "ContactT.h"

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

/* constructor */
ContactT::ContactT(FEManagerT& fe_manager, int numfacetnodes):
	ElementBaseT(fe_manager),
	fNumFacetNodes(numfacetnodes)
{
	/* override base class parameters */
	fNumElemNodes = fNumFacetNodes + 1; // facet nodes + 1 striker for each
	      fNumDOF = fNumSD;             // contact interaction
	fNumElemEqnos = fNumElemNodes*fNumDOF;
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

/* solution calls */
void ContactT::AddNodalForce(int node, dArrayT& force)
{
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
void ContactT::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)
	/* contact statistics */
	ostream& out = fFEManager.Output();
	out << "\n Contact tracking: group " << fFEManager.ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << fFEManager.Time() << '\n';
	out << " Active strikers                = " << fActiveStrikers.Length() << '\n';
	if (fActiveStrikers.Length() > 0)
	{
		out << setw(kIntWidth) << "striker";
		out << setw(kIntWidth) << "surface";
		out << setw(kIntWidth) << "facet";
		out << setw(fNumFacetNodes*kIntWidth) << "facet nodes" << '\n';
		for (int i = 0; i < fActiveStrikers.Length(); i++)
		{
			out << setw(kIntWidth) << fActiveStrikers[i] + 1;
			out << setw(kIntWidth) << fHitSurface[i] + 1; int hit_facet = fHitFacets[i];
			out << setw(kIntWidth) << hit_facet + 1;
			for (int j = 0; j < fNumFacetNodes; j++)
//out << setw(kIntWidth) << fSurfaces[fHitSurface[i]](hit_facet, j) + 1;
			out << '\n';
		}
		out << endl;
	}
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

/* echo contact surfaces */
void ContactT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces;
	in >> num_surfaces;
	out << " Number of contact surfaces. . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) throw eBadInputValue;

	/* read contact surfaces */
	fSurfaces.Allocate(num_surfaces); // need to call Initialize
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		int spec_mode;
		in >> spec_mode;
		switch (spec_mode)
		{
			case kSideSets:
				fSurfaces[i].InputSideSets
				    (ElementBaseT::FEManager(),in, out);
				break;
			
			default:
				cout << "\n ContactT::EchoSurfaceData:"
                                     << " unknown surface specification\n";
				cout <<   "     mode " << spec_mode 
                                     << " for surface " << i+1 << '\n';
				throw eBadInputValue;
		}
	}
	
	throw ; // HACK
	/* allocate contact nodes coords */ // MOVE THIS
	fStrikerCoords.Allocate(fStrikerTags.Length(),fNumSD);
}

void ContactT::SetWorkSpace(void)
{
	/* allocate map to active strikers data */
	fActiveMap.Allocate(fStrikerCoords.MajorDim());
	fActiveMap = -1;

	/* set the managed array - can only be set once */
	fConnectivities_man.SetWard(0, fConnectivities, fNumElemNodes);
	fEqnos_man.SetWard(0, fEqnos, fNumElemEqnos);

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
bool ContactT::SetContactConfiguration(void)
{
	bool contact_changed = 1;
//bool contact_changed = SetActiveInteractions();
	if (contact_changed)
	{
		/* resize */
		int num_active = fActiveStrikers.Length();
		fConnectivities_man.SetMajorDimension(num_active, false);
		fEqnos_man.SetMajorDimension(num_active, false);

		/* generate connectivities */
//SetConnectivities();	
	}
	
	return contact_changed;
}

/***********************************************************************
* Private
***********************************************************************/
