/* $Id: ElementT.cpp,v 1.7 2002-06-08 20:20:13 paklein Exp $ */

#include "ElementT.h"

#include <iostream.h>
#include "ExceptionCodes.h"

/* stream extraction operator */ 
istream& operator>>(istream& in, ElementT::TypeT& type)
{
	int i_type;
	in >> i_type;
	switch (i_type)
	{
		case ElementT::kRod:
			type = ElementT::kRod;
			break;
		case ElementT::kElastic:
			type = ElementT::kElastic;
			break;
		case ElementT::kHyperElastic:
			type = ElementT::kHyperElastic;
			break;
		case ElementT::kLocalizing:
			type = ElementT::kLocalizing;
			break;
		case ElementT::kVariTri:
		{
			cout << "\n operator>>ElementT::TypeT: element type is not longer\n"
			     <<   "     supported. Support for changing number of elements is being re-\n"
			     <<   "     written: " << i_type << endl;
			throw eBadInputValue;
		}
		case ElementT::kSWDiamond:
			type = ElementT::kSWDiamond;
			break;
		case ElementT::kMixedSWDiamond:
			type = ElementT::kMixedSWDiamond;
			break;
		case ElementT::kUnConnectedRod:
			type = ElementT::kUnConnectedRod;
			break;
		case ElementT::kVirtualRod:
			type = ElementT::kVirtualRod;
			break;
		case ElementT::kVirtualSWDC:
			type = ElementT::kVirtualSWDC;
			break;
		case ElementT::kCohesiveSurface:
			type = ElementT::kCohesiveSurface;
			break;
		case ElementT::kThermalSurface:
			type = ElementT::kThermalSurface;
			break;
		case ElementT::kPenaltyContact:
			type = ElementT::kPenaltyContact;
			break;
		case ElementT::kBEMelement:
			type = ElementT::kBEMelement;
			break;
		case ElementT::kAugLagContact2D:
			type = ElementT::kAugLagContact2D;
			break;
		case ElementT::kTotLagHyperElastic:
			type = ElementT::kTotLagHyperElastic;
			break;
		case ElementT::kMeshFreeElastic:
			type = ElementT::kMeshFreeElastic;
			break;
		case ElementT::kMeshFreeFDElastic:
			type = ElementT::kMeshFreeFDElastic;
			break;
		case ElementT::kD2MeshFreeFDElastic:
			type = ElementT::kD2MeshFreeFDElastic;
			break;
		case ElementT::kLinearDiffusion:
			type = ElementT::kLinearDiffusion;
			break;
		case ElementT::kMFCohesiveSurface:
			type = ElementT::kMFCohesiveSurface;
			break;
		case ElementT::kACME_Contact:
			type = ElementT::kACME_Contact;
			break;
		case ElementT::kMultiplierContact3D:
			type = ElementT::kMultiplierContact3D;
			break;
		case ElementT::kPenaltyContactElement2D:
			type = ElementT::kPenaltyContactElement2D;
			break;
		case ElementT::kTotLagrExternalField:
			type = ElementT::kTotLagrExternalField;
			break;
		case ElementT::kNonsingularContinuum:
			type = ElementT::kNonsingularContinuum;
			break;
		case ElementT::kMultiplierContactElement2D:
			type = ElementT::kMultiplierContactElement2D;
			break;
		case ElementT::kSimoFiniteStrain:
			type = ElementT::kSimoFiniteStrain;
			break;
		case ElementT::kMultiScale:
			type = ElementT::kMultiScale;
			break;
		case ElementT::kCoarseScale:
			type = ElementT::kCoarseScale;
			break;
		case ElementT::kFinePhest:
			type = ElementT::kFinePhest;
			break;
		default:
			cout << "\n operator>>ElementT::TypeT: unknown type: "
			<< i_type<< endl;
			throw eBadInputValue;	
	}
	return in;
}
