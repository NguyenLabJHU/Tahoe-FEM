/* $Id: ElementT.cpp,v 1.15 2002-11-25 07:24:51 paklein Exp $ */
#include "ElementT.h"

#include <iostream.h>
#include "ExceptionT.h"

namespace Tahoe {

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
			throw ExceptionT::kBadInputValue;
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
		case ElementT::kPenaltyContactElement3D:
			type = ElementT::kPenaltyContactElement3D;
			break;
		case ElementT::kTotLagrExternalField:
			type = ElementT::kTotLagrExternalField;
			break;
		case ElementT::kNonsingularContinuum:
		{
			cout << "\n operator>>ElementT::TypeT: type is not longer supported: " << i_type << endl;
			throw ExceptionT::kBadInputValue;
		}
		case ElementT::kMultiplierContactElement2D:
			type = ElementT::kMultiplierContactElement2D;
			break;
		case ElementT::kSimoFiniteStrain:
			type = ElementT::kSimoFiniteStrain;
			break;
		case ElementT::kStaggeredMultiScale:
			type = ElementT::kStaggeredMultiScale;
			break;
		case ElementT::kBridgingScale:
			type = ElementT::kBridgingScale;
			break;	
		case ElementT::kSimoQ1P0:
			type = ElementT::kSimoQ1P0;
			break;	
		case ElementT::kAdhesion:
			type = ElementT::kAdhesion;
			break;	
		case ElementT::kParticlePair:
			type = ElementT::kParticlePair;
			break;	
		default:
			cout << "\n operator>>ElementT::TypeT: unknown type: "
			<< i_type<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}

}
