/* $Id: ElementListT.cpp,v 1.34 2002-12-01 19:50:48 paklein Exp $ */
/* created: paklein (04/20/1998) */
#include "ElementListT.h"

#include <iostream.h>
#include "fstreamT.h"
#include "StringT.h"
#include "ElementT.h"
#include "ElementSupportT.h"
#include "GeometryT.h"

/* elements */
#include "ElementBaseT.h"
#include "RodT.h"
#include "SmallStrainT.h"
#include "UpdatedLagrangianT.h"
#include "TotalLagrangianT.h"
#include "LocalizerT.h"
#include "SWDiamondT.h"
#include "MixedSWDiamondT.h"
#include "UnConnectedRodT.h"
#include "VirtualRodT.h"
#include "VirtualSWDC.h"
#include "BEMelement.h"
#include "CSEIsoT.h"
#include "CSEAnisoT.h"
#include "ThermalSurfaceT.h"
#include "SimoFiniteStrainT.h"
#include "StaggeredMultiScaleT.h"
#include "BridgingScaleT.h"
#include "SimoQ1P0.h"
#include "AdhesionT.h"

/* contact */
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "ACME_Contact3DT.h"
#include "MultiplierContact3DT.h"
#include "MultiplierContactElement2DT.h"
#include "PenaltyContactElement2DT.h"
#include "PenaltyContactElement3DT.h"

//TEMP
#include "MeshFreeElasticT.h"
#include "MeshFreeFDElasticT.h"
#include "D2MeshFreeFDElasticT.h"

/* linear diffusion */
#include "DiffusionT.h"

/* meshfree cohesive surface elements */
#include "MeshFreeCSEAnisoT.h"

/* class to read external field from file */
#include "UpLagr_ExternalFieldT.h"

/* particle classes */
#include "ParticlePairT.h"

using namespace Tahoe;

/* constructors */
ElementListT::ElementListT(void)
{

}

/* initialization functions */
void ElementListT::EchoElementData(ifstreamT& in, ostream& out, FEManagerT& fe)
{
	/* initialize element support */
	fSupport.SetFEManager(&fe);

	/* print header */
	out << "\n E l e m e n t   G r o u p   D a t a :\n\n";
	out << " Number of element groups. . . . . . . . . . . . = ";
	out << Length() << "\n";

	/* construct element groups */
	for (int i = 0; i < Length(); i++)
	{
		int	group;
		in >> group;
		group--;

		/* read code */
		ElementT::TypeT code;
		in >> code;

		/* check */
		if (group < 0 || group >= Length())
		{
			cout << "\n ElementListT::EchoElementData: Element group number is out of\n";
			cout <<   "     range: " << group + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}

		/* no over-writing existing groups */
		if (fArray[group]) {
			cout << "\n ElementListT::EchoElementData: group already exists" << group + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}

		/* no predefined field names */
		const FieldT* field = NULL;
		if (fSupport.Analysis() == GlobalT::kMultiField)
		{
			StringT name;
			in >> name;
			field = fSupport.Field(name);
		}
		else /* legacy - field set by analysis code */
		{
			switch (fSupport.Analysis())
			{
				case GlobalT::kPML:
				{
					cout << "\n ElementListT::EchoElementData: PML not fully implemented" << endl;
					throw ExceptionT::kGeneralFail;
				}
				case GlobalT::kLinStaticHeat:
				case GlobalT::kLinTransHeat:
				{
					field = fSupport.Field("temperature");
					break;
				}
				case GlobalT::kLinStatic:
				case GlobalT::kLinDynamic:
				case GlobalT::kLinExpDynamic:
				case GlobalT::kNLExpDynamic:
				case GlobalT::kNLExpDynKfield:
				case GlobalT::kDR:
				case GlobalT::kNLStatic:
				case GlobalT::kNLDynamic:
				case GlobalT::kNLStaticKfield:
				{
					field = fSupport.Field("displacement");
					break;
				}
				default:
					cout << "\n ElementListT::EchoElementData: recognized analysis type: "
					     << fSupport.Analysis() << endl;
					throw ExceptionT::kBadInputValue;
			}
		}
		
		/* check field */
		if (!field) {
			cout << "\n ElementListT::EchoElementData: could not resolve field" << endl;
			throw ExceptionT::kGeneralFail;
		}

		/* write parameters */
		out << "\n Group number. . . . . . . . . . . . . . . . . . = " << group + 1 << '\n';
		out <<   " Element type code . . . . . . . . . . . . . . . = " <<      code << '\n';
		out << "    eq. " << ElementT::kRod                << ", rod\n";
		out << "    eq. " << ElementT::kElastic            << ", elastic\n";
		out << "    eq. " << ElementT::kHyperElastic       << ", hyperelastic\n";
		out << "    eq. " << ElementT::kLocalizing         << ", hyperelastic with localization\n";
		out << "    eq. " << ElementT::kSWDiamond          << ", diamond cubic lattice\n";   	
		out << "    eq. " << ElementT::kMixedSWDiamond     << ", diamond cubic lattice with evolving params\n";   	
		out << "    eq. " << ElementT::kUnConnectedRod     << ", self-connecting rods\n";   	
		out << "    eq. " << ElementT::kVirtualRod         << ", self-connecting rods with periodic BC's\n";   	
		out << "    eq. " << ElementT::kVirtualSWDC        << ", diamond cubic lattice with periodic BC's\n";   	
		out << "    eq. " << ElementT::kCohesiveSurface    << ", cohesive surface element\n";   	
		out << "    eq. " << ElementT::kThermalSurface    << ", thermal surface element\n";   	
		out << "    eq. " << ElementT::kPenaltyContact     << ", penalty contact\n";
		out << "    eq. " << ElementT::kAugLagContact2D    << ", augmented Lagrangian contact\n";
		out << "    eq. " << ElementT::kTotLagHyperElastic << ", hyperelastic (total Lagrangian)\n";
		out << "    eq. " << ElementT::kMeshFreeElastic    << ", elastic with MLS displacements\n";
		out << "    eq. " << ElementT::kMeshFreeFDElastic  << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << ElementT::kD2MeshFreeFDElastic << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << ElementT::kLinearDiffusion    << ", linear diffusion element\n";
		out << "    eq. " << ElementT::kMFCohesiveSurface  << ", meshfree cohesive surface element\n";
		out << "    eq. " << ElementT::kStaggeredMultiScale << ", Staggered MultiScale Element (for VMS) \n";
		out << "    eq. " << ElementT::kACME_Contact       << ", 3D contact using ACME\n";
		out << "    eq. " << ElementT::kMultiplierContact3D       << ", 3D contact using Lagrange multipliers\n";
		out << "    eq. " << ElementT::kMultiplierContactElement2D       << ", 2D Lagrange multiplier contact elements\n";
		out << "    eq. " << ElementT::kPenaltyContactElement2D       << ", 2D penalty contact elements\n";
		out << "    eq. " << ElementT::kPenaltyContactElement3D       << ", 3D penalty contact elements\n";
		out << "    eq. " << ElementT::kBridgingScale      << ", Bridging Scale\n";
		out << "    eq. " << ElementT::kSimoQ1P0           << ", Q1P0 mixed element\n";
		out << "    eq. " << ElementT::kAdhesion           << ", surface adhesion\n";

		/* create new element group - read control parameters
		   and allocate space in the contructors */
		switch (code)
		{
			case ElementT::kRod:
			{
				fArray[group] = new RodT(fSupport, *field);
				break;
			}
			case ElementT::kElastic:
				fArray[group] = new SmallStrainT(fSupport, *field);
				break;

			case ElementT::kMeshFreeElastic:
				fArray[group] = new MeshFreeElasticT(fSupport, *field);
				break;

			case ElementT::kHyperElastic:
				fArray[group] = new UpdatedLagrangianT(fSupport, *field);
				break;

			case ElementT::kTotLagHyperElastic:
				fArray[group] = new TotalLagrangianT(fSupport, *field);
				break;

			case ElementT::kSimoFiniteStrain:
				fArray[group] = new SimoFiniteStrainT(fSupport, *field);
				break;

			case ElementT::kSimoQ1P0:
				fArray[group] = new SimoQ1P0(fSupport, *field);
				break;

			case ElementT::kStaggeredMultiScale:
			{
				/* must be using multi-field solver */
				if (fSupport.Analysis() != GlobalT::kMultiField)				
					ExceptionT::BadInputValue("ElementListT::EchoElementData", "multi field required");
			
				/* coarse scale field read above */
				const FieldT* coarse_scale = field;

				/* fine scale field */				
				StringT fine_field_name;
				in >> fine_field_name;
				const FieldT* fine_scale = fSupport.Field(fine_field_name);
				if (!coarse_scale || !fine_scale)
					ExceptionT::BadInputValue("ElementListT::EchoElementData", "error resolving field names");
			
				fArray[group] = new StaggeredMultiScaleT(fSupport, *coarse_scale, *fine_scale);
				break;
			}
			case ElementT::kMeshFreeFDElastic:
				fArray[group] = new MeshFreeFDElasticT(fSupport, *field);
				break;

			case ElementT::kD2MeshFreeFDElastic:
				fArray[group] = new D2MeshFreeFDElasticT(fSupport, *field);
				break;

			case ElementT::kLocalizing:
				fArray[group] = new LocalizerT(fSupport, *field);
				break;

			case ElementT::kSWDiamond:
				fArray[group] = new SWDiamondT(fSupport, *field);
				break;
				
			case ElementT::kMixedSWDiamond:
				fArray[group] = new MixedSWDiamondT(fSupport, *field);
				break;

			case ElementT::kUnConnectedRod:
				fArray[group] = new UnConnectedRodT(fSupport, *field);
				break;
			
			case ElementT::kVirtualRod:
				fArray[group] = new VirtualRodT(fSupport, *field);
				break;

			case ElementT::kVirtualSWDC:
				fArray[group] = new VirtualSWDC(fSupport, *field);
				break;

			case ElementT::kCohesiveSurface:
			{
				int CSEcode;
				in >> CSEcode;
				out << " Cohesive surface formulation code . . . . . . . = " << CSEcode << '\n';
				out << "    eq. " << CSEBaseT::Isotropic   << ", isotropic\n";
				out << "    eq. " << CSEBaseT::Anisotropic << ", anisotropic\n";
				out << "    eq. " << CSEBaseT::NoRotateAnisotropic << ", fixed-frame anisotropic\n";

				if (CSEcode == CSEBaseT::Isotropic)
					fArray[group] = new CSEIsoT(fSupport, *field);	
				else if (CSEcode == CSEBaseT::Anisotropic)
					fArray[group] = new CSEAnisoT(fSupport, *field, true);
				else if (CSEcode == CSEBaseT::NoRotateAnisotropic)
					fArray[group] = new CSEAnisoT(fSupport, *field, false);
				else
				{
					cout << "\n ElementListT::EchoElementData: unknown CSE formulation: ";
					cout << CSEcode << '\n';
					throw ExceptionT::kBadInputValue;
				}
				break;
			}

			case ElementT::kThermalSurface:
				fArray[group] = new ThermalSurfaceT(fSupport, *field);	
				break;

			case ElementT::kPenaltyContact:
			{
				int nsd = fSupport.NumSD();
				if (nsd == 2)
					fArray[group] = new PenaltyContact2DT(fSupport, *field);
				else
					fArray[group] = new PenaltyContact3DT(fSupport, *field);
					
				break;
			}
			case ElementT::kAugLagContact2D:
			{
				fArray[group] = new AugLagContact2DT(fSupport, *field);	
				break;
			}
			case ElementT::kBEMelement:
			{
				StringT BEMfilename;
				in >> BEMfilename;
			
				fArray[group] = new BEMelement(fSupport, *field, BEMfilename);	
				break;
			}
			case ElementT::kLinearDiffusion:
				fArray[group] = new DiffusionT(fSupport, *field);
				break;

			case ElementT::kMFCohesiveSurface:
				fArray[group] = new MeshFreeCSEAnisoT(fSupport, *field);
				break;

			case ElementT::kTotLagrExternalField:
				fArray[group] = new UpLagr_ExternalFieldT(fSupport, *field);
				break;
			case ElementT::kACME_Contact:
#ifdef __ACME__
				fArray[group] = new ACME_Contact3DT(fSupport, *field);
#else
				cout << "\n ElementListT::EchoElementData: ACME not installed.";
				throw ExceptionT::kGeneralFail;					
#endif /* __ACME__ */			
				break;

			case ElementT::kMultiplierContact3D:
				fArray[group] = new MultiplierContact3DT(fSupport, *field);
				break;

			case ElementT::kMultiplierContactElement2D:
			{
				fArray[group] = new MultiplierContactElement2DT(fSupport, *field);
				break;
			}
		case ElementT::kPenaltyContactElement2D:
		  {
		    fArray[group] = new PenaltyContactElement2DT(fSupport, *field);
		    break;
		  }
		case ElementT::kPenaltyContactElement3D:
		  {
		    fArray[group] = new PenaltyContactElement3DT(fSupport, *field);
		    break;
		  }
		case ElementT::kBridgingScale:
		{
			/* associated group numbers */
			int particle_group = -99;
			int solid_group = -99;
			in >> particle_group >> solid_group;
			const RodT* particle = dynamic_cast<const RodT*>(&(fSupport.ElementGroup(--particle_group)));
			if (!particle) {
				cout << "\n ElementListT::EchoElementData: unable to cast pointer to group " << particle_group+1 << '\n'
				     <<   "     to type RodT" << endl;
				throw ExceptionT::kBadInputValue;
			}
			const ElasticT* solid = dynamic_cast<const ElasticT*>(&(fSupport.ElementGroup(--solid_group)));
			if (!solid) {
				cout << "\n ElementListT::EchoElementData: unable to cast pointer to group " << solid_group+1 << '\n'
				     <<   "     to type ElasticT" << endl;
				throw ExceptionT::kBadInputValue;
			}
			fArray[group] = new BridgingScaleT(fSupport, *field, *particle, *solid);
		    break;
		}
		case ElementT::kAdhesion:
		{
			fArray[group] = new AdhesionT(fSupport, *field);
			break;
		}
		case ElementT::kParticlePair:
		{
			fArray[group] = new ParticlePairT(fSupport, *field);
			break;
		}
		default:
		  
		  cout << "\n ElementListT::EchoElementData: unknown element type:" << code << endl;
		  throw ExceptionT::kBadInputValue;
		}
		
		if (!fArray[group]) throw ExceptionT::kOutOfMemory;
		fArray[group]->Initialize();
	}
}

/* returns true of ALL element groups have interpolant DOF's */
bool ElementListT::InterpolantDOFs(void) const
{
	bool are_interp = true;
	for (int i = 0; i < Length(); i++)
		if (!fArray[i]->InterpolantDOFs()) are_interp = false;

	return are_interp;
}

/* returns true if contact group present */
bool ElementListT::HasContact(void) const
{
#ifdef __NO_RTTI__
	cout << "\n ElementListT::HasContact: needs RTTI: returning false" << endl;
	return false; // safe default??
#else
	for (int i = 0; i < Length(); i++)
	{
		ContactT* contact_elem = dynamic_cast<ContactT*>(fArray[i]);
		if (contact_elem) return true;
	}
	return false;
#endif
}
