/* $Id: ElementListT.cpp,v 1.36 2003-01-28 22:35:25 paklein Exp $ */
/* created: paklein (04/20/1998) */
#include "ElementListT.h"
#include "ElementsConfig.h"

#include <iostream.h>
#include "fstreamT.h"
#include "StringT.h"
#include "ElementT.h"
#include "ElementSupportT.h"
#include "GeometryT.h"

/* elements */
#include "ElementBaseT.h"
#include "RodT.h"
#include "UnConnectedRodT.h"
#include "VirtualRodT.h"
#include "BEMelement.h"
#include "AdhesionT.h"

#ifdef MULTISCALE_ELEMENT
#include "StaggeredMultiScaleT.h"
#endif

#ifdef COHESIVE_SURFACE_ELEMENT
#include "CSEIsoT.h"
#include "CSEAnisoT.h"
#include "MeshFreeCSEAnisoT.h"
#include "ThermalSurfaceT.h"
#endif

#ifdef CONTINUUM_ELEMENT
#include "SmallStrainT.h"
#include "UpdatedLagrangianT.h"
#include "TotalLagrangianT.h"
#include "LocalizerT.h"
#include "SimoFiniteStrainT.h"
#include "SimoQ1P0.h"
#include "DiffusionT.h"
#include "MeshFreeElasticT.h"
#include "MeshFreeFDElasticT.h"
#include "D2MeshFreeFDElasticT.h"
#include "UpLagr_ExternalFieldT.h"
#include "BridgingScaleT.h"
#endif

#ifdef DEV_CONTACT_ELEMENT
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "ACME_Contact3DT.h"
#include "MultiplierContact3DT.h"
#include "MultiplierContactElement2DT.h"
#include "PenaltyContactElement2DT.h"
#include "PenaltyContactElement3DT.h"
#endif

#ifdef PARTICLE_ELEMENT
#include "ParticlePairT.h"
#endif

#ifdef LATTICE_ELEMENT
#include "SWDiamondT.h"
#include "MixedSWDiamondT.h"
#include "VirtualSWDC.h"
#endif

using namespace Tahoe;

/* constructors */
ElementListT::ElementListT(void)
{

}

/* initialization functions */
void ElementListT::EchoElementData(ifstreamT& in, ostream& out, FEManagerT& fe)
{
	const char caller[] = "ElementListT::EchoElementData";

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
			ExceptionT::BadInputValue(caller, "element group %d is out of range", group+1);

		/* no over-writing existing groups */
		if (fArray[group])
			ExceptionT::BadInputValue(caller, "group %d already exists", group+1);

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
					ExceptionT::BadInputValue(caller, "PML not fully implemented");
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
					ExceptionT::BadInputValue(caller, "unrecognized analysis type: %d", fSupport.Analysis());
			}
		}
		
		/* check field */
		if (!field)
			ExceptionT::BadInputValue(caller, "could not resolve field for group %d", group+1);

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
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new SmallStrainT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kMeshFreeElastic:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new MeshFreeElasticT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kHyperElastic:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new UpdatedLagrangianT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kTotLagHyperElastic:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new TotalLagrangianT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kSimoFiniteStrain:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new SimoFiniteStrainT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kSimoQ1P0:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new SimoQ1P0(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kStaggeredMultiScale:
			{
#ifdef MULTISCALE_ELEMENT
				/* must be using multi-field solver */
				if (fSupport.Analysis() != GlobalT::kMultiField)				
					ExceptionT::BadInputValue(caller, "multi field required");
			
				/* coarse scale field read above */
				const FieldT* coarse_scale = field;

				/* fine scale field */				
				StringT fine_field_name;
				in >> fine_field_name;
				const FieldT* fine_scale = fSupport.Field(fine_field_name);
				if (!coarse_scale || !fine_scale)
					ExceptionT::BadInputValue(caller, "error resolving field names");
			
				fArray[group] = new StaggeredMultiScaleT(fSupport, *coarse_scale, *fine_scale);
				break;
#else
				ExceptionT::BadInputValue(caller, "MULTISCALE_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kMeshFreeFDElastic:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new MeshFreeFDElasticT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kD2MeshFreeFDElastic:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new D2MeshFreeFDElasticT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kLocalizing:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new LocalizerT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kSWDiamond:
			{
#ifdef LATTICE_ELEMENT
				fArray[group] = new SWDiamondT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "LATTICE_ELEMENT not enabled: %d", code);
#endif
			}	
			case ElementT::kMixedSWDiamond:
			{
#ifdef LATTICE_ELEMENT
				fArray[group] = new MixedSWDiamondT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "LATTICE_ELEMENT not enabled: %d", code);
#endif
			}	
			case ElementT::kUnConnectedRod:
				fArray[group] = new UnConnectedRodT(fSupport, *field);
				break;
			
			case ElementT::kVirtualRod:
				fArray[group] = new VirtualRodT(fSupport, *field);
				break;

			case ElementT::kVirtualSWDC:
			{
#ifdef LATTICE_ELEMENT
				fArray[group] = new VirtualSWDC(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "LATTICE_ELEMENT not enabled: %d", code);
#endif
			}	
			case ElementT::kCohesiveSurface:
			{
#ifdef COHESIVE_SURFACE_ELEMENT
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
					ExceptionT::BadInputValue(caller, "unknown CSE formulation: %d", CSEcode);
				}
				break;
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kThermalSurface:
			{
#ifdef COHESIVE_SURFACE_ELEMENT
				fArray[group] = new ThermalSurfaceT(fSupport, *field);	
				break;
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kPenaltyContact:
			{
#ifdef DEV_CONTACT_ELEMENT
				int nsd = fSupport.NumSD();
				if (nsd == 2)
					fArray[group] = new PenaltyContact2DT(fSupport, *field);
				else
					fArray[group] = new PenaltyContact3DT(fSupport, *field);
					
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kAugLagContact2D:
			{
#ifdef DEV_CONTACT_ELEMENT
				fArray[group] = new AugLagContact2DT(fSupport, *field);	
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kBEMelement:
			{
				StringT BEMfilename;
				in >> BEMfilename;
			
				fArray[group] = new BEMelement(fSupport, *field, BEMfilename);	
				break;
			}
			case ElementT::kLinearDiffusion:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new DiffusionT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kMFCohesiveSurface:
			{
#ifdef COHESIVE_SURFACE_ELEMENT
				fArray[group] = new MeshFreeCSEAnisoT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT not enabled: %d", code);
#endif
			}
			case ElementT::kTotLagrExternalField:
			{
#ifdef CONTINUUM_ELEMENT
				fArray[group] = new UpLagr_ExternalFieldT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif
			}				
			case ElementT::kACME_Contact:
			{
#ifdef DEV_CONTACT_ELEMENT
#ifdef __ACME__
				fArray[group] = new ACME_Contact3DT(fSupport, *field);
#else
				ExceptionT::GeneralFail(caller, "ACME not installed");
#endif /* __ACME__ */			
				break;
#else /* DEV_CONTACT_ELEMENT */
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif				
			}
			case ElementT::kMultiplierContact3D:
			{
#ifdef DEV_CONTACT_ELEMENT
				fArray[group] = new MultiplierContact3DT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif				
			}
			case ElementT::kMultiplierContactElement2D:
			{
#ifdef DEV_CONTACT_ELEMENT
				fArray[group] = new MultiplierContactElement2DT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif				
			}
			case ElementT::kPenaltyContactElement2D:
			{
#ifdef DEV_CONTACT_ELEMENT
				fArray[group] = new PenaltyContactElement2DT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif				
			}
			case ElementT::kPenaltyContactElement3D:
			{
#ifdef DEV_CONTACT_ELEMENT
				fArray[group] = new PenaltyContactElement3DT(fSupport, *field);
				break;
#else
				ExceptionT::BadInputValue(caller, "DEV_CONTACT_ELEMENT not enabled: %d", code);
#endif				
			}
		case ElementT::kBridgingScale:
		{
#ifdef CONTINUUM_ELEMENT
			/* associated group numbers */
			int particle_group = -99;
			int solid_group = -99;
			in >> particle_group >> solid_group;

			const RodT* particle = dynamic_cast<const RodT*>(&(fSupport.ElementGroup(--particle_group)));
			if (!particle)
				ExceptionT::BadInputValue(caller, "unable to cast pointer to group %d to type RodT", particle_group+1);

			const ElasticT* solid = dynamic_cast<const ElasticT*>(&(fSupport.ElementGroup(--solid_group)));
			if (!solid)
				ExceptionT::BadInputValue(caller, "unable to cast pointer to group %d to type ElasticT", solid_group+1);

			fArray[group] = new BridgingScaleT(fSupport, *field, *particle, *solid);
		    break;
#else
				ExceptionT::BadInputValue(caller, "CONTINUUM_ELEMENT not enabled: %d", code);
#endif				
		}
		case ElementT::kAdhesion:
		{
			fArray[group] = new AdhesionT(fSupport, *field);
			break;
		}
		case ElementT::kParticlePair:
		{
#ifdef PARTICLE_ELEMENT
			fArray[group] = new ParticlePairT(fSupport, *field);
			break;
#else
			ExceptionT::BadInputValue(caller, "PARTICLE_ELEMENT not enabled: %d", code);
#endif				
		}
		default:
			ExceptionT::BadInputValue(caller, "unknown element type: %d", code);
		}
		
		if (!fArray[group]) ExceptionT::OutOfMemory();
		fArray[group]->Initialize();
	}

#ifdef PARTICLE_ELEMENT
	/* can only have one particle group */
	int count = 0;
	for (int i = 0; i < Length(); i++)
	{
		ElementBaseT* e_group = fArray[i];
		ParticleT* particle = dynamic_cast<ParticleT*>(e_group);
		if (particle) count++;
	}
	if (count > 1) ExceptionT::BadInputValue(caller, "only one particle group allowed: %d", count);
#endif
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
#else /* __NO_RTTI__ */
#ifdef DEV_CONTACT_ELEMENT
	for (int i = 0; i < Length(); i++)
	{
		ContactT* contact_elem = dynamic_cast<ContactT*>(fArray[i]);
		if (contact_elem) return true;
	}
	return false;
#else  /* DEV_CONTACT_ELEMENT */
	return false;
#endif /* DEV_CONTACT_ELEMENT */
#endif /* __NO_RTTI__ */
}
