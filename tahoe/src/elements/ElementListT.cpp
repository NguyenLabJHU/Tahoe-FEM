/* $Id: ElementListT.cpp,v 1.20.2.3 2002-04-30 08:21:59 paklein Exp $ */
/* created: paklein (04/20/1998) */

#include "ElementListT.h"
#include <iostream.h>
#include "fstreamT.h"
//#include "FEManagerT.h"
#include "StringT.h"
#include "ElementT.h"
#include "ElementSupportT.h"

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
#include "GeometryT.h"
#include "SimoFiniteStrainT.h"
#include "MultiScaleT.h"
#include "CoarseScaleT.h"
#include "FinePhestT.h"

/* contact */
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "ACME_Contact3DT.h"
#include "MultiplierContact3DT.h"
#include "MultiplierContactElement2DT.h"
#include "PenaltyContactElement2DT.h"

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
		ElementT::TypeT code;
		in >> group >> code;
		group--;
		
		/* no predefined field names */
		const FieldT* field = NULL;
		if (fSupport.Analysis() == GlobalT::kMultiField)
		{
			StringT name;
			in >> name;
			field = fSupport.Field(name);
			break;
		}
		else /* legacy - field set by analysis code */
		{
			switch (fSupport.Analysis())
			{
				case GlobalT::kPML:
				{
					cout << "\n ElementListT::EchoElementData: PML not fully implemented" << endl;
					throw eGeneralFail;
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
					throw eBadInputValue;
			}
		}
		
		/* check field */
		if (!field) {
			cout << "\n ElementListT::EchoElementData: could not resolve field" << endl;
			throw eGeneralFail;
		}

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
		out << "    eq. " << ElementT::kPenaltyContact     << ", penalty contact\n";
		out << "    eq. " << ElementT::kAugLagContact2D    << ", augmented Lagrangian contact\n";
		out << "    eq. " << ElementT::kTotLagHyperElastic << ", hyperelastic (total Lagrangian)\n";
		out << "    eq. " << ElementT::kMeshFreeElastic    << ", elastic with MLS displacements\n";
		out << "    eq. " << ElementT::kMeshFreeFDElastic  << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << ElementT::kD2MeshFreeFDElastic << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << ElementT::kLinearDiffusion    << ", linear diffusion element\n";
		out << "    eq. " << ElementT::kMFCohesiveSurface  << ", meshfree cohesive surface element\n";
		out << "    eq. " << ElementT::kMultiScale              << ", Variational Multi-Scale (VMS) Element \n";
		out << "    eq. " << ElementT::kCoarseScale             << ", Coarse Scale Element (for VMS) \n";
		out << "    eq. " << ElementT::kFinePhest               << ", Fine Sclale Phenomenological Str. Grad\n";

		out << "    eq. " << ElementT::kACME_Contact       << ", 3D contact using ACME\n";
		out << "    eq. " << ElementT::kMultiplierContact3D       << ", 3D contact using Lagrange multipliers\n";
		out << "    eq. " << ElementT::kMultiplierContactElement2D       << ", 2D Lagrange multiplier contact elements\n";
		out << "    eq. " << ElementT::kPenaltyContactElement2D       << ", 2D penalty contact elements\n";
		
		/* check */
		if (group < 0 || group >= Length())
		{
			cout << "\n ElementListT::EchoElementData: Element group number is out of\n";
			cout <<   "     range: " << group + 1 << endl;
			throw eBadInputValue;
		}

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

			case ElementT::kMultiScale:
				fArray[group] = new MultiScaleT(fSupport, *field);
				break;

			case ElementT::kCoarseScale:
				fArray[group] = new CoarseScaleT(fSupport, *field);
				break;

			case ElementT::kFinePhest:
				fArray[group] = new FinePhestT(fSupport, *field);
				break;

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
					throw eBadInputValue;
				}
				break;
			}
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
				throw eGeneralFail;					
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
			default:
			
				cout << "\n ElementListT::EchoElementData: unknown element type:" << code << endl;
				throw eBadInputValue;
		}
		
		if (!fArray[group]) throw eOutOfMemory;
//		fArray[group]->SetController(e_controller); //TEMP: this is dangerous. should pass
		fArray[group]->Initialize();                //      controller in during construction
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
