/* $Id: ElementListT.cpp,v 1.19 2002-03-21 22:36:32 creigh Exp $ */
/* created: paklein (04/20/1998) */

#include "ElementListT.h"
#include <iostream.h>
#include "fstreamT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "StringT.h"
#include "ElementT.h"

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
#include "CourseScaleT.h"
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
ElementListT::ElementListT(FEManagerT& fe_manager):
	fFEManager(fe_manager)
{

}

/* initialization functions */
void ElementListT::EchoElementData(ifstreamT& in, ostream& out,
	eControllerT* e_controller)
{
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
		out << "    eq. " << ElementT::kCourseScale             << ", Course Scale Element (for VMS) \n";
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

				fArray[group] = new RodT(fFEManager);
				break;

			case ElementT::kElastic:
				fArray[group] = new SmallStrainT(fFEManager);
				break;

			case ElementT::kMeshFreeElastic:
				fArray[group] = new MeshFreeElasticT(fFEManager);
				break;

			case ElementT::kHyperElastic:
				fArray[group] = new UpdatedLagrangianT(fFEManager);
				break;

			case ElementT::kTotLagHyperElastic:
				fArray[group] = new TotalLagrangianT(fFEManager);
				break;

			case ElementT::kSimoFiniteStrain:
				fArray[group] = new SimoFiniteStrainT(fFEManager);
				break;

			case ElementT::kMultiScale:
				fArray[group] = new MultiScaleT(fFEManager);
				break;

			case ElementT::kCourseScale:
				fArray[group] = new CourseScaleT(fFEManager);
				break;

			case ElementT::kFinePhest:
				fArray[group] = new FinePhestT(fFEManager);
				break;

			case ElementT::kMeshFreeFDElastic:
				fArray[group] = new MeshFreeFDElasticT(fFEManager);
				break;

			case ElementT::kD2MeshFreeFDElastic:
				fArray[group] = new D2MeshFreeFDElasticT(fFEManager);
				break;

			case ElementT::kLocalizing:
				fArray[group] = new LocalizerT(fFEManager);
				break;

			case ElementT::kSWDiamond:
				fArray[group] = new SWDiamondT(fFEManager);
				break;
				
			case ElementT::kMixedSWDiamond:
				fArray[group] = new MixedSWDiamondT(fFEManager);
				break;

			case ElementT::kUnConnectedRod:
				fArray[group] = new UnConnectedRodT(fFEManager);
				break;
			
			case ElementT::kVirtualRod:
				fArray[group] = new VirtualRodT(fFEManager);
				break;

			case ElementT::kVirtualSWDC:
				fArray[group] = new VirtualSWDC(fFEManager);
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
					fArray[group] = new CSEIsoT(fFEManager);	
				else if (CSEcode == CSEBaseT::Anisotropic)
					fArray[group] = new CSEAnisoT(fFEManager, true);
				else if (CSEcode == CSEBaseT::NoRotateAnisotropic)
					fArray[group] = new CSEAnisoT(fFEManager, false);
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
				int nsd = (fFEManager.NodeManager())->NumSD();
				if (nsd == 2)
					fArray[group] = new PenaltyContact2DT(fFEManager);
				else
					fArray[group] = new PenaltyContact3DT(fFEManager);
					
				break;
			}
			case ElementT::kAugLagContact2D:
			{
				fArray[group] = new AugLagContact2DT(fFEManager);	
				break;
			}
			case ElementT::kBEMelement:
			{
				StringT BEMfilename;
				in >> BEMfilename;
			
				fArray[group] = new BEMelement(fFEManager, BEMfilename);	
				break;
			}
			case ElementT::kLinearDiffusion:
				fArray[group] = new DiffusionT(fFEManager);
				break;

			case ElementT::kMFCohesiveSurface:
				fArray[group] = new MeshFreeCSEAnisoT(fFEManager);
				break;

			case ElementT::kTotLagrExternalField:
				fArray[group] = new UpLagr_ExternalFieldT(fFEManager);
				break;

#if 0
			case ElementT::kNonsingularContinuum:
				fArray[group] = new NonsingularContinuumT(fFEManager);
				break;
#endif

			case ElementT::kACME_Contact:
#ifdef __ACME__
				fArray[group] = new ACME_Contact3DT(fFEManager);
#else
				cout << "\n ElementListT::EchoElementData: ACME not installed.";
				throw eGeneralFail;					
#endif /* __ACME__ */			
				break;

			case ElementT::kMultiplierContact3D:
				fArray[group] = new MultiplierContact3DT(fFEManager);
				break;

			case ElementT::kMultiplierContactElement2D:
			{
#ifdef __NO_RTTI__
				if (fFEManager.Analysis() != GlobalT::kAugLagStatic) throw eGeneralFail;
				XDOF_ManagerT* XDOF_man	= (XDOF_ManagerT*) fFEManager.NodeManager();
#else
				XDOF_ManagerT* XDOF_man = dynamic_cast<XDOF_ManagerT*>(fFEManager.NodeManager());
				if (!XDOF_man)
				{
				cout<< "\n ElementListT::EchoElementData: "
				      << "failed to cast node manager to "
			              << "XDOF_ManagerT\n"
				      << "     as needed with analysis code: " 
				      << ElementT::kMultiplierContactElement2D << endl;
				throw eBadInputValue;
				}
#endif /* __NO_RTTI__ */

				fArray[group]
				= new MultiplierContactElement2DT
				          (fFEManager,XDOF_man);
				break;
			}

            case ElementT::kPenaltyContactElement2D:
                fArray[group] = new PenaltyContactElement2DT(fFEManager);
                break;


			default:
			
				cout << "\n ElementListT::EchoElementData: unknown element type:" << code << endl;
				throw eBadInputValue;
		}
		
		if (!fArray[group]) throw eOutOfMemory;
		fArray[group]->SetController(e_controller); //TEMP: this is dangerous. should pass
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
