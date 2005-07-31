/* $Id: ElementListT.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (04/20/1998)                                          */

#include "ElementListT.h"
#include <iostream.h>
#include "fstreamT.h"
#include "FEManagerT.h"
#include "XDOF_FDNodesT.h"
#include "StringT.h"

/* elements */
#include "ElementBaseT.h"
#include "RodT.h"
#include "ElasticT.h"
#include "UpLag_FDElasticT.h"
#include "LocalizerT.h"
#include "SWDiamondT.h"
#include "MixedSWDiamondT.h"
#include "UnConnectedRodT.h"
#include "VirtualRodT.h"
#include "VirtualSWDC.h"
#include "BEMelement.h"
#include "VariTriT.h" //TEMP
#include "TotLag_FDElasticT.h"
#include "CSEIsoT.h"
#include "CSEAnisoT.h"
#include "GeometryT.h"

/* contact */
#include "PenaltyContact2DT.h"
#include "PenaltyContact3DT.h"
#include "AugLagContact2DT.h"
#include "ACME_Contact3DT.h"

//TEMP
#include "MeshFreeElasticT.h"
#include "MeshFreeFDElasticT.h"
#include "D2MeshFreeFDElasticT.h"

/* linear diffusion */
#include "DiffusionT.h"

/* meshfree cohesive surface elements */
#include "MeshFreeCSEAnisoT.h"

/* Element Types */
const int kRod                = 1;
const int kElastic            = 2;
const int kHyperElastic       = 3;
const int kLocalizing         = 4;
const int kVariTri            = 5; //TEMP

const int kSWDiamond          = 6;
const int kMixedSWDiamond     = 7;
const int kUnConnectedRod     = 8;
const int kVirtualRod         = 9;
const int kVirtualSWDC        = 10;

const int kCohesiveSurface    = 11;
const int kPenaltyContact     = 14;
const int kBEMelement         = 15;
const int kAugLagContact2D    = 16;
const int kTotLagHyperElastic = 17;

const int kMeshFreeElastic    = 18;
const int kMeshFreeFDElastic  = 19;
const int kD2MeshFreeFDElastic = 20;

const int kLinearDiffusion    = 21;
const int kMFCohesiveSurface  = 22;

const int kACME_Contact       = 23;

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
		int	group, code;
		in >> group >> code;
		group--;

		out << "\n Group number. . . . . . . . . . . . . . . . . . = " << group + 1 << '\n';
		out <<   " Element type code . . . . . . . . . . . . . . . = " <<      code << '\n';
		out << "    eq. " << kRod                << ", rod\n";
		out << "    eq. " << kElastic            << ", elastic\n";
		out << "    eq. " << kHyperElastic       << ", hyperelastic\n";
		out << "    eq. " << kLocalizing         << ", hyperelastic with localization\n";
		out << "    eq. " << kSWDiamond          << ", diamond cubic lattice\n";   	
		out << "    eq. " << kMixedSWDiamond     << ", diamond cubic lattice with evolving params\n";   	
		out << "    eq. " << kUnConnectedRod     << ", self-connecting rods\n";   	
		out << "    eq. " << kVirtualRod         << ", self-connecting rods with periodic BC's\n";   	
		out << "    eq. " << kVirtualSWDC        << ", diamond cubic lattice with periodic BC's\n";   	
		out << "    eq. " << kCohesiveSurface    << ", cohesive surface element\n";   	
		out << "    eq. " << kPenaltyContact     << ", penalty contact\n";
		out << "    eq. " << kAugLagContact2D    << ", augmented Lagrangian contact\n";
		out << "    eq. " << kTotLagHyperElastic << ", hyperelastic (total Lagrangian)\n";
		out << "    eq. " << kMeshFreeElastic    << ", elastic with MLS displacements\n";
		out << "    eq. " << kMeshFreeFDElastic  << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << kD2MeshFreeFDElastic << ", hyperelastic MLS (total Lagrangian)\n";
		out << "    eq. " << kLinearDiffusion    << ", linear diffusion element\n";
		out << "    eq. " << kMFCohesiveSurface  << ", meshfree cohesive surface element\n";

		out << "    eq. " << kACME_Contact       << ", 3D contact using ACME\n";
		
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
			case kRod:

				fArray[group] = new RodT(fFEManager);
				break;

			case kElastic:
				fArray[group] = new ElasticT(fFEManager);
				break;

			case kMeshFreeElastic:
				fArray[group] = new MeshFreeElasticT(fFEManager);
				break;

			case kHyperElastic:
				fArray[group] = new UpLag_FDElasticT(fFEManager);
				break;

			case kTotLagHyperElastic:
				fArray[group] = new TotLag_FDElasticT(fFEManager);
				break;

			case kMeshFreeFDElastic:
				fArray[group] = new MeshFreeFDElasticT(fFEManager);
				break;

			case kD2MeshFreeFDElastic:
				fArray[group] = new D2MeshFreeFDElasticT(fFEManager);
				break;

			case kLocalizing:
				fArray[group] = new LocalizerT(fFEManager);
				break;

			case kVariTri:
				fArray[group] = new VariTriT(fFEManager);
				break;

			case kSWDiamond:
				fArray[group] = new SWDiamondT(fFEManager);
				break;
				
			case kMixedSWDiamond:
				fArray[group] = new MixedSWDiamondT(fFEManager);
				break;

			case kUnConnectedRod:
				fArray[group] = new UnConnectedRodT(fFEManager);
				break;
			
			case kVirtualRod:
				fArray[group] = new VirtualRodT(fFEManager);
				break;

			case kVirtualSWDC:
				fArray[group] = new VirtualSWDC(fFEManager);
				break;

			case kCohesiveSurface:
			{
				int CSEcode;
				in >> CSEcode;
				out << " Cohesive surface formulation code . . . . . . . = " << CSEcode << '\n';
				out << "    eq. " << CSEBaseT::Isotropic   << ", isotropic\n";
				out << "    eq. " << CSEBaseT::Anisotropic << ", anisotropic\n";

				if (CSEcode == CSEBaseT::Isotropic)
					fArray[group] = new CSEIsoT(fFEManager);	
				else if (CSEcode == CSEBaseT::Anisotropic)
					fArray[group] = new CSEAnisoT(fFEManager);
				else
				{
					cout << "\n ElementListT::EchoElementData: unknown CSE formulation: ";
					cout << CSEcode << '\n';
					throw eBadInputValue;
				}
				break;
			}
			case kPenaltyContact:
			{
				int nsd = (fFEManager.NodeManager())->NumSD();
				if (nsd == 2)
					fArray[group] = new PenaltyContact2DT(fFEManager);
				else
					fArray[group] = new PenaltyContact3DT(fFEManager);
					
				break;
			}
			case kAugLagContact2D:
			{
#ifdef __NO_RTTI__
				if (fFEManager.Analysis() != GlobalT::kAugLagStatic) throw eGeneralFail;
				XDOF_FDNodesT* XDOF_man = (XDOF_FDNodesT*) fFEManager.NodeManager();
#else
				XDOF_FDNodesT* XDOF_man = dynamic_cast<XDOF_FDNodesT*>(fFEManager.NodeManager());
				if (!XDOF_man)
				{
					cout << "\n ElementListT::EchoElementData: failed to cast node manager to XDOF_FDNodesT\n"
					     <<   "     as needed with analysis code: " << kAugLagContact2D << endl;
					throw eBadInputValue;
				}
#endif /* __NO_RTTI__ */
			
				fArray[group] = new AugLagContact2DT(fFEManager, XDOF_man);	
				break;
			}
			case kBEMelement:
			{
				StringT BEMfilename;
				in >> BEMfilename;
			
				fArray[group] = new BEMelement(fFEManager, BEMfilename);	
				break;
			}
			case kLinearDiffusion:
				fArray[group] = new DiffusionT(fFEManager);
				break;

			case kMFCohesiveSurface:
				fArray[group] = new MeshFreeCSEAnisoT(fFEManager);
				break;

			case kACME_Contact:
#ifdef __ACME__
				fArray[group] = new ACME_Contact3DT(fFEManager);
#else
				cout << "\n ElementListT::EchoElementData: ACME not installed.";
				throw eGeneralFail;					
#endif /* __ACME__ */			
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
