/* $Id: MeshFreeCSEAnisoT.cpp,v 1.23 2004-06-17 07:13:20 paklein Exp $ */
/* created: paklein (06/08/2000) */

#include "MeshFreeCSEAnisoT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "toolboxConstants.h"
#include "ElementSupportT.h"
#include "eIntegratorT.h"
#include "MeshFreeSurfaceShapeT.h"

/* potential functions */
#include "SurfacePotentialT.h"
#include "XuNeedleman2DT.h"
#include "XuNeedleman3DT.h"
#include "TvergHutch2DT.h"
#include "LinearDamageT.h"
#include "Tijssens2DT.h"
#include "RateDep2DT.h"

/* meshfree domain element types */
#include "MeshFreeFractureSupportT.h"

/* array behavior */

using namespace Tahoe;

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<MeshFreeCSEAnisoT::StatusFlagT>::fByteCopy = true;
} /* namespace Tahoe */

/* parameters */
const int kHeadRoom = 0;

/* constructor */
MeshFreeCSEAnisoT::MeshFreeCSEAnisoT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fMFSurfaceShape(NULL),
	fSurfacePotential(NULL),
	fLocDisp(LocalArrayT::kDisp),
	fFractureArea(0.0),
	fQ(NumSD()),
	fdelta(NumSD()),
	fT(NumSD()),
	fddU_l(NumSD()), fddU_g(NumSD()),
	fdQ(NumSD()),
	fElemEqnosEX(kHeadRoom),
	fActiveFlag(kHeadRoom, true),

	/* dynamic work space managers */
	fLocGroup(kHeadRoom),
	fNEEArray(kHeadRoom, true),
	fNEEMatrix(kHeadRoom, true),
	fMatrixManager(kHeadRoom, true)
{
	const char caller[] = "MeshFreeCSEAnisoT::MeshFreeCSEAnisoT";

	/* set format of element stiffness matrix */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);

	/* read control parameters */
	ifstreamT& in = ElementSupport().Input();
	in >> fGeometryCode;
	in >> fNumIntPts;
	in >> fOutputArea;
	in >> fMFElementGroup;

	/* checks */
	if (NumSD() == 2 && fGeometryCode != GeometryT::kLine)
		ExceptionT::BadInputValue(caller, "expecting geometry code %d for 2D: %d", 
			GeometryT::kLine, fGeometryCode);
	else if (NumSD() == 3 &&
	         fGeometryCode != GeometryT::kQuadrilateral &&
	         fGeometryCode != GeometryT::kTriangle)
		ExceptionT::BadInputValue(caller, "expecting geometry code %d or %d for 3D: %d", 
			GeometryT::kQuadrilateral, GeometryT::kTriangle , fGeometryCode);
	if (fOutputArea != 0 && fOutputArea != 1) throw ExceptionT::kBadInputValue;

	/* check element group */
	fMFElementGroup--;
	ElementBaseT* element_group = &(ElementSupport().ElementGroup(fMFElementGroup));
	
	/* check cast to meshfree group */
	fMFFractureSupport = TB_DYNAMIC_CAST(MeshFreeFractureSupportT*, element_group);
	if (!fMFFractureSupport)
		ExceptionT::BadInputValue(caller, "domain element group %d is not meshfree", fMFElementGroup + 1);

	/* RTTI is required */
#ifdef __NO_RTTI__
	ExceptionT::GeneralFail(caller, "requires RTTI");
#endif
}

/* destructor */
MeshFreeCSEAnisoT::~MeshFreeCSEAnisoT(void)
{
	delete fMFSurfaceShape;
	fMFSurfaceShape = NULL;
	
	delete fSurfacePotential;
	fSurfacePotential = NULL;
}

/* form of tangent matrix */
GlobalT::SystemTypeT MeshFreeCSEAnisoT::TangentType(void) const
{
	/* tangent matrix is not symmetric */
	return GlobalT::kNonSymmetric;
}

/* allocates space and reads connectivity data */
void MeshFreeCSEAnisoT::Initialize(void)
{
	/* override all inherited */
	ElementBaseT::Initialize();

	/* streams */
	ifstreamT& in = ElementSupport().Input();
	ostream&   out = ElementSupport().Output();
		
	/* initialize local arrays */
	fLocDisp.Dimension(0, NumDOF()); // set minor dimension
	Field().RegisterLocal(fLocDisp);
	fLocGroup.Register(fLocDisp);

	/* check */
	if (fMFFractureSupport->NumFacetNodes() < 0)
	{
		cout << "\n MeshFreeCSEAnisoT::Initialize: facets not dimensioned: meshfree\n"
		     <<   "     domain from element group " << fMFElementGroup+1
		     << " must have at least one cutting\n"
		     <<   "     facet or sampling surface" << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* construct surface shape functions */
	MeshFreeSupportT& mf_support = fMFFractureSupport->MeshFreeSupport();
	fMFSurfaceShape = new MeshFreeSurfaceShapeT(fGeometryCode, fNumIntPts,
		mf_support, fLocDisp, fMFFractureSupport->Facets(),
		fMFFractureSupport->NumFacetNodes(), true);	
	if (!fMFSurfaceShape) throw ExceptionT::kOutOfMemory;
	
	/* set up initial cutting facets */
	fMFSurfaceShape->Initialize();

	/* work space */
	fNEEMatrix.Register(fLHS);
	fNEEArray.Register(fRHS);
	fNEEMatrix.Register(fNEEmat);
	fNEEArray.Register(fNEEvec);
	fMatrixManager.Register(fnsd_nee_1);
	fMatrixManager.Register(fnsd_nee_2);
	for (int k = 0; k < NumSD(); k++)
		fMatrixManager.Register(fdQ[k]);

	/* output stream */
	if (fOutputArea == 1)
	{
		/* generate file name */
		StringT name = (ElementSupport().Input()).filename();
		name.Root();
		name.Append(".grp", ElementSupport().ElementGroupNumber(this) + 1);
		name.Append(".fracture");
		
		/* initialize file */
		ofstream out(name);
	}

	/* construct surface potential */
	int code;
	in >> code;
	switch (code)
	{
		case SurfacePotentialT::kXuNeedleman:
		{			
			if (NumDOF() == 2)
				fSurfacePotential = new XuNeedleman2DT(in);
			else
				fSurfacePotential = new XuNeedleman3DT(in);
			break;
		}
		case SurfacePotentialT::kTvergaardHutchinson:
		{
			if (NumDOF() == 2)
				fSurfacePotential = new TvergHutch2DT(in);
			else
			{
				cout << "\n MeshFreeCSEAnisoT::Initialize: Tvergaard-Hutchinson potential not\n"
				     <<   "     implemented for 3D: " << code << endl; 				
				throw ExceptionT::kBadInputValue;
			}
			break;
		}
		case SurfacePotentialT::kLinearDamage:
		{
			fInitTraction.Dimension(NumDOF());
			LinearDamageT* lin_damage = new LinearDamageT(in, fInitTraction);
			if (!lin_damage) throw ExceptionT::kOutOfMemory;
			
			/* cast down */
			fSurfacePotential = lin_damage;
			break;
		}
		case SurfacePotentialT::kTijssens:
		{	
		       if (NumDOF() == 2)
			 fSurfacePotential = new Tijssens2DT(in, ElementSupport().TimeStep());
		       else
		       {
			  cout << "MeshFreeCSEAnisoT::Initialize potential not implemented for 3D: " << code << endl;

	                 throw ExceptionT::kBadInputValue;
		       }
		       break;
		}
	        case SurfacePotentialT::kRateDep:
		{
		       if (NumDOF() == 2)
			 fSurfacePotential = new RateDep2DT(in, ElementSupport().TimeStep());
		       else
		       {
			 cout << "\n MeshFreeCSEAnisoT::Initialize potential not implemented for 3D: " << code << endl;

			 throw ExceptionT::kBadInputValue;
		       }
		       break;
		}
		default:
			cout << "\n MeshFreeCSEAnisoT::Initialize: unknown potential code: " << code << endl;
			throw ExceptionT::kBadInputValue;
	}
	if (!fSurfacePotential) throw ExceptionT::kOutOfMemory;

	/* configure material storage */
	int num_state = fSurfacePotential->NumStateVariables();
	fd_Storage_man.SetWard(0, fd_Storage, num_state);
	fd_Storage_last_man.SetWard(0, fd_Storage_last, num_state);

	/* write */
	out << "\n Cohesive surface potential:\n";
	out << " Potential name:\n";
	fSurfacePotential->PrintName(out);
	fSurfacePotential->Print(out);

	/* allocate flags array */
	const dArray2DT& facets = fMFFractureSupport->Facets();
	fActiveFlag.Dimension(facets.MajorDim());
	fActiveFlag = kON;

	/* initialize decohesion laws */
	InitializeNewFacets();
}

/* start of new time sequence */
void MeshFreeCSEAnisoT::InitialCondition(void)
{
	/* inherited */
	ElementBaseT::InitialCondition();
	
	/* activate all facets */
	fActiveFlag = kON;
}

/* finalize time increment */
void MeshFreeCSEAnisoT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* state variable storage */
	if (fd_Storage_last.MajorDim() != fd_Storage.MajorDim())
		fd_Storage_last_man.SetMajorDimension(fd_Storage.MajorDim(), false);
	fd_Storage_last = fd_Storage;

	/* update state flags */
	for (int i = 0; i < fActiveFlag.Length(); i++)
	{
		/* deactivate facet */
		StatusFlagT& flag = fActiveFlag[i];
		flag = (flag == kMarked) ? kOFF : flag;
	}
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT MeshFreeCSEAnisoT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

	/* mismatch could occur with misuse of managers */
	if (fd_Storage.MajorDim() != fd_Storage_last.MajorDim())
		ExceptionT::GeneralFail("MeshFreeCSEAnisoT::ResetStep", "state variable storage mismatch");
	
	/* restore last state */
	fd_Storage = fd_Storage_last;

	/* unset marks */
	for (int i = 0; i < fActiveFlag.Length(); i++)
	{
		StatusFlagT& flag = fActiveFlag[i];
		flag = (flag == kMarked) ? kON : flag;
	}

	return relax;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT MeshFreeCSEAnisoT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* check for facets to reset */
	const ArrayT<int>& reset_facets = fMFFractureSupport->ResetFacets();
	if (reset_facets.Length() > 0)
	{
		/* signal shape functions */
		//TEMP
		//fMFSurfaceShape->ResetFacets(&reset_facets);
		fMFSurfaceShape->ResetFacets(); // recompute for all facets rather than
		                                // estimating affected field points.
	
		/* resize active flags - fill with kON */
		const dArray2DT& facets = fMFFractureSupport->Facets();
		fActiveFlag.Resize(facets.MajorDim(), kON);

		/* initialize cohesive laws */
		InitializeNewFacets();

		/* override flag */
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQRelax);
	}
	return relax;
}

/* solution calls */
void MeshFreeCSEAnisoT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
//TEMP - not implemented
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* returns the energy as defined by the derived class types */
double MeshFreeCSEAnisoT::InternalEnergy(void) { return 0.0; } //not implemented

/* append element equations numbers to the list */
void MeshFreeCSEAnisoT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* configure equation array */
	ArrayT<int> counts;
	fMFSurfaceShape->NeighborCounts(counts);
	fElemEqnosEX.Configure(counts, NumDOF());

//TEMP - assume all cutting facets are configured.
//       is this the best place to do this???
	fActiveFlag.Dimension(fElemEqnosEX.MajorDim());

	/* get facet neighbors data */
	RaggedArray2DT<int> neighbors;
	neighbors.Configure(counts);
	fMFSurfaceShape->Neighbors(neighbors);

	/* get local equations numbers */
	Field().SetLocalEqnos(neighbors, fElemEqnosEX);

	/* add to list */
	eq_2.Append(&fElemEqnosEX);
}

/* writing output */
void MeshFreeCSEAnisoT::RegisterOutput(void)
{
//TEMP - no output for now
}

/* write integration point data to the output stream */
void MeshFreeCSEAnisoT::WriteOutput(void)
{
	if (fOutputArea)
	{
		/* generate file name */
		StringT name = (ElementSupport().Input()).filename();
		name.Root();
		name.Append(".grp", ElementSupport().ElementGroupNumber(this) + 1);
		name.Append(".fracture");
		
		/* open output file */
		ofstreamT out(name, ios::app);

		/* facet data */
		const dArray2DT& facets = fMFFractureSupport->Facets();

		/* count active facets */
		int count = 0;
		for (int ii = 0; ii < fActiveFlag.Length(); ii++)
			if (fActiveFlag[ii] != kOFF)
				count++;

		/* header */
		out << "\n time = " << setw(kDoubleWidth) << ElementSupport().Time() << '\n';
		out << " fracture area = " << setw(kDoubleWidth) << fFractureArea << '\n';
		out << " number of facets = " << count << '\n';

		/* quick exit */
		if (facets.MajorDim() == 0) return;	
		
		/* header */
		out << " surface tractions (global frame):\n";
		int d_width = out.precision() + kDoubleExtra;
		const char* coord_labels[] = {"x[1]", "x[2]", "x[3]"};
		const char* gap_labels[] = {"du[1]", "du[2]", "du[3]"};
		const char* traction_labels[] = {"t[1]", "t[2]", "t[3]"};
		out << setw(kIntWidth) << "facet";
		out << setw(kIntWidth) << "ip";
		for (int j1 = 0; j1 < NumSD() && j1 < 3; j1++)
			out << setw(d_width) << coord_labels[j1];
		for (int j2 = 0; j2 < NumDOF() && j2 < 3; j2++)
			out << setw(d_width) << gap_labels[j2];
		for (int j3 = 0; j3 < NumSD() && j3 < 3; j3++)
			out << setw(d_width) << traction_labels[j3];
		out << '\n';

		/* loop over cutting facets */
		dArrayT state;
		for (int i = 0; i < facets.MajorDim(); i++)
			if (fActiveFlag[i] != kOFF)
			{  			
				/* set facet dimensions */
				int nnd = fMFSurfaceShape->SetFacet(i);
				SetNumberOfNodes(nnd);
		
				/* get local arrays */
				fLocDisp.SetLocal(fMFSurfaceShape->NodesOnFacet());

				/* loop over integration points */
				fMFSurfaceShape->TopIP();
				while (fMFSurfaceShape->NextIP())
				{
					/* write facet and ip number */
					out << setw(kIntWidth) << i+1
					    << setw(kIntWidth) << fMFSurfaceShape->CurrIP() + 1;
				
					/* fetch state variables */
					fd_Storage.RowAlias(i*fNumIntPts + fMFSurfaceShape->CurrIP(), state);

					/* coordinate transformations */
					double j0, j;
					fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
				
					/* check */
					if (j0 <= 0.0 || j <= 0.0)
					{
						cout << "\n MeshFreeCSEAnisoT::WriteOutput: jacobian error" << endl;
						throw ExceptionT::kBadJacobianDet;
					}
	
					/* gap vector (from side 1 to 2) */
					const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);

					dArrayT tensorIP(3);
		
					/* gap -> traction, in/out of local frame */
					fQ.MultTx(delta, fdelta);
					fQ.Multx(fSurfacePotential->Traction(fdelta, state, tensorIP, false), fT);

					/* coordinates */
					out << fMFSurfaceShape->IPCoords().no_wrap();
					
					/* gap */
					out << fdelta.no_wrap();
					
					/* tractions */
					out << fT.no_wrap() << '\n';
				}									
			}	
	}
}

/* compute specified output parameter and send for smoothing */
void MeshFreeCSEAnisoT::SendOutput(int kincode)
{
//TEMP - no output for now
#pragma unused(kincode)
}

/* appends group connectivities to the array (X -> geometry, U -> field) */
void MeshFreeCSEAnisoT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused(connects)
//Nothing to send
}

void MeshFreeCSEAnisoT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
	
	/* connectivities for surface facets */
	connects_2.Append(&(fMFSurfaceShape->Neighbors(0)));
	connects_2.Append(&(fMFSurfaceShape->Neighbors(1)));
}	

/* returns 1 if DOF's are interpolants of the nodal values */
int MeshFreeCSEAnisoT::InterpolantDOFs(void) const { return 0; }

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void MeshFreeCSEAnisoT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	/* control parameters */
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode << '\n';
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIntPts << '\n';
	out << " Output fracture surface area. . . . . . . . . . = " << fOutputArea << '\n';
	out << " Meshfree domain element group . . . . . . . . . = " << fMFElementGroup + 1 << '\n';
}

/* element data */
void MeshFreeCSEAnisoT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
#pragma unused(in)
#pragma unused(out)
}

void MeshFreeCSEAnisoT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time-integration parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	/* facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();

	/* loop over cutting facets */
	dArrayT state;
	iArrayT eqnos;
	for (int i = 0; i < facets.MajorDim(); i++)
		if (fActiveFlag[i] != kOFF)
		{
		/* set facet dimensions */
		int nnd = fMFSurfaceShape->SetFacet(i);
		SetNumberOfNodes(nnd);
	
		/* get local arrays */
		fLocDisp.SetLocal(fMFSurfaceShape->NodesOnFacet());

		/* initialize */
		fLHS = 0.0;

		/* loop over integration points */
		fMFSurfaceShape->TopIP();
		while (fMFSurfaceShape->NextIP())
		{
			/* load history data */
			fd_Storage.RowAlias(i*fNumIntPts + fMFSurfaceShape->CurrIP(), state);

			/* coordinate transformations */
			double w = fMFSurfaceShape->IPWeight();		
			double j0, j;
			fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
			if (j0 <= 0.0 || j <= 0.0) throw ExceptionT::kBadJacobianDet;
		
			/* gap vector (from side 1 to 2) */
			const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);

			/* gap -> {traction, stiffness} in local frame */
			fQ.MultTx(delta, fdelta);

			dArrayT tensorIP(3);
			const dMatrixT& K = fSurfacePotential->Stiffness(fdelta, state,tensorIP);
			fddU_l.SetToScaled(j0*w*constK, K);
			fddU_g.MultQBQT(fQ, K);
			fddU_g *= j0*w*constK;
		
			const dArrayT& T = fSurfacePotential->Traction(fdelta, state, tensorIP, false);
			fT.SetToScaled(j0*w*constK, T);

			/* shape function table */
			const dMatrixT& d_delta = fMFSurfaceShape->Grad_d();

			/* 1st term */
			Q_ijk__u_j(fdQ, fT, fnsd_nee_1);
			fNEEmat.MultATB(d_delta, fnsd_nee_1);
			fLHS += fNEEmat;
	
			/* 2st term */
			fnsd_nee_1.MultATB(fQ, d_delta);
			fnsd_nee_2.MultATB(fddU_l, fnsd_nee_1);
			u_i__Q_ijk(delta, fdQ, fNEEvec, fnsd_nee_1);
			fNEEmat.MultATB(fnsd_nee_2, fnsd_nee_1);
			fLHS += fNEEmat;

			/* 3rd term */
			fLHS.MultQTBQ(d_delta, fddU_g, dMatrixT::kWhole,
				dMatrixT::kAccumulate);			
		}

		/* assemble */
		fElemEqnosEX.RowAlias(i, eqnos);
		ElementSupport().AssembleLHS(Group(), fLHS, eqnos);
	}
}

void MeshFreeCSEAnisoT::RHSDriver(void)
{
	/* time-integration parameters */
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* set state to start of current step */
	fd_Storage = fd_Storage_last;

	/* fracture surface area */
	fFractureArea = 0.0;

	/* facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();

	/* loop over cutting facets */
	dArrayT state;
	iArrayT eqnos;
	for (int i = 0; i < facets.MajorDim(); i++)
	{
		/* set facet dimensions */
		int nnd = fMFSurfaceShape->SetFacet(i);
		SetNumberOfNodes(nnd);
	
		/* get local arrays */
		fLocDisp.SetLocal(fMFSurfaceShape->NodesOnFacet());

		if (fActiveFlag[i] != kOFF)
		{
	  		/* initialize */
	  		fRHS = 0.0;
			
			/* loop over integration points */
			int all_failed = 1;
			fMFSurfaceShape->TopIP();
			while (fMFSurfaceShape->NextIP())
			{
				/* load history data */
				fd_Storage.RowAlias(i*fNumIntPts + fMFSurfaceShape->CurrIP(), state);

				/* coordinate transformations */
				double w = fMFSurfaceShape->IPWeight();		
				double j0, j;
				fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
				
				/* check */
				if (j0 <= 0.0 || j <= 0.0)
				{
					cout << "\n MeshFreeCSEAnisoT::RHSDriver: jacobian error" << endl;
					throw ExceptionT::kBadJacobianDet;
				}
	
				/* gap vector (from side 1 to 2) */
				const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);

				dArrayT tensorIP(3);	
				/* gap -> traction, in/out of local frame */
				fQ.MultTx(delta, fdelta);
				fQ.Multx(fSurfacePotential->Traction(fdelta, state, tensorIP, true), fT);

				/* expand */
				fMFSurfaceShape->Grad_d().MultTx(fT, fNEEvec);

				/* accumulate */
				fRHS.AddScaled(-j0*w*constKd, fNEEvec);
				
				/* check status */
				SurfacePotentialT::StatusT status = fSurfacePotential->Status(fdelta, state);
				if (status != SurfacePotentialT::Failed) all_failed = 0;
				
				/* fracture area */
				if (fOutputArea && status != SurfacePotentialT::Precritical)
					fFractureArea += j0*w;
			}
									
			/* assemble */
			fElemEqnosEX.RowAlias(i, eqnos);
			ElementSupport().AssembleRHS(Group(), fRHS, eqnos);

			/* mark elements */
			if (all_failed)
			{
				StatusFlagT& flag = fActiveFlag[i];
				if (flag == kON) flag = kMarked;
			}
		}
		else if (fOutputArea)
		{
			/* integrate fracture area */
			fMFSurfaceShape->TopIP();
			while (fMFSurfaceShape->NextIP())
			{
				/* area */
				double j0, j;
				fMFSurfaceShape->Jacobian(j0, j);
			
				/* accumulate */
				fFractureArea += j0*(fMFSurfaceShape->IPWeight());
			}
		}
	}	
}

/* initialize facets in the reset list */
void MeshFreeCSEAnisoT::InitializeNewFacets(void)
{
	/* list of new facets */
	const ArrayT<int>& reset_facets = fMFFractureSupport->ResetFacets();

	/* no initialization needed */
	if (reset_facets.Length() == 0) return;

	/* new facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();
	const dArray2DT& init_tractions = fMFFractureSupport->InitTractions();
	
	/* resize storage space */
	int size = facets.MajorDim()*fNumIntPts;
	fd_Storage_man.SetMajorDimension(size, true);
	
	/* new facets are initialized at construction and during runtime
	 * during RelaxSystem. RelaxSystem is called before CloseStep. Only
	 * database for the new facets is set here. The updated values in
	 * fd_Storage will be copied in during CloseStep. */
	fd_Storage_last_man.SetMajorDimension(size, true);
		
	/* loop over all [facets] x [ip] */
	dArrayT state;
	for (int i = 0; i < reset_facets.Length(); i++)
	{
		/* indices */
		int facet_dex = reset_facets[i];
		int point_dex = facet_dex*fNumIntPts;
	
		/* init traction constant over facet */
		init_tractions.RowAlias(i, fInitTraction);
		for (int j = 0; j < fNumIntPts; j++)
		{
			/* set facet data */
			fd_Storage.RowAlias(point_dex, state);

			/* initialize */
			fSurfacePotential->InitStateVariables(state);

			/* set facet data */
			fd_Storage_last.SetRow(point_dex, state);
		}
		
		/* check activation */
		if (fInitTraction.Magnitude() < kSmall)
			fActiveFlag[facet_dex] = kOFF;
		else
			fActiveFlag[facet_dex] = kON;
	}
}

/* write all current element information to the stream */
void MeshFreeCSEAnisoT::CurrElementInfo(ostream& out) const
{
#pragma unused(out)

	/* override all */
	out << "\n MeshFreeCSEAnisoT::CurrElementInfo: NOT AVAILABLE" << endl;
	
//TEMP - write decohesion T data???
}

/***********************************************************************
* Private
***********************************************************************/

/* set element work space dimensions */
void MeshFreeCSEAnisoT::SetNumberOfNodes(int nnd)
{
	fLocGroup.SetNumberOfNodes(nnd);
	
	int nee = nnd*NumDOF();

	fNEEArray.Dimension(nee, false);
	fNEEMatrix.Dimension(nee, nee);
	fMatrixManager.Dimension(NumSD(), nee);
}

/* operations with pseudo rank 3 (list in j) matrices */
void MeshFreeCSEAnisoT::u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
	dArrayT& nee_vec, dMatrixT& Qu)
{
	for (int i = 0; i < u.Length(); i++)
	{	
		Q[i].MultTx(u, nee_vec);
		Qu.SetRow(i, nee_vec);
	}
}

void MeshFreeCSEAnisoT::Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
	dMatrixT& Qu)
{
	if (Q.Length() == 2)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1]);
	else if (Q.Length() == 3)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1], u[2], Q[2]);
	else
		throw ExceptionT::kGeneralFail;
}
