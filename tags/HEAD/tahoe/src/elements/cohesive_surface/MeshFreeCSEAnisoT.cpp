/* $Id: MeshFreeCSEAnisoT.cpp,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (06/08/2000)                                          */

#include "MeshFreeCSEAnisoT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "SurfacePotentialT.h"
#include "eControllerT.h"
#include "MeshFreeSurfaceShapeT.h"

/* potential functions */
#include "XuNeedleman2DT.h"
#include "XuNeedleman3DT.h"
#include "TvergHutch2DT.h"
#include "LinearDamageT.h"

/* meshfree domain element types */
#include "MeshFreeFractureSupportT.h"

/* array behavior */
const bool ArrayT<MeshFreeCSEAnisoT::StatusFlagT>::fByteCopy = true;

/* parameters */
const int kHeadRoom = 0;

/* constructor */
MeshFreeCSEAnisoT::MeshFreeCSEAnisoT(FEManagerT& fe_manager):
	ElementBaseT(fe_manager),
	fMFSurfaceShape(NULL),
	fSurfacePotential(NULL),
	fIsDecohesion(false),
	fLocDisp(LocalArrayT::kDisp),
	fFractureArea(0.0),
	fQ(fNumSD),
	fdelta(fNumSD),
	fT(fNumSD),
	fddU_l(fNumSD), fddU_g(fNumSD),
	fdQ(fNumSD),
	fElemEqnosEX(kHeadRoom),
	fActiveFlag(kHeadRoom, true),

	/* dynamic work space managers */
	fLocGroup(kHeadRoom),
	fNEEArray(kHeadRoom),
	fNEEMatrix(kHeadRoom),
	fMatrixManager(kHeadRoom)
{
	/* set format of element stiffness matrix */
	fLHS.SetFormat(ElementMatrixT::kNonSymmetric);

	/* read control parameters */
	ifstreamT& in = fFEManager.Input();
	in >> fGeometryCode;
	in >> fNumIntPts;
	in >> fOutputArea;
	in >> fMFElementGroup;

	/* checks */
	if (fNumSD == 2 && fGeometryCode != GeometryT::kLine)
	{
		cout << "\n MeshFreeCSEAnisoT::MeshFreeCSEAnisoT: expecting geometry code "
		     << GeometryT::kLine<< " for 2D: " << fGeometryCode << endl;
		throw eBadInputValue;
	}
	else if (fNumSD == 3 &&
	         fGeometryCode != GeometryT::kQuadrilateral &&
	         fGeometryCode != GeometryT::kTriangle)
	{
		cout << "\n MeshFreeCSEAnisoT::MeshFreeCSEAnisoT: expecting geometry code "
<< GeometryT::kQuadrilateral
		     << " or\n" <<   "     " << GeometryT::kTriangle << " for 3D: "
		     << fGeometryCode << endl;
		throw eBadInputValue;
	}
	if (fOutputArea != 0 && fOutputArea != 1) throw eBadInputValue;

	/* check element group */
	fMFElementGroup--;
	ElementBaseT* element_group = fFEManager.ElementGroup(fMFElementGroup);
	if (!element_group)
	{
		cout << "\n MeshFreeCSEAnisoT::MeshFreeCSEAnisoT: domain element group\n"
		     <<   "     " << fMFElementGroup + 1 << " not found" << endl;
		throw eBadInputValue;
	}
	
	/* check cast to meshfree group */
#ifdef __NO_RTTI__
	cout << "\n MeshFreeCSEAnisoT::MeshFreeCSEAnisoT: NO RTTI: Domain element\n"
	     <<   "     group " << fMFElementGroup + 1
	     << " cannot be verified as meshfree" << endl;
	fMFFractureSupport = (MeshFreeFractureSupportT*) element_group;
#else
	fMFFractureSupport = dynamic_cast<MeshFreeFractureSupportT*>(element_group);
	if (!fMFFractureSupport)
	{
		cout << "\n MeshFreeCSEAnisoT::MeshFreeCSEAnisoT: domain element group\n"
		     <<   "    " << fMFElementGroup + 1 << " is not meshfree" << endl;
		throw eBadInputValue;
	}
#endif
}

/* destructor */
MeshFreeCSEAnisoT::~MeshFreeCSEAnisoT(void)
{
	delete fMFSurfaceShape;
	fMFSurfaceShape = NULL;
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
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();
		
	/* initialize local arrays */
	fLocDisp.Allocate(0, fNumDOF); // set minor dimension
	fFEManager.RegisterLocal(fLocDisp);
	fLocGroup.Register(fLocDisp);

	/* check */
	if (fMFFractureSupport->NumFacetNodes() < 0)
	{
		cout << "\n MeshFreeCSEAnisoT::Initialize: facets not dimensioned: meshfree\n"
		     <<   "     domain from element group " << fMFElementGroup+1
		     << " must have at least one cutting\n"
		     <<   "     facet or sampling surface" << endl;
		throw eBadInputValue;
	}

	/* construct surface shape functions */
	MeshFreeSupportT& mf_support = fMFFractureSupport->MeshFreeSupport();
	fMFSurfaceShape = new MeshFreeSurfaceShapeT(fGeometryCode, fNumIntPts,
		mf_support, fLocDisp, fMFFractureSupport->Facets(),
		fMFFractureSupport->NumFacetNodes(), true);	
	if (!fMFSurfaceShape) throw eOutOfMemory;
	
	/* set up initial cutting facets */
	fMFSurfaceShape->Initialize();

	/* work space */
	fNEEMatrix.Register(fLHS);
	fNEEArray.Register(fRHS);
	fNEEMatrix.Register(fNEEmat);
	fNEEArray.Register(fNEEvec);
	fMatrixManager.Register(fnsd_nee_1);
	fMatrixManager.Register(fnsd_nee_2);
	for (int k = 0; k < fNumSD; k++)
		fMatrixManager.Register(fdQ[k]);

	/* output stream */
	if (fOutputArea == 1)
	{
		/* generate file name */
		StringT name = (fFEManager.Input()).filename();
		name.Root();
		name.Append(".grp", fFEManager.ElementGroupNumber(this) + 1);
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
			if (fNumDOF == 2)
				fSurfacePotential = new XuNeedleman2DT(in);
			else
				fSurfacePotential = new XuNeedleman3DT(in);
			break;
		}
		case SurfacePotentialT::kTvergaardHutchinson:
		{
			if (fNumDOF == 2)
				fSurfacePotential = new TvergHutch2DT(in);
			else
			{
				cout << "\n MeshFreeCSEAnisoT::Initialize: Tvergaard-Hutchinson potential not\n"
				     <<   "     implemented for 3D: " << code << endl; 				
				throw eBadInputValue;
			}
			break;
		}
		case SurfacePotentialT::kLinearDamage:
		{
			/* set flag */
			fIsDecohesion = true; // use flag instead of RTTI
		
			fInitTraction.Allocate(fNumDOF);
			LinearDamageT* lin_damage = new LinearDamageT(in, fInitTraction, fi_vec, fd_vec);
			if (!lin_damage) throw eOutOfMemory;
		
			/* configure material storage */
			fi_StorageMan.SetWard(0, fi_Storage, lin_damage->IntegerStorage());
			fd_StorageMan.SetWard(0, fd_Storage, lin_damage->DoubleStorage());
			
			/* cast down */
			fSurfacePotential = lin_damage;
			break;
		}
		default:
			cout << "\n MeshFreeCSEAnisoT::Initialize: unknown potential code: " << code << endl;
			throw eBadInputValue;
	}
	if (!fSurfacePotential) throw eOutOfMemory;

	/* write */
	out << "\n Cohesive surface potential:\n";
	out << " Potential name:\n";
	fSurfacePotential->PrintName(out);
	fSurfacePotential->Print(out);

	/* allocate flags array */
	const dArray2DT& facets = fMFFractureSupport->Facets();
	fActiveFlag.Allocate(facets.MajorDim());
	fActiveFlag = kON;

	/* initialize decohesion laws */
	if (fIsDecohesion) InitializeNewFacets();
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

	if (fIsDecohesion)
	{
		/* cohesive law */
		DecohesionLawT* cohesive = (DecohesionLawT*) fSurfacePotential;
	
		for (int i = 0; i < fActiveFlag.Length(); i++)
		{
			StatusFlagT& flag = fActiveFlag[i];
			if (flag != kOFF)
			{
				for (int j = 0; j < fNumIntPts; j++)
				{
					/* load ip data */
					int row = i*fNumIntPts + j;
					fi_Storage.RowAlias(row, fi_vec);
					fd_Storage.RowAlias(row, fd_vec);
					
					/* set internal variables */
					cohesive->UpdateHistory();
				}
			}
			
			/* deactivate facet */
			flag = (flag == kMarked) ? kOFF : flag;
		}
	}
	else
	{
		for (int i = 0; i < fActiveFlag.Length(); i++)
		{
			/* deactivate facet */
			StatusFlagT& flag = fActiveFlag[i];
			flag = (flag == kMarked) ? kOFF : flag;
		}
	}
}

/* resets to the last converged solution */
void MeshFreeCSEAnisoT::ResetStep(void)
{
	/* inherited */
	ElementBaseT::ResetStep();

	if (fIsDecohesion)
	{
		/* cohesive law */
		DecohesionLawT* cohesive = (DecohesionLawT*) fSurfacePotential;
	
		for (int i = 0; i < fActiveFlag.Length(); i++)
		{
			StatusFlagT& flag = fActiveFlag[i];
			if (flag != kOFF)
			{
				for (int j = 0; j < fNumIntPts; j++)
				{
					/* load ip data */
					int row = i*fNumIntPts + j;
					fi_Storage.RowAlias(row, fi_vec);
					fd_Storage.RowAlias(row, fd_vec);
					
					/* set internal variables */
					cohesive->ResetHistory();
				}
			}
			
			/* reset facet */
			flag = (flag == kMarked) ? kON : flag;
		}
	}
	else
	{
		/* unset marks */
		for (int i = 0; i < fActiveFlag.Length(); i++)
		{
			StatusFlagT& flag = fActiveFlag[i];
			flag = (flag == kMarked) ? kON : flag;
		}
	}
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

		/* initialize decohesion laws */
		if (fIsDecohesion) InitializeNewFacets();

		/* override flag */
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQRelax);
	}
	return relax;
}

/* solution calls */
void MeshFreeCSEAnisoT::AddNodalForce(int node, dArrayT& force)
{
//TEMP - not implemented
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
	fElemEqnosEX.Configure(counts, fNumDOF);

//TEMP - assume all cutting facets are configured.
//       is this the best place to do this???
	fActiveFlag.Allocate(fElemEqnosEX.MajorDim());

	/* get facet neighbors data */
	RaggedArray2DT<int> neighbors;
	neighbors.Configure(counts);
	fMFSurfaceShape->Neighbors(neighbors);

	/* get local equations numbers */
	fNodes->SetLocalEqnos(neighbors, fElemEqnosEX);

	/* add to list */
	eq_2.Append(&fElemEqnosEX);
}

/* writing output */
void MeshFreeCSEAnisoT::RegisterOutput(void)
{
//TEMP - no output for now
}

/* write integration point data to the output stream */
void MeshFreeCSEAnisoT::WriteOutput(IOBaseT::OutputModeT mode)
{
	if (fOutputArea && mode == IOBaseT::kAtInc)
	{
		/* generate file name */
		StringT name = (fFEManager.Input()).filename();
		name.Root();
		name.Append(".grp", fFEManager.ElementGroupNumber(this) + 1);
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
		out << "\n time = " << setw(kDoubleWidth) << fFEManager.Time() << '\n';
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
		for (int j1 = 0; j1 < fNumSD && j1 < 3; j1++)
			out << setw(d_width) << coord_labels[j1];
		for (int j2 = 0; j2 < fNumDOF && j2 < 3; j2++)
			out << setw(d_width) << gap_labels[j2];
		for (int j3 = 0; j3 < fNumSD && j3 < 3; j3++)
			out << setw(d_width) << traction_labels[j3];
		out << '\n';

		/* loop over cutting facets */
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
				
					/* load history data */
					if (fIsDecohesion)
					{
						int row = i*fNumIntPts + fMFSurfaceShape->CurrIP();
						fi_Storage.RowAlias(row, fi_vec);
						fd_Storage.RowAlias(row, fd_vec);
					}

					/* coordinate transformations */
					double j0, j;
					fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
				
					/* check */
					if (j0 <= 0.0 || j <= 0.0)
					{
						cout << "\n MeshFreeCSEAnisoT::WriteOutput: jacobian error" << endl;
						throw eBadJacobianDet;
					}
	
					/* gap vector (from side 1 to 2) */
					const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);
		
					/* gap -> traction, in/out of local frame */
					fQ.MultTx(delta, fdelta);
					fQ.Multx(fSurfacePotential->Traction(fdelta), fT);

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

void MeshFreeCSEAnisoT::LHSDriver(void)
{
	/* time-integration parameters */
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

	/* facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();

	/* loop over cutting facets */
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
			if (fIsDecohesion)
			{
				int row = i*fNumIntPts + fMFSurfaceShape->CurrIP();
				fi_Storage.RowAlias(row, fi_vec);
				fd_Storage.RowAlias(row, fd_vec);
			}

			/* coordinate transformations */
			double w = fMFSurfaceShape->IPWeight();		
			double j0, j;
			fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
			if (j0 <= 0.0 || j <= 0.0) throw eBadJacobianDet;
		
			/* gap vector (from side 1 to 2) */
			const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);

			/* gap -> {traction, stiffness} in local frame */
			fQ.MultTx(delta, fdelta);

			const dMatrixT& K = fSurfacePotential->Stiffness(fdelta);
			fddU_l.SetToScaled(j0*w*constK, K);
			fddU_g.MultQBQT(fQ, K);
			fddU_g *= j0*w*constK;
			
			const dArrayT& T = fSurfacePotential->Traction(fdelta);
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
		fFEManager.AssembleLHS(fLHS, eqnos);
	}
}

void MeshFreeCSEAnisoT::RHSDriver(void)
{
	/* time-integration parameters */
	double constKd = 0.0;
	int formKd = fController->FormKd(constKd);
	if (!formKd) return;

	/* fracture surface area */
	fFractureArea = 0.0;

	/* facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();

	/* loop over cutting facets */
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
				if (fIsDecohesion)
				{
					int row = i*fNumIntPts + fMFSurfaceShape->CurrIP();
					fi_Storage.RowAlias(row, fi_vec);
					fd_Storage.RowAlias(row, fd_vec);
				}

				/* coordinate transformations */
				double w = fMFSurfaceShape->IPWeight();		
				double j0, j;
				fMFSurfaceShape->Jacobian(j0, j, fQ, fdQ);
				
				/* check */
				if (j0 <= 0.0 || j <= 0.0)
				{
					cout << "\n MeshFreeCSEAnisoT::RHSDriver: jacobian error" << endl;
					throw eBadJacobianDet;
				}
	
				/* gap vector (from side 1 to 2) */
				const dArrayT& delta = fMFSurfaceShape->InterpolateJumpU(fLocDisp);
	
				/* gap -> traction, in/out of local frame */
				fQ.MultTx(delta, fdelta);
				fQ.Multx(fSurfacePotential->Traction(fdelta), fT);

				/* expand */
				fMFSurfaceShape->Grad_d().MultTx(fT, fNEEvec);

				/* accumulate */
				fRHS.AddScaled(-j0*w*constKd, fNEEvec);
				
				/* check status */
				SurfacePotentialT::StatusT status = fSurfacePotential->Status(fdelta);
				if (status != SurfacePotentialT::Failed) all_failed = 0;
				
				/* fracture area */
				if (fOutputArea && status != SurfacePotentialT::Precritical)
					fFractureArea += j0*w;
			}
									
			/* assemble */
			fElemEqnosEX.RowAlias(i, eqnos);
			fFEManager.AssembleRHS(fRHS, eqnos);

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
	if (!fIsDecohesion || reset_facets.Length() == 0) return;

	/* new facet data */
	const dArray2DT& facets = fMFFractureSupport->Facets();
	const dArray2DT& init_tractions = fMFFractureSupport->InitTractions();
	
	/* resize storage space */
	int size = facets.MajorDim()*fNumIntPts;
	fi_StorageMan.SetMajorDimension(size, true);
	fd_StorageMan.SetMajorDimension(size, true);
		
	/* loop over all [facets] x [ip] */
	DecohesionLawT* cohesive = (DecohesionLawT*) fSurfacePotential;
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
			fi_Storage.RowAlias(point_dex, fi_vec);
			fd_Storage.RowAlias(point_dex, fd_vec);
			point_dex++;

			/* initialize */
			cohesive->InitializeFacet();
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
	
	int nee = nnd*fNumDOF;

	fNEEArray.Dimension(nee, false);
	fNEEMatrix.Dimension(nee, nee);
	fMatrixManager.Dimension(fNumSD, nee);
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
		throw eGeneralFail;
}
