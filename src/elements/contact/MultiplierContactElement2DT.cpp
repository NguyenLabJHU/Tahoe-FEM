/* $Id $ */
#include "MultiplierContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>


#include "MultiplierContactElement2DT.h"
#include "ContactNodeT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"

/* vector functions */
#include "vector2D.h"

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
MultiplierContactElement2DT::MultiplierContactElement2DT
(FEManagerT& fe_manager, XDOF_ManagerT* xdof_nodes):
	ContactElementT(fe_manager, kNumEnfParameters, xdof_nodes)
{
	fNumMultipliers = 1;
}

/* print/compute element output quantities */
void MultiplierContactElement2DT::WriteOutput(IOBaseT::OutputModeT mode)
{
	/* call base class */
	ContactElementT::WriteOutput(mode);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* print element group data */
void MultiplierContactElement2DT::PrintControlData(ostream& out) const
{
	ContactElementT::PrintControlData(out);
}

/* called before LHSDriver during iteration process */
void MultiplierContactElement2DT::SetStatus(void)
{ 
  int opp_surf_tag;
  for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	surface.PrintMultipliers(cout);
	ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	for(int n = 0; n < nodes.Length(); n++){	
		ContactNodeT* node = nodes[n];
		if (node->HasMultiplier()){
			double& pre = node->Pressure(); 
			if (node->HasProjection()){
				opp_surf_tag = node->OpposingFace()->Surface().Tag();
				int ipass = PassType(surf_tag,opp_surf_tag);
				dArrayT& parameters = 
				fEnforcementParameters(surf_tag,opp_surf_tag);
				double tolP = parameters[kTolP];
				if (ipass == kSecondary) { /* pressure collocation */
					node->EnforcementStatus() = kPJump;	
				}
				else if(pre > -tolP) { /* contact */
					node->EnforcementStatus() = kGapZero;	
				}
				else { /* no contact */
					node->EnforcementStatus() = kPZero;	
				}
			}
			else { /* no contact */
					node->EnforcementStatus() = kPZero;	
			}
		}
		else {
			node->EnforcementStatus() = kNoP;	
		}
	}
  }
}

/* called before LHSDriver during iteration process */
void MultiplierContactElement2DT::RHSDriver(void)
{ /* form RESIDUAL */ 
  /* update kinematic data */
  UpdateContactConfiguration();

  /* set status of all surface nodes */
  SetStatus();

  bool elem_in_contact = 0;
  int opp_surf_tag=-1, status=-1;
  ContactNodeT* node;
  double gap, pen, pre, opp_pre=0.0;

  for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	RHS_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	xRHS_man.SetLength(surface.NumNodesPerFace()*fNumMultipliers,false);
	tmp_RHS_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	N1_man.SetDimensions(surface.NumNodesPerFace()*fNumSD,fNumSD);
	P1_man.SetDimensions(surface.NumNodesPerFace(),fNumMultipliers);
	weights_man.SetLength(surface.NumNodesPerFace(),false);
	eqnums_man.SetMajorDimension(surface.NumNodesPerFace(),false);
	xeqnums_man.SetMajorDimension(surface.NumNodesPerFace(),false);
		
	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
		const FaceT* face = faces[f];
		face->Quadrature(points,weights);
		/* primary face */
		const iArrayT& conn = face->GlobalConnectivity();
#if 0
		const iArrayT& xconn = face->XConnectivity();
#endif
        iArrayT xconn;

		RHS = 0.0;
		elem_in_contact = 0;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			status = node->EnforcementStatus();
			if (status > kNoP )  {
				elem_in_contact = 1;
				opp_surf_tag = node->OpposingFace()->Surface().Tag();
				dArrayT& parameters = 
				fEnforcementParameters(surf_tag,opp_surf_tag);
				if (status == kGapZero || status == kPJump){
					/* pressure */
					pre = node->Pressure() ;
					/* gap */
					gap = node->Gap();
					/* pressure =  penalty + Lagrange multiplier */
					if (status == kGapZero) 
						{ pre += parameters[kPenalty]*gap;}

					/* BLM */
					face->ComputeShapeFunctions(points(i),N1);
					//n1.Set(fNumSD,(double*) node->Normal());
					for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
					N1.Multx(n1, tmp_RHS);
					tmp_RHS.SetToScaled(pre*weights[i], tmp_RHS);
					RHS += tmp_RHS;
				}

				/* Constraint */
				face->ComputeShapeFunctions(points(i),P1);
				if (status == kGapZero) {
					P1.SetToScaled(gap*parameters[kGScale], P1);
				}
				else if (status == kPJump) {
					/* calculate pressure on opposing face */
#if 0
					opp_pre = node->OpposingFace()->
						Interpolate(node->OpposingLocalCoordinates(),pvalues);
#endif
					P1.SetToScaled((pre-opp_pre)*parameters[kPScale], P1);
				}
				else if (status == kPZero) {
					P1.SetToScaled(pre*parameters[kPScale], P1);
				}
				xRHS += P1;
			}
		} 
		/* assemble */
		if (elem_in_contact) {
			ElementBaseT::fNodes-> SetLocalEqnos(conn, eqnums);
			fFEManager.AssembleRHS(RHS, eqnums);
			ElementBaseT::fNodes-> SetLocalEqnos(xconn, xeqnums);
			fFEManager.AssembleRHS(xRHS, xeqnums);
		}
	}
  }
}

void MultiplierContactElement2DT::LHSDriver(void)
{ /* form STIFFNESS */
  /* del g =  (n1.N2 * del u2 -  n1.N1 * del u1) */
  /* primary (X) primary block */
  /* primary (X) secondary block */

  bool in_contact;
  int consistent;
  int opp_num_nodes;
  ContactNodeT* node;
  double pen, pre, dpre_dg, gap;
  dArrayT l1;
  l1.Allocate(fNumSD);
  double lm2[3];
  dArrayT n1alphal1;
  n1alphal1.Allocate(fNumSD);

  /* for consistent stiffness */    
  dArrayT N1nl;
  VariArrayT<double> N1nl_man(kMaxNumFaceDOF,N1nl);
  dMatrixT T1;
  nVariMatrixT<double> T1_man(kMaxNumFaceDOF,T1);
  dArrayT T1n;
  VariArrayT<double> T1n_man(kMaxNumFaceDOF,T1n);
  dMatrixT Perm(fNumSD);
  Perm(0,0) = 0.0 ; Perm(0,1) = -1.0;
  Perm(1,0) = 1.0 ; Perm(1,1) =  0.0;
  double alpha;


  for(int s = 0; s < fSurfaces.Length(); s++) {
	ContactSurfaceT& surface = fSurfaces[s];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	LHS_man.SetDimensions(surface.NumNodesPerFace()*fNumSD);
	tmp_LHS_man.SetDimensions(surface.NumNodesPerFace()*fNumSD);
	N1_man.SetDimensions(surface.NumNodesPerFace()*fNumSD,fNumSD);
	N1n_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	//if consistent
	N1nl_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	T1_man.SetDimensions(surface.NumNodesPerFace()*fNumSD,fNumSD);
	T1n_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	weights_man.SetLength(surface.NumNodesPerFace(),false);
	eqnums_man.SetMajorDimension(surface.NumNodesPerFace(),false);

	/* form stiffness */
	for (int f = 0;  f < faces.Length(); f++) {
		const FaceT* face = faces[f];
		face->Quadrature(points,weights);
		/* primary face */
		const iArrayT& conn = face->GlobalConnectivity();
		/* get equation numbers */
		ElementBaseT::fNodes-> SetLocalEqnos(conn, eqnums);
		LHS = 0.0;
		in_contact = 0;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			if (node->OpposingFace()
			&& node->Status() != ContactNodeT::kNoProjection
			&& fEnforcementParameters(s,
			node->OpposingFace()->Surface().Tag()).Length() !=0) {
				in_contact = 1;
				gap =  node->Gap();
				face->ComputeShapeFunctions(points(i),N1);
				dArrayT& parameters =
				fEnforcementParameters(s,
				node->OpposingFace()->Surface().Tag());

				pen = parameters[kPenalty];
				pre = pen*gap; 
				dpre_dg = pen;

				consistent = (int) parameters[kConsistentTangent];
				const FaceT* opp_face = node->OpposingFace();
				opp_num_nodes = opp_face->NumNodes();
				N2_man.SetDimensions(opp_num_nodes*fNumSD,fNumSD);
				N2n_man.SetLength(opp_num_nodes*fNumSD,false);
				opp_LHS_man.SetDimensions(opp_num_nodes*fNumSD);
				opp_eqnums_man.SetMajorDimension(opp_num_nodes,false);

				const iArrayT& opp_conn = 
					opp_face->GlobalConnectivity();
				opp_face->ComputeShapeFunctions
					(node->OpposingLocalCoordinates(),N2);
				// n1.Set(fNumSD,node->Normal());
				for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
				N1.Multx(n1, N1n);
				if (consistent) {
					//l1.Set(fNumSD,node->Tangent1());
					for (int j =0; j < fNumSD; j++) {l1[j] = node->Tangent1()[j];}
					node->OpposingFace()->ComputeTangent1(points(i),lm2);
					alpha = Dot(node->Normal(),lm2) / Dot(l1.Pointer(),  lm2);
					for (int j =0; j < fNumSD; j++) 
						{n1alphal1[j] = n1[j] - alpha*l1[j];}
				}
				N2.Multx(n1, N2n); 

				/* N1n (x) D g */
				/* Part:  dx1 (x) dx2 */
				opp_LHS.Outer(N1n, N2n);
				opp_LHS.SetToScaled(-dpre_dg*weights[i], opp_LHS);

				/* get equation numbers */
				ElementBaseT::fNodes-> SetLocalEqnos(opp_conn, opp_eqnums);
				/* assemble primary-secondary face stiffness */
				fFEManager.AssembleLHS(opp_LHS, eqnums,opp_eqnums);

				/* Part:  dx1 (x) dx1 */
				if (consistent) {
					N1.Multx(n1alphal1, N1nl);
					tmp_LHS.Outer(N1n, N1nl);
					tmp_LHS.SetToScaled(dpre_dg*weights[i], tmp_LHS);
					LHS += tmp_LHS;

					face->ComputeShapeFunctionDerivatives(points(i),T1);
					double jac = face->ComputeJacobian(points(i));

					T1.Multx(n1, T1n);
					tmp_LHS.Outer(N1n, T1n);
					tmp_LHS.SetToScaled(pre*alpha*weights[i]/jac, tmp_LHS);
					LHS += tmp_LHS;

					T1.MultAB(T1,Perm);
					tmp_LHS.MultABT(N1, T1);
					tmp_LHS.SetToScaled(-pre*weights[i], tmp_LHS);
					LHS += tmp_LHS;

				} else {
					tmp_LHS.Outer(N1n, N1n);
					tmp_LHS.SetToScaled(dpre_dg*weights[i], tmp_LHS);
					LHS += tmp_LHS;
				}
			}
		}
		/* assemble primary-primary face stiffness */
		if (in_contact) fFEManager.AssembleLHS(LHS, eqnums);


// MORE ....
	}
  }
}
