/* $Id: MultiplierContactElement2DT.cpp,v 1.4 2002-04-01 19:04:29 rjones Exp $ */
// created by : rjones 2001
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
	/* warning : req. solver that pivots */
	// CHECK SOLVER TYPE

    /* write out search parameter matrix */
    out << " Interaction parameters ............................\n";
    int num_surfaces = fSearchParameters.Rows();
    for (int i = 0; i < num_surfaces ; i++)
    {
        for (int j = i ; j < num_surfaces ; j++)
        {
            dArrayT& search_parameters = fSearchParameters(i,j);
            dArrayT& enf_parameters = fEnforcementParameters(i,j);
            /* only print allocated parameter arrays */
            if (search_parameters.Length() == kSearchNumParameters) {
              out << "  surface pair: ("  << i << "," << j << ")\n" ;
              out << "  gap tolerance:      "
                    << search_parameters[kGapTol] << '\n';
              out << "  xi tolerance :      "
                    << search_parameters[kXiTol] << '\n';
			  out << "  pass flag    :      "
                    << (int) search_parameters[kPass] << '\n';
              out << "  consistent tangent: "
                    << (int) enf_parameters[kConsistentTangent] << '\n';
              out << "  penalty :           "
                    << enf_parameters[kPenalty] << '\n';
              out << "  gap scale :         "
                    << enf_parameters[kGScale] << '\n';
              out << "  pressure scale:     "
                    << enf_parameters[kPScale] << '\n';
              out << "  pressure tolerance: "
                    << enf_parameters[kTolP] << '\n';
			}
		}
	}

}

/* called before LHSDriver during iteration process */
void MultiplierContactElement2DT::SetStatus(void)
{ 
  int opp_surf_tag;
  for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	for(int n = 0; n < nodes.Length(); n++){	
		ContactNodeT* node = nodes[n];
		node->EnforcementStatus() = kNoP;//initialize
		if (node->HasMultiplier()){
			double& pre = node->Pressure(); 
			if (node->HasProjection() ){
				opp_surf_tag = node->OpposingSurface()->Tag();
				dArrayT& parameters 
					= fEnforcementParameters(surf_tag,opp_surf_tag);
				int ipass = PassType(surf_tag,opp_surf_tag);
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
	}
// DEBUG
	surface.PrintGaps(cout);
	surface.PrintMultipliers(cout);
	surface.PrintStatus(cout);
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
  int num_nodes;
  ContactNodeT* node;
  double gap, pen, pre, opp_pre=0.0;

  for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	RHS_man.SetLength(num_nodes*fNumSD,false);
	xRHS_man.SetLength(num_nodes*fNumMultipliers,false);
	tmp_RHS_man.SetLength(num_nodes*fNumSD,false);
	N1_man.SetDimensions(num_nodes*fNumSD,fNumSD);
	P1_man.SetDimensions(num_nodes,fNumMultipliers);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	xeqnums1_man.SetMajorDimension(num_nodes,false);
	xconn1_man.SetLength(num_nodes,false);
		
	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
		const FaceT* face = faces[f];
		face->Quadrature(points,weights);
		/* primary face */
		const iArrayT& conn1 = face->GlobalConnectivity();
		surface.MultiplierTags(face->Connectivity(),xconn1);

		RHS = 0.0;
		xRHS = 0.0;
		elem_in_contact = 0;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			status = node->EnforcementStatus();
			if (status > kNoP )  {
				elem_in_contact = 1;
				const ContactSurfaceT* opp_surf = node->OpposingSurface();
				opp_surf_tag = opp_surf->Tag();
				dArrayT& parameters = 
				fEnforcementParameters(surf_tag,opp_surf_tag);

/* BLM: U dof on primary surface  -------------------------------------*/
				if (status == kGapZero || status == kPJump){
					/* pressure */
					pre = node->Pressure() ;
					/* gap */
					gap = node->Gap();
					/* pressure =  Lagrange multiplier + penalty */
					if (status == kGapZero) 
						{ pre += -parameters[kPenalty]*gap;}
					face->ComputeShapeFunctions(points(i),N1);
					for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
					N1.Multx(n1, tmp_RHS);
					tmp_RHS.SetToScaled(pre*weights[i], tmp_RHS);
					RHS += tmp_RHS;
				}

/* Constraint : X dof on primary surface ----------------------------- */
				face->ComputeShapeFunctions(points(i),P1);
				if (status == kGapZero) {
					gap = node->Gap();
					P1.SetToScaled(-parameters[kGScale]*gap, P1);
				}
				else if (status == kPJump) {
					/* calculate pressure on opposing face */
					const double* opp_xi  = node->OpposingLocalCoordinates();
					P2values_man.SetLength(opp_surf->NumNodesPerFace()
						*fNumMultipliers,false);
					opp_surf->MultiplierValues(
						node->OpposingFace()->Connectivity(),P2values);	
					opp_pre = node->OpposingFace()->
						Interpolate(opp_xi,P2values);
					double pj = opp_pre - node->Pressure() ;
					P1.SetToScaled(-parameters[kPScale]*pj, P1);
				}
				else if (status == kPZero) {
					double pj = - node->Pressure() ;
					P1.SetToScaled(-parameters[kPScale]*pj, P1);
				}
				xRHS += P1;
			}
		} 
		/* assemble */
		if (elem_in_contact) {
			ElementBaseT::fNodes-> SetLocalEqnos(conn1, eqnums1);
			fFEManager.AssembleRHS(RHS, eqnums1);
			ElementBaseT::fNodes-> XDOF_SetLocalEqnos(xconn1, xeqnums1);
			fFEManager.AssembleRHS(xRHS, xeqnums1);
		}
	}
  }
}

void MultiplierContactElement2DT::LHSDriver(void)
{ /* form STIFFNESS */
  bool elem_in_contact = 0;
  int opp_surf_tag=-1, status=-1;
  int consistent, num_nodes, opp_num_nodes;
  ContactNodeT* node;
  double sfac, gap;
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


  for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	LHS_man.SetDimensions(num_nodes*fNumSD);
	N1_man.SetDimensions(num_nodes*fNumSD,fNumSD);
	N1n_man.SetLength(num_nodes*fNumSD,false);
	//if consistent
	N1nl_man.SetLength(num_nodes*fNumSD,false);
	T1_man.SetDimensions(num_nodes*fNumSD,fNumSD);
	T1n_man.SetLength(num_nodes*fNumSD,false);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	xeqnums1_man.SetMajorDimension(num_nodes,false);

	/* form stiffness */
	for (int f = 0;  f < faces.Length(); f++) {
		/* primary face */
		const FaceT* face = faces[f];
		const iArrayT& conn1 = face->GlobalConnectivity();
		surface.MultiplierTags(face->Connectivity(),xconn1);
		ElementBaseT::fNodes-> XDOF_SetLocalEqnos(xconn1, xeqnums1);
		face->Quadrature(points,weights);
		/* get equation numbers */
		ElementBaseT::fNodes-> SetLocalEqnos(conn1, eqnums1);
		LHS = 0.0;
		elem_in_contact = 0;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			status = node->EnforcementStatus();
			if (status > kNoP )  {
				elem_in_contact = 1;
				/* set-up */
				gap = node->Gap();
				const FaceT* opp_face = node->OpposingFace();
				opp_num_nodes = opp_face->NumNodes();
                const double* opp_xi  = node->OpposingLocalCoordinates();
				const ContactSurfaceT* opp_surf = node->OpposingSurface();
				opp_surf_tag = opp_surf->Tag();
				dArrayT& parameters =
                	fEnforcementParameters(surf_tag,opp_surf_tag);
				consistent = (int) parameters[kConsistentTangent];
				sfac = parameters[kPenalty];

				/* pressure shape function matrix */
				face->ComputeShapeFunctions(points(i),P1);

/* BLM ============================================================= */
/* K =  N1n (x) { P1 + pen (n1.N2 * D u2 -  n1.N1 * D u1) }  */
                if (status == kGapZero || status == kPJump){
					face->ComputeShapeFunctions(points(i),N1);				
					for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
					N1.Multx(n1, N1n);

/* primary U (X) primary   P block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumSD,num_nodes*fNumMultipliers);
					tmp_LHS.Outer(N1n, P1);
					fFEManager.AssembleLHS(tmp_LHS, eqnums1,xeqnums1);

                	if (status == kGapZero){
/* primary U (X) primary   U block */
						if (consistent) { 
							/* slip */
							tmp_LHS_man.SetDimensions
								(num_nodes*fNumSD,num_nodes*fNumSD);
							for (int j =0; j < fNumSD; j++) 
								{l1[j] = node->Tangent1()[j];}
							opp_face->ComputeTangent1(points(i),lm2);
							alpha = Dot(node->Normal(),lm2) 
							      / Dot(l1.Pointer(),lm2);
							for (int j =0; j < fNumSD; j++) 
								{n1alphal1[j] = n1[j] - alpha*l1[j];}
							N1.Multx(n1alphal1, N1nl);
							tmp_LHS.Outer(N1n, N1nl);
							face->ComputeShapeFunctionDerivatives(points(i),T1);
							double jac = face->ComputeJacobian(points(i));
							T1.Multx(n1, T1n);
							T1n.AddCombination(sfac*weights[i], N1nl,
								sfac*gap*alpha*weights[i]/jac,T1n);
							tmp_LHS.Outer(N1n, T1n);
							LHS += tmp_LHS;
							/* jacobian */
							T1.MultAB(T1,Perm);
							tmp_LHS.MultABT(N1, T1);
							tmp_LHS.SetToScaled
								(node->Pressure()-sfac*gap*weights[i], tmp_LHS);
							LHS += tmp_LHS;

						} else {
							LHS.Outer(N1n, N1n);
							LHS.SetToScaled(sfac*weights[i], LHS);
						}
						fFEManager.AssembleLHS(LHS, eqnums1);

/* primary U (X) secondary U block */
						/* get connectivity */
						const iArrayT& conn2 = opp_face->GlobalConnectivity();

						N2_man.SetDimensions(opp_num_nodes*fNumSD,fNumSD);
						N2n_man.SetLength(opp_num_nodes*fNumSD,false);
						eqnums2_man.SetMajorDimension(opp_num_nodes,false);
						opp_face->ComputeShapeFunctions (opp_xi,N2);
						N2.Multx(n1, N2n); 

						tmp_LHS_man.SetDimensions
							(num_nodes*fNumSD,opp_num_nodes*fNumSD);
						tmp_LHS.Outer(N1n, N2n);
						tmp_LHS.SetToScaled(-sfac*weights[i], tmp_LHS);

						/* get equation numbers */
						ElementBaseT::fNodes-> SetLocalEqnos(conn2, eqnums2);
						/* assemble primary-secondary face stiffness */
						fFEManager.AssembleLHS(tmp_LHS, eqnums1,eqnums2);

					}
				} 
/* Constraint ====================================================== */
                if (status == kGapZero) {
					double gfac = parameters[kGScale];
/* primary P (X) primary   U block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,num_nodes*fNumSD);
					if (consistent) {
						tmp_LHS.Outer(P1, T1n);
					} else {
						tmp_LHS.Outer(P1, N1n);
					}
					tmp_LHS.SetToScaled(gfac, tmp_LHS);
					fFEManager.AssembleLHS(tmp_LHS, xeqnums1,eqnums1);
/* primary P (X) secondary U block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,opp_num_nodes*fNumSD);
					tmp_LHS.Outer(P1, N2n);
					tmp_LHS.SetToScaled(-gfac, tmp_LHS);
					fFEManager.AssembleLHS(tmp_LHS, xeqnums1,eqnums2);
                }
                else if (status == kPJump || status == kPZero) {
					double pfac = parameters[kPScale];
/* primary P (X) primary   P block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,num_nodes*fNumMultipliers);
					tmp_LHS.Outer(P1, P1);
					tmp_LHS.SetToScaled(pfac, tmp_LHS);
					fFEManager.AssembleLHS(tmp_LHS, xeqnums1,xeqnums1);
                	if (status == kPJump) {
/* primary P (X) secondary P block */		
					    xconn2_man.SetLength(opp_num_nodes,false);
						opp_surf->MultiplierTags
							(opp_face->Connectivity(),xconn2);
					    xeqnums2_man.SetMajorDimension(opp_num_nodes,false);
						ElementBaseT::fNodes-> 
							XDOF_SetLocalEqnos(xconn2, xeqnums2);
						P2_man.SetDimensions
							(opp_num_nodes*fNumMultipliers,fNumMultipliers);
						opp_face->ComputeShapeFunctions(opp_xi,P2);
						tmp_LHS_man.SetDimensions (num_nodes*fNumMultipliers,
							opp_num_nodes*fNumMultipliers);
						tmp_LHS.Outer(P1, P2);
						tmp_LHS.SetToScaled(-pfac, tmp_LHS);
						fFEManager.AssembleLHS(tmp_LHS, xeqnums1,xeqnums2);
					}
                }
			}
		}
	}
 }
}
