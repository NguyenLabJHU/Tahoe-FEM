/* $Id: PenaltyContactElement2DT.cpp,v 1.6 2002-02-01 20:03:39 dzeigle Exp $ */

#include "PenaltyContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ContactNodeT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "GreenwoodWilliamson.h"

/* vector functions */
#include "vector2D.h"

/* constants */
const double PI = 2.0*acos(0.0);

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
PenaltyContactElement2DT::PenaltyContactElement2DT(FEManagerT& fe_manager):
	ContactElementT(fe_manager, kNumEnfParameters)
{
}

/* print/compute element output quantities */
void PenaltyContactElement2DT::WriteOutput(IOBaseT::OutputModeT mode)
{
	/* call base class */
	ContactElementT::WriteOutput(mode);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* print element group data */
void PenaltyContactElement2DT::PrintControlData(ostream& out) const
{
	ContactElementT::PrintControlData(out);
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
				out << "  gap tolerance:    "
					<< search_parameters[kGapTol] << '\n';
				out << "  xi tolerance :    "
					<< search_parameters[kXiTol] << '\n';
				out << "  penalty :         "
					<< enf_parameters[kPenalty] << '\n';
				out << "  Average asperity height : "
					<< enf_parameters[kAsperityHeightMean] << '\n';
				out << "  Asperity height standard deviation : "
					<< enf_parameters[kAsperityHeightStandardDeviation] << '\n';
				out << "  (work of adhesion) : " 
					<< enf_parameters[kAsperityHeightMean]
				*enf_parameters[kAsperityHeightStandardDeviation]*enf_parameters[kAsperityHeightStandardDeviation] << '\n';
			}
		}
	}
	out <<'\n';
}

/* called before LHSDriver during iteration process */
void PenaltyContactElement2DT::RHSDriver(void)
{ /* form RESIDUAL */ 
  /* update kinematic data */
  UpdateContactConfiguration();

  bool in_contact = 0;
  ContactNodeT* node;
  double gap, pen, pre;
  double gw_m,gw_s;
  double material_coeff=kAsperityDensity*kHertzianModulus*sqrt(kAsperityTipRadius)/(6.0*sqrt(PI));

  /* residual */
  for(int s = 0; s < fSurfaces.Length(); s++) {
	ContactSurfaceT& surface = fSurfaces[s];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	RHS_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	tmp_RHS_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
	N1_man.SetDimensions(surface.NumNodesPerFace()*fNumSD,fNumSD);
	weights_man.SetLength(surface.NumNodesPerFace(),false);
	eqnums_man.SetMajorDimension(surface.NumNodesPerFace(),false);
		
	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
	  const FaceT* face = faces[f];
	  face->Quadrature(points,weights);
	  /* primary face */
	  const iArrayT& conn = face->GlobalConnectivity();
	  RHS = 0.0;
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
		  dArrayT& parameters = 
			fEnforcementParameters(s,
                        node->OpposingFace()->Surface().Tag());

                  /* parameters for Smith-Ferrante Potential */
                  gw_m = parameters[kAsperityHeightMean];
                  gw_s = parameters[kAsperityHeightStandardDeviation];
		  
                  GreenwoodWilliamson GW(gw_m,gw_s);
		  /* First derivative of Greenwood-Williamson represents force */
                  pre  = material_coeff*GW.DFunction(gap);

		  face->ComputeShapeFunctions(points(i),N1);
		  for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
		  N1.Multx(n1, tmp_RHS);
		  /* pressure =  penalty + Lagrange multiplier */
		  tmp_RHS.SetToScaled(pre*weights[i], tmp_RHS);
		  RHS += tmp_RHS;
		}
	  } 
          /* get equation numbers */
          ElementBaseT::fNodes-> SetLocalEqnos(conn, eqnums);
          /* assemble */
          if (in_contact) fFEManager.AssembleRHS(RHS, eqnums);
	}
  }
}

void PenaltyContactElement2DT::LHSDriver(void)
{ /* form STIFFNESS */
  /* del g =  (n1.N2 * del u2 -  n1.N1 * del u1) */
  /* primary (X) primary block */
  /* primary (X) secondary block */

  bool in_contact;
  int consistent;
  int opp_num_nodes;
  ContactNodeT* node;
  double pen, pre, dpre_dg;
  double gap,gw_m,gw_s;
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
   double material_coeff=kAsperityDensity*kHertzianModulus*sqrt(kAsperityTipRadius)/(6.0*sqrt(PI));


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
                  gw_m = parameters[kAsperityHeightMean];
                  gw_s = parameters[kAsperityHeightStandardDeviation];

                  GreenwoodWilliamson GW(gw_m,gw_s);
                  pre  = material_coeff*GW.DFunction(gap);
                  dpre_dg = material_coeff*GW.DDFunction(gap);
		  
                  consistent = (int) parameters[kConsistentTangent];
		  const FaceT* opp_face = node->OpposingFace();
		  opp_num_nodes = opp_face->NumNodes();
                  N2_man.SetDimensions(opp_num_nodes*fNumSD,fNumSD);
                  N2n_man.SetLength(opp_num_nodes*fNumSD,false);
        	  opp_LHS_man.SetDimensions(opp_num_nodes*fNumSD);
                  opp_eqnums_man.SetMajorDimension(opp_num_nodes,false);
		  const iArrayT& opp_conn = opp_face->GlobalConnectivity();
                  opp_face->ComputeShapeFunctions
			(node->OpposingLocalCoordinates(),N2);
		  const double* nm1 = node->Normal();
		  for (int j =0; j < fNumSD; j++) {n1[j] = nm1[j];}
                  N1.Multx(n1, N1n);
		  if (consistent) {
		   for (int j =0; j < fNumSD; j++) {l1[j] =node->Tangent1()[j];}
		   node->OpposingFace()->ComputeTangent1(points(i),lm2);
		   alpha = Dot(node->Normal(),lm2)
			 / Dot(l1.Pointer(),  lm2);
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
        }
  }
}

/***********************************************************************
 * Private
 ***********************************************************************/

