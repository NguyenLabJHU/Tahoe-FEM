/* $Id: MultiplierContactElement3DT.cpp,v 1.3 2004-07-15 08:28:08 paklein Exp $ */
#include "MultiplierContactElement3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ContactNodeT.h"

/* vector functions */
#include "vector3D.h"

using namespace Tahoe;

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
MultiplierContactElement3DT::MultiplierContactElement3DT
(const ElementSupportT& support, const FieldT& field):
	ContactElementT(support, field, kNumEnfParameters)
{
	fNumMultipliers = 1;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

#if 0
/* print element group data */
void MultiplierContactElement3DT::PrintControlData(ostream& out) const
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
            const dArrayT& search_parameters = fSearchParameters(i,j);
            const dArrayT& enf_parameters = fEnforcementParameters(i,j);
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
#endif

/* called before LHSDriver during iteration process */
void MultiplierContactElement3DT::SetContactStatus(void)
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
#if PRINT_DEBUG
	cout << "\n\n";
	surface.PrintGaps(cout);
  	surface.PrintNormals(cout);
	surface.PrintMultipliers(cout);
  	surface.PrintStatus(cout);
#endif
  }
}

/* called before LHSDriver during iteration process */
void MultiplierContactElement3DT::RHSDriver(void)
{ /* form RESIDUAL */ 
#if 0
  /* update kinematic data */
  UpdateContactConfiguration();

  ContactNodeT* node;
  double gap, pen;
  dArrayT n1;
  n1.Dimension(fNumSD);

  /* residual */
  dArrayT RHS;
  VariArrayT<double> RHS_man(kMaxNumFaceNodes,RHS);
  dArrayT tmp_RHS;
  VariArrayT<double> tmp_RHS_man(kMaxNumFaceNodes,tmp_RHS);
  /* shape functions */
  dMatrixT N1;
  nVariMatrixT<double> N1_man(kMaxNumFaceDOF,N1);
  /* integration weights */
  dArrayT weights;
  VariArrayT<double> weights_man(kMaxNumFaceNodes,weights);
  /* integration points */
  dArray2DT points; // Maybe should be a pointer and const
  /* equation numbers */
  iArray2DT eqnums;
  nVariArray2DT<int> eqnums_man(kMaxNumFaceDOF,eqnums,fNumSD);

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
		

                  dArrayT& parameters =
                        fEnforcementParameters(s,
                        node->OpposingFace()->Surface().Tag());

                  /* linear penalty function */
                  pen = parameters[kPenalty];

	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
	  const FaceT* face = faces[f];
	  face->Quadrature(points,weights);
	  /* primary face */
	  const iArrayT& conn = face->GlobalConnectivity();
	  RHS = 0.0;
	  /*loop over (nodal) quadrature points */
	  /*NOTE: these CORRESPOND to local node numbers */
	  for (int i = 0 ; i < weights.Length() ; i++) {
		node = nodes[face->Node(i)];
		if (node->OpposingFace()) {
		  gap =  node->Gap();
		  face->ComputeShapeFunctions(points(i),N1);
//  		  n1 = node->Normal();
//		  n1.Set(fNumSD,node->Normal());
		  for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
		  N1.Multx(n1, tmp_RHS);
		  /* pressure =  penalty + Lagrange multiplier */
		  tmp_RHS.SetToScaled(pen*gap*weights[i], tmp_RHS);
		  RHS += tmp_RHS;
		}
	  } 
          /* get equation numbers */
          ElementBaseT::fNodes-> SetLocalEqnos(conn, eqnums);
          /* assemble */
	  cout << "RHS\n" << RHS <<'\n';
          fFEManager.AssembleRHS(RHS, eqnums);
	}
  }
#endif
}

void MultiplierContactElement3DT::LHSDriver(GlobalT::SystemTypeT)
{ /* form STIFFNESS */
#if 0
  /* del g =  (n1.N2 * del u2 -  n1.N1 * del u1) */
  /* primary (X) primary block */
  /* primary (X) secondary block */
  int opp_num_nodes;
  ContactNodeT* node;
  double gap, pen;
  dArrayT n1;
  n1.Dimension(fNumSD);

  /* stiffness */    
  ElementMatrixT LHS(ElementMatrixT::kNonSymmetric); //should be using fLHS
  nVariMatrixT <double> LHS_man(kMaxNumFaceNodes,LHS);
  ElementMatrixT tmp_LHS(ElementMatrixT::kNonSymmetric); //should be using fLHS
  nVariMatrixT <double> tmp_LHS_man(kMaxNumFaceNodes,tmp_LHS);
  ElementMatrixT opp_LHS(ElementMatrixT::kNonSymmetric);
  nVariMatrixT <double> opp_LHS_man(kMaxNumFaceNodes,opp_LHS);
  /* shape functions */
  dMatrixT N1;
  nVariMatrixT<double> N1_man(kMaxNumFaceDOF,N1);
  dMatrixT N2;
  nVariMatrixT<double> N2_man(kMaxNumFaceDOF,N2);
  dArrayT N1n;
  VariArrayT<double> N1n_man(kMaxNumFaceDOF,N1n);
  dArrayT N2n;
  VariArrayT<double> N2n_man(kMaxNumFaceDOF,N2n);
  /* integration weights */
  dArrayT weights;
  VariArrayT<double> weights_man(kMaxNumFaceNodes,weights);
  /* integration points */
  dArray2DT points; // NOTE these are copied right now
  /* equation numbers */
  iArray2DT eqnums;
  nVariArray2DT<int> eqnums_man(kMaxNumFaceDOF,eqnums,fNumSD);
  iArray2DT opp_eqnums;
  nVariArray2DT<int> opp_eqnums_man(kMaxNumFaceDOF,opp_eqnums,fNumSD);


  for(int s = 0; s < fSurfaces.Length(); s++) {
        ContactSurfaceT& surface = fSurfaces[s];
        /* all faces on a surface the same size */
        const ArrayT<FaceT*>& faces = surface.Faces();
        const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
        LHS_man.SetDimensions(surface.NumNodesPerFace()*fNumSD);
        tmp_LHS_man.SetDimensions(surface.NumNodesPerFace()*fNumSD);
        N1_man.SetDimensions(surface.NumNodesPerFace()*fNumSD,fNumSD);
        N1n_man.SetLength(surface.NumNodesPerFace()*fNumSD,false);
        weights_man.SetLength(surface.NumNodesPerFace(),false);
        eqnums_man.SetMajorDimension(surface.NumNodesPerFace(),false);

                  dArrayT& parameters =
                        fEnforcementParameters(s,
                        node->OpposingFace()->Surface().Tag());

                  /* linear penalty function */
                  pen = parameters[kPenalty];

        /* form stiffness */
        for (int f = 0;  f < faces.Length(); f++) {
          const FaceT* face = faces[f];
          face->Quadrature(points,weights);
          /* primary face */
	  const iArrayT& conn = face->GlobalConnectivity();
          /* get equation numbers */
          ElementBaseT::fNodes-> SetLocalEqnos(conn, eqnums);
          LHS = 0.0;
          /*loop over (nodal) quadrature points */
	  /*NOTE: these CORRESPOND to local node numbers */
          for (int i = 0 ; i < weights.Length() ; i++) {
		node = nodes[face->Node(i)];
		if (const FaceT* opp_face = node->OpposingFace()) {
		  opp_num_nodes = face->NumNodes();
                  N2_man.SetDimensions(opp_num_nodes*fNumSD,fNumSD);
                  N2n_man.SetLength(opp_num_nodes*fNumSD,false);
        	  opp_LHS_man.SetDimensions(opp_num_nodes*fNumSD);
                  opp_eqnums_man.SetMajorDimension(opp_num_nodes,false);
                  face->ComputeShapeFunctions(points(i),N1);
		  const iArrayT& opp_conn = opp_face->GlobalConnectivity();
                  opp_face->ComputeShapeFunctions
			(node->OpposingLocalCoordinates(),N2);
		  for (int j =0; j < fNumSD; j++) {n1[j] = node->Normal()[j];}
                  N1.Multx(n1, N1n);
                  N2.Multx(n1, N2n); 
		  /* N1n (x) D g */
		  tmp_LHS.Outer(N1n, N1n);
                  tmp_LHS.SetToScaled(pen*weights[i], tmp_LHS);
		  LHS += tmp_LHS;
		  opp_LHS.Outer(N1n, N2n);
                  opp_LHS.SetToScaled(-pen*weights[i], opp_LHS);
          	  /* get equation numbers */
          	  ElementBaseT::fNodes-> SetLocalEqnos(opp_conn, opp_eqnums);
                  /* assemble primary-secondary face stiffness */
	          cout << "p-s LHS \n"<< opp_LHS << '\n';
                  fFEManager.AssembleLHS(opp_LHS, eqnums,opp_eqnums);
		}
          }
          /* assemble primary-primary face stiffness */
	  cout << "p-p LHS \n"<< LHS << '\n';
          fFEManager.AssembleLHS(LHS, eqnums);
        }
  }
#endif
}


