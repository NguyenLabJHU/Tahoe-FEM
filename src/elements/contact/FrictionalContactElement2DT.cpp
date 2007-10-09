/* $Id: FrictionalContactElement2DT.cpp,v 1.5 2007-10-09 23:24:46 rjones Exp $ */
// created by : rjones 2003
#include "FrictionalContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "FrictionalContactElement2DT.h"
#include "ContactNodeT.h"
#include "ElementSupportT.h"
#include "XDOF_ManagerT.h"
#include "ofstreamT.h"

/* vector functions */
#include "vector2D.h"

#undef  PRINT_DEBUG
#define PRINT_DEBUG 1
#undef  HACK
#define HACK 0

using namespace Tahoe;

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
FrictionalContactElement2DT::FrictionalContactElement2DT(const ElementSupportT& support):
	ContactElementT(support)
{
	SetName("Jones_frictional_contact_2D");
	fNumMultipliers = 2;
}

/* accept parameter list */
void FrictionalContactElement2DT::TakeParameterList(const ParameterListT& list)
{
	/* set pointer to the XDOF manager */
	fXDOF_Nodes = &(ElementSupport().XDOF_Manager());

	/* inherited */
	ContactElementT::TakeParameterList(list);

    /* write out search parameter matrix */
    ofstreamT& out = ElementSupport().Output();
    out << " Interaction parameters ............................\n";
    int num_surfaces = fSearchParameters.Rows();
    for (int i = 0; i < num_surfaces ; i++)
    {
        for (int j = i ; j < num_surfaces ; j++)
        {
            const dArrayT& search_parameters = fSearchParameters(i,j);
            const dArrayT& enf_parameters = fEnforcementParameters(i,j);
			if (enf_parameters.Length() != kNumEnfParameters)
				ExceptionT::GeneralFail("FrictionalContactElement2DT::TakeParameterList",
					"expecting %d enforcement parameters not %d",
					kNumEnfParameters, enf_parameters.Length());

            /* only print allocated parameter arrays */
            if (search_parameters.Length() == kNumSearchParameters) {
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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called before LHSDriver during iteration process */
void FrictionalContactElement2DT::SetContactStatus(void)
{ 
  int iteration = ElementSupport().IterationNumber();
  cout << "ITERATION : " << iteration << "\n";
  /* update kinematic data */
  UpdateContactConfiguration();

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
#if PRINT_DEBUG
		cout << "surface: " << surf_tag << ", node: " << node->Tag();
		if      (node->EnforcementStatus() == kNoP)      cout << " NoP\n";
		else if (node->EnforcementStatus() == kPJump)    cout << " [[P]]\n";
		else if (node->EnforcementStatus() == kGapZero)  cout << " [[x]]\n";
		else if (node->EnforcementStatus() == kPZero)    cout << " P=0\n";
		else                                             cout << " unk\n";
#endif
	}
#if PRINT_DEBUG && 0
	cout << "\n";
	surface.PrintGaps(cout);
	cout << "\n";
  	surface.PrintNormals(cout);
	cout << "\n";
	surface.PrintMultipliers(cout);
	cout << "\n";
  	surface.PrintStatus(cout);
#endif
  }
}

/* called before LHSDriver during iteration process */
void FrictionalContactElement2DT::RHSDriver(void)
{ /* form RESIDUAL */ 

/* set status of all surface nodes */
SetContactStatus();

bool elem_in_contact = 0;
int opp_surf_tag=-1, status=-1;
int num_nodes;
ContactNodeT* node=0;
double pre=0.0, opp_pre=0.0, tau=0.0, opp_tau=0.0;
double* xi=0;
int len_eqnum1=0, len_xeqnum1=0, len_eqnum2=0, len_xeqnum2=0;

dArrayT traction(NumSD());

iArrayT node_tag(1);
iArrayT xc(1);
iArray2DT xe(1,fNumMultipliers);
dArrayT xR(fNumMultipliers);
dArrayT xv(fNumMultipliers);
xv = 0.0;

int nsd = NumSD();
for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	len_eqnum1  = num_nodes*nsd;
	len_xeqnum1 = num_nodes*fNumMultipliers;
	RHS_man.SetLength(len_eqnum1,false);
	tmp_RHS_man.SetLength(len_eqnum1,false);
	xRHS_man.SetLength(len_xeqnum1,false);
	tmp_xRHS_man.SetLength(len_xeqnum1,false);
	N1_man.SetDimensions(len_eqnum1, nsd);
	P1_man.SetDimensions(len_xeqnum1,fNumMultipliers);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	xeqnums1_man.SetMajorDimension(num_nodes,false);
	xconn1_man.SetLength(num_nodes,false);
		
	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
		const FaceT* face = faces[f];
		bool has_projections = 1, has_multipliers = 1;
		for (int i = 0 ; i < num_nodes ; i++) {
			node = nodes[face->Node(i)];
			if (!(node->HasProjection())) has_projections = 0;
			if (!(node->HasMultiplier())) has_multipliers = 0;
		}
		if ( ! has_projections || ! has_multipliers ) {
			double pScale = fEnforcementParameters(0,1)[kPScale];
			for (int i = 0 ; i < num_nodes ; i++) {
				node = nodes[face->Node(i)];
				if (node->EnforcementStatus() == kPZero) {
					node_tag[0] = node->Tag();
					surface.MultiplierTags(node_tag,xc);
					xR[0] = pScale* ( - node->Pressure() ) ;
					xR[1] = pScale* ( - node->TangTraction() ) ;
					const ElementSupportT& support = ElementSupport();
					support.XDOF_Manager().XDOF_SetLocalEqnos(Group(), xc, xe);
					support.AssembleRHS(Group(), xR, xe);
				}
			} 
		}
		else {
		face->Quadrature(points,weights);
		/* primary face */
		const iArrayT& conn1 = face->GlobalConnectivity();
		surface.MultiplierTags(face->Connectivity(),xconn1);

		RHS = 0.0;
		xRHS = 0.0;
		elem_in_contact = 1;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			xi = points(i);
			status = node->EnforcementStatus();
			if (status > kPZero )  {
				const ContactSurfaceT* opp_surf = node->OpposingSurface();
				opp_surf_tag = opp_surf->Tag();
				dArrayT& parameters = 
					fEnforcementParameters(surf_tag,opp_surf_tag);
/* BLM: U dof on primary surface  -------------------------------------*/
				if (status == kGapZero || status == kPJump){
					/* traction */
					pre = node->Pressure() ;
					tau = node->TangTraction() ;
					/* traction =  Lagrange multiplier + penalty */
					if (status == kGapZero) { 
						double pen = parameters[kPenalty];
						pre += -pen * node->Gap();
						tau += -pen * node->ComputeSlip();
					}
					for (int j =0; j < nsd; j++) 
						{traction[j] = (-pre*node->Normal()[j] 
								        -tau*node->Tangent()[j] )*weights[i];}
#if PRINT_DEBUG
cout << node->Tag() << ", pre: " << pre << ", tau: " << tau 
		<< ", traction = (" << traction[0]/weights[i] << "," << traction[1]/weights[i] << ")\n";
#endif
					face->ComputeShapeFunctions(xi,N1);
					N1.Multx(traction, tmp_RHS);
					RHS += tmp_RHS;
				}
/* Constraint : X dof on primary surface ----------------------------- */
				if (status == kGapZero) {
					double gfac = parameters[kGScale];
					xv[0] = gfac*node->Gap();
					xv[1] = gfac*node->ComputeSlip();
#if PRINT_DEBUG
cout << node->Tag() << ", gap: " << xv[0]/gfac << ", slip: " << xv[1]/gfac << "\n"; 
#endif
				}
				else if (status == kPJump) {
					/* calculate pressure on opposing face */
					const double* opp_xi  = node->OpposingLocalCoordinates();
					P2values_man.Dimension
						(opp_surf->NumNodesPerFace(),fNumMultipliers);
					opp_surf->MultiplierValues(
						node->OpposingFace()->Connectivity(),P2values);	
					node->OpposingFace()->
						InterpolateVector(opp_xi,P2values,traction.Pointer());
					double pfac = parameters[kPScale];
					xv[0]=pfac*(traction[0]-node->Pressure());
					xv[1]=pfac*(traction[1]-node->TangTraction());
#if PRINT_DEBUG
cout << node->Tag() << ", [[p]]: " << xv[0]/pfac << ", [[t]]: " << xv[1]/pfac << "\n"; 
#endif
				}
				face->ComputeShapeFunctions(xi,P1);
				P1.Multx(xv, tmp_xRHS);
				xRHS += tmp_xRHS;
			}
			else if (status == kPZero) {
				double pScale = fEnforcementParameters(0,1)[kPScale];
				xv[0] = pScale*( - node->Pressure() );
				xv[1] = pScale*( - node->TangTraction() );
				face->ComputeShapeFunctions(xi,P1);
				P1.Multx(xv, tmp_xRHS);
				xRHS += tmp_xRHS;
			} 
			else {
				elem_in_contact = 0;
			}
		} 
		/* assemble */
		if (elem_in_contact) {
			const ElementSupportT& support = ElementSupport();
		
			Field().SetLocalEqnos(conn1, eqnums1);
			support.AssembleRHS(Group(), RHS, eqnums1);
			support.XDOF_Manager().XDOF_SetLocalEqnos(Group(),xconn1,xeqnums1);
			support.AssembleRHS(Group(), xRHS, xeqnums1);
		}
		}
	}
}
}

void FrictionalContactElement2DT::LHSDriver(GlobalT::SystemTypeT)
{ /* form STIFFNESS */
bool elem_in_contact = 0;
int opp_surf_tag=-1, status=-1;
int consistent, num_nodes, opp_num_nodes;
ContactNodeT* node;
double sfac, gap;
double l1[2],lm2[2];
dArrayT n1alphal1;
n1alphal1.Dimension(NumSD());

int len_eqnum1=0, len_xeqnum1=0, len_eqnum2=0, len_xeqnum2=0;
double* xi=0;

/* for consistent stiffness */    
dArrayT N1nl;
VariArrayT<double> N1nl_man(kMaxNumFaceDOF,N1nl);
dMatrixT T1;
nVariMatrixT<double> T1_man(kMaxNumFaceDOF,T1);
dArrayT T1n;
VariArrayT<double> T1n_man(kMaxNumFaceDOF,T1n);
dMatrixT Perm(NumSD());
Perm(0,0) = 0.0 ; Perm(0,1) = -1.0;
Perm(1,0) = 1.0 ; Perm(1,1) =  0.0;
double alpha;

iArrayT node_tag(1);
iArrayT xc(1);
iArray2DT xe(1,fNumMultipliers);
ElementMatrixT xL(fNumMultipliers,ElementMatrixT::kNonSymmetric);
dMatrixT nXn(NumSD());

dArrayT e1(fNumMultipliers);
e1[0] = 1; e1[1] = 0;
dArrayT e2(fNumMultipliers);
e2[0] = 0; e2[1] = 1;

int nsd = NumSD();
for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	len_eqnum1  = num_nodes*nsd;
	len_xeqnum1 = num_nodes*fNumMultipliers;
	LHS_man.SetDimensions(len_eqnum1);
	N1_man.SetDimensions(len_eqnum1, nsd);
	P1_man.SetDimensions(len_xeqnum1,fNumMultipliers);
	T1_man.SetDimensions(len_eqnum1, nsd);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	xeqnums1_man.SetMajorDimension(num_nodes,false);

	/* form stiffness */
	for (int f = 0;  f < faces.Length(); f++) {
		/* primary face */
		const FaceT* face = faces[f];
		bool has_projections = 1, has_multipliers = 1;
		for (int i = 0 ; i < num_nodes ; i++) {
			node = nodes[face->Node(i)];
			if (!(node->HasProjection())) has_projections = 0;
			if (!(node->HasMultiplier())) has_multipliers = 0;
		}
		if ( ! has_projections || ! has_multipliers ) {
			for (int i = 0 ; i < num_nodes ; i++) {
				node = nodes[face->Node(i)];
				if (node->EnforcementStatus() == kPZero) {
					double pScale = fEnforcementParameters(0,1)[kPScale];
					xL(0,0) = pScale;
					xL(1,1) = pScale;
					node_tag[0] = node->Tag();
#if PRINT_DEBUG
					cout << "Zeroing multipliers for node : " << node_tag[0] << "\n";
#endif
					surface.MultiplierTags(node_tag,xc);
					const ElementSupportT& support = ElementSupport();
					support.XDOF_Manager().XDOF_SetLocalEqnos
						(Group(), xc, xe);
					support.AssembleLHS(Group(), xL, xe,xe);
				}
			} 
		}
		else {
		/* get equation numbers */
		const iArrayT& conn1 = face->GlobalConnectivity();
		Field().SetLocalEqnos(conn1, eqnums1);
		surface.MultiplierTags(face->Connectivity(),xconn1);
		ElementSupport().XDOF_Manager()
				.XDOF_SetLocalEqnos(Group(), xconn1, xeqnums1);

		LHS = 0.0;
		elem_in_contact = 0;
		face->Quadrature(points,weights);
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			xi = points(i);
			status = node->EnforcementStatus();

			/* pressure shape function matrix */
			face->ComputeShapeFunctions(points(i),P1);

			if (status > kPZero )  {
				elem_in_contact = 1;
				/* set-up */
				const FaceT* opp_face = node->OpposingFace();
				opp_num_nodes = opp_face->NumNodes();
				len_eqnum2 = opp_num_nodes*nsd;
				len_xeqnum2 = opp_num_nodes*fNumMultipliers;
                const double* opp_xi  = node->OpposingLocalCoordinates();
				const ContactSurfaceT* opp_surf = node->OpposingSurface();
				opp_surf_tag = opp_surf->Tag();
				dArrayT& parameters =
                	fEnforcementParameters(surf_tag,opp_surf_tag);
				consistent = (int) parameters[kConsistentTangent];

/* BLM ============================================================= */
/* K =  N1n (x) { e1.P1 + pen (n1.N2 * D u2 -  n1.N1 * D u1) }  */
                if (status == kGapZero || status == kPJump){
					face->ComputeShapeFunctions(points(i),N1);				
					for (int j =0; j < NumSD(); j++) {
						n1[j] = node->Normal()[j];
					    l1[j] = node->Tangent()[j];
					}	
/* primary U (X) primary   P block */
					tmp_LHS_man.SetDimensions(len_eqnum1,len_xeqnum1);
					nXn(0,0) = n1[0]; nXn(0,1) = l1[0];
					nXn(1,0) = n1[1]; nXn(1,1) = l1[1];
					tmp_LHS.MultABCT(N1,nXn,P1);
					tmp_LHS.SetToScaled(weights[i], tmp_LHS);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, eqnums1,xeqnums1);
                	if (status == kGapZero){
/* primary U (X) primary   U block */
						tmp_LHS_man.SetDimensions(len_eqnum1,len_eqnum1);
						sfac = parameters[kPenalty];
						if (0) { // was consistent
						} else {
							nXn(0,0) = n1[0]*n1[0]; nXn(0,1) = n1[0]*n1[1];
							nXn(1,0) = n1[0]*n1[1]; nXn(1,1) = n1[1]*n1[1];
							nXn(0,0)+= l1[0]*l1[0]; nXn(0,1)+= l1[0]*l1[1];
							nXn(1,0)+= l1[0]*l1[1]; nXn(1,1)+= l1[1]*l1[1];
							LHS.MultABCT(N1,nXn,N1);
							LHS.SetToScaled(sfac*weights[i], LHS);
						}
						ElementSupport().AssembleLHS(Group(), LHS, eqnums1);

/* primary U (X) secondary U block */
						tmp_LHS_man.SetDimensions(len_eqnum1,len_eqnum2);
						N2_man.SetDimensions(len_eqnum2,nsd);
						opp_face->ComputeShapeFunctions (opp_xi,N2);

						/* gap */
						nXn.Outer(n1,n1);
						tmp_LHS.MultABCT(N1,nXn, N2);
						tmp_LHS.SetToScaled(-sfac*weights[i], tmp_LHS);

						/* get equation numbers and assemble */
						const iArrayT& conn2 = opp_face->GlobalConnectivity();
						eqnums2_man.SetMajorDimension(opp_num_nodes,false);
						Field().SetLocalEqnos(conn2, eqnums2);
						ElementSupport().AssembleLHS
								(Group(), tmp_LHS, eqnums1,eqnums2);

						/* slip */
                        node->OriginalOpposingFace()->ComputeShapeFunctions 
							(node->OriginalLocalCoordinates(),N2);
                        nXn.Outer(l1,l1);
                        tmp_LHS.MultABCT(N1,nXn, N2);
                        tmp_LHS.SetToScaled(-sfac*weights[i], tmp_LHS);

                        /* get equation numbers and assemble */
                        const iArrayT& conn2o = node->
							OriginalOpposingFace()->GlobalConnectivity();
                        Field().SetLocalEqnos(conn2o, eqnums2);
                        ElementSupport().AssembleLHS
                                (Group(), tmp_LHS, eqnums1,eqnums2);


					}
				} 
/* Constraint ====================================================== */
                if (status == kGapZero) {
					double gfac = parameters[kGScale];
/* primary P (X) primary   U block */
					tmp_LHS_man.SetDimensions(len_xeqnum1,len_eqnum1);
					if (0) {
					} else {
						nXn(0,0) = n1[0]; nXn(0,1) = n1[1];
						nXn(1,0) = l1[0]; nXn(1,1) = l1[1];
						tmp_LHS.MultABCT(P1,nXn,N1);
					}
					tmp_LHS.SetToScaled(gfac, tmp_LHS);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, xeqnums1,eqnums1);
/* primary P (X) secondary U block */
					tmp_LHS_man.SetDimensions(len_xeqnum1,len_eqnum2);

					/* gap */
					nXn(0,0) = n1[0]; nXn(0,1) = n1[1];
					nXn(1,0) = 0.0;   nXn(1,1) = 0.0;
					opp_face->ComputeShapeFunctions (opp_xi,N2);
					tmp_LHS.MultABCT(P1,nXn,N2);
					tmp_LHS.SetToScaled(-gfac, tmp_LHS);

					/* get equation numbers and assemble */
					const iArrayT& conn2 = opp_face->GlobalConnectivity();
					eqnums2_man.SetMajorDimension(opp_num_nodes,false);
					Field().SetLocalEqnos(conn2, eqnums2);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, xeqnums1,eqnums2);

					/* slip */
					nXn(0,0) = 0.0;   nXn(0,1) = 0.0;  
					nXn(1,0) = l1[0]; nXn(1,1) = l1[1];
					node->OriginalOpposingFace()->ComputeShapeFunctions 
						(node->OriginalLocalCoordinates(),N2);
					tmp_LHS.MultABCT(P1,nXn,N2);
					tmp_LHS.SetToScaled(-gfac, tmp_LHS);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, xeqnums1,eqnums2);
                }
                else if (status == kPJump) {
					double pfac = parameters[kPScale];
/* primary P (X) primary   P block */
					tmp_LHS_man.SetDimensions(len_xeqnum1, len_xeqnum1);
					tmp_LHS.MultABT(P1, P1);
					tmp_LHS.SetToScaled(pfac, tmp_LHS);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, xeqnums1,xeqnums1);

/* primary P (X) secondary P block */		
					tmp_LHS_man.SetDimensions (len_xeqnum1,len_xeqnum2);
					xconn2_man.SetLength(opp_num_nodes,false);
					P2_man.SetDimensions(len_xeqnum1,fNumMultipliers);
					opp_face->ComputeShapeFunctions(opp_xi,P2);
					tmp_LHS.MultABT(P1, P2);
					tmp_LHS.SetToScaled(-pfac, tmp_LHS);

					/* get equation numbers and assemble */
					opp_surf->MultiplierTags
							(opp_face->Connectivity(),xconn2);
					xeqnums2_man.SetMajorDimension(opp_num_nodes,false);
					ElementSupport().XDOF_Manager()
							.XDOF_SetLocalEqnos(Group(), xconn2, xeqnums2);
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, xeqnums1,xeqnums2);
                }
			}
            else if (status == kPZero) {

/* primary P (X) primary   P block */
				tmp_LHS_man.SetDimensions (len_xeqnum1,len_xeqnum1);
				tmp_LHS.MultABT(P1, P1);
				double pScale = fEnforcementParameters(0,1)[kPScale];
				tmp_LHS.SetToScaled(pScale, tmp_LHS);
				ElementSupport().AssembleLHS
						(Group(), tmp_LHS, xeqnums1,xeqnums1);
			}
		}
	}
	}
}
}
