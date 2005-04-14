/* $Id: MultiplierContactElement2DT.cpp,v 1.24 2005-04-14 01:18:53 paklein Exp $ */
// created by : rjones 2001
#include "MultiplierContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "MultiplierContactElement2DT.h"
#include "ContactNodeT.h"
#include "ElementSupportT.h"
#include "XDOF_ManagerT.h"

/* vector functions */
#include "vector2D.h"

#undef  PRINT_DEBUG
#define PRINT_DEBUG 0
#undef  HACK
#define HACK 0

using namespace Tahoe;

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

/* constructor */
MultiplierContactElement2DT::MultiplierContactElement2DT(const ElementSupportT& support):
	ContactElementT(support)
{
	SetName("Jones_multiplier_contact_2D");
	fNumMultipliers = 1;
}

/* accept parameter list */
void MultiplierContactElement2DT::TakeParameterList(const ParameterListT& list)
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
				ExceptionT::GeneralFail("MultiplierContactElement2DT::TakeParameterList",
					"expecting %d enforcement parameters not %d",
					kNumEnfParameters, enf_parameters.Length());

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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called before LHSDriver during iteration process */
void MultiplierContactElement2DT::SetContactStatus(void)
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
void MultiplierContactElement2DT::RHSDriver(void)
{ /* form RESIDUAL */ 
/* update kinematic data */
UpdateContactConfiguration();

/* set status of all surface nodes */
SetContactStatus();

bool elem_in_contact = 0;
int opp_surf_tag=-1, status=-1;
int num_nodes;
ContactNodeT* node;
double gap, /*pen,*/ pre, opp_pre=0.0;

iArrayT node_tag(1);
iArrayT xc(1);
iArray2DT xe(1,fNumMultipliers);
dArrayT xR(1);

int nsd = NumSD();
for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	RHS_man.SetLength(num_nodes*nsd,false);
	xRHS_man.SetLength(num_nodes*fNumMultipliers,false);
	tmp_RHS_man.SetLength(num_nodes*nsd,false);
	N1_man.SetDimensions(num_nodes*nsd, nsd);
	P1_man.SetDimensions(num_nodes,fNumMultipliers);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	xeqnums1_man.SetMajorDimension(num_nodes,false);
	xconn1_man.SetLength(num_nodes,false);
		
//	surface.PrintStatus(cout);
//	surface.PrintMultipliers(cout);

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
			for (int i = 0 ; i < num_nodes ; i++) {
				node = nodes[face->Node(i)];
				if (node->EnforcementStatus() == kPZero) {
					node_tag[0] = node->Tag();
					surface.MultiplierTags(node_tag,xc);
					double pScale = fEnforcementParameters(0,1)[kPScale];
					double pj = - node->Pressure() ;
					xR[0] = pScale*pj;
					const ElementSupportT& support = ElementSupport();
					support.XDOF_Manager().XDOF_SetLocalEqnos
						(Group(), xc, xe);
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
			status = node->EnforcementStatus();
			if (status > kPZero )  {
				const ContactSurfaceT* opp_surf = node->OpposingSurface();
				opp_surf_tag = opp_surf->Tag();
				dArrayT& parameters = 
				fEnforcementParameters(surf_tag,opp_surf_tag);

/* BLM: U dof on primary surface  -------------------------------------*/
				if (status == kGapZero || status == kPJump){
					/* pressure */
					pre = node->Pressure() ;
#if PRINT_DEBUG
					cout << "node: " << node->Tag() << ", pressure: "
							<< pre << "\n"; 
#endif
					/* gap */
					gap = node->Gap();
					/* pressure =  Lagrange multiplier + penalty */
					if (status == kGapZero) 
						{ pre += -parameters[kPenalty]*gap;}
					face->ComputeShapeFunctions(points(i),N1);
					for (int j =0; j < nsd; j++) {n1[j] = node->Normal()[j];}
					N1.Multx(n1, tmp_RHS);
					tmp_RHS.SetToScaled(-pre*weights[i], tmp_RHS);
					RHS += tmp_RHS;
				}

/* Constraint : X dof on primary surface ----------------------------- */
				face->ComputeShapeFunctions(points(i),P1);
				if (status == kGapZero) {
					gap = node->Gap();
					P1.SetToScaled(parameters[kGScale]*gap, P1);
				}
				else if (status == kPJump) {
					/* calculate pressure on opposing face */
					const double* opp_xi  = node->OpposingLocalCoordinates();
					P2values_man.Dimension
						(opp_surf->NumNodesPerFace(),fNumMultipliers);
					opp_surf->MultiplierValues(
						node->OpposingFace()->Connectivity(),P2values);	
					opp_pre = node->OpposingFace()->
						Interpolate(opp_xi,(dArrayT&) P2values);
					double pj = opp_pre - node->Pressure() ;
					P1.SetToScaled(parameters[kPScale]*pj, P1);
				}
				xRHS += P1;
			}
			else if (status == kPZero) {
				double pScale = fEnforcementParameters(0,1)[kPScale];
				double pj = - node->Pressure() ;
				face->ComputeShapeFunctions(points(i),P1);
				P1.SetToScaled(pScale*pj, P1);
				xRHS += P1;
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
#if PRINT_DEBUG
cout << " RHS : " << eqnums1[0] << ", " << eqnums1[1] 
<< ", " << eqnums1[2] << ", " << eqnums1[3] << "\n";
cout << " RHS : " << RHS[0] << ", " << RHS[1] 
<< ", " << RHS[2] << ", " << RHS[3] << "\n";
#endif
			support.XDOF_Manager().XDOF_SetLocalEqnos(Group(), xconn1, xeqnums1);
			support.AssembleRHS(Group(), xRHS, xeqnums1);
#if PRINT_DEBUG
cout << " xRHS : " << xeqnums1[0] << ", " << xeqnums1[1] << "\n";
cout << " xRHS : " << xRHS[0] << ", " << xRHS[1] << "\n";
#endif
		}
		}
	}
}
}

void MultiplierContactElement2DT::LHSDriver(GlobalT::SystemTypeT)
{ /* form STIFFNESS */
bool elem_in_contact = 0;
int opp_surf_tag=-1, status=-1;
int consistent, num_nodes, opp_num_nodes;
ContactNodeT* node;
double sfac, gap;
double l1[3],lm2[3];
dArrayT n1alphal1;
n1alphal1.Dimension(NumSD());

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

int nsd = NumSD();
for(int surf_tag = 0; surf_tag < fSurfaces.Length(); surf_tag++) {
	ContactSurfaceT& surface = fSurfaces[surf_tag];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	LHS_man.SetDimensions(num_nodes*nsd);
	N1_man.SetDimensions(num_nodes*nsd, nsd);
	N1n_man.SetLength(num_nodes*nsd,false);
	//if consistent
	N1nl_man.SetLength(num_nodes*nsd,false);
	T1_man.SetDimensions(num_nodes*nsd, nsd);
	T1n_man.SetLength(num_nodes*nsd,false);
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
					node_tag[0] = node->Tag();
					surface.MultiplierTags(node_tag,xc);
					double pScale = fEnforcementParameters(0,1)[kPScale];
					xL(0,0) = pScale;
					const ElementSupportT& support = ElementSupport();
					support.XDOF_Manager().XDOF_SetLocalEqnos
						(Group(), xc, xe);
					support.AssembleLHS(Group(), xL, xe,xe);
				}
			} 
		}
		else {
		const iArrayT& conn1 = face->GlobalConnectivity();
		surface.MultiplierTags(face->Connectivity(),xconn1);
		ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), xconn1, xeqnums1);
		face->Quadrature(points,weights);
		/* get equation numbers */
		Field().SetLocalEqnos(conn1, eqnums1);
		LHS = 0.0;
		elem_in_contact = 0;
		/*loop over (nodal) quadrature points */
		/*NOTE: these CORRESPOND to local node numbers */
		for (int i = 0 ; i < weights.Length() ; i++) {
			node = nodes[face->Node(i)];
			status = node->EnforcementStatus();
			if (status > kPZero )  {
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
					for (int j =0; j < NumSD(); j++) {n1[j] = node->Normal()[j];}
					N1.Multx(n1, N1n);

/* primary U (X) primary   P block */
					tmp_LHS_man.SetDimensions
						(num_nodes*NumSD(),num_nodes*fNumMultipliers);
					tmp_LHS.Outer(N1n, P1);
					tmp_LHS.SetToScaled(weights[i], tmp_LHS);
					ElementSupport().AssembleLHS(Group(), tmp_LHS, eqnums1,xeqnums1);
#if PRINT_DEBUG
cout << "N1n (x) P1 "  << xeqnums1[0] << ", " << xeqnums1[1] << "\n"
<< eqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) << "\n"
<< eqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) << "\n"
<< eqnums1[2] << " | " << tmp_LHS(2,0) << ", " << tmp_LHS(2,1) << "\n"
<< eqnums1[3] << " | " << tmp_LHS(3,0) << ", " << tmp_LHS(3,1) << "\n";
#endif


                	if (status == kGapZero){
/* primary U (X) primary   U block */
						if (consistent) { 
							/* slip */
							tmp_LHS_man.SetDimensions
								(num_nodes*NumSD(),num_nodes*NumSD());
							for (int j =0; j < NumSD(); j++) 
								{l1[j] = node->Tangent1()[j];}
							opp_face->ComputeTangent1(
							  node->OpposingLocalCoordinates(),lm2);
							alpha = Dot(node->Normal(),lm2) 
							      / Dot(l1,lm2);
							for (int j =0; j < NumSD(); j++) 
								{n1alphal1[j] = n1[j] - alpha*l1[j];}
							N1.Multx(n1alphal1, N1nl);
							face->ComputeShapeFunctionDerivatives(points(i),T1);
							double jac = face->ComputeJacobian(points(i));
							T1.Multx(n1, T1n);
							T1n.AddCombination(sfac*weights[i], N1nl,
								sfac*gap*alpha*weights[i]/jac,T1n);
							LHS.Outer(N1n, T1n);
							/* jacobian */
							T1.MultAB(T1,Perm);
							tmp_LHS.MultABT(N1, T1);
							tmp_LHS.SetToScaled
							 ((node->Pressure()-sfac*gap)*weights[i], tmp_LHS);
							LHS += tmp_LHS;
						} else {
							LHS.Outer(N1n, N1n);
							LHS.SetToScaled(sfac*weights[i], LHS);
#if PRINT_DEBUG
if ( sfac > 0.0 ) {
cout << "N1n (x) N1n " << eqnums1[0] << ", " << eqnums1[1] << ", " << eqnums1[2] << ", " << eqnums1[3] << "\n"
<< eqnums1[0] << " | " << LHS(0,0) << ", " << LHS(0,1) << ", " << LHS(0,2) << ", " << LHS(0,3)  << "\n"
<< eqnums1[1] << " | " << LHS(1,0) << ", " << LHS(1,1) << ", " << LHS(1,2) << ", " << LHS(1,3)  << "\n"
<< eqnums1[2] << " | " << LHS(2,0) << ", " << LHS(2,1) << ", " << LHS(2,2) << ", " << LHS(2,3)  << "\n"
<< eqnums1[3] << " | " << LHS(3,0) << ", " << LHS(3,1) << ", " << LHS(3,2) << ", " << LHS(3,3)  << "\n";
}
#endif
						}
						ElementSupport().AssembleLHS(Group(), LHS, eqnums1);

/* primary U (X) secondary U block */
						/* get connectivity */
						const iArrayT& conn2 = opp_face->GlobalConnectivity();

						N2_man.SetDimensions(opp_num_nodes*NumSD(), NumSD());
						N2n_man.SetLength(opp_num_nodes*NumSD(),false);
						eqnums2_man.SetMajorDimension(opp_num_nodes,false);
						opp_face->ComputeShapeFunctions (opp_xi,N2);
						N2.Multx(n1, N2n); 

						tmp_LHS_man.SetDimensions
							(num_nodes*NumSD(),opp_num_nodes*NumSD());
						tmp_LHS.Outer(N1n, N2n);
						tmp_LHS.SetToScaled(-sfac*weights[i], tmp_LHS);
#if PRINT_DEBUG
if ( sfac > 0.0 ) {
cout << "N1n (x) N2n " << eqnums1[0] << ", " << eqnums1[1] << ", " << eqnums1[2] << ", " << eqnums1[3] << "\n"
<< eqnums2[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) << ", " << tmp_LHS(0,2) << ", " << tmp_LHS(0,3)  << "\n"
<< eqnums2[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) << ", " << tmp_LHS(1,2) << ", " << tmp_LHS(1,3)  << "\n"
<< eqnums2[2] << " | " << tmp_LHS(2,0) << ", " << tmp_LHS(2,1) << ", " << tmp_LHS(2,2) << ", " << tmp_LHS(2,3)  << "\n"
<< eqnums2[3] << " | " << tmp_LHS(3,0) << ", " << tmp_LHS(3,1) << ", " << tmp_LHS(3,2) << ", " << tmp_LHS(3,3)  << "\n";
}
#endif

						/* get equation numbers */
						Field().SetLocalEqnos(conn2, eqnums2);
						/* assemble primary-secondary face stiffness */
						ElementSupport().AssembleLHS(Group(), tmp_LHS, eqnums1,eqnums2);

					}
				} 
/* Constraint ====================================================== */
                if (status == kGapZero) {
					double gfac = parameters[kGScale];
/* primary P (X) primary   U block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,num_nodes*NumSD());
					if (consistent) {
						/* T1n is the consistent D_u1 g(u1,u2) */
						tmp_LHS.Outer(P1, T1n);
					} else {
						tmp_LHS.Outer(P1, N1n);
					}
					tmp_LHS.SetToScaled(gfac, tmp_LHS);
					ElementSupport().AssembleLHS(Group(), tmp_LHS, xeqnums1,eqnums1);
#if PRINT_DEBUG
cout << "P1 (x) N1n "  << eqnums1[0] << ", " << eqnums1[1] << ", " << eqnums1[2] << ", " << eqnums1[3] << "\n"
<< xeqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) << ", " << tmp_LHS(0,2) << ", " << tmp_LHS(0,3)  << "\n"
<< xeqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) << ", " << tmp_LHS(1,2) << ", " << tmp_LHS(1,3)  << "\n";
#endif
/* primary P (X) secondary U block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,opp_num_nodes*NumSD());
					tmp_LHS.Outer(P1, N2n);
					tmp_LHS.SetToScaled(-gfac, tmp_LHS);
					ElementSupport().AssembleLHS(Group(), tmp_LHS, xeqnums1,eqnums2);
#if PRINT_DEBUG
cout << "P1 (x) N2n " << eqnums2[0] << ", " << eqnums2[1] << ", " << eqnums2[2] << ", " << eqnums2[3] << "\n"
<< xeqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) << ", " << tmp_LHS(0,2) << ", " << tmp_LHS(0,3)  << "\n"
<< xeqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) << ", " << tmp_LHS(1,2) << ", " << tmp_LHS(1,3)  << "\n";
#endif
                }
                else if (status == kPJump) {
					double pfac = parameters[kPScale];
/* primary P (X) primary   P block */
					tmp_LHS_man.SetDimensions
						(num_nodes*fNumMultipliers,num_nodes*fNumMultipliers);
					tmp_LHS.Outer(P1, P1);
					tmp_LHS.SetToScaled(pfac, tmp_LHS);
					ElementSupport().AssembleLHS(Group(), tmp_LHS, xeqnums1,xeqnums1);
#if PRINT_DEBUG
cout << "P1 (x) P1 " << xeqnums1[0] << ", " << xeqnums1[1] << "\n"
<< xeqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) <<  "\n"
<< xeqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) <<  "\n";
#endif
/* primary P (X) secondary P block */		
					xconn2_man.SetLength(opp_num_nodes,false);
					opp_surf->MultiplierTags
							(opp_face->Connectivity(),xconn2);
					xeqnums2_man.SetMajorDimension(opp_num_nodes,false);
					ElementSupport().XDOF_Manager().XDOF_SetLocalEqnos(Group(), xconn2, xeqnums2);
					P2_man.SetDimensions
							(opp_num_nodes*fNumMultipliers,fNumMultipliers);
					opp_face->ComputeShapeFunctions(opp_xi,P2);
					tmp_LHS_man.SetDimensions (num_nodes*fNumMultipliers,
							opp_num_nodes*fNumMultipliers);
					tmp_LHS.Outer(P1, P2);
					tmp_LHS.SetToScaled(-pfac, tmp_LHS);
					ElementSupport().AssembleLHS(Group(), tmp_LHS, xeqnums1,xeqnums2);
#if PRINT_DEBUG
cout << "P1 (x) P2 " << xeqnums2[0] << ", " << xeqnums2[1] << "\n"
<< xeqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) <<  "\n"
<< xeqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) <<  "\n";
#endif
                }
			}
            else if (status == kPZero) {
				double pScale = fEnforcementParameters(0,1)[kPScale];
/* primary P (X) primary   P block */
				tmp_LHS_man.SetDimensions
					(num_nodes*fNumMultipliers,num_nodes*fNumMultipliers);
				tmp_LHS.Outer(P1, P1);
				tmp_LHS.SetToScaled(pScale, tmp_LHS);
				ElementSupport().AssembleLHS(Group(), tmp_LHS, xeqnums1,xeqnums1);
#if PRINT_DEBUG
cout << "P1 (x) P1 " << xeqnums1[0] << ", " << xeqnums1[1] << "\n"
<< xeqnums1[0] << " | " << tmp_LHS(0,0) << ", " << tmp_LHS(0,1) <<  "\n"
<< xeqnums1[1] << " | " << tmp_LHS(1,0) << ", " << tmp_LHS(1,1) <<  "\n";
#endif
			}
		}
	}
	}
}
}
