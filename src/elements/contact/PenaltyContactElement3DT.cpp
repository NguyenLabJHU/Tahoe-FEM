/* $Id: PenaltyContactElement3DT.cpp,v 1.15 2005-04-14 01:18:53 paklein Exp $ */
#include "PenaltyContactElement3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ContactNodeT.h"
#include "ParabolaT.h"
#include "ModSmithFerrante.h"
#include "GreenwoodWilliamson.h"
#include "ofstreamT.h"

/* vector functions */
#include "vector3D.h"

/* constants */
const double PI = 2.0*acos(0.0);

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

using namespace Tahoe;

/* constructor */
PenaltyContactElement3DT::PenaltyContactElement3DT(const ElementSupportT& support):
	ContactElementT(support)
{
	SetName("Jones_penalty_contact_3D");
}

/* accept parameter list */
void PenaltyContactElement3DT::TakeParameterList(const ParameterListT& list)
{
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
				ExceptionT::GeneralFail("PenaltyContactElement3DT::TakeParameterList",
					"expecting %d enforcement parameters not %d",
					kNumEnfParameters, enf_parameters.Length());
			
			/* only print allocated parameter arrays */
			if (search_parameters.Length() == kSearchNumParameters) {
		  	  out << "  surface pair: ("  << i << "," << j << ")\n" ;
			  out << "  gap tolerance:    "
					<< search_parameters[kGapTol] << '\n';
			  out << "  xi tolerance :    "
					<< search_parameters[kXiTol] << '\n';
			  out << "  pass flag    :    "
					<< (int) search_parameters[kPass] << '\n';
			  out << "  consis. tangent:  "
                    << (int) enf_parameters[kConsistentTangent] << '\n';
			  out << "  penalty :         "
					<< enf_parameters[kPenalty] << '\n';
			  out << "  penalty types:\n"
				  << "     Default (linear)    " 
 				  << PenaltyContactElement3DT::kDefault << "\n"
				  << "     ModSmithFerrante    " 
				  << PenaltyContactElement3DT::kModSmithFerrante << "\n"
				  << "     GreenwoodWilliamson " 
			      << PenaltyContactElement3DT::kGreenwoodWilliamson << "\n";
			  out << "  penalty Type :         "
					<< (int) enf_parameters[kMaterialType] << '\n';
			  switch ((int) enf_parameters[kMaterialType]) 
			  {
			  case kDefault: // no other parameters
			    out << "  <no parameters> \n";
				break;	
			  case kModSmithFerrante:
				out << "  Smith-Ferrante A : "
					<< enf_parameters[kSmithFerranteA] << '\n';
				out << "  Smith-Ferrante B : "
					<< enf_parameters[kSmithFerranteB] << '\n';
				break;	
			  case kGreenwoodWilliamson:
				out << "  Average asperity height            : "
					<< enf_parameters[kAsperityHeightMean] << '\n';
				out << "  Asperity height standard deviation : "
					<< enf_parameters[kAsperityHeightStandardDeviation] << '\n';
				out << "  Asperity density                   : "
					<< enf_parameters[kAsperityDensity] << '\n';
				out << "  Asperity Radius                    : "
					<< enf_parameters[kAsperityTipRadius] << '\n';
				out << "  Hertzian Modulus                   : "
					<< enf_parameters[kHertzianModulus] << '\n';
				break;	
			  default:
				throw ExceptionT::kBadInputValue;
		  	  }
			}
		}
	}
	out <<'\n';

	/* set up Penalty functions */
	fPenaltyFunctions.Dimension(num_surfaces*(num_surfaces-1));
    for (int i = 0; i < num_surfaces ; i++)
    {
        for (int j = 0 ; j < num_surfaces ; j++)
        {
          dArrayT& enf_parameters = fEnforcementParameters(i,j);
		  if (enf_parameters.Length()) {
            dArrayT& parameters = fEnforcementParameters(i,j);
			switch ((int) enf_parameters[kMaterialType]) 
			{
			case PenaltyContactElement3DT::kDefault:
				// Macauley bracket:  <-x> ???
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] = new ParabolaT(1.0);
				break;
			case PenaltyContactElement3DT::kModSmithFerrante:
				{
                double A = parameters[kSmithFerranteA];
                double B = parameters[kSmithFerranteB];
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] = new ModSmithFerrante(A,B);
				//parameters[kPenalty] = 1.0; // overwrite pen value
				}
				break;
			case PenaltyContactElement3DT::kGreenwoodWilliamson:
				{
                /* parameters for Greenwood-Williamson load formulation */
                double gw_m = parameters[kAsperityHeightMean];
                double gw_s = parameters[kAsperityHeightStandardDeviation];
                double gw_dens = parameters[kAsperityDensity];
                double gw_mod = parameters[kHertzianModulus];
                double gw_rad = parameters[kAsperityTipRadius];
                double material_coeff=(4.0/3.0)*gw_dens*gw_mod*sqrt(gw_rad);
          		double area_coeff = PI*gw_dens*gw_rad;
				parameters[kPenalty] *= material_coeff; // overwrite pen value
				fPenaltyFunctions[LookUp(i,j,num_surfaces)]
                                        = new GreenwoodWilliamson(1.5,gw_m,gw_s);
				}
				break;
			default:
				throw ExceptionT::kBadInputValue;
			}
		  }
		}
	}

	/* subsidary data for GW models */
    fRealArea.Dimension(fSurfaces.Length());
	fRealArea = 0.0;
}

/* print/compute element output quantities */
void PenaltyContactElement3DT::WriteOutput(void)
{
	/* call base class */
	ContactElementT::WriteOutput();
	
	if (fOutputFlags[kArea] ) 
//  (parameters[kMaterialType]==PenaltyContactElement3DT::kGreenwoodWilliamson))
	{
		cout << "\n";
		for (int i=0; i<fSurfaces.Length(); i++)
			if ( fRealArea[i] > 0.0) cout << "Surface" << i 
			     << ": real contact area = " << fRealArea[i] << "\n";
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called before LHSdriver during iteration process */
void PenaltyContactElement3DT::RHSDriver(void)
{ /* form RESIDUAL */ 
  /* update kinematic data */
  UpdateContactConfiguration();

  bool in_contact = 0;
  ContactNodeT* node;
  double gap, pre;
  int num_nodes,s2;


  int num_surfaces = fSurfaces.Length(); 
  int nsd = NumSD();
  /* residual */
  for(int s = 0; s < num_surfaces; s++) {
	ContactSurfaceT& surface = fSurfaces[s];
	/* all faces on a surface the same size */
	const ArrayT<FaceT*>& faces = surface.Faces();
	const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	num_nodes = surface.NumNodesPerFace();
	RHS_man.SetLength(num_nodes*nsd,false);
	tmp_RHS_man.SetLength(num_nodes*nsd,false);
	N1_man.SetDimensions(num_nodes*nsd, nsd);
	weights_man.SetLength(num_nodes,false);
	eqnums1_man.SetMajorDimension(num_nodes,false);
	fRealArea[s] = 0;
		
	/*form residual for this surface */
	for (int f = 0;  f < faces.Length(); f++) {
	  const FaceT* face = faces[f];
	  face->Quadrature(points,weights);
	  /* primary face */
	  const iArrayT& conn1 = face->GlobalConnectivity();
	  RHS = 0.0;
	  in_contact = 0;
	  /*loop over (nodal) quadrature points */
	  /*NOTE: these CORRESPOND to local node numbers */
	  for (int i = 0 ; i < weights.Length() ; i++) {
		node = nodes[face->Node(i)];
		if (node->Status() > ContactNodeT::kNoProjection) {
		  in_contact = 1;
		  s2 = node->OpposingFace()->Surface().Tag();
		  gap =  node->Gap();
		  dArrayT& parameters = fEnforcementParameters(s,s2);
		  /* First derivative represents force */
          pre  = -parameters[kPenalty] * fPenaltyFunctions[LookUp(s,s2,num_surfaces)]->DFunction(gap);
		  node->nPressure() = pre; // store value on ContactNode for output
		  face->ComputeShapeFunctions(points(i),N1);
		  for (int j =0; j < nsd; j++) {n1[j] = node->Normal()[j];}
		  N1.Multx(n1, tmp_RHS);
		  /* pressure = -e <g> and t = - p n   */
		  tmp_RHS.SetToScaled(-pre*weights[i], tmp_RHS);
		  RHS += tmp_RHS;


		  /* real area computation */
		  if (parameters[kMaterialType] 
				== PenaltyContactElement3DT::kGreenwoodWilliamson) {
			double gw_m = parameters[kAsperityHeightMean];
			double gw_s = parameters[kAsperityHeightStandardDeviation];
			double gw_dens = parameters[kAsperityDensity];
			double gw_mod = parameters[kHertzianModulus];
			double gw_rad = parameters[kAsperityTipRadius];

			GreenwoodWilliamson GWArea(1.0,gw_m,gw_s); 	
  		  	double area_coeff = PI*gw_dens*gw_rad;
		  	fRealArea[s] += (area_coeff*GWArea.Function(gap)*weights[i]);
		  }
		}
	  } 
          /* get equation numbers */
          Field().SetLocalEqnos(conn1, eqnums1);
          /* assemble */
          if (in_contact) ElementSupport().AssembleRHS(Group(), RHS, eqnums1);
	}
  }
}

void PenaltyContactElement3DT::LHSDriver(GlobalT::SystemTypeT)
{ /* form STIFFNESS */
  /* del g =  (n1.N2 * del u2 -  n1.N1 * del u1) */
  /* primary (X) primary block */
  /* primary (X) secondary block */

  bool in_contact;
  int consistent;
  int num_nodes,opp_num_nodes;
  int s2;
  ContactNodeT* node;
  double pre=0.0, dpre_dg=0.0;
  double gap=0.0;
  double l11[3], l12[3], l21[3], l22[3];
  dArrayT n1alphal1;
  n1alphal1.Dimension(NumSD());

  /* for consistent stiffness */    
  dArrayT N1nl;
  VariArrayT<double> N1nl_man(kMaxNumFaceDOF,N1nl);
  dMatrixT T1;
  nVariMatrixT<double> T1_man(kMaxNumFaceDOF,T1);
  dMatrixT T2;
  nVariMatrixT<double> T2_man(kMaxNumFaceDOF,T2);
  dMatrixT Z1;
  nVariMatrixT<double> Z1_man(kMaxNumFaceDOF,Z1);
  dMatrixT Z2;
  nVariMatrixT<double> Z2_man(kMaxNumFaceDOF,Z2);
//dArrayT T1n;
//VariArrayT<double> T1n_man(kMaxNumFaceDOF,T1n);

  dMatrixT W1(NumSD()), W2(NumSD());
  dMatrixT LL(2), LLi(2);
  double deti, b[2];


	int nsd = NumSD();
    int num_surfaces = fSurfaces.Length(); 
	for(int s = 0; s < num_surfaces; s++) {
		ContactSurfaceT& surface = fSurfaces[s];
        /* all faces on a surface the same size */
        const ArrayT<FaceT*>& faces = surface.Faces();
        const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
		num_nodes = surface.NumNodesPerFace();	
        LHS_man.SetDimensions(num_nodes*nsd);
        tmp_LHS_man.SetDimensions(num_nodes*nsd);
        N1_man.SetDimensions(num_nodes*nsd, nsd);
        N1n_man.SetLength(num_nodes*nsd,false);
        N1nl_man.SetLength(num_nodes*nsd,false);
  		/* for consistent stiffness */    
        T1_man.SetDimensions(num_nodes*nsd, nsd);
        T2_man.SetDimensions(num_nodes*nsd, nsd);
        Z1_man.SetDimensions(num_nodes*nsd, nsd);
        Z2_man.SetDimensions(num_nodes*nsd, nsd);
//      T1n_man.SetLength(num_nodes*nsd,false);
        weights_man.SetLength(num_nodes,false);
        eqnums1_man.SetMajorDimension(num_nodes,false);
        /* form stiffness */
        for (int f = 0;  f < faces.Length(); f++) {
			const FaceT* face = faces[f];
			face->Quadrature(points,weights);
			/* primary face */
			const iArrayT& conn1 = face->GlobalConnectivity();
			/* get equation numbers */
			Field().SetLocalEqnos(conn1, eqnums1);
			LHS = 0.0;
			in_contact = 0;
			/*loop over (nodal) quadrature points */
			/*NOTE: these CORRESPOND to local node numbers */
			for (int i = 0 ; i < weights.Length() ; i++) {
				node = nodes[face->Node(i)];
                if (node->Status() > ContactNodeT::kNoProjection) {
					in_contact = 1;
					s2 = node->OpposingFace()->Surface().Tag();
					gap =  node->Gap();
					face->ComputeShapeFunctions(points(i),N1);
					dArrayT& parameters = fEnforcementParameters(s,s2);		
					pre  = 
					 parameters[kPenalty]*fPenaltyFunctions[LookUp(s,s2,num_surfaces)]->DFunction(gap);
					dpre_dg = 
					 parameters[kPenalty]*fPenaltyFunctions[LookUp(s,s2,num_surfaces)]->DDFunction(gap);
					consistent = (int) parameters[kConsistentTangent];
					const FaceT* opp_face = node->OpposingFace();
					opp_num_nodes = opp_face->NumNodes();
					N2_man.SetDimensions(opp_num_nodes*nsd, nsd);
					N2n_man.SetLength(opp_num_nodes*nsd,false);
					eqnums2_man.SetMajorDimension(opp_num_nodes,false);
					const iArrayT& conn2 = opp_face->GlobalConnectivity();
					opp_face->ComputeShapeFunctions
						(node->OpposingLocalCoordinates(),N2);
					const double* nm1 = node->Normal();
					for (int j =0; j < nsd; j++) {n1[j] = nm1[j];}
					N1.Multx(n1, N1nl);
					if (consistent)  {
						/* D xi_a,  assuming g approx 0 */
						/* unique tangents of node on surface 1*/
						for (int j =0; j < nsd; j++) { 
							l11[j] =node->Tangent1()[j];
							l12[j] =node->Tangent2()[j];
						}
						/* face tangents on surface 2 */
						opp_face->
						  ComputeTangent1(node->OpposingLocalCoordinates(),l21);
						opp_face->
						  ComputeTangent2(node->OpposingLocalCoordinates(),l22);
						LL(0,0) = Dot(l11,l21);
						LL(1,0) = Dot(l12,l21);
						LL(0,1) = Dot(l11,l22);
						LL(1,1) = Dot(l12,l22);
						deti = 1/(LL(0,0)*LL(1,1)-LL(0,1)*LL(1,0));
						LLi(0,0) = deti * LL(1,1);
						LLi(1,1) = deti * LL(0,0);
						LLi(0,1) = -deti * LL(0,1);
						LLi(1,0) = -deti * LL(1,0);
						
						b[0] = Dot(nm1,l21);
						b[1] = Dot(nm1,l22);
#if 0
cout << "LL" << "\n";
cout << LL(0,0) << " " << LL(0,1) << "\n";
cout << LL(1,0) << " " << LL(1,1) << "\n";
cout << "LLi" << "\n";
cout << LLi(0,0) << " " << LLi(0,1) << "\n";
cout << LLi(1,0) << " " << LLi(1,1) << "\n";
cout << "b" << "\n";
cout << b[0] << " " << b[1] << "\n";

#endif
						 
						for (int j =0; j < nsd; j++) {
							n1[j] -= (b[0]*LLi(0,0) + b[1]*LLi(1,0))*l11[j]
							       + (b[0]*LLi(0,1) + b[1]*LLi(1,1))*l12[j];
						}
					}
					N1.Multx(n1, N1n);
					N2.Multx(n1, N2n); 

					/* Part:  dx1 (x) dx2 */
					tmp_LHS_man.SetDimensions(opp_num_nodes*nsd);
					tmp_LHS.Outer(N1nl, N2n);
					tmp_LHS.SetToScaled(-dpre_dg*weights[i], tmp_LHS);
					/* get equation numbers */
					Field().SetLocalEqnos(conn2, eqnums2);
					/* assemble primary-secondary face stiffness */
					ElementSupport().AssembleLHS(Group(), tmp_LHS, eqnums1,eqnums2);

					/* Part:  dx1 (x) dx1 */
					tmp_LHS.Outer(N1nl, N1n);
					tmp_LHS.SetToScaled(dpre_dg*weights[i], tmp_LHS);
					LHS += tmp_LHS;
					if (consistent) { /* D jn [dx] */
						/* D j */
						face->ComputeShapeFunctionDerivatives(points(i),T1,T2);
						/* face tangents */
						face->ComputeTangent1(points(i),l11);
						face->ComputeTangent2(points(i),l12);
						double jac = face->ComputeJacobian(points(i));
						W2 = 0;
						W2(1,2) = -l12[0]; // -(1),[2,3]
						W2(2,1) =  l12[0]; // -(1),[3,2]
						W2(0,2) =  l12[1]; // -(2),[1,3]
						W2(2,0) = -l12[1]; // -(2),[3,1]
						W2(0,1) = -l12[2]; // -(3),[1,2]
						W2(1,0) =  l12[2]; // -(3),[2,1]
						W1 = 0;
						W1(1,2) = -l11[0];
						W1(2,1) =  l11[0];
						W1(0,2) =  l11[1];
						W1(2,0) = -l11[1];
						W1(0,1) = -l11[2];
						W1(1,0) =  l11[2];
						Z1.MultABT(T1,W2);
						Z2.MultABT(T2,W1);
						Z2.AddScaled(-1.0,Z1);
						tmp_LHS.MultABT(N1, Z2);
						tmp_LHS.SetToScaled(-pre*weights[i]/jac, tmp_LHS);
						LHS += tmp_LHS;
					} 
				}
			}
			/* assemble primary-primary face stiffness */
			if (in_contact) ElementSupport().AssembleLHS(Group(), LHS, eqnums1);
		}
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

