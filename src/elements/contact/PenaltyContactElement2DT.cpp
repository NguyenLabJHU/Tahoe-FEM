/* $Id: PenaltyContactElement2DT.cpp,v 1.37 2003-05-22 22:56:17 rjones Exp $ */
#include "PenaltyContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ContactNodeT.h"
#include "ParabolaT.h"
#include "ModSmithFerrante.h"
#include "GreenwoodWilliamson.h"
#include "MajumdarBhushan.h"

/* vector functions */
#include "vector2D.h"

/* constants */
const double PI = 2.0*acos(0.0);

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

using namespace Tahoe;

/* constructor */
PenaltyContactElement2DT::PenaltyContactElement2DT(const ElementSupportT& support, const FieldT& field):
	ContactElementT(support, field, kNumEnfParameters)
{
}

void PenaltyContactElement2DT::Initialize(void)
{
    ContactElementT::Initialize();

	/* set up Penalty functions */
	int num_surfaces = fSearchParameters.Rows();
	fPenaltyFunctions.Dimension(num_surfaces*(num_surfaces-1));
    for (int i = 0; i < num_surfaces ; i++)
    {
        for (int j = 0 ; j < num_surfaces ; j++)
        {
          dArrayT& parameters = fEnforcementParameters(i,j);
		  if (parameters.Length()) {
			switch ((int) parameters[kPenaltyType]) 
			{
			case PenaltyContactElement2DT::kLinear:
				// Macauley bracket:  <-x> ???
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] = new ParabolaT(1.0);
				break;
			case PenaltyContactElement2DT::kModSmithFerrante:
				{
                double A = parameters[kSmithFerranteA];
                double B = parameters[kSmithFerranteB];
				fPenaltyFunctions[LookUp(i,j,num_surfaces)] = new ModSmithFerrante(A,B);
				//parameters[kPenalty] = 1.0; // overwrite pen value
				}
				break;
			case PenaltyContactElement2DT::kGreenwoodWilliamson:
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
			case PenaltyContactElement2DT::kMajumdarBhushan:
				{
                /* parameters for Majumdar-Bhushan load formulation */
                double mb_s = parameters[kSigma];
                double mb_g = parameters[kRoughnessScale];
                double mb_mod = parameters[kEPrime];
                double mb_f = parameters[kFractalDimension];
                double mb_c = parameters[kAreaFraction];
                double material_coeff=(4.0/3.0)*sqrt(2.0*PI)*mb_mod*pow(mb_g,mb_f-1);
				parameters[kPenalty] *= material_coeff; // overwrite pen value
				fPenaltyFunctions[LookUp(i,j,num_surfaces)]
                                        = new MajumdarBhushan(mb_f,mb_s,mb_c);
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
	fRealArea = 0;
}

/* print/compute element output quantities */
void PenaltyContactElement2DT::WriteOutput(void)
{
	/* call base class */
	ContactElementT::WriteOutput();
	
    if (fOutputFlags[kArea] )
//  (parameters[kPenaltyType]==PenaltyContactElement2DT::kGreenwoodWilliamson))
 	{
		cout << "\n";
		for (int i=0; i<fSurfaces.Length(); i++)
			cout << "real contact area = " << fRealArea[i] << "\n";
	}
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
			  out << "  pass flag    :    "
					<< (int) search_parameters[kPass] << '\n';
			  out << "  consis. tangent:  "
                    << (int) enf_parameters[kConsistentTangent] << '\n';
			  out << "  penalty :         "
					<< enf_parameters[kPenalty] << '\n';
			  out << "  penalty types:\n"
				  << "     Linear              " 
 				  << PenaltyContactElement2DT::kLinear << "\n"
				  << "     ModSmithFerrante    " 
				  << PenaltyContactElement2DT::kModSmithFerrante << "\n"
				  << "     GreenwoodWilliamson " 
			      << PenaltyContactElement2DT::kGreenwoodWilliamson << "\n";
			  out << "  penalty Type :         "
					<< (int) enf_parameters[kPenaltyType] << '\n';
			  switch ((int) enf_parameters[kPenaltyType]) 
			  {
			  case kLinear: // no other parameters
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
			  case kMajumdarBhushan:
			  	out << " Asperity height standard deviation : "
			  		<< enf_parameters[kSigma] << '\n';
			  	out << "  Asperity roughness scale : "
			  		<< enf_parameters[kRoughnessScale] << '\n';
			  	out << "  Fractal dimension : "
			  		<< enf_parameters[kFractalDimension] << '\n';
			  	out << "  Hertzian Modulus                   : "
					<< enf_parameters[kEPrime] << '\n';
				out << "  Area Fraction                   : "
					<< enf_parameters[kAreaFraction] << '\n';
				break;	
			  default:
				throw ExceptionT::kBadInputValue;
		  	  }
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
		  if (parameters[kPenaltyType] == PenaltyContactElement2DT::kMajumdarBhushan)
		  {
			double frdim = parameters[kFractalDimension];
			double jac = face->ComputeJacobian(points(i));
			pre *= pow(jac,0.5*(1.0-frdim));
		  }
		  tmp_RHS.SetToScaled(-pre*weights[i], tmp_RHS);
		  RHS += tmp_RHS;


		  /* real area computation */
		  if (parameters[kPenaltyType] 
				== PenaltyContactElement2DT::kGreenwoodWilliamson) {
			double gw_m = parameters[kAsperityHeightMean];
			double gw_s = parameters[kAsperityHeightStandardDeviation];
			double gw_dens = parameters[kAsperityDensity];
			double gw_mod = parameters[kHertzianModulus];
			double gw_rad = parameters[kAsperityTipRadius];

			GreenwoodWilliamson GWArea(1.0,gw_m,gw_s); 	
  		  	double area_coeff = PI*gw_dens*gw_rad;
		  	fRealArea[s] += (area_coeff*GWArea.Function(gap)*weights[i]);
		  }
		  
		  if (parameters[kPenaltyType] 
				== PenaltyContactElement2DT::kMajumdarBhushan) {
			double mb_s = parameters[kAsperityHeightStandardDeviation];
			double mb_f = parameters[kFractalDimension];
			double mb_c = parameters[kAreaFraction];

			MajumdarBhushan MBArea(mb_f,mb_s,mb_c); 	
		  	fRealArea[s] += (MBArea.Function(gap)*weights[i]);
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

void PenaltyContactElement2DT::LHSDriver(GlobalT::SystemTypeT)
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
  double l1[3],l2[3];

  /* for consistent stiffness */    
   dArrayT N1nl;
   VariArrayT<double> N1nl_man(kMaxNumFaceDOF,N1nl);
   dMatrixT T1;
   nVariMatrixT<double> T1_man(kMaxNumFaceDOF,T1);
   dArrayT T1n;
   VariArrayT<double> T1n_man(kMaxNumFaceDOF,T1n);
   dMatrixT Z1;
   nVariMatrixT<double> Z1_man(kMaxNumFaceDOF,Z1);

   dMatrixT Perm(NumSD());
   // exoII CCW coordinate
   Perm(0,0) =  0.0 ; Perm(0,1) =  1.0;
   Perm(1,0) = -1.0 ; Perm(1,1) =  0.0;
   double alpha;


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
        T1_man.SetDimensions(num_nodes*nsd, nsd);
        Z1_man.SetDimensions(num_nodes*nsd, nsd);
        T1n_man.SetLength(num_nodes*nsd,false);
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
					if (consistent) {
						/* average tangent */
						for (int j =0;j<nsd; j++) {l1[j] =node->Tangent1()[j];}
						node->OpposingFace()->
						  ComputeTangent1(node->OpposingLocalCoordinates(),l2);
						alpha = Dot(nm1,l2)/Dot(l1,l2);
						for (int j =0;j< nsd; j++) {n1[j] -= alpha*l1[j];}

						/* missing component that scales with g */
					}
					N1.Multx(n1, N1n);
					N2.Multx(n1, N2n); 

					/* N1n (x) D g */
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
					if (consistent) {
						face->ComputeShapeFunctionDerivatives(points(i),T1);
						double jac = face->ComputeJacobian(points(i));

						Z1.MultABT(T1,Perm);
						tmp_LHS.MultABT(N1, Z1);
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

