/* $Id: PenaltyContactElement2DT.cpp,v 1.54 2007-10-09 23:24:47 rjones Exp $ */
#include "PenaltyContactElement2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ofstreamT.h"
#include "ifstreamT.h"
#include "ContactNodeT.h"
#include "ParabolaT.h"
#include "ModSmithFerrante.h"
#include "GreenwoodWilliamson.h"
#include "MajumdarBhushan.h"
#include "GWPlastic.h"
#include "ofstreamT.h"

/* vector functions */
#include "vector2D.h"

/* constants */
const double PI = 2.0*acos(0.0);

/* parameters */
static const int kMaxNumFaceNodes = 4;
static const int kMaxNumFaceDOF   = 12;

using namespace Tahoe;

/* constructor */
PenaltyContactElement2DT::PenaltyContactElement2DT(const ElementSupportT& support):
	ContactElementT(support)
{
	SetName("Jones_penalty_contact_2D");
}

/* accept parameter list */
void PenaltyContactElement2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContactElementT::TakeParameterList(list);

  fContactArea.Dimension(fSurfaces.Length());
	fContactArea = 0.0;
	/* subsidary data for GW models */
  fRealArea.Dimension(fSurfaces.Length());
	fRealArea = 0.0;
  fPlasticArea.Dimension(fSurfaces.Length());
	fPlasticArea = 0.0;
}

/* print/compute element output quantities */
void PenaltyContactElement2DT::WriteOutput(void)
{
	/* call base class */
	ContactElementT::WriteOutput();

	int step_num = ElementSupport().StepNumber();

	StringT filename;
	filename.Root(ElementSupport().InputFile());
	filename.Append(".contact_data");
	ofstreamT data_file;
	if (fFirstPass) {
		data_file.open(filename);
        data_file << "# 1.step 2.surf 3.node 4.x 5.y 6.gap 7.gmin 8.f 9.df \n";
		fFirstPass = false;
	}
	else
		data_file.open_append(filename);

    if (fOutputFlags[kArea] )
 	{
		cout << "\n";
		for (int i=0; i<fSurfaces.Length(); i++) {
			cout << "Surface : " << i << 
					"         contact area = " << fContactArea[i] << "\n";
			cout << "Surface : " << i << 
					"    real contact area = " << fRealArea[i] << "\n";
			cout << "Surface : " << i << 
					" plastic contact area = " << fPlasticArea[i] << "\n";
		}
	}

    if (fOutputFlags[kGaps] ) {
    for(int s = 0; s < fSurfaces.Length(); s++) {
        ContactSurfaceT& surface = fSurfaces[s];
	    const ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
        for (int j = 0; j < nodes.Length(); j++) {
            ContactNodeT* node = nodes[j];
			if (node->Status() > ContactNodeT::kNoProjection) {
				int s2 = node->OpposingFace()->Surface().Tag();
				dArrayT& enf_parameters = fEnforcementParameters(s,s2);
				C1FunctionT*  pen_function 
				  = fPenaltyFunctions[LookUp(s,s2,fSurfaces.Length())];
				double gap = node->Gap();
				double f = -enf_parameters[kPenalty]*pen_function->DFunction(gap);
				double df = -enf_parameters[kPenalty]*pen_function->DDFunction(gap);
            	double gmin = node->MinGap();
				int node_num = surface.NodeNumber(j);
            	data_file <<  step_num << " " << s << " " << ++node_num << " "   // node->Tag()
					<< node->Position()[0] << " " << node->Position()[1] << "       "
					<< gap << " " << gmin << " " << f << " " << df 
					<< "\n";
			}

        }
	}
	}

}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called before LHSDriver during iteration process */
void PenaltyContactElement2DT::RHSDriver(void)
{ /* form RESIDUAL */ 
  /* update kinematic data */
  UpdateContactConfiguration();

  bool in_contact = 0;
  ContactNodeT* node;
  double gap, gmin=0.0, pre;
  int num_nodes,s2;


  int num_surfaces = fSurfaces.Length(); 
  int nsd = NumSD();
  fContactArea = 0.0;
  fRealArea = 0.0;
  fPlasticArea = 0.0;

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
		  dArrayT& enf_parameters = fEnforcementParameters(s,s2);
		  dArrayT& mat_parameters = fMaterialParameters(s,s2);
		  C1FunctionT*  pen_function 
				  = fPenaltyFunctions[LookUp(s,s2,num_surfaces)];
		  if (enf_parameters[kMaterialType] 
						  == PenaltyContactElement2DT::kGWPlastic)
		  {
			gmin = node->MinGap(); 
			((GWPlastic*) pen_function)->ResetParameters(gmin);
		  }
		  /* First derivative represents force */
          pre  = -enf_parameters[kPenalty] * pen_function->DFunction(gap);
		  node->nPressure() = pre; // store value on ContactNode for output
		  face->ComputeShapeFunctions(points(i),N1);
		  for (int j =0; j < nsd; j++) {n1[j] = node->Normal()[j];}
		  N1.Multx(n1, tmp_RHS);
		  /* pressure = -e <g> and t = - p n   */
		  if (enf_parameters[kMaterialType] 
						  == PenaltyContactElement2DT::kMajumdarBhushan)
		  {
			double frdim = mat_parameters[kFractalDimension];
			double jac = face->ComputeJacobian(points(i));
			pre *= pow(jac,0.5*(1.0-frdim));
		  }
		  tmp_RHS.SetToScaled(-pre*weights[i], tmp_RHS);
		  RHS += tmp_RHS;


		  fContactArea[s] += weights[i];

		  /* real area computation */
		  if (enf_parameters[kMaterialType] 
				== PenaltyContactElement2DT::kGreenwoodWilliamson) {
			double gw_m = mat_parameters[kAsperityHeightMean];
			double gw_s = mat_parameters[kAsperityHeightStandardDeviation];
			double gw_dens = mat_parameters[kAsperityDensity];
			double gw_mod = mat_parameters[kHertzianModulus];
			double gw_rad = mat_parameters[kAsperityTipRadius];

			GreenwoodWilliamson GWArea(1.0,gw_m,gw_s); 	
  		  	double area_coeff = PI*gw_dens*gw_rad;
		  	fRealArea[s] += (area_coeff*GWArea.Function(gap)*weights[i]);
		  } 
		  else if (enf_parameters[kMaterialType] 
				== PenaltyContactElement2DT::kMajumdarBhushan) {
			double mb_s = mat_parameters[kAsperityHeightStandardDeviation];
			double mb_f = mat_parameters[kFractalDimension];
			double mb_c = mat_parameters[kAreaFraction];

			MajumdarBhushan MBArea(mb_f,mb_s,mb_c); 
			double area_coeff = 0.5*mb_f/(mb_f-2.0); //0.5/mb_f;
		  	fRealArea[s] += (area_coeff*MBArea.Function(gap)*weights[i]);
		  }
		  else if (enf_parameters[kMaterialType] 
				== PenaltyContactElement2DT::kGWPlastic) {
			double gp_dens = mat_parameters[kDensity];
			double area_coeff = gp_dens;
			double real_area = ((GWPlastic*) pen_function)->Function(gap);
		  	fRealArea[s] += (area_coeff*real_area*weights[i]);
			double plastic_area = 0.0;
			if (real_area > 0.0) {
				gmin = node->MinGap(); 
				if (gap < gmin) gmin = gap;
				plastic_area = ((GWPlastic*) pen_function)->PlasticArea(gmin);
			}
		  	fPlasticArea[s] += (area_coeff*plastic_area*weights[i]);
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
					dArrayT& enf_parameters = fEnforcementParameters(s,s2);		
					dArrayT& mat_parameters = fMaterialParameters(s,s2);		
					C1FunctionT* pen_function 
							= fPenaltyFunctions[LookUp(s,s2,num_surfaces)];
					if (enf_parameters[kMaterialType]
							== PenaltyContactElement2DT::kGWPlastic)
					{
						double gmin = node->MinGap(); 
						((GWPlastic*) pen_function)->ResetParameters(gmin);
					}
					pre  = enf_parameters[kPenalty]
					     * pen_function->DFunction(gap);
					dpre_dg = enf_parameters[kPenalty]
							* pen_function->DDFunction(gap);
					double jac = face->ComputeJacobian(points(i));
					if (enf_parameters[kMaterialType] 
							== PenaltyContactElement2DT::kMajumdarBhushan)
		  			{
							double frdim = mat_parameters[kFractalDimension];
							pre *= pow(jac,0.5*(1.0-frdim));
							dpre_dg *= pow(jac,0.5*(1.0-frdim));
		  			}

					consistent = (int) enf_parameters[kConsistentTangent];
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
					ElementSupport().AssembleLHS
							(Group(), tmp_LHS, eqnums1,eqnums2);

					/* Part:  dx1 (x) dx1 */
					tmp_LHS.Outer(N1nl, N1n);
					tmp_LHS.SetToScaled(dpre_dg*weights[i], tmp_LHS);
					LHS += tmp_LHS;
					if (consistent) {
						face->ComputeShapeFunctionDerivatives(points(i),T1);

						Z1.MultABT(T1,Perm);
						tmp_LHS.MultABT(N1, Z1);
						
						if (enf_parameters[kMaterialType] 
							== PenaltyContactElement2DT::kMajumdarBhushan)
						{
							double frdim = mat_parameters[kFractalDimension];
							double mbexp = 0.5*(1.0-frdim);
							// this just adding d pre /d jac
							tmp_LHS.SetToScaled
								(-pre*weights[i]*(1.0+mbexp/jac)/jac, tmp_LHS);
						}
		  				else
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
