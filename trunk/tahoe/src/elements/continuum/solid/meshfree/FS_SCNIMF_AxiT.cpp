/* $Id: FS_SCNIMF_AxiT.cpp,v 1.20 2005-01-27 22:04:20 cjkimme Exp $ */
#include "FS_SCNIMF_AxiT.h"

//#define VERIFY_B

#include "ArrayT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "fstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "ElementSupportT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "BasicFieldT.h"
#include "ParentDomainT.h"
#include "ParameterContainerT.h"


#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "FSSolidMatT.h"
#include "FSSolidMatList3DT.h"
#include "SolidMatSupportT.h"

/* materials list */
#include "FSSolidMatList2DT.h"

using namespace Tahoe;

const double Pi = acos(-1.0);
const double twoPi = 2.*acos(-1.0);
const int kRadialDirection = 0; /* the x direction is radial */
const int kNSD = 2;
const double kZeroTol = 0.0000001;

/* constructor */
FS_SCNIMF_AxiT::FS_SCNIMF_AxiT(const ElementSupportT& support):
	SCNIMFT(support),
	fFSMatSupport(NULL)
{
	SetName("fd_mfparticle_axi");
}

/* destructor */
FS_SCNIMF_AxiT::~FS_SCNIMF_AxiT(void)
{
	delete fFSMatSupport;
}

/* generate labels for output data */
void FS_SCNIMF_AxiT::GenerateOutputLabels(ArrayT<StringT>& labels)
{

  	/* Reference Configuration */
	const char* ref[2] = {"R", "Z"};

	/* displacement labels */
	int num_labels = 4; // displacements
	int num_stress=4;

	const char* stress[4];
	const char* strain[4];

	stress[0] = "srr";
	stress[1] = "szz";
	stress[2] = "srz";
	stress[3] = "stt";
	strain[0] = "Err";
	strain[1] = "Ezz";
	strain[2] = "Erz";
	strain[3] = "Ett";

	num_labels += 2 * num_stress + 1; 
	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < fSD; dex++)
		labels[dex] = ref[dex];

	const ArrayT<StringT>& disp_labels = Field().Labels();
	for (int ns = 0 ; ns < fSD; ns++)
	  	labels[dex++] = disp_labels[ns];

	labels[dex++] = "mass/rho";
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = strain[ns];
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = stress[ns];
}

void FS_SCNIMF_AxiT::WriteOutput(void)
{
	const char caller[] = "FS_SCNIMF_AxiT::WriteOutput";

	/* dimensions */
	int ndof = NumDOF();
	int num_stress = 4;
	int num_output = 2*ndof + 1 + 2 * num_stress; /* ref. coords, displacements, mass, e, s */

	/* number of output nodes */
	int non = fNodes.Length();

	/* output arrays length number of active nodes */
	dArray2DT n_values(non, num_output), e_values;
	n_values = 0.0;

	/* global coordinates */
	const dArray2DT& coords = ElementSupport().InitialCoordinates();

	/* the field */
	const FieldT& field = Field();
	const dArray2DT* velocities = NULL;
	if (field.Order() > 0) velocities = &(field[1]);

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
		ExceptionT::GeneralFail("FS_SCNIMF_AxiT::RHSDriver","Cannot get material\n");	

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	dMatrixT& F3D = fF_list[0];
	dMatrixT F2D(2);
	dMatrixT BJ(fSD*fSD, fSD);
	dMatrixT E3D(3);

	/* displacements */
	const dArray2DT& u = Field()(0,0);

	/* collect displacements */
	dArrayT vec, values_i;
	double F_33;
	for (int i = 0; i < non; i++) 
	{
		/* set current element */
		fElementCards.Current(i);

		int tag_i = fNodes[i];

		/* values for particle i */
		n_values.RowAlias(i, values_i);

		/* copy in */
		vec.Set(ndof, values_i.Pointer());
		coords.RowCopy(tag_i, vec);
		vec.Set(ndof, values_i.Pointer() + ndof);
		vec = 0.;

		int* nodal_supp = fNodalSupports(i);
		double* phi_i = fNodalPhi(i);
		for (int j = 0; j < fNodalPhi.MinorDim(i); j++)
			vec.AddScaled(*phi_i++, u(*nodal_supp++));

		// Convert initial coordinates to current
		for (int j = 0; j < ndof; j++)
			values_i[j] += vec[j];

		// support size 
		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed deformation gradient
		F2D = 0.0; F_33 = 0.;
		dArrayT* bVec_i = bVectorArray(i);
		double* b_33 = circumferential_B(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			F_33 += u(*supp_i,0) * *b_33++; 
			BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
		}
		F2D.PlusIdentity(); // convert to F
		F3D.Rank2ExpandFrom2D(F2D);
		F_33 += 1.;
		F3D(2,2) = F_33;

		E3D.MultATB(F3D, F3D);
		E3D.PlusIdentity(-1.0);
		E3D *= 0.5;

		/* last deformation gradient */
		if (fF_last_list.Length() > 0) 
		{
			/* last displacement */
			const dArray2DT& u = Field()(-1,0);

			/* destination */
			dMatrixT& F3D = fF_list[0];

			F2D = 0.0; F_33 = 0.;
			dArrayT* bVec_i = bVectorArray(i);
			double* b_33 = circumferential_B(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				bVec_i++;
				F_33 += u(*supp_i,0) * *b_33++; 
				BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
			}
			F2D.PlusIdentity(); // convert to F
			F3D.Rank2ExpandFrom2D(F2D);
			F_33 += 1.;
			F3D(2,2) = F_33;
		}

		const double* stress = fCurrMaterial->s_ij().Pointer();
		double* inp_val = values_i.Pointer() + 2*ndof;

		/* mass */		
		*inp_val++ = fCellVolumes[i] * 2. * Pi * fCellCentroids(i,0);;

		/* strain */
		*inp_val++ = E3D(0,0); //rr
		*inp_val++ = E3D(1,1); //zz
		*inp_val++ = E3D(0,1); //rz
		*inp_val++ = E3D(2,2); //tt

		/* stress */
		*inp_val++ = stress[0]; //rr
		*inp_val++ = stress[1]; //zz
		*inp_val++ = stress[5]; //rz
		*inp_val = stress[2]; //tt
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);

}

/* compute specified output parameter(s) */
void FS_SCNIMF_AxiT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT FS_SCNIMF_AxiT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();
	
	return relax;
}

/* write restart data to the output stream */
void FS_SCNIMF_AxiT::WriteRestart(ostream& out) const
{
	ElementBaseT::WriteRestart(out);
}

/* read restart data to the output stream */
void FS_SCNIMF_AxiT::ReadRestart(istream& in)
{
	ElementBaseT::ReadRestart(in);
}

/* information about subordinate parameter lists */
void FS_SCNIMF_AxiT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SCNIMFT::DefineSubs(sub_list);

	/* element blocks for underlying connectivity -- TEMP */
	sub_list.AddSub("fd_scni_axi_element_block");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FS_SCNIMF_AxiT::NewSub(const StringT& name) const
{
   if  (name == "fd_scni_axi_element_block") {
	  ParameterContainerT* block = new ParameterContainerT(name);

	  block->AddSub("block_ID_list",ParameterListT::Once);
	  block->AddSub("large_strain_material_3D");

	  block->SetSubSource(this);

	  return block;

  }
  else /* inherited */
	return SCNIMFT::NewSub(name);

}

/* accept parameter list */
void FS_SCNIMF_AxiT::TakeParameterList(const ParameterListT& list)
{
	/* we are axisymmetric */
	qIsAxisymmetric = true;

	/* inherited */
	SCNIMFT::TakeParameterList(list);

	/* deformation gradients */
	fF_list.Dimension(1);
	fF_list[0].Dimension(3);

	/* casts are safe since class contructs materials list - just one material */
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	FSSolidMatT* mat = (FSSolidMatT*) pcont_mat;
	if (mat->Need_F_last()) {
		fF_last_list.Dimension(1);
		fF_last_list[0].Dimension(3);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void FS_SCNIMF_AxiT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "FS_SCNIMFT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();
	mat_params.SetName("large_strain_material_3D");	

	/* collect material parameters */
	int num_blocks = all_params.NumLists("fd_scni_axi_element_block");
	for (int i = 0; i < num_blocks; i++) {

		const ParameterListT& block = all_params.GetList("fd_scni_axi_element_block",i);

		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/* form group contribution to the stiffness matrix */
void FS_SCNIMF_AxiT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	/* time integration parameters */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formM = fIntegrator->FormM(constM);

	/* quick exit */
	if ((formM == 0 && formK == 0) ||
	    (fabs(constM) < kSmall &&
	     fabs(constK) < kSmall)) return;

	/* assemble particle mass */
	if (formM) {
	
		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!fCurrMaterial)
		{
			ExceptionT::GeneralFail("FS_SCNIMF_AxiT::LHSDriver","Cannot get material\n");
		}

		AssembleParticleMass(fCurrMaterial->Density());
	}

	if (formK)
	{
		/* hold the smoothed strain */
		dMatrixT& F3D = fF_list[0];
		dMatrixT F2D(2);

		/* displacements */
		const dArray2DT& u = Field()(0,0);

		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!fCurrMaterial)
		{
			ExceptionT::GeneralFail("FS_SCNIMF_AxiT::LHSDriver","Cannot get material\n");
		}

		int nNodes = fNodes.Length();

		/* assembly information */
		int group = Group();
		int ndof = NumDOF();
		fLHS.Dimension(ndof);
		const ElementSupportT& support = ElementSupport();

		const iArray2DT& field_eqnos = Field().Equations();
		iArrayT row_eqnos(ndof); 
		iArrayT col_eqnos(ndof);
		dMatrixT BJ(4, ndof), BK(4, ndof), FTs(2), fStress3D(3), fStress2D(2), fCauchy(3);
		dMatrixT Tijkl(fSD*fSD), BJTCijkl(fSD, fSD*fSD), K_JK, mod2D(3), Finverse(3);
		dArrayT c_theta(4), BKdotc_t(2), BJdotc_t(2);
		K_JK.Alias(fLHS);
		LinkedListT<dArrayT> bVectors_j;
		LinkedListT<int> nodeSupport_j;
		double F_33, c_theta_theta, J;
		for (int i = 0; i < nNodes; i++)
		{	
			/* set current element */
			fElementCards.Current(i);

			double w_i = fCellVolumes[i]*constK*twoPi*fNodalCoordinates(i,0); // integration weights

			int n_supp = nodalCellSupports.MinorDim(i);

			// Compute smoothed deformation gradient
			F2D = 0.0; F_33 = 0.;
			dArrayT* bVec_i = bVectorArray(i);
			double* b_33 = circumferential_B(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				bVec_i++;
				F_33 += u(*supp_i,0) * *b_33++; 
				BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
			}
			F2D.PlusIdentity(); // convert to F
			F3D.Rank2ExpandFrom2D(F2D);
			F_33 += 1.;
			F3D(2,2) = F_33;
			J = F3D.Det();
			if (J <= 0.0)
				ExceptionT::BadJacobianDet("FS_SCNIMF_AxiT::FormKd");
//			fFSMatSupport->SetDeformationGradient(&Flist);

			/* last deformation gradient */
			if (fF_last_list.Length() > 0) 
			{
				/* last displacement */
				const dArray2DT& u = Field()(-1,0);

				/* destination */
				dMatrixT& F3D = fF_list[0];

				F2D = 0.0; F_33 = 0.;
				dArrayT* bVec_i = bVectorArray(i);
				double* b_33 = circumferential_B(i);
				int* supp_i = nodalCellSupports(i);

				for (int j = 0; j < n_supp; j++) { 
					bVectorToMatrix(bVec_i->Pointer(), BJ);
					bVec_i++;
					F_33 += u(*supp_i,0) * *b_33++; 
					BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
				}
				F2D.PlusIdentity(); // convert to F
				F3D.Rank2ExpandFrom2D(F2D);
				F_33 += 1.;
				F3D(2,2) = F_33;
			}
			
			const dMatrixT& cijkl = fCurrMaterial->C_IJKL();
			mod2D.Rank4ReduceFrom3D(cijkl);
			fCurrMaterial->s_ij().ToMatrix(fCauchy);
			Finverse.Inverse(F3D);
			Finverse *= J;   // compute J F^-1
			fStress3D.MultABT(fCauchy, Finverse); // compute PK1
			fStress2D.Rank2ReduceFrom3D(fStress3D);

			// FTs = F^T sigma
			FTs.MultATB(F2D, fStress2D);
			// T_11 = T_22 = FTs_11 T_13 = T_24 = FTs_12
			// T_31 = T_42 = FTs_21 T_33 = T_44 = FTs_22
			Tijkl.Expand(FTs, 2, dMatrixT::kOverwrite);
			// Tijkl = 0.; used to debug material stiffness
			// use brute force, low-level until I find optimal way
			Tijkl[0] += mod2D[0]; // += C_11
			Tijkl[1] += mod2D[6]; // += C_13
			Tijkl[2] += mod2D[6]; // += C_13
			Tijkl[3] += mod2D[3]; // += C_12
			Tijkl[4] += mod2D[6]; // += C_13
			Tijkl[5] += mod2D[8]; // += C_33
			Tijkl[6] += mod2D[8]; // += C_33
			Tijkl[7] += mod2D[7]; // += C_23
			Tijkl[8] += mod2D[6]; // += C_13
			Tijkl[9] += mod2D[8]; // += C_33
			Tijkl[10] += mod2D[8]; // += C_33
			Tijkl[11] += mod2D[7]; // += C_23
			Tijkl[12] += mod2D[3]; // += C_12
			Tijkl[13] += mod2D[7]; // += C_23
			Tijkl[14] += mod2D[7]; // += C_23
			Tijkl.Last() += mod2D[4]; // += C_22

			// c_theta is the last row of mod2D with the last entry dropped
			c_theta[0] = cijkl[2];
			c_theta[1] = cijkl[8];
			c_theta[2] = cijkl[8];
			c_theta[3] = cijkl[32];
			c_theta_theta = cijkl[14];
			
			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			b_33 = circumferential_B(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++)
			{
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				//BJTCijkl.MultATB(BJ, cijkl, 0);
				BJ.MultTx(c_theta, BJdotc_t, *b_33);

				/* simultanesouly compute material and stress stiffnesses */
				BJTCijkl.MultATB(BJ, Tijkl, dMatrixT::kWhole); // accumulate stress stiffness

				col_eqnos.Copy(field_eqnos(*supp_i));

				dArrayT* bVec_j = bVectorArray(i);
				int* supp_j = nodalCellSupports(i);
				double* b_33_k = circumferential_B(i);
				for (int k = 0; k < n_supp; k++)
				{
					bVectorToMatrix(bVec_j->Pointer(), BK);
					bVec_j++;

					// K_JK = BT_J x Cijkl x B_K 
					K_JK.MultAB(BJTCijkl, BK, dMatrixT::kWhole);
					BK.MultTx(c_theta, BKdotc_t, *b_33);
					K_JK[0] += BKdotc_t[0];
					K_JK[2] += BKdotc_t[1]; 
					K_JK[0] += (BJdotc_t[0] + *b_33 * c_theta_theta)* *b_33_k;

					K_JK *= w_i*constK; 

					/* assemble */
					row_eqnos.Copy(field_eqnos(*supp_j++));
					support.AssembleLHS(group, fLHS, col_eqnos, row_eqnos);
					
					b_33_k++;
				}
				b_33++;
			}	
		}
	}
}

/* assemble particle mass matrix into LHS of global equation system */
void FS_SCNIMF_AxiT::AssembleParticleMass(const double rho)
{

  fForce = 0.0;
  int* nodes = fNodes.Pointer();
  double* volume = fCellVolumes.Pointer();
  for (int i = 0; i < fNodes.Length(); i++) {

    double* m = fForce(fNodes[i]);

    for (int j = 0; j < fSD; j++)
      *m++ = *volume * 2. * Pi * fCellCentroids(i,0);
    volume++;
  }

  fForce *= rho;
  
  /* assemble all */
  ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}

void FS_SCNIMF_AxiT::RHSDriver(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver2D";

	/* contribution from natural boundary conditions */
	SCNIMFT::RHSDriver();

	/* check 2D */
	if (NumDOF() != 2) ExceptionT::GeneralFail(caller, "2D only: %d", NumDOF());

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
	{
		ExceptionT::GeneralFail("FS_SCNIMF_AxiT::RHSDriver","Cannot get material\n");
	}	

	int nNodes = fNodes.Length();

	if (formMa) 
	{
		if (Field().Order() < 2)
			ExceptionT::GeneralFail(caller,"Field's Order does not have accelerations\n");
	
		fLHS = 0.0; // fLHS.Length() = nNodes * fSD;
		const dArray2DT& a = Field()(0,2); // accelerations
		double* ma = fLHS.Pointer();
		const double* acc;

		int* nodes = fNodes.Pointer();
		double* volume = fCellVolumes.Pointer();
		for (int i = 0; i < nNodes; i++)
		{
			acc = a(*nodes++);
			for (int j = 0; j < fSD; j++)
				*ma++ = *volume * 2. * Pi * fCellCentroids(i,0) * *acc++;
			volume++;
		}
		fLHS *= fCurrMaterial->Density();
	}

	fForce = 0.0;
	dMatrixT& F3D = fF_list[0];
	dMatrixT BJ(4, 2), fStress3D(3), fStress2D(2), fCauchy(3), Finverse(3);
	double F_33, S_33, J;
	dMatrixT F2D(2);

	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++)
	{
		/* set current element */
		fElementCards.Current(i);

		double w_i = fCellVolumes[i]*twoPi*fNodalCoordinates(i,0); // integration weight

		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed deformation gradient
		F2D = 0.0; F_33 = 0.;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		double* b_33 = circumferential_B(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			F_33 += u(*supp_i,0) * *b_33++; 
			BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
		}
		F2D.PlusIdentity(); // convert to F
		F_33 += 1.; 	
		F3D.Rank2ExpandFrom2D(F2D);
		F3D(2,2) = F_33;
		J = F3D.Det();
		if (J <= 0.0)
			ExceptionT::BadJacobianDet("FS_SCNIMF_AxiT::FormKd");

		/* last deformation gradient */
		if (fF_last_list.Length() > 0) 
		{
			/* last displacement */
			const dArray2DT& u = Field()(-1,0);

			/* destination */
			dMatrixT& F3D = fF_last_list[0];

			F2D = 0.0; F_33 = 0.;
			dArrayT* bVec_i = bVectorArray(i);
			double* b_33 = circumferential_B(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				bVec_i++;
				F_33 += u(*supp_i,0) * *b_33++; 
				BJ.Multx(u(*supp_i++), F2D.Pointer(), 1.0, dMatrixT::kAccumulate);
			}
			F2D.PlusIdentity(); // convert to F
			F3D.Rank2ExpandFrom2D(F2D);
			F_33 += 1.;
			F3D(2,2) = F_33;
		}

		fCurrMaterial->s_ij().ToMatrix(fCauchy);
		Finverse.Inverse(F3D);
		Finverse *= J;   // compute J F^-1
		fStress3D.MultABT(fCauchy, Finverse); // compute PK1
		fStress2D.Rank2ReduceFrom3D(fStress3D);
		S_33 = fStress3D(2,2);

		supp_i = nodalCellSupports(i);
		bVec_i = bVectorArray(i); b_33 = circumferential_B(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			double* fint = fForce(*supp_i++);
			BJ.MultTx(fStress2D.Pointer(), fint, w_i, dMatrixT::kAccumulate);
			*fint += w_i * S_33 * *b_33++;
		}
	}

	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

void FS_SCNIMF_AxiT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
#if __option(extended_errorcheck)
	if (BJ.Rows() != fSD*2) 
		ExceptionT::SizeMismatch("FS_SCNIMF_AxiT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != fSD) 
		ExceptionT::SizeMismatch("FS_SCNIMF_AxiT::bVectorToMatrix","Matrix has bad minorDim");
#endif

	double* Bptr = BJ.Pointer();

	BJ = 0.;
	Bptr[0] = *bVector;
	if (fSD == 2)
	{
		Bptr[5] = *bVector++;
		Bptr[2] = *bVector;
		Bptr[7] = *bVector;
	}
	else
		ExceptionT::GeneralFail("FS_SCNIMF_AxiT::bVectorToMatrix","Two-D please\n");
}

/* return a pointer to a new material list */
MaterialListT* FS_SCNIMF_AxiT::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	int nsd = -1;
	if (name == "large_strain_material_3D")
		nsd = 3;
		
	/* no match */
	if (nsd == -1) return NULL;

	if (size > 0) {
		 /* material support */
		 if (!fFSMatSupport) {
		 	fFSMatSupport = new FSMatSupportT(nsd, 1);      
		 	if (!fFSMatSupport)
		 		ExceptionT::GeneralFail("FS_SCNIMFT::NewMaterialList","Could not instantiate material support\n");

			/* initializations */
			const FEManagerT& fe_man = ElementSupport().FEManager();
			fFSMatSupport->SetFEManager(&fe_man);
			fFSMatSupport->SetDeformationGradient(&fF_list);
			fFSMatSupport->SetDeformationGradient_last(&fF_last_list);
			fFSMatSupport->SetElementCards(&fElementCards);
			fFSMatSupport->SetGroup(Group());
		 }	 

		return new FSSolidMatList3DT(size, *fFSMatSupport);
	} else 
	  return new FSSolidMatList3DT;

	/* no match */
	return NULL;
}

