/* $Id: FS_SCNIMF_AxiT.cpp,v 1.26 2005-03-01 08:26:29 paklein Exp $ */
#include "FS_SCNIMF_AxiT.h"

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
const double twoPi = 2.0*Pi;

/* constructor */
FS_SCNIMF_AxiT::FS_SCNIMF_AxiT(const ElementSupportT& support):
	FS_SCNIMFT(support)
{
	SetName("fd_mfparticle_axi");
}

/* generate labels for output data */
void FS_SCNIMF_AxiT::GenerateOutputLabels(ArrayT<StringT>& labels)
{
	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

	/* initialize */
	labels.Dimension(num_output);

	/* coordinates */
	if (fOutputFlags[kCoordinates]) {
		const char* ref[2] = {"R", "Z"};
		int index = offsets[kCoordinates];
		for (int i = 0; i < 2; i++)
			labels[index++] = ref[i];
	}

	/* displacements */
	if (fOutputFlags[kDisplacement]) {

		/* labels from the field */
		const ArrayT<StringT>& field_labels = Field().Labels();

		int index = offsets[kDisplacement];
		for (int i = 0; i < NumDOF(); i++)
			labels[index++] = field_labels[i];
	}

	/* mass */
	if (fOutputFlags[kMass])
		labels[offsets[kMass]] = "volume";

	/* strain */
	if (fOutputFlags[kStrain]) {
		const char* elabels[4] = {"err", "ezz", "erz", "ett"};
		int index = offsets[kStrain];
		for (int i = 0; i < 4; i++)
			labels[index++] = elabels[i];
	}

	/* stress */
	if (fOutputFlags[kStress]) {
		const char* slabels[4] = {"srr", "szz", "srz", "stt"};
		int index = offsets[kStress];
		for (int i = 0; i < 4; i++)
			labels[index++] = slabels[i];
	}

	/* material output labels */
	if (fOutputFlags[kMaterialOutput])
	{
		ArrayT<StringT> mat_labels;
		(*fMaterialList)[0]->OutputLabels(mat_labels);	
		
		int index = offsets[kMaterialOutput];
		for (int i = 0; i < mat_labels.Length(); i++)
			labels[index++] = mat_labels[i];
	}
}

void FS_SCNIMF_AxiT::WriteOutput(void)
{
	const char caller[] = "FS_SCNIMF_AxiT::WriteOutput";

	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	int nstrs = 4;
	int non = fNodes.Length();

	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

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
	SolidMaterialT* curr_material = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!curr_material) ExceptionT::GeneralFail(caller, "cannot get material");	

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	dMatrixT& F3D = fF_list[0];
	dMatrixT F2D(2);
	dMatrixT BJ(nsd*nsd, nsd);
	dSymMatrixT E3D(3);
	dSymMatrixT s_axi(dSymMatrixT::k3D_plane), E_axi(dSymMatrixT::k3D_plane);

	/* displacements */
	const dArray2DT& u = Field()(0,0);
	const dArray2DT& u_last = Field()(-1,0);

	/* material outputs */
	dArrayT mat_output(counts[kMaterialOutput]);

	/* collect displacements */
	dArrayT vec, values_i;
	for (int i = 0; i < non; i++) 
	{
		/* set current element */
		fElementCards.Current(i);

		/* global ID */
		int tag_i = fNodes[i];

		/* values for particle i */
		n_values.RowAlias(i, values_i);

		/* coordinates */
		if (fOutputFlags[kCoordinates]) {
			vec.Alias(nsd, values_i.Pointer(offsets[kCoordinates]));
			coords.RowCopy(tag_i, vec);
		}

		/* displacement */
		if (fOutputFlags[kDisplacement]) {
			vec.Alias(ndof, values_i.Pointer(offsets[kDisplacement]));

			const int* nodal_supp = fNodalSupports(i);
			const double* phi_i = fNodalPhi(i);
			for (int j = 0; j < fNodalPhi.MinorDim(i); j++)
				vec.AddScaled(*phi_i++, u(*nodal_supp++));
		}

		/* mass */		
		if (fOutputFlags[kMass])
			values_i[offsets[kMass]] = fCellVolumes[i] * 2. * Pi * fCellCentroids(i,0);

		/* compute smoothed deformation gradient */
		if (fOutputFlags[kStrain] || fOutputFlags[kStress] || fOutputFlags[kMaterialOutput])
		{
			/* support size */
			int n_supp = nodalCellSupports.MinorDim(i);

			/* destination */
			dMatrixT& F3D = fF_last_list[0];
			
			/* compute smoothed deformation gradient */
			double F_33 = 0.0;
			F2D = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			const double* b_33 = circumferential_B(i);
			const int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				F_33 += u(*supp_i,0) * *b_33++; 
				F2D.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++;
			}
			F3D.Rank2ExpandFrom2D(F2D);
			F3D.Last() = F_33;
			F3D.PlusIdentity(); // convert to F	

			/* Green-Lagrangian strain */
			curr_material->Strain(E3D);
			E_axi.Translate(E3D);

			/* last deformation gradient */
			if (fF_last_list.Length() > 0) 
			{
				/* destination */
				dMatrixT& F3D = fF_last_list[0];

				double F_33 = 0.0;
				F2D = 0.0;
				dArrayT* bVec_i = bVectorArray(i);
				double* b_33 = circumferential_B(i);
				int* supp_i = nodalCellSupports(i);
				for (int j = 0; j < n_supp; j++) { 
					F_33 += u_last(*supp_i,0) * *b_33++; 
					F2D.Outer(u_last(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
					bVec_i++;
				}
				F3D.Rank2ExpandFrom2D(F2D);
				F3D.Last() = F_33;
				F3D.PlusIdentity(); // convert to F
			}
		}

		/* strain */
		if (fOutputFlags[kStrain]) {
			int index = offsets[kStrain];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = E_axi[j];
		}

		/* stress */
		if (fOutputFlags[kStress])
		{
			s_axi.Translate(curr_material->s_ij());
			int index = offsets[kStress];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = s_axi[j];		
		}

		/* material output parameters */
		if (fOutputFlags[kMaterialOutput])
		{
			/* update stress */
			if (!fOutputFlags[kStress]) curr_material->s_ij();

			/* compute material output */
			curr_material->ComputeOutput(mat_output);
			int index = offsets[kMaterialOutput];
			for (int j = 0; j < mat_output.Length(); j++)
				values_i[index++] = mat_output[j];		
		}
	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT FS_SCNIMF_AxiT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();
	return relax;
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

/* return number of values for each output variable */
void FS_SCNIMF_AxiT::SetOutputCount(const iArrayT& flags, iArrayT& counts) const
{
	/* inherited */
	FS_SCNIMFT::SetOutputCount(flags, counts);

	/* redimension stress and strain */
	if (flags[kStrain]) counts[kStrain] = 4;
	if (flags[kStress]) counts[kStress] = 4;
}

void FS_SCNIMF_AxiT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "FS_SCNIMF_AxiT::CollectMaterialInfo";

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
	const char caller[] = "FS_SCNIMF_AxiT::LHSDriver";

#pragma unused(sys_type)
	/* time integration parameters */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formM = fIntegrator->FormM(constM);
	int nsd = NumSD();

	/* quick exit */
	if ((formM == 0 && formK == 0) ||
	    (fabs(constM) < kSmall &&
	     fabs(constK) < kSmall)) return;

	/* assemble particle mass */
	if (formM) {
	
		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* curr_material = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!curr_material) ExceptionT::GeneralFail(caller, "cannot get material");

		AssembleParticleMass(curr_material->Density());
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
		SolidMaterialT* curr_material = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!curr_material) ExceptionT::GeneralFail(caller, "cannot get material");

		int nNodes = fNodes.Length();

		/* assembly information */
		int group = Group();
		int ndof = NumDOF();
		fLHS.Dimension(ndof);
		const ElementSupportT& support = ElementSupport();

		const iArray2DT& field_eqnos = Field().Equations();
		iArrayT row_eqnos(ndof); 
		iArrayT col_eqnos(ndof);
		dMatrixT BJ(4, ndof), BK(4, ndof), fStress3D(3), fStress2D(2);
		dMatrixT Tijkl(nsd*nsd), BJTCijkl(nsd, nsd*nsd), K_JK, mod2D(3), Finverse(3);
		dMatrixT Cijklsigma(nsd*nsd);
		K_JK.Alias(fLHS);
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
				F_33 += u(*supp_i,0) * *b_33++;
				F2D.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++; 
			}
			F2D.PlusIdentity(); // convert to F
			F3D.Rank2ExpandFrom2D(F2D);
			F_33 += 1.;
			F3D(2,2) = F_33;
			J = F3D.Det();
			if (J <= 0.0)
				ExceptionT::BadJacobianDet(caller);

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
					F_33 += u(*supp_i,0) * *b_33++;
					F2D.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
					bVec_i++; 
				}
				F2D.PlusIdentity(); // convert to F
				F3D.Rank2ExpandFrom2D(F2D);
				F_33 += 1.;
				F3D(2,2) = F_33;
			}
			
			curr_material->s_ij().ToMatrix(fStress3D);
			fStress2D.Rank2ReduceFrom3D(fStress3D);

			// stress stiffness
			// T_11 = T_22 = FTs_11 T_13 = T_24 = FTs_12
			// T_31 = T_42 = FTs_21 T_33 = T_44 = FTs_22
			Tijkl.Expand(fStress2D, 2, dMatrixT::kOverwrite);
			// Tijkl = 0.; used to debug material stiffness
			
			// material stiffness
			const dMatrixT& moduli = curr_material->C_IJKL();
			mod2D.Rank4ReduceFrom3D(moduli);
			Tijkl += TransformModuli(moduli, F2D, Cijklsigma);

			// c_theta is the last row of mod2D with the last entry dropped
			c_theta_theta = moduli[14] * F_33 * F_33;
			
			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			b_33 = circumferential_B(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++)
			{
				bVectorToMatrix(bVec_i->Pointer(), BJ);

				BJTCijkl.MultATB(BJ, Tijkl, dMatrixT::kWhole); 

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
					K_JK[0] += *b_33 * c_theta_theta * *b_33_k;

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
	int nsd = NumSD();
  fForce = 0.0;
  int* nodes = fNodes.Pointer();
  double* volume = fCellVolumes.Pointer();
  for (int i = 0; i < fNodes.Length(); i++) {

    double* m = fForce(fNodes[i]);

    for (int j = 0; j < nsd; j++)
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
	int nsd = NumSD();

	if (formMa) 
	{
		if (Field().Order() < 2)
			ExceptionT::GeneralFail(caller,"Field's Order does not have accelerations\n");
	
		fLHS = 0.0; // fLHS.Length() = nNodes * nsd;
		const dArray2DT& a = Field()(0,2); // accelerations
		double* ma = fLHS.Pointer();
		const double* acc;

		int* nodes = fNodes.Pointer();
		double* volume = fCellVolumes.Pointer();
		for (int i = 0; i < nNodes; i++)
		{
			acc = a(*nodes++);
			for (int j = 0; j < nsd; j++)
				*ma++ = *volume * 2. * Pi * fCellCentroids(i,0) * *acc++;
			volume++;
		}
		fLHS *= fCurrMaterial->Density();
	}

	fForce = 0.0;
	dMatrixT& F3D = fF_list[0];
	dMatrixT fStress3D(3), fStress2D(2), fCauchy(3), Finverse(3);
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
			F_33 += u(*supp_i,0) * *b_33++;
			F2D.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
			bVec_i++; 
		}
		F3D.Rank2ExpandFrom2D(F2D);
		F3D.Last() = F_33;
		F3D.PlusIdentity(); // convert to F
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
				F_33 += u(*supp_i,0) * *b_33++;
				F2D.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++; 
			}
			F3D.Rank2ExpandFrom2D(F2D);
			F3D.Last() = F_33;
			F3D.PlusIdentity(); // convert to F
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
			double* fint = fForce(*supp_i++);
			*fint += w_i * S_33 * *b_33++;
			fStress2D.Multx(bVec_i->Pointer(), fint, w_i, dMatrixT::kAccumulate);
			bVec_i++;
		}
	}

	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

