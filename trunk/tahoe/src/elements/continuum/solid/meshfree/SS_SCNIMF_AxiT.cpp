/* $Id: SS_SCNIMF_AxiT.cpp,v 1.2 2005-04-11 17:39:44 cjkimme Exp $ */

#include "SS_SCNIMF_AxiT.h"

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
#include "LinkedListT.h"
#include "ParameterContainerT.h"

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SSSolidMatT.h"
#include "SSMatSupportT.h"

/* materials lists */
#include "SSSolidMatList3DT.h"

using namespace Tahoe;

const double twoPi = 2.0*acos(-1.0);

/* constructors */
SS_SCNIMF_AxiT::SS_SCNIMF_AxiT(const ElementSupportT& support, const FieldT& field):
	SS_SCNIMFT(support, field)
{
	SetName("ss_mfparticle_axi");
}

SS_SCNIMF_AxiT::SS_SCNIMF_AxiT(const ElementSupportT& support):
	SS_SCNIMFT(support)
{
	SetName("ss_mfparticle_axi");
}

void SS_SCNIMF_AxiT::WriteOutput(void)
{
	const char caller[] = "SS_SCNIMF_AxiT::WriteOutput";

	/* dimensions */
	int nsd  = NumSD();
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

	dSymMatrixT& strain3D = fStrain_list[0];
	dSymMatrixT s_axi(dSymMatrixT::k3D_plane), e_axi(dSymMatrixT::k3D_plane);
	dMatrixT asym(nsd);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	
	/* material outputs */
	dArrayT mat_output(counts[kMaterialOutput]);
	
	/* collect output values */
	dArrayT vec, values_i;
	for (int i = 0; i < non; i++) {
	
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
			values_i[offsets[kMass]] = fCellVolumes[i];

		/* compute smoothed strain */
		if (fOutputFlags[kStrain] || fOutputFlags[kStress] || fOutputFlags[kMaterialOutput]) {
			asym = 0.0;
			double e_33 = 0.;
			int n_supp = nodalCellSupports.MinorDim(i);
			dSymMatrixT& strain3D = fStrain_list[0];
			dArrayT* bVec_i = bVectorArray(i);
			double* b_33 = circumferential_B(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++, bVec_i++) {
				e_33 += u(*supp_i,0) * *b_33++;
				asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
			}
			strain3D[0] = asym[0];
			strain3D[1] = asym[3];
			strain3D[2] = e_33;
			strain3D[5] = (asym[1] + asym[2])/2.;	
			e_axi.Translate(strain3D);
			
			/* last strain tensor */
			if (fStrain_last_list.Length() > 0) 
			{
				/* last displacement */
				const dArray2DT& u = Field()(-1,0);

				/* destination */
				dSymMatrixT& strain3D = fStrain_last_list[0];

				asym = 0.0; e_33 = 0.;
				strain3D = 0.0;
				bVec_i = bVectorArray(i);
				b_33 = circumferential_B(i);
				supp_i = nodalCellSupports(i);
				for (int j = 0; j < n_supp; j++, bVec_i++) { 
					e_33 += u(*supp_i,0) * *b_33++;
					asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
				}
				strain3D[0] = asym[0];
				strain3D[1] = asym[3];
				strain3D[2] = e_33;
				strain3D[5] = (asym[1] + asym[2])/2.;
			}
		}

		/* strain */
		if (fOutputFlags[kStrain]) {
			int index = offsets[kStrain];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = e_axi[j];
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
GlobalT::RelaxCodeT SS_SCNIMF_AxiT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();
	return relax;
}

/***********************************************************************
 * Protected
 ***********************************************************************/
 
/* return number of values for each output variable */
void SS_SCNIMF_AxiT::SetOutputCount(const iArrayT& flags, iArrayT& counts) const
{
	/* inherited */
	SS_SCNIMFT::SetOutputCount(flags, counts);

	/* redimension stress and strain */
	if (flags[kStrain]) counts[kStrain] = 4;
	if (flags[kStress]) counts[kStress] = 4;
}

/* generate labels for output data */
void SS_SCNIMF_AxiT::GenerateOutputLabels(ArrayT<StringT>& labels)
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

/* assemble particle mass matrix into LHS of global equation system */
void SS_SCNIMF_AxiT::AssembleParticleMass(const double rho)
{
	int nsd = NumSD();
  	fForce = 0.0;
  	int* nodes = fNodes.Pointer();
  	double* volume = fCellVolumes.Pointer();
  	for (int i = 0; i < fNodes.Length(); i++) {

    	double* m = fForce(fNodes[i]);

    	for (int j = 0; j < nsd; j++)
      		*m++ = *volume * twoPi * fCellCentroids(i,0);

    	volume++;
  	}

  	fForce *= rho;
  
  	/* assemble all */
  	ElementSupport().AssembleLHS(Group(), fForce, Field().Equations());
}


/* form group contribution to the stiffness matrix */
void SS_SCNIMF_AxiT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	const char caller[] = "SS_SCNIMF_AxiT::LHSDriver";

	int nsd = NumSD();

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
			ExceptionT::GeneralFail(caller, "Cannot get material");

		AssembleParticleMass(fCurrMaterial->Density());
	}

	if (formK) {
		/* hold the smoothed strain */
		dSymMatrixT& strain3D = fStrain_list[0];

		/* displacements */
		const dArray2DT& u = Field()(0,0);

		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);

		if (!fCurrMaterial)
			ExceptionT::GeneralFail(caller, "Cannot get material");

		int nNodes = fNodes.Length();

		/* assembly information */
		int group = Group();
		int ndof = NumDOF();

		fLHS.Dimension(ndof);
		const ElementSupportT& support = ElementSupport();

		const iArray2DT& field_eqnos = Field().Equations();
		iArrayT row_eqnos(ndof); 
		iArrayT col_eqnos(ndof);
		dMatrixT BJ(4, ndof), BK(4, ndof), K_JK;
		dMatrixT BJTCijkl(4, nsd);
		dMatrixT stress3D(3), moduli2D(4);
		dMatrixT asym(nsd);
		double e_33;
		K_JK.Alias(fLHS);
		for (int i = 0; i < nNodes; i++) {	
		
			/* set current element */
			fElementCards.Current(i);
		
			double w_i = fCellVolumes[i]*constK*twoPi*fNodalCoordinates(i,0); // integration weight
			int n_supp = nodalCellSupports.MinorDim(i);
			
			// Compute smoothed strain
			asym = 0.0; e_33 = 0.;
			strain3D = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			double* b_33 = circumferential_B(i);
			for (int j = 0; j < n_supp; j++, bVec_i++) {
				e_33 += u(*supp_i,0) * *b_33++;
				asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
			}
			strain3D[0] = asym[0];
			strain3D[1] = asym[3];
			strain3D[2] = e_33;
			strain3D[5] = (asym[1] + asym[2])/2;

			/* last strain tensor */
			if (fStrain_last_list.Length() > 0) 
			{
				/* last displacement */
				const dArray2DT& u = Field()(-1,0);

				/* destination */
				dSymMatrixT& strain3D = fStrain_last_list[0];

				asym = 0.0; e_33 = 0.;
				strain3D = 0.0;
				bVec_i = bVectorArray(i);
				b_33 = circumferential_B(i);
				supp_i = nodalCellSupports(i);
				for (int j = 0; j < n_supp; j++, bVec_i++) { 
					e_33 += u(*supp_i,0) * *b_33++;
					asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
				}
				strain3D[0] = asym[0];
				strain3D[1] = asym[3];
				strain3D[2] = e_33;
				strain3D[5] = (asym[1] + asym[2])/2.;
			}

			const dMatrixT& cijkl = fCurrMaterial->c_ijkl();
			moduli2D.Rank4ReduceFrom3D(cijkl);

			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			b_33 = circumferential_B(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++) {
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				BJ[3] = *b_33++;
				BJTCijkl.MultAB(moduli2D, BJ, 0);
				col_eqnos.Copy(field_eqnos(*supp_i));

				dArrayT* bVec_j = bVectorArray(i);
				int* supp_j = nodalCellSupports(i);
				double* b_33_k = circumferential_B(i);
				for (int k = 0; k < n_supp; k++) {
					bVectorToMatrix(bVec_j->Pointer(), BK);
					BK[3] = *b_33_k++;
					bVec_j++;

					// K_JK = BT_K x Cijkl x B_J 
					K_JK.MultATB(BK, BJTCijkl);
					K_JK *= w_i;

					/* assemble */
					row_eqnos.Copy(field_eqnos(*supp_j++));			
					support.AssembleLHS(group, fLHS, row_eqnos, col_eqnos);
				}
			}	
		}
	}
}


void SS_SCNIMF_AxiT::RHSDriver(void)
{
	/* function name */
	const char caller[] = "SS_SCNIMF_AxiT::RHSDriver2D";

	/* contribution from natural boundary conditions */
	SCNIMFT::RHSDriver();

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
		ExceptionT::GeneralFail(caller, "cannot get material");

	int nNodes = fNodes.Length();
	int nsd = NumSD();
	
	fForce = 0.0;
	dMatrixT stress3D(3), stress_vec(2);
	dSymMatrixT& strain3D = fStrain_list[0];
	dMatrixT BJ(4, nsd);
	dMatrixT asym(nsd);
	double e_33;

	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++) {
	
		/* set current element */
		fElementCards.Current(i);
	
		double w_i = fCellVolumes[i]*twoPi*fNodalCoordinates(i,0); // integration weight

		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed strain
		asym = 0.0; e_33 = 0.;
		strain3D = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		double* b_33 = circumferential_B(i);
		for (int j = 0; j < n_supp; j++, bVec_i++) {
			e_33 += u(*supp_i,0) * *b_33++;
			asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
		}
		strain3D[0] = asym[0];
		strain3D[1] = asym[3];
		strain3D[2] = e_33;
		strain3D[5] = (asym[1] + asym[2])/2.;

		/* last strain tensor */
		if (fStrain_last_list.Length() > 0) 
		{
			/* last displacement */
			const dArray2DT& u = Field()(-1,0);

			/* destination */
			dSymMatrixT& strain3D = fStrain_last_list[0];

			asym = 0.0; e_33 = 0.;
			strain3D = 0.0;
			bVec_i = bVectorArray(i);
			b_33 = circumferential_B(i);
			supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++, bVec_i++) { 
				e_33 += u(*supp_i,0) * *b_33++;
				asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
			}
			strain3D[0] = asym[0];
			strain3D[1] = asym[3];
			strain3D[2] = e_33;
			strain3D[5] = (asym[1] + asym[2])/2.;
		}

		fCurrMaterial->s_ij().ToMatrix(stress3D);
		stress_vec[0] = stress3D[0];
		stress_vec[1] = stress3D[4];
		stress_vec[2] = stress3D[3];
		stress_vec[3] = stress3D[8];

		supp_i = nodalCellSupports(i);
		bVec_i = bVectorArray(i);
		b_33 = circumferential_B(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			BJ[3] = *b_33++;
			bVec_i++;
			double* fint = fForce(*supp_i++);
			BJ.MultTx(stress_vec.Pointer(), fint, w_i, dMatrixT::kAccumulate);
		}

	}

	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}


void SS_SCNIMF_AxiT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "SS_SCNIMF_AxiT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();
	mat_params.SetName("small_strain_material_3D");

	int num_blocks = all_params.NumLists("ss_scni_axi_element_block");
	for (int i = 0; i < num_blocks; i++) {
	  const ParameterListT& block = all_params.GetList("ss_scni_axi_element_block",i);

	  /* collect material parameters */
	  const ParameterListT& mat_list = block.GetList(mat_params.Name());
	  const ArrayT<ParameterListT>& mat = mat_list.Lists();
	  mat_params.AddList(mat[0]);
	}
}

// XML stuff below

/* accept parameter list */
void SS_SCNIMF_AxiT::TakeParameterList(const ParameterListT& list)
{
	/* we are axisymmetric */
	qIsAxisymmetric = true;

	/* inherited */
	SCNIMFT::TakeParameterList(list);
	
	/* strains */
	int nsd = NumSD();
	fStrain_list.Dimension(1);
	fStrain_list[0].Dimension(3);

	/* casts are safe since class contructs materials list - just one material */
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;
	if (mat->Need_Strain_last()) {
		fStrain_last_list.Dimension(1);
		fStrain_last_list[0].Dimension(3);
	}

}

/* information about subordinate parameter lists */
void SS_SCNIMF_AxiT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SCNIMFT::DefineSubs(sub_list);

	/* element blocks for underlying connectivity -- TEMP */
	sub_list.AddSub("ss_scni_axi_element_block");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SS_SCNIMF_AxiT::NewSub(const StringT& name) const
{
  if  (name == "ss_scni_axi_element_block") {
	  ParameterContainerT* block = new ParameterContainerT(name);
	  block->AddSub("block_ID_list",ParameterListT::Once);
	  block->AddSub("small_strain_material_3D");
	  block->SetSubSource(this);
	  return block;
  }
  else /* inherited */
    return SCNIMFT::NewSub(name);
}

