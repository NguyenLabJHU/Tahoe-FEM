/* $Id: SS_SCNIMFT.cpp,v 1.17 2005-03-01 08:26:29 paklein Exp $ */

#include "SS_SCNIMFT.h"

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

//#define  VERIFY_INTEGRATION_CONSTRAINT

/* materials lists */
#include "SSSolidMatList1DT.h"
#include "SSSolidMatList2DT.h"
#include "SSSolidMatList3DT.h"

using namespace Tahoe;

/* constructors */
SS_SCNIMFT::SS_SCNIMFT(const ElementSupportT& support, const FieldT& field):
	SCNIMFT(support, field),
	fSSMatSupport(NULL)
{
	SetName("ss_mfparticle");
}

SS_SCNIMFT::SS_SCNIMFT(const ElementSupportT& support):
	SCNIMFT(support),
	fSSMatSupport(NULL)
{
	SetName("ss_mfparticle");
}

/* destructor */
SS_SCNIMFT::~SS_SCNIMFT(void)
{
	delete fSSMatSupport;
}

void SS_SCNIMFT::WriteOutput(void)
{
	const char caller[] = "SS_SCNIMFT::WriteOutput";

//TEMP
if (NumDOF() != 2) ExceptionT::GeneralFail(caller, "2D only: %d", NumDOF());

	/* dimensions */
	int nsd  = NumSD();
	int ndof = NumDOF();
	int nstrs = dSymMatrixT::NumValues(nsd);
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

	ArrayT<dSymMatrixT> strainList(1);
	strainList[0].Dimension(nsd);
	dSymMatrixT& strain = strainList[0];
	dMatrixT asym(nsd);

	/* displacements */
	const dArray2DT& u = Field()(0,0);
	
	/* material outputs */
	dArrayT mat_output(counts[kMaterialOutput]);
	
	/* collect output values */
	dArrayT vec, values_i;
	for (int i = 0; i < non; i++) {
	
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
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < nodalCellSupports.MinorDim(i); j++, bVec_i++) 
				asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
			strain.Symmetrize(asym);
			fSSMatSupport->SetLinearStrain(&strainList);
		}

		/* strain */
		if (fOutputFlags[kStrain]) {
			int index = offsets[kStrain];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = strain[j];
		}

		/* stress */
		if (fOutputFlags[kStress])
		{
			const dSymMatrixT& stress = curr_material->s_ij();
			int index = offsets[kStress];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = stress[j];		
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
GlobalT::RelaxCodeT SS_SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();
	return relax;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void SS_SCNIMFT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)
	const char caller[] = "SS_SCNIMFT::LHSDriver";

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
		ArrayT<dSymMatrixT> strainList(1);
		strainList[0].Dimension(nsd);
		dSymMatrixT& strain = strainList[0];

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
		dMatrixT BJ(nsd == 2 ? 3 : 6, ndof), BK(nsd == 2 ? 3 : 6, ndof), K_JK;
		dMatrixT BJTCijkl(nsd == 2 ? 3 : 6, nsd);
		dMatrixT asym(nsd);
		K_JK.Alias(fLHS);

		LinkedListT<dArrayT> bVectors_j;
		LinkedListT<int> nodeSupport_j;

		for (int i = 0; i < nNodes; i++) {	
			double w_i = fCellVolumes[i]*constK; // integration weight
			int n_supp = nodalCellSupports.MinorDim(i);
			
			// Compute smoothed strain 
			asym = 0.0;

			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++, bVec_i++) 
				asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
			strain.Symmetrize(asym);
			
			fSSMatSupport->SetLinearStrain(&strainList);
			const dMatrixT& cijkl = fCurrMaterial->c_ijkl();

			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);

			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++) {
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				BJTCijkl.MultAB(cijkl, BJ, 0);
				col_eqnos.Copy(field_eqnos(*supp_i));

				dArrayT* bVec_j = bVectorArray(i);
				int* supp_j = nodalCellSupports(i);
				for (int k = 0; k < n_supp; k++) {
					bVectorToMatrix(bVec_j->Pointer(), BK);
					bVec_j++;

					if (nsd == 2) {
						BK(2,0) *= 2.; // I either have to do this here or on the RHS 
						BK(2,1) *= 2.; // It's the C_1212 e_12 + C_1221 e_21 factor of 2
					} else {
						BK(3,1) *= 2.;
						BK(3,2) *= 2.;
						BK(4,0) *= 2.;
						BK(4,2) *= 2.;
						BK(5,0) *= 2.;
						BK(5,1) *= 2.;
					}

					// K_JK = BT_K x Cijkl x B_J 
					K_JK.MultATB(BK,BJTCijkl, 0);
					K_JK *= w_i;

					/* assemble */
					row_eqnos.Copy(field_eqnos(*supp_j++));			
					support.AssembleLHS(group, fLHS, row_eqnos, col_eqnos);
				}
			}	
		}
	}
}


void SS_SCNIMFT::RHSDriver(void)
{
	/* function name */
	const char caller[] = "SS_SCNIMFT::RHSDriver2D";

//TEMP
if (fMaterialList->HasHistoryMaterials())
	ExceptionT::GeneralFail(caller, "SS_SCNIMFT does not support materials with history");

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

	if (formMa) {

		if (Field().Order() < 2)
			ExceptionT::GeneralFail(caller, "Field's Order does not have accelerations");

		fLHS = 0.0; // fLHS.Length() = nNodes * nsd;
		const dArray2DT& a = Field()(0,2); // accelerations
		double* ma = fLHS.Pointer();
		const double* acc;
		int* nodes = fNodes.Pointer();
		double* volume = fCellVolumes.Pointer();

		for (int i = 0; i < nNodes; i++) {
			acc = a(*nodes++);
			for (int j = 0; j < nsd; j++)
				*ma++ = *volume * *acc++;
			volume++;
		}
		fLHS *= fCurrMaterial->Density();
	}

	fForce = 0.0;
	ArrayT<dSymMatrixT> strainList(1);
	strainList[0].Dimension(nsd);
	dSymMatrixT& strain = strainList[0];
	dMatrixT BJ(nsd == 2 ? 3 : 6, nsd);
	dMatrixT asym(nsd);

#ifdef VERIFY_INTEGRATION_CONSTRAINT
	// TEMP -- verify that \sum_L \mathbf{B}_{Ii} = 0
	dArray2DT test_sum(nNodes,nsd);
	test_sum = 0.;
#endif

	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++) {
		double w_i = fCellVolumes[i]; // integration weight

		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed strain
		asym = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < n_supp; j++, bVec_i++) {
#ifdef VERIFY_INTEGRATION_CONSTRAINT
			test_sum.SetRow(*supp_i, bVec_i->Pointer());
			test_sum.ScaleRow(*supp_i, w_i);
#endif
			asym.Outer(bVec_i->Pointer(), u(*supp_i++), 1.0, dMatrixT::kAccumulate);
		}	
		strain.Symmetrize(asym);	

		fSSMatSupport->SetLinearStrain(&strainList);
		const double* stress = fCurrMaterial->s_ij().Pointer();

		supp_i = nodalCellSupports(i);
		bVec_i = bVectorArray(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			double* fint = fForce(*supp_i++);
			BJ.MultTx(stress, fint, w_i, dMatrixT::kAccumulate);
		}

	}

#ifdef VERIFY_INTEGRATION_CONSTRAINT
	static int firstTime = 0;
	if (!firstTime) {
		firstTime++;
		for (int i = 0; i < nNodes; i++) { 
			cout << " i = " 	<< i << " ts "; 
			for (int j = 0; j < nsd; j++) 
				cout << test_sum(i,j) << " ";
			cout << "\n";
		}
	}	
#endif

	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}



void SS_SCNIMFT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
	int nsd = NumSD();
	
#if __option(extended_errorcheck)
	if (BJ.Rows() != nsd*(nsd+1)/2) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != nsd) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad minorDim");
#endif

	double* Bptr = BJ.Pointer();

	BJ = 0.;
	Bptr[0] = *bVector;
	if (nsd == 2) {
		Bptr[5] = .5 * *bVector++;
		Bptr[4] = *bVector;
		Bptr[2] = .5 * *bVector;
	} else { // nsd == 3
		Bptr[11] = Bptr[16] = 0.5 * *bVector++;
		Bptr[7] = *bVector;
		Bptr[5] = Bptr[15] = 0.5 * *bVector++;
		Bptr[14] = *bVector;
		Bptr[4] = Bptr[9] = 0.5 * *bVector;
	}
}

void SS_SCNIMFT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "SS_SCNIMFT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();

	int num_blocks = all_params.NumLists("ss_connectivity_element_block");
	for (int i = 0; i < num_blocks; i++) {
	  const ParameterListT& block = all_params.GetList("ss_connectivity_element_block",i);

	  if (i == 0) {
	    const ParameterListT& mat_list_params = block.GetListChoice(*this, "small_strain_material_choice");
	    mat_params.SetName(mat_list_params.Name());
	  }

	  /* collect material parameters */
	  const ParameterListT& mat_list = block.GetList(mat_params.Name());
	  const ArrayT<ParameterListT>& mat = mat_list.Lists();
	  mat_params.AddList(mat[0]);
	}
}

/* return a pointer to a new material list */
MaterialListT* SS_SCNIMFT::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	int nsd = -1;

	if (name == "small_strain_material_2D")
		nsd = 2;
	else if (name == "small_strain_material_3D")
		nsd = 3;

	/* no match */
	if (nsd == -1) return NULL;

	if (size > 0) {
		/* material support */
		 if (!fSSMatSupport) {
		 	fSSMatSupport = new SSMatSupportT(nsd, 1);      
		 	if (!fSSMatSupport)
		 		ExceptionT::GeneralFail("SS_SCNIMFT::NewMaterialList","Could not instantiate material support\n");

			fSSMatSupport->SetFEManager(&ElementSupport().FEManager());
		 }

		if (nsd == 2)
			return new SSSolidMatList2DT(size, *fSSMatSupport);
		else if (nsd == 3)
			return new SSSolidMatList3DT(size, *fSSMatSupport);
	} else {
		if (nsd == 2)
			return new SSSolidMatList2DT;
		else if (nsd == 3)
			return new SSSolidMatList3DT;
	}

	/* no match */
	return NULL;
}

// XML stuff below

/* describe the parameters needed by the interface */
void SS_SCNIMFT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SCNIMFT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void SS_SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SCNIMFT::DefineSubs(sub_list);

	/* element blocks for underlying connectivity -- TEMP */
	sub_list.AddSub("ss_connectivity_element_block");
}

/* return the description of the given inline subordinate parameter list */
void SS_SCNIMFT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "small_strain_material_choice") {
		order = ParameterListT::Choice;

		/* list of choices */
		sub_lists.AddSub("small_strain_material_2D");
		sub_lists.AddSub("small_strain_material_3D");
	} else /* inherited */
		SCNIMFT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SS_SCNIMFT::NewSub(const StringT& name) const
{
  if  (name == "ss_connectivity_element_block") {
	  ParameterContainerT* block = new ParameterContainerT(name);
	  block->AddSub("block_ID_list",ParameterListT::Once);
	  block->AddSub("small_strain_material_choice", ParameterListT::Once, true);
	  block->SetSubSource(this);
	  return block;
  }
  else /* inherited */
    return SCNIMFT::NewSub(name);
}

