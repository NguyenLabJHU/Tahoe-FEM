/* $Id: FS_SCNIMFT.cpp,v 1.24 2005-03-01 17:46:50 cjkimme Exp $ */
#include "FS_SCNIMFT.h"

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
#include "ParameterContainerT.h"

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"
#include "TensorTransformT.h"

/* materials lists */
#include "FSSolidMatList1DT.h"
#include "FSSolidMatList2DT.h"
#include "FSSolidMatList3DT.h"

using namespace Tahoe;

/* constructors */
FS_SCNIMFT::FS_SCNIMFT(const ElementSupportT& support, const FieldT& field):
	SCNIMFT(support, field),
	fFSMatSupport(NULL)
{
	SetName("fd_mfparticle");
}

FS_SCNIMFT::FS_SCNIMFT(const ElementSupportT& support):
	SCNIMFT(support),
	fFSMatSupport(NULL)
{
	SetName("fd_mfparticle");
}

/* destructor */
FS_SCNIMFT::~FS_SCNIMFT(void)
{
	delete fFSMatSupport;
}

void FS_SCNIMFT::WriteOutput(void)
{
	const char caller[] = "FS_SCNIMFT::WriteOutput";

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
	FSSolidMatT* curr_material = TB_DYNAMIC_CAST(FSSolidMatT*,mat);
	if (!curr_material) ExceptionT::GeneralFail(caller , "cannot get material");

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	dMatrixT BJ(nsd*nsd, nsd), Finverse(nsd);
	dSymMatrixT E(nsd);
	double J;
	
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
			values_i[offsets[kMass]] = fCellVolumes[i];
		
		/* compute smoothed deformation gradient */
		if (fOutputFlags[kStrain] || fOutputFlags[kStress] || fOutputFlags[kMaterialOutput]) {		
			dMatrixT& Fdef = fF_list[0];
			Fdef = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < nodalCellSupports.MinorDim(i); j++) {
				Fdef.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++;
			}			
			Fdef.PlusIdentity();
			J = Fdef.Det();
			
			/* Green-Lagrangian strain */
			curr_material->Strain(E);
	
			/* compute smoothed deformation gradient from the previous time increment */
			if (fF_last_list.Length() > 0) {
				dMatrixT& Fdef_last = fF_last_list[0];
				Fdef_last = 0.0;
				dArrayT* bVec_i = bVectorArray(i);
				int* supp_i = nodalCellSupports(i);
				for (int j = 0; j < nodalCellSupports.MinorDim(i); j++) { 
					Fdef_last.Outer(u_last(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
					bVec_i++;
				}
				Fdef_last.PlusIdentity();
			}
		}

		/* strain */
		if (fOutputFlags[kStrain]) {
			int index = offsets[kStrain];
			for (int j = 0; j < nstrs; j++)
				values_i[index++] = E[j];
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
GlobalT::RelaxCodeT FS_SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	return relax;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void FS_SCNIMFT::LHSDriver(GlobalT::SystemTypeT sys_type)
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
		FSSolidMatT* curr_material = TB_DYNAMIC_CAST(FSSolidMatT*,mat);
		if (!curr_material)
			ExceptionT::GeneralFail("FS_SCNIMFT::LHSDriver","Cannot get material\n");
	
		AssembleParticleMass(curr_material->Density());
	}

	if (formK) {
		/* hold the smoothed strain */
		dMatrixT& Fdef = fF_list[0];
		
		/* displacements */
		const dArray2DT& u = Field()(0,0);
		const dArray2DT& u_last = Field()(-1,0);
	
		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		FSSolidMatT* curr_material = TB_DYNAMIC_CAST(FSSolidMatT*,mat);
		if (!curr_material)
			ExceptionT::GeneralFail("FS_SCNIMFT::LHSDriver","Cannot get material\n");

		int nNodes = fNodes.Length();

		/* assembly information */
		int group = Group();
		int ndof = NumDOF();
		int nsd  = NumSD();
		fLHS.Dimension(ndof);
		const ElementSupportT& support = ElementSupport();
		
		const iArray2DT& field_eqnos = Field().Equations();
		iArrayT row_eqnos(ndof); 
		iArrayT col_eqnos(ndof);
		dMatrixT ctmpin(nsd), ctmpout;
		dMatrixT BJ(nsd*nsd, ndof), BK(nsd*nsd, ndof), FTs(nsd), fStress(nsd), fCauchy(nsd);
		dMatrixT Tijkl(nsd*nsd), BJTCijkl(nsd, nsd*nsd), K_JK, Finverse(nsd);
		dMatrixT Cijklsigma(nsd*nsd);
		K_JK.Alias(fLHS);
		LinkedListT<dArrayT> bVectors_j;
		LinkedListT<int> nodeSupport_j;
		for (int i = 0; i < nNodes; i++) 
		{	
			/* set current element */
			fElementCards.Current(i);

			double w_i = fCellVolumes[i]*constK; // integration weights
			
			int n_supp = nodalCellSupports.MinorDim(i);
			
			// Compute smoothed deformation gradient
			Fdef = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				Fdef.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++;
			}		
			Fdef.PlusIdentity(); // convert to F

			/* need deformation gradient from the previous time increment */
			if (fF_last_list.Length() > 0) {
				dMatrixT& Fdef_last = fF_last_list[0];	
				Fdef_last = 0.0;
				dArrayT* bVec_i = bVectorArray(i);
				int* supp_i = nodalCellSupports(i);
				for (int j = 0; j < n_supp; j++) { 
					Fdef_last.Outer(u_last(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
					bVec_i++;
				}
				Fdef_last.PlusIdentity(); // convert to F
			}
		
			const dMatrixT& cijkl = curr_material->C_IJKL();
			curr_material->s_ij().ToMatrix(fCauchy);
			Finverse.Inverse(Fdef);
			Finverse *= Fdef.Det(); // compute J F^-1
			fStress.MultABT(fCauchy, Finverse); // compute PK1
		
			// FTs = F^T sigma
			FTs.MultATB(Fdef, fStress);

			curr_material->S_IJ().ToMatrix(fStress);
			Tijkl.Expand(fStress, nsd, dMatrixT::kOverwrite);
			// Tijkl = 0.; //used to debug material stiffness
			
			// material stiffness
			const dMatrixT& moduli = curr_material->C_IJKL();
			Tijkl += TransformModuli(moduli, Fdef, Cijklsigma);
			
			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++) {
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				
				/* simultanesouly compute material and stress stiffnesses */
				BJTCijkl.MultATB(BJ, Tijkl, dMatrixT::kWhole); // accumulate stress stiffness
				
				col_eqnos.Copy(field_eqnos(*supp_i));
				
				dArrayT* bVec_j = bVectorArray(i);
				int* supp_j = nodalCellSupports(i);
				for (int k = 0; k < n_supp; k++) {
					bVectorToMatrix(bVec_j->Pointer(), BK);
					bVec_j++;
					
					// K_JK = BT_J x Cijkl x B_K 
					K_JK.MultAB(BJTCijkl, BK, dMatrixT::kWhole);
					
					K_JK *= w_i*constK;
					
					/* assemble */
					row_eqnos.Copy(field_eqnos(*supp_j++));
					support.AssembleLHS(group, fLHS, col_eqnos, row_eqnos);
				}
			}	
		}
	}
}

void FS_SCNIMFT::RHSDriver(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver2D";
	
	/* contribution from natural boundary conditions */
	SCNIMFT::RHSDriver();

	/* time integration parameters */
	double constMa = 0.0;
	double constKd = 0.0;
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	FSSolidMatT* curr_material = TB_DYNAMIC_CAST(FSSolidMatT*,mat);
	if (!curr_material)
		ExceptionT::GeneralFail("FS_SCNIMFT::RHSDriver","Cannot get material\n");
	
	int nNodes = fNodes.Length();
	int nsd = NumSD();

	if (formMa) {
		if (Field().Order() < 2)
			ExceptionT::GeneralFail(caller,"Field's Order does not have accelerations\n");
	
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
		fLHS *= curr_material->Density();
	}

	fForce = 0.0;
	dMatrixT& Fdef = fF_list[0];
	dMatrixT BJ(nsd*nsd, nsd), fCauchy(nsd), Finverse(nsd), fStress(nsd);
	double J;
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	const dArray2DT& u_last = Field()(-1,0);
	for (int i = 0; i < nNodes; i++) 
	{
		/* set current element */
		fElementCards.Current(i);
	
		double w_i = fCellVolumes[i]; // integration weight
		
		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed deformation gradient
		Fdef = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < n_supp; j++) { 
			Fdef.Outer(u(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
			bVec_i++;
		}
		Fdef.PlusIdentity(); // convert to F 	
		J = Fdef.Det();
		if (J <= 0.0)
			ExceptionT::BadJacobianDet("FS_SCNIMFT::FormKd");

		/* compute deformation gradient from the previous time increment */
		if (fF_last_list.Length() > 0) {
			dMatrixT& Fdef_last = fF_last_list[0];
			Fdef_last = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				Fdef_last.Outer(u_last(*supp_i++), bVec_i->Pointer(), 1.0, dMatrixT::kAccumulate);
				bVec_i++;
			}
			Fdef_last.PlusIdentity(); // convert to F 			
		}
		
		curr_material->s_ij().ToMatrix(fCauchy);
		Finverse.Inverse(Fdef);
		Finverse *= J; // compute J F^-1
		fStress.MultABT(fCauchy, Finverse); // compute PK1
				
		supp_i = nodalCellSupports(i);
		bVec_i = bVectorArray(i);
		for (int j = 0; j < n_supp; j++) { 
			double* fint = fForce(*supp_i++);
			fStress.Multx(bVec_i->Pointer(), fint, w_i, dMatrixT::kAccumulate);
			bVec_i++;
		}
	}
	
	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

void FS_SCNIMFT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
	int nsd = NumSD();
#if __option(extended_errorcheck)
	if (BJ.Rows() != nsd*nsd) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != nsd) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad minorDim");
#endif

	double* Bptr = BJ.Pointer();
	BJ = 0.;
	Bptr[0] = *bVector;
	if (nsd == 2) {
		Bptr[5] = *bVector++;
		Bptr[2] = *bVector;
		Bptr[7] = *bVector;
	} else { // nsd == 3
		Bptr[10] = Bptr[19] = *bVector++;
		Bptr[3] = *bVector;
		Bptr[13] = Bptr[23] = *bVector++;
		Bptr[6] = *bVector;
		Bptr[16] = Bptr[26] = *bVector;
	}
}

dMatrixT& FS_SCNIMFT::TransformModuli(const dMatrixT& moduli, const dMatrixT& F, dMatrixT& Csig) {
		
		Csig = 0.;
		
		int nsd = F.Rows();
		dMatrixT mtmp(nsd*nsd);
		
		mtmp = 0.;
		
		if (nsd == 2) {
			// use brute force, low-level until I find optimal way
			// numbering scheme for CIJKL is standard C_AB, I need matrix T_AB with 
			// indexing 1 <-> 11, 2 <-> 21, 3 <-> 12, 4 <-> 22 in 2D
			
			// Column 1
			mtmp[0] = moduli[0]; // C_11
			mtmp[1] = mtmp[2] = moduli[6]; // C_13
			mtmp[3] = moduli[3]; // C_12
			
			// Column 2
			mtmp[4] = moduli[6]; // C_13
			mtmp[5] = mtmp[6] = moduli[8]; // C_33
			mtmp[7] = moduli[7]; // C_23
			
			// Column 3
			mtmp[8] = moduli[6]; // C_13
			mtmp[9] = mtmp[10] = moduli[8]; // C_33
			mtmp[11] = moduli[7]; // C_23
			
			// Column 4
			mtmp[12] = moduli[3]; // C_12
			mtmp[13] = mtmp[14] = moduli[7]; // C_23
			mtmp.Last() = moduli[4]; // C_22
			
		} else {
			// use brute force, low-level until I find optimal way
			// numbering scheme for CIJKL is standard C_AB, I need matrix T_AB with DIFFERENT
			// indexing 1 <-> 11, 2 <-> 21, 3 <-> 31, 4 <-> 12 5 <-> 22 6 <-> 32 
			//          7 <-> 13  8 <-> 23  9 <-> 33
			// in 3D
			
			// Column 1
			int offs1 = 0; 
			int offs2 = 0;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_11
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_61
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_51
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_21
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_41
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_31
			
			// Column 2
			offs1 = 9; offs2 = 30;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_16
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_66
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_56
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_26
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_46
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_36
			
			// Column 3
			offs1 = 18; offs2 = 24;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_15
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_65
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_55
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_25
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_45
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_35
			
			// Column 4
			offs1 = 27; offs2 = 30;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_16
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_66
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_56
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_26
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_46
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_36
			
			// Column 5
			offs1 = 36; offs2 = 6;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_12
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_62
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_52
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_22
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_42
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_32
			
			// Column 6
			offs1 = 45; offs2 = 18;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_14
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_64
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_54
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_24
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_44
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_34
			
			// Column 7 
			offs1 = 54; offs2 = 24;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_15
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_65
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_55
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_25
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_45
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_35
			
			// Column 8
			offs1 = 63; offs2 = 18;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_14
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_64
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_54
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_24
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_44
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_34
			
			// Column 9
			offs1 = 72; offs2 = 12;
			mtmp[offs1 + 0] = moduli[offs2 + 0]; // C_13
			mtmp[offs1 + 1] = mtmp[offs1 + 3] = moduli[offs2 + 5]; // C_63
			mtmp[offs1 + 2] = mtmp[offs1 + 6] = moduli[offs2 + 4]; // C_53
			mtmp[offs1 + 4] = moduli[offs2 + 1]; // C_23
			mtmp[offs1 + 5] = mtmp[offs1 + 7] = moduli[offs2 + 3]; // C_43
			mtmp[offs1 + 8] = moduli[offs2 + 2]; // C_33
		}
				
		for (int i = 0; i < nsd*nsd; i++) {
			dMatrixT col_in(nsd, nsd, mtmp.Pointer(i*nsd*nsd));
			dMatrixT col_out(nsd, nsd, Csig.Pointer(i*nsd*nsd));
			col_out.MultAB(F, col_in);
		}
		
		Csig.Transpose(Csig);
		mtmp = Csig;
		
		for (int i = 0; i < nsd*nsd; i++) {
			dMatrixT col_in(nsd, nsd, mtmp.Pointer(i*nsd*nsd));
			dMatrixT col_out(nsd, nsd, Csig.Pointer(i*nsd*nsd));
			col_out.MultAB(F, col_in);
		}
		
		return Csig;
		
}

void FS_SCNIMFT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "FS_SCNIMFT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();

    int num_blocks = all_params.NumLists("fd_connectivity_element_block");
	for (int i = 0; i < num_blocks; i++) {
	  
	  const ParameterListT& block = all_params.GetList("fd_connectivity_element_block",i);

	  if (i == 0) {
	    const ParameterListT& mat_list_params = block.GetListChoice(*this, "large_strain_material_choice");
	    mat_params.SetName(mat_list_params.Name());
	  }

	  /* collect material parameters */
	  const ParameterListT& mat_list = block.GetList(mat_params.Name());
	  const ArrayT<ParameterListT>& mat = mat_list.Lists();
	  mat_params.AddList(mat[0]);
	}
}

/* return a pointer to a new material list */
MaterialListT* FS_SCNIMFT::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	int nsd = -1;
	if (name == "large_strain_material_2D")
		nsd = 2;
	else if (name == "large_strain_material_3D")
		nsd = 3;
	
	/* no match */
	if (nsd == -1) return NULL;
	
	if (qIsAxisymmetric && nsd != 3) // need 3D material
		ExceptionT::GeneralFail("FS_SCNIMFT::NewMaterialList","Axisymmetric formulation needs 3D material\n");

	if (size > 0) {
		 /* material support */
		 if (!fFSMatSupport) {
		 	fFSMatSupport = new FSMatSupportT(nsd, 1);      
			if (qIsAxisymmetric)
				fFSMatSupport->SetNumSD(3);
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

		if (nsd == 2)
			return new FSSolidMatList2DT(size, *fFSMatSupport);
		else if (nsd == 3)
			return new FSSolidMatList3DT(size, *fFSMatSupport);
	} else {
	 	if (nsd == 2)
			return new FSSolidMatList2DT;
		else if (nsd == 3)
			return new FSSolidMatList3DT;
	}
	
	/* no match */
	return NULL;
}
	
// XML stuff below

/* describe the parameters needed by the interface */
void FS_SCNIMFT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SCNIMFT::DefineParameters(list);

}

/* information about subordinate parameter lists */
void FS_SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SCNIMFT::DefineSubs(sub_list);

	/* element blocks for underyling connectivity -- TEMP */
	sub_list.AddSub("fd_connectivity_element_block");
}

/* return the description of the given inline subordinate parameter list */
void FS_SCNIMFT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "large_strain_material_choice") {
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("large_strain_material_2D");
		sub_lists.AddSub("large_strain_material_3D");
	} else /* inherited */
		SCNIMFT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FS_SCNIMFT::NewSub(const StringT& name) const
{
   if  (name == "fd_connectivity_element_block") {
    
	  ParameterContainerT* block = new ParameterContainerT(name);

	  block->AddSub("block_ID_list",ParameterListT::Once);

	  block->AddSub("large_strain_material_choice", ParameterListT::Once, true);

	  block->SetSubSource(this);

	  return block;

  }
  else /* inherited */
	return SCNIMFT::NewSub(name);
}

/* accept parameter list */
void FS_SCNIMFT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SCNIMFT::TakeParameterList(list);

	/* deformation gradients */
	int nsd = NumSD();
	fF_list.Dimension(1);
	fF_list[0].Dimension(nsd);

	/* casts are safe since class contructs materials list - just one material */
	ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
	FSSolidMatT* mat = (FSSolidMatT*) pcont_mat;
	if (mat->Need_F_last()) {
		fF_last_list.Dimension(1);
		fF_last_list[0].Dimension(nsd);
	}
}
