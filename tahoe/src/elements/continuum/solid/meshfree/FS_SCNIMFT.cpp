/* $Id: FS_SCNIMFT.cpp,v 1.14 2005-01-25 02:23:58 cjkimme Exp $ */
#include "FS_SCNIMFT.h"

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
#include "ParameterContainerT.h"

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"

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

/* generate labels for output data */
void FS_SCNIMFT::GenerateOutputLabels(ArrayT<StringT>& labels)
{
  	/* Reference Configuration */
	const char* ref[3] = {"X", "Y", "Z"};

	/* displacement labels */
	const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	int num_labels = 2*fSD; // displacements
	int num_stress=0;

	const char* stress[6];
	const char* strain[6];
	
	stress[0] = "s11";
	stress[1] = "s22";
	strain[0] = "E11";
	strain[1] = "E22";
	num_stress = 3;
	
	if (fSD == 3) {
		num_stress=6;
		stress[2] = "s33";
	  	stress[3] = "s23";
	 	stress[4] = "s13";
	  	stress[5] = "s12";
	  	strain[2] = "E33";
	  	strain[3] = "E23";
	  	strain[4] = "E13";
	  	strain[5] = "E12";
	} 
	
	if (fSD == 2) {
	  	stress[2] = "s12";
	  	strain[2] = "E12";
	}
		
	num_labels += 2 * num_stress + 1; 
	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < fSD; dex++)
		labels[dex] = ref[dex];
	for (int ns = 0 ; ns < fSD; ns++)
	  	labels[dex++] = disp[ns];
	labels[dex++] = "mass";
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = strain[ns];
	for (int ns = 0 ; ns < num_stress; ns++)
		labels[dex++] = stress[ns];
}

void FS_SCNIMFT::WriteOutput(void)
{
	ofstreamT& out = ElementSupport().Output();

	const char caller[] = "FS_SCNIMFT::WriteOutput";
	
	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
	int rank = ElementSupport().Rank();

	/* dimensions */
	int ndof = NumDOF();
	int num_stress = fSD == 2 ? 3 : 6;
	int num_output = 2*ndof + 1 + 2 * num_stress; /* ref. coords, displacements, mass, e, s */

	/* number of nodes */
	const ArrayT<int>* partition_nodes = comm_manager.PartitionNodes();
	int non = (partition_nodes) ? partition_nodes->Length() : 
								ElementSupport().NumNodes();

	/* map from partition node index */
	const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();

	/* output arrays length number of active nodes */
	dArray2DT n_values(non, num_output), e_values;
	n_values = 0.0;

	/* global coordinates */
	const dArray2DT& coords = ElementSupport().InitialCoordinates();

	/* the field */
	const FieldT& field = Field();
	const dArray2DT* velocities = NULL;
	if (field.Order() > 0) velocities = &(field[1]);

	/* check 2D */
	if (NumDOF() != 2) ExceptionT::GeneralFail(caller, "2D only: %d", NumDOF());

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
		ExceptionT::GeneralFail("FS_SCNIMFT::RHSDriver","Cannot get material\n");
	
	int nNodes = fNodes.Length();

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	ArrayT<dMatrixT> Flist(1);
	Flist[0].Dimension(fSD);
	dMatrixT& Fdef = Flist[0];
	dMatrixT BJ(fSD*fSD, fSD), Finverse(fSD);
	dMatrixT E(fSD), fStress(fSD), fCauchy(fSD);
	double J;
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);

	/* collect displacements */
	dArrayT vec, values_i;
	for (int i = 0; i < nNodes; i++) {
		int   tag_i = (partition_nodes) ? (*partition_nodes)[fNodes[i]] : fNodes[i];
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

		/* values for particle i */
		n_values.RowAlias(local_i, values_i);

		/* copy in */
		vec.Set(ndof, values_i.Pointer());
		coords.RowCopy(tag_i, vec);
		vec.Set(ndof, values_i.Pointer() + ndof);
		vec = 0.;
			
		int* nodal_supp = fNodalSupports(i);
		double* phi_i = fNodalPhi(i);
		for (int j = 0; j < fNodalPhi.MinorDim(i); j++)
			vec.AddScaled(*phi_i++, u(*nodal_supp++));

		// Convert initial coordinates to current coordinates
		for (int j = 0; j < ndof; j++) 
			values_i[j] += vec[j];
		
		// Compute smoothed deformation gradient
		Fdef = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < nodalCellSupports.MinorDim(i); j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			BJ.Multx(u(*supp_i++), Fdef.Pointer(), 1.0, dMatrixT::kAccumulate);
		}
		E = 0.;
		E.MultATB(Fdef, Fdef, dMatrixT::kWhole);
		E += Fdef;
		for (int rows = 0; rows < fSD; rows++)
			for (int cols = 0; cols < fSD; cols++)
				E(rows,cols) += Fdef(cols,rows);
		Fdef.PlusIdentity();
		fFSMatSupport->SetDeformationGradient(&Flist);
		J = Fdef.Det();
		
		fCurrMaterial->s_ij().ToMatrix(fCauchy);
		Finverse.Inverse(Fdef);
		Finverse *= J; // compute J F^-1
		fStress.MultABT(fCauchy, Finverse); // compute PK1
		double* stress = fStress.Pointer();

		double* inp_val = values_i.Pointer() + 2*ndof;
		
		*inp_val++ = fCellVolumes[i];
		
		E *= 0.5;
		//for (int j = 0; j < num_stress; j++)
		//	*inp_val++ = E[j];
		*inp_val++ = E[0];
		*inp_val++ = E[fSD+1];
		if (fSD == 3) {
			*inp_val++ = E.Last();
			*inp_val++ = E[fSD];
			*inp_val++ = E[2*fSD];
			*inp_val++ = E[2*fSD+1];
		} else
			*inp_val++ = E[fSD];	
		for (int j = 0; j < num_stress; j++)
			*inp_val++ = stress[j]; 

	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);

}

/* compute specified output parameter(s) */
void FS_SCNIMFT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT FS_SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	return relax;
}

/* write restart data to the output stream */
void FS_SCNIMFT::WriteRestart(ostream& out) const
{
	ElementBaseT::WriteRestart(out);
	
}

/* read restart data to the output stream */
void FS_SCNIMFT::ReadRestart(istream& in)
{
	ElementBaseT::ReadRestart(in);
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
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!fCurrMaterial)
			ExceptionT::GeneralFail("FS_SCNIMFT::LHSDriver","Cannot get material\n");
	
		AssembleParticleMass(fCurrMaterial->Density());
	}
	
	if (formK) {
		/* hold the smoothed strain */
		ArrayT<dMatrixT> Flist(1);
		Flist[0].Dimension(fSD);
		dMatrixT& Fdef = Flist[0];
		
		/* displacements */
		const dArray2DT& u = Field()(0,0);
	
		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!fCurrMaterial)
			ExceptionT::GeneralFail("FS_SCNIMFT::LHSDriver","Cannot get material\n");

		int nNodes = fNodes.Length();

		/* assembly information */
		int group = Group();
		int ndof = NumDOF();
		fLHS.Dimension(ndof);
		const ElementSupportT& support = ElementSupport();
		
		const iArray2DT& field_eqnos = Field().Equations();
		iArrayT row_eqnos(ndof); 
		iArrayT col_eqnos(ndof);
		dMatrixT BJ(fSD*fSD, ndof), BK(fSD*fSD, ndof), FTs(fSD), fStress(fSD), fCauchy(fSD);
		dMatrixT Tijkl(fSD*fSD), BJTCijkl(fSD, fSD*fSD), K_JK, Finverse(fSD);
		K_JK.Alias(fLHS);
		LinkedListT<dArrayT> bVectors_j;
		LinkedListT<int> nodeSupport_j;
		for (int i = 0; i < nNodes; i++) {	
			double w_i = fCellVolumes[i]*constK; // integration weights
			
			int n_supp = nodalCellSupports.MinorDim(i);
			
			// Compute smoothed deformation gradient
			Fdef = 0.0;
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				bVec_i++;
				BJ.Multx(u(*supp_i++), Fdef.Pointer(), 1.0, dMatrixT::kAccumulate);
			}
			Fdef.PlusIdentity(); // convert to F
			fFSMatSupport->SetDeformationGradient(&Flist);
		
			const dMatrixT& cijkl = fCurrMaterial->C_IJKL();
			fCurrMaterial->s_ij().ToMatrix(fCauchy);
			Finverse.Inverse(Fdef);
			Finverse *= Fdef.Det(); // compute J F^-1
			fStress.MultABT(fCauchy, Finverse); // compute PK1
		
			// FTs = F^T sigma
			FTs.MultATB(Fdef, fStress);
			// T_11 = T_22 = FTs_11 T_13 = T_24 = FTs_12
			// T_31 = T_42 = FTs_21 T_33 = T_44 = FTs_22
			Tijkl.Expand(FTs, fSD, dMatrixT::kOverwrite);
			// Tijkl = 0.; used to debug material stiffness
			// use brute force, low-level until I find optimal way
			Tijkl[0] += cijkl[0]; // += C_11
			Tijkl[1] += cijkl[6]; // += C_13
			Tijkl[2] += cijkl[6]; // += C_13
			Tijkl[3] += cijkl[3]; // += C_12
			Tijkl[4] += cijkl[6]; // += C_13
			Tijkl[5] += cijkl[8]; // += C_33
			Tijkl[6] += cijkl[8]; // += C_33
			Tijkl[7] += cijkl[7]; // += C_23
			Tijkl[8] += cijkl[6]; // += C_13
			Tijkl[9] += cijkl[8]; // += C_33
			Tijkl[10] += cijkl[8]; // += C_33
			Tijkl[11] += cijkl[7]; // += C_23
			Tijkl[12] += cijkl[3]; // += C_12
			Tijkl[13] += cijkl[7]; // += C_23
			Tijkl[14] += cijkl[7]; // += C_23
			Tijkl.Last() += cijkl[4]; // += C_22
			
			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++) {
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				//BJTCijkl.MultATB(BJ, cijkl, 0);
				
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
		ExceptionT::GeneralFail("FS_SCNIMFT::RHSDriver","Cannot get material\n");
	
	int nNodes = fNodes.Length();

	if (formMa) {
		if (Field().Order() < 2)
			ExceptionT::GeneralFail(caller,"Field's Order does not have accelerations\n");
	
		fLHS = 0.0; // fLHS.Length() = nNodes * fSD;
		const dArray2DT& a = Field()(0,2); // accelerations
		double* ma = fLHS.Pointer();
		const double* acc;
		int* nodes = fNodes.Pointer();
		double* volume = fCellVolumes.Pointer();
		for (int i = 0; i < nNodes; i++) {
			acc = a(*nodes++);
			for (int j = 0; j < fSD; j++)
				*ma++ = *volume * *acc++;
			volume++;
		}
		fLHS *= fCurrMaterial->Density();
	}

	fForce = 0.0;
	ArrayT<dMatrixT> Flist(1);
	Flist[0].Dimension(fSD);
	dMatrixT& Fdef = Flist[0];
	dMatrixT BJ(fSD*fSD, fSD), fCauchy(fSD), Finverse(fSD), fStress(fSD);
	double J;
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++) {
		double w_i = fCellVolumes[i]; // integration weight
		
		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed deformation gradient
		Fdef = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			BJ.Multx(u(*supp_i++), Fdef.Pointer(), 1.0, dMatrixT::kAccumulate);
		}
		Fdef.PlusIdentity(); // convert to F 	
		fFSMatSupport->SetDeformationGradient(&Flist);
		J = Fdef.Det();
		if (J <= 0.0)
			ExceptionT::BadJacobianDet("FS_SCNIMFT::FormKd");
		
		fCurrMaterial->s_ij().ToMatrix(fCauchy);
		Finverse.Inverse(Fdef);
		Finverse *= J; // compute J F^-1
		fStress.MultABT(fCauchy, Finverse); // compute PK1
		
		supp_i = nodalCellSupports(i);
		bVec_i = bVectorArray(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			double* fint = fForce(*supp_i++);
			BJ.MultTx(fStress.Pointer(), fint, w_i, dMatrixT::kAccumulate);
		}
	}
	
	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

void FS_SCNIMFT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
#if __option(extended_errorcheck)
	if (BJ.Rows() != fSD*2) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != fSD) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad minorDim");
#endif

	double* Bptr = BJ.Pointer();
	
	BJ = 0.;
	Bptr[0] = *bVector;
	if (fSD == 2) {
		Bptr[5] = *bVector++;
		Bptr[2] = *bVector;
		Bptr[7] = *bVector;
	} else { // fSD == 3
	    // I haven't changed this yet
		Bptr[11] = Bptr[16] = *bVector++;
		Bptr[7] = *bVector;
		Bptr[5] = Bptr[15] = *bVector++;
		Bptr[14] = *bVector;
		Bptr[4] = Bptr[9] = *bVector;
	}
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

	if (size > 0) {
		 /* material support */
		 if (!fFSMatSupport) {
		 	fFSMatSupport = new FSMatSupportT(nsd, 1);      
		 	if (!fFSMatSupport)
		 		ExceptionT::GeneralFail("FS_SCNIMFT::NewMaterialList","Could not instantiate material support\n");
			fFSMatSupport->SetFEManager(&ElementSupport().FEManager());
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
