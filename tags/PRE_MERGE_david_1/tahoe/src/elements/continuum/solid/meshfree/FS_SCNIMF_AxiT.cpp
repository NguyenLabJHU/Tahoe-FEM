/* $Id: FS_SCNIMF_AxiT.cpp,v 1.5 2004-09-24 23:44:25 cjkimme Exp $ */
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

const double twoPi = 2.*acos(-1.0);
const int kRadialDirection = 0; /* the x direction is radial */
const int kNSD = 2;
const double kZeroTol = 0.0000001;

/* constructors */
FS_SCNIMF_AxiT::FS_SCNIMF_AxiT(const ElementSupportT& support, const FieldT& field):
	SCNIMFT(support, field),
	fFSMatSupport(NULL),
	circumferential_B()
{
	SetName("fd_mfparticle_axi");
}

FS_SCNIMF_AxiT::FS_SCNIMF_AxiT(const ElementSupportT& support):
	SCNIMFT(support),
	fFSMatSupport(NULL),
	circumferential_B()
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
	const char* disp[2] = {"D_R", "D_Z"};
	int num_labels = 4; // displacements
	int num_stress=0;

	const char* stress[6];
	const char* strain[6];
	
	stress[0] = "s11";
	stress[1] = "s22";
	strain[0] = "E11";
	strain[1] = "E22";
	num_stress=6;
	stress[2] = "s33";
	stress[3] = "s23";
	stress[4] = "s13";
	stress[5] = "s12";
	strain[2] = "E33";
	strain[3] = "E23";
	strain[4] = "E13";
	strain[5] = "E12"; 
	
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

void FS_SCNIMF_AxiT::WriteOutput(void)
{
	ofstreamT& out = ElementSupport().Output();

	const char caller[] = "FS_SCNIMF_AxiT::WriteOutput";
	
	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
	int rank = ElementSupport().Rank();

	/* dimensions */
	int ndof = NumDOF();
	int num_stress = 6;
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

	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
	{
		ExceptionT::GeneralFail("FS_SCNIMF_AxiT::RHSDriver","Cannot get material\n");
	}
	
	int nNodes = fNodes.Length();

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	ArrayT<dMatrixT> Flist(1);
	Flist[0].Dimension(fSD);
	dMatrixT& Fdef = Flist[0];
	dMatrixT BJ(fSD*fSD, fSD);
	dMatrixT E(fSD);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);

	/* collect displacements */
	dArrayT vec, values_i;
	double F_33;
	for (int i = 0; i < nNodes; i++) 
	{
		int   tag_i = (partition_nodes) ? (*partition_nodes)[i] : i;
		int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

		/* values for particle i */
		n_values.RowAlias(local_i, values_i);

		/* copy in */
		vec.Set(ndof, values_i.Pointer());
		coords.RowCopy(tag_i, vec);
		vec.Set(ndof, values_i.Pointer() + ndof);
		vec = 0.;
			
		LinkedListT<int>& nodal_supp = fNodalSupports[i];
		LinkedListT<double>& phi_i = fNodalPhi[i];
		nodal_supp.Top(); phi_i.Top();
		while (nodal_supp.Next() && phi_i.Next()) 
		{
			vec.AddScaled(*(phi_i.CurrentValue()), u(*(nodal_supp.CurrentValue())));
		}

		// Convert initial coordinates to current
	        for (int j = 0; j < ndof; j++) 
		  values_i[j] += vec[j];
		
		// Compute smoothed deformation gradient
		Fdef = 0.0; F_33 = 0.;
		dArrayT* bVec_i = bVectorArray(i);
		double* b_33 = circumferential_B(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < nodalCellSupports.MinorDim(i); j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			F_33 += u(*supp_i,0) * *b_33++; 
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
		
		const double* stress = fCurrMaterial->s_ij().Pointer();

		double* inp_val = values_i.Pointer() + 2*ndof;
		
		*inp_val++ = fVoronoiCellVolumes[i];
		
		E *= 0.5;
		*inp_val++ = E[0];
		*inp_val++ = E[fSD+1];
		*inp_val++ = E.Last();
		*inp_val++ = E[fSD];
		*inp_val++ = E[2*fSD];
		*inp_val++ = E[2*fSD+1];
		
		for (int j = 0; j < num_stress; j++)
			*inp_val++ = stress[j]; 

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

/* contribution to the nodal residual forces */
//const dArray2DT& FS_SCNIMF_AxiT::InternalForce(int group)
//{
	/* check */
//	if (group != Group())
//		ExceptionT::GeneralFail("FS_SCNIMF_AxiT::InternalForce", 
//			"expecting solver group %d not %d", Group(), group);
//	return fForce;
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

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
		ArrayT<dMatrixT> Flist(1);
		Flist[0].Dimension(3);
		dMatrixT& F3D = Flist[0];
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
			double w_i = fVoronoiCellVolumes[i]*constK*twoPi*fDeloneVertices(i,0); // integration weights
			
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
			fFSMatSupport->SetDeformationGradient(&Flist);
		
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

void FS_SCNIMF_AxiT::RHSDriver(void)
{
	/* function name */
	const char caller[] = "ParticlePairT::RHSDriver2D";

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
		double* volume = fVoronoiCellVolumes.Pointer();
		for (int i = 0; i < nNodes; i++)
		{
			acc = a(*nodes++);
			for (int j = 0; j < fSD; j++)
				*ma++ = *volume * *acc++;
			volume++;
		}
		fLHS *= fCurrMaterial->Density();
	}

	fForce = 0.0;
	ArrayT<dMatrixT> Flist(1);
	Flist[0].Dimension(3);
	dMatrixT& F3D = Flist[0];
	dMatrixT BJ(4, 2), fStress3D(3), fStress2D(2), fCauchy(3), Finverse(3);
	double F_33, S_33, J;
	dMatrixT F2D(2);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++)
	{
		double w_i = fVoronoiCellVolumes[i]*twoPi*fDeloneVertices(i,0); // integration weight
		
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
		fFSMatSupport->SetDeformationGradient(&Flist);
		
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
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != fSD) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad minorDim");
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
	else // fSD == 3
	{   // I haven't changed this yet
		Bptr[11] = Bptr[16] = *bVector++;
		Bptr[7] = *bVector;
		Bptr[5] = Bptr[15] = *bVector++;
		Bptr[14] = *bVector;
		Bptr[4] = Bptr[9] = *bVector;
	}
}

void FS_SCNIMF_AxiT::CollectMaterialInfo(const ParameterListT& all_params,
				  ParameterListT& mat_params) const
{
	const char caller[] = "FS_SCNIMFT::CollectMaterialInfo";

	/* initialize */
	mat_params.Clear();

        int num_blocks = all_params.NumLists("large_strain_element_block");
	for (int i = 0; i < num_blocks; i++) {
	  
	  const ParameterListT& block = all_params.GetList("large_strain_element_block",i);

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
		 }
	 
		return new FSSolidMatList3DT(size, *fFSMatSupport);
	} else 
	  return new FSSolidMatList3DT;
	
	/* no match */
	return NULL;
}

void FS_SCNIMF_AxiT::ComputeBMatrices(void)
{
	/* possible best implementation is to loop over all Delone edges
	 * and compute all the necessary values only once per Voronoi
	 * facet. This approach minimizes number of times that the support of
	 * an arbitrary point in space (the Voronoi facet centroid) has to be
	 * found.
	 */

	/* Here for the Axisymmetric case, I need to add a row to each B that is
	 * {Psi_I(X_L)/R)L,0.} I have to redimension each B-vector. 
	 */

	const char caller[] = "FS_SCNIMF_AxiT::ComputeBMatrices";

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	int nNodes = fNodes.Length();
	nodeWorkSpace.Dimension(nNodes);
	facetWorkSpace.Dimension(nNodes);
	circumferentialWorkSpace.Dimension(nNodes);
	
	dArrayT zeroFacet(3);
	zeroFacet = 0.0;
	double zeroSingle = 0.;
	for (int i = 0; i < nNodes; i++)
	{
		int l_supp_i = nodeSupport.MinorDim(i);
		iArrayT supp_i(l_supp_i);
		supp_i.Copy(nodeSupport(i));
		supp_i.SortAscending();
		nodeWorkSpace[i].AppendArray(l_supp_i, supp_i.Pointer());
		facetWorkSpace[i].AppendArray(l_supp_i, zeroFacet);
		circumferentialWorkSpace[i].AppendArray(l_supp_i, zeroSingle);
	}
	
	/* integration */
	int nfn = 2;
	int nsd = 2;
	ParentDomainT domain(GeometryT::kLine, fNumIP, nfn);
	domain.Initialize();
	LocalArrayT facet_coords(LocalArrayT::kInitCoords, nfn, nsd);
	facet_coords.SetGlobal(fVoronoiVertices);
	iArrayT keys;
	dArrayT ip_coords(nsd);
	dMatrixT jacobian(nsd, 1);
	const double* ip_weight = domain.Weight();

	dArrayT facetNormal(fSD), facetIntegral(fSD);
	double* currentB, *currentI;
	int n_0, n_1;
	bool traverseQ_0, traverseQ_1;
	int *next_0, *next_1;
	for (int i = 0; i < fDeloneEdges.MajorDim(); i++)
	{
		n_0 = fDeloneEdges(i,0);
		n_1 = fDeloneEdges(i,1);
		facetNormal.DiffOf(fDeloneVertices(n_1), fDeloneVertices(n_0));
		facetNormal.UnitVector();
		
		fDualFacets.RowAlias(i,keys);
		facet_coords.SetLocal(keys);
		for (int ii = 0; ii < fNumIP; ii++) {

		  /* jacobian of the coordinate transformation */
		  domain.DomainJacobian(facet_coords, ii, jacobian);
		  double jw = ip_weight[ii]*domain.SurfaceJacobian(jacobian);
		  
		  /* integration point coordinates */
		  domain.Interpolate(facet_coords, ip_coords, ii);

		  if (!fNodalShapes->SetFieldAt(ip_coords, NULL)) // shift = 0 or not ?
		    ExceptionT::GeneralFail(caller,"Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
		
		iArrayT ip_cover(fNodalShapes->Neighbors());	
		int n_ip_cover = ip_cover.Length();	
		iArrayT ip_cover_key(n_ip_cover);
		ip_cover_key.SetValueToPosition();
		ip_cover_key.SortAscending(ip_cover);

		LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
		LinkedListT<int>& supp_1 = nodeWorkSpace[n_1];
		LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
		LinkedListT< dArrayT >& bVectors_1 = facetWorkSpace[n_1];
		LinkedListT< double >& circumf_0 = circumferentialWorkSpace[n_0];
		LinkedListT< double >& circumf_1 = circumferentialWorkSpace[n_1];
		int s_0 = -1;
		int s_1 = -1;
		/* Simultaneously loop over support of the two nodes that are endpoints of the
		 * current Delone edge and the nodes in the support of the midpoint of this
		 * edge. If a node covering the centroid is not in the support of n_0 or n_1,
		 * insert that covering node into the sorted list.
		 */
		 
		int* c = ip_cover.Pointer();
		int* c_j = ip_cover_key.Pointer();
		
		supp_0.Top(); bVectors_0.Top(); circumf_0.Top();
		supp_1.Top(); bVectors_1.Top(); circumf_1.Top();
		next_0 = supp_0.CurrentValue();
		next_1 = supp_1.CurrentValue();
		for (int j = 0; j < n_ip_cover; j++, c++, c_j++)
		{
			facetIntegral = facetNormal;
			facetIntegral *= fDualAreas[i]*phiValues[*c_j];	
			
			if (next_0)
				traverseQ_0 = *next_0 <= *c;
			else
				traverseQ_0 = false;
					
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next())
			{
				next_0 = supp_0.PeekAhead(); 
				if (!next_0)
					traverseQ_0 = false;
				else
					if (*next_0 > *c)
						traverseQ_0 = false;
			}
				
			if (s_0 != *c) // means we're not at the end of the linked list
			{
				supp_0.InsertAtCurrent(*c);
				bVectors_0.InsertAtCurrent(zeroFacet);
				circumf_0.InsertAtCurrent(0.);
				s_0 = *c;
				if (supp_0.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_0.Next(); 
					bVectors_0.Next();
					circumf_0.Next();
				}
			}
				
			currentI = facetIntegral.Pointer();
			currentB = bVectors_0.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ += *currentI++;
				
			if (next_1)
				traverseQ_1 = *next_1 <= *c;
			else
				traverseQ_1 = false;
				
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_1 && supp_1.Next(s_1) && bVectors_1.Next())
			{
				next_1 = supp_1.PeekAhead(); 
				if (!next_1)
					traverseQ_1 = false;
				else
					if (*next_1 > *c)
						traverseQ_1 = false;
			}		
								
			if (s_1 != *c)
			{
				supp_1.InsertAtCurrent(*c);
				bVectors_1.InsertAtCurrent(zeroFacet);
				circumf_1.InsertAtCurrent(0.);
				s_1 = *c;
				if (supp_1.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_1.Next(); 
					bVectors_1.Next();
					circumf_1.Next();
				}
			}
				 
			currentI = facetIntegral.Pointer();
			currentB =  bVectors_1.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ -= *currentI++; //NB change in sign; facet normal is inverted!
		}
		}
	}
	
	/** Loop over remaining edges */
	for (int i = 0; i < fNonDeloneEdges.Length(); i++)
	{
		n_0 = fNonDeloneEdges[i];
		facetNormal.Set(fSD, fNonDeloneNormals(i));
		facetNormal.UnitVector();

		/* copy face coordinates with local ordering */
		fSelfDualFacets.RowAlias(i,keys);
		facet_coords.SetLocal(keys);
		for (int ii = 0; ii < fNumIP; ii++) {
		  
		  /* jacobian of the coordinate transformation */
		  domain.DomainJacobian(facet_coords, ii, jacobian);
		  double jw = ip_weight[ii]*domain.SurfaceJacobian(jacobian);

		  /* integration point coordinates */
		  domain.Interpolate(facet_coords, ip_coords, ii);
		
		if (!fNodalShapes->SetFieldAt(ip_coords, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail(caller,"Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
				
		iArrayT ip_cover(fNodalShapes->Neighbors());	
		int n_ip_cover = ip_cover.Length();	
		iArrayT ip_cover_key(n_ip_cover);
		ip_cover_key.SetValueToPosition();
		ip_cover_key.SortAscending(ip_cover);
		
		LinkedListT<int>& supp_0 = nodeWorkSpace[n_0];
		LinkedListT< dArrayT >& bVectors_0 = facetWorkSpace[n_0];
		LinkedListT< double >& circumf_0 = circumferentialWorkSpace[n_0];
		int s_0;
		
		/* Merge support of the boundary node with covering of integration point
		 */
		int* c = ip_cover.Pointer();
		int* c_j = ip_cover_key.Pointer();
		
		supp_0.Top(); bVectors_0.Top(); circumf_0.Top();
		next_0 = supp_0.CurrentValue();
		for (int j = 0; j < n_ip_cover; j++, c++, c_j++)
		{
			facetIntegral = facetNormal;
			facetIntegral *= fBoundaryIntegrationWeights[i]*phiValues[*c_j];		
		
			if (next_0)
				traverseQ_0 = *next_0 <= *c;
			else
				traverseQ_0 = false;
					
			// advance supp_0 and supp_1 until they are greater than or equal to current node
			while (traverseQ_0 && supp_0.Next(s_0) && bVectors_0.Next())
			{
				next_0 = supp_0.PeekAhead(); 
				if (!next_0)
					traverseQ_0 = false;
				else
					if (*next_0 > *c)
						traverseQ_0 = false;
			}
				
			if (s_0 != *c) // means we're not at the end of the linked list
			{
				supp_0.InsertAtCurrent(*c);
				bVectors_0.InsertAtCurrent(zeroFacet);
				circumf_0.InsertAtCurrent(0.);
				s_0 = *c;
				if (supp_0.AtTop()) // if we're inserting at the front, LinkedListT's behavior requires more work
				{
					supp_0.Next(); 
					bVectors_0.Next();
					circumf_0.Next();
				}
			}
				
			currentI = facetIntegral.Pointer();
			currentB =  bVectors_0.CurrentValue()->Pointer();
			for (int k = 0; k < fSD; k++)
				*currentB++ += *currentI++;
		}
		}
	}
	
	// scale integrals by volumes of Voronoi cells
	dArrayT* currFacetIntegral;
	for (int i = 0; i < nNodes; i++)
	{
		LinkedListT<dArrayT>& bVectors_i = facetWorkSpace[i];
		LinkedListT<int>& nodes_i = nodeWorkSpace[i];
		bVectors_i.Top(); nodes_i.Top();
		while ((currFacetIntegral = bVectors_i.Next()))
			*currFacetIntegral *= 1./fVoronoiCellVolumes[i];
	}
	
	// calculate Psi/R terms. These are evaluated nodally, so this additional loop
	// is required
	dArrayT phis, nodal_init_coords;
	for (int i = 0; i < fDeloneVertices.MajorDim(); i++) {
		nodal_init_coords.Set(fSD, fDeloneVertices(i)); // This is the nodal coordinate.
		double R_i = nodal_init_coords[0];
		
		if (R_i >! kZeroTol)
		{
			if (!fNodalShapes->SetFieldAt(nodal_init_coords, NULL)) 
				ExceptionT::GeneralFail(caller,"Shape Function evaluation"
					"failed at node %d\n",i);
					
			const dArrayT& phiValues = fNodalShapes->FieldAt();	
			
			phis.Dimension(phiValues.Length());
			phis = phiValues;	
			phis /= R_i;
		}
		else
		{
			if (!fNodalShapes->SetDerivativesAt(nodal_init_coords))
				ExceptionT::GeneralFail(caller,"Shape Function derivate evaluation"
					"failed at node %d\n",i);
			
			const dArray2DT& DphiValues = fNodalShapes->DFieldAt();

			phis.Dimension(DphiValues.MajorDim());
			//Copy the first column of DphiValues, i.e. d phi / d R = lim_{R -> 0} phi/R
			phis.Copy(DphiValues.Pointer());
		} 		
		
		
		iArrayT node_cover(fNodalShapes->Neighbors());	
		int n_node_cover = node_cover.Length();	
		iArrayT node_cover_key(node_cover);
		node_cover_key.SetValueToPosition();
		node_cover_key.SortAscending(node_cover);

		LinkedListT<int>& supp_0 = nodeWorkSpace[i];
		LinkedListT< double >& circumf_0 = circumferentialWorkSpace[i];
		int s_0 = -1;
		/* Simultaneously loop over support of the two nodes that are endpoints of the
		 * current Delone edge and the nodes in the support of the midpoint of this
		 * edge. If a node covering the centroid is not in the support of n_0 or n_1,
		 * insert that covering node into the sorted list.
		 */
		 
		int* n = node_cover.Pointer();
		int* n_j = node_cover_key.Pointer();
		
		supp_0.Top(); circumf_0.Top();
		next_0 = supp_0.CurrentValue();
		for (int j = 0; j < n_node_cover; j++, n++, n_j++)
		{			
			if (next_0)
				traverseQ_0 = *next_0 <= *n;
			else
				ExceptionT::GeneralFail(caller,"Support list does not exist\n");
					
			// advance supp_0 until it is equal to current node
			// What we're really doing is skipping nodes in support of facet 
			// centroids (our smoothed strain integration points) that are 
			// not in support of the node
			while (traverseQ_0 && supp_0.Next(s_0) && circumf_0.Next())
			{
				next_0 = supp_0.PeekAhead(); 
				if (!next_0)
					traverseQ_0 = false;
				else
					if (*next_0 > *n)
						traverseQ_0 = false;
			}
				
			if (s_0 != *n) 
				ExceptionT::GeneralFail(caller,"Node %d in support of node %d but not in data\n",s_0,*n);
			
			currentI = facetIntegral.Pointer();
			*(circumf_0.CurrentValue()) = phis[*n_j];
				
		}
	}
	
	// move into more efficient storage for computation
	nodalCellSupports.Configure(nodeWorkSpace);
	bVectorArray.Configure(facetWorkSpace);
	circumferential_B.Configure(circumferentialWorkSpace);
	
	for (int i = 0; i < nodalCellSupports.MajorDim(); i++) {
		int* irow_i = nodalCellSupports(i);
		dArrayT* drow_i = bVectorArray(i);
		double* crow_i = circumferential_B(i);
		LinkedListT<int>& ilist = nodeWorkSpace[i];
		LinkedListT<dArrayT>& dlist = facetWorkSpace[i];
		LinkedListT< double >& clist = circumferentialWorkSpace[i];
		ilist.Top(); dlist.Top(); clist.Top();
		while (ilist.Next() && dlist.Next()) {
			*irow_i++ = *(ilist.CurrentValue());
			*drow_i++ = *(dlist.CurrentValue());
			*crow_i++ = *(clist.CurrentValue());
		}
	}
}
	
// XML stuff below

/* describe the parameters needed by the interface */
void FS_SCNIMF_AxiT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	SCNIMFT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void FS_SCNIMF_AxiT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SCNIMFT::DefineSubs(sub_list);
	
	/* element blocks for underlying connectivity -- TEMP */
	sub_list.AddSub("fd_axi_connectivity_block");

}

/* return the description of the given inline subordinate parameter list */
void FS_SCNIMF_AxiT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "large_strain_material_choice") {
		order = ParameterListT::Choice;
		
		/* list of choices */
		sub_lists.AddSub("large_strain_material_3D");
	} else /* inherited */
		SCNIMFT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FS_SCNIMF_AxiT::NewSub(const StringT& name) const
{
   if  (name == "fd_axi_connectivity_block") {
    
	  ParameterContainerT* block = new ParameterContainerT(name);

	  block->AddSub("mfparticle_block_ID_list",ParameterListT::Once);

	  block->AddSub("large_strain_material_choice", ParameterListT::Once, true);

	  block->SetSubSource(this);

	  return block;

  }
  else /* inherited */
	return SCNIMFT::NewSub(name);
}
