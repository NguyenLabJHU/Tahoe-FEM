/* $Id: SS_SCNIMFT.cpp,v 1.5 2004-07-15 08:29:39 paklein Exp $ */
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

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SSSolidMatT.h"
#include "SolidMatSupportT.h"

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
	SCNIMFT(support)
{
	SetName("ss_mfparticle");
}

/* destructor */
SS_SCNIMFT::~SS_SCNIMFT(void)
{
	delete fSSMatSupport;
}

/* generate labels for output data */
void SS_SCNIMFT::GenerateOutputLabels(ArrayT<StringT>& labels)
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
	strain[0] = "e11";
	strain[1] = "e22";
	num_stress = 3;
	
	if (fSD == 3)
	{
		num_stress=6;
		stress[2] = "s33";
	  	stress[3] = "s23";
	 	stress[4] = "s13";
	  	stress[5] = "s12";
	  	strain[2] = "e33";
	  	strain[3] = "e23";
	  	strain[4] = "e13";
	  	strain[5] = "e12";
	} 
	
	if (fSD == 2) 
	{
	  	stress[2] = "s12";
	  	strain[2] = "e12";
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

void SS_SCNIMFT::WriteOutput(void)
{
	ofstreamT& out = ElementSupport().Output();

	const char caller[] = "SS_SCNIMFT::WriteOutput";
	
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
	const dArray2DT& coords = ElementSupport().CurrentCoordinates();

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
	{
		ExceptionT::GeneralFail("SS_SCNIMFT::RHSDriver","Cannot get material\n");
	}
	
	int nNodes = fNodes.Length();

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	
	ArrayT<dSymMatrixT> strainList(1);
	strainList[0].Dimension(fSD);
	dSymMatrixT& strain = strainList[0];
	dMatrixT BJ(fSD == 2 ? 3 : 6, fSD);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);

	/* collect displacements */
	dArrayT vec, values_i;
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
			vec.AddScaled(*(phi_i.CurrentValue()), u(*(nodal_supp.CurrentValue())));
		
		// Compute smoothed strain
		strain = 0.0;
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < nodalCellSupports.MinorDim(i); j++) {
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			BJ.Multx(u(*supp_i++), strain.Pointer(), 1.0, dMatrixT::kAccumulate);
		}	
		
		fSSMatSupport->SetLinearStrain(&strainList);
		
		const double* stress = fCurrMaterial->s_ij().Pointer();

		double* inp_val = values_i.Pointer() + 2*ndof;
		
		*inp_val++ = fVoronoiCellVolumes[i];

		for (int j = 0; j < num_stress; j++)
			*inp_val++ = strain[j];
		for (int j = 0; j < num_stress; j++)
			*inp_val++ = stress[j]; 

	}

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);

}

/* compute specified output parameter(s) */
void SS_SCNIMFT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT SS_SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	return relax;
}

/* write restart data to the output stream */
void SS_SCNIMFT::WriteRestart(ostream& out) const
{
	ElementBaseT::WriteRestart(out);
	
}

/* read restart data to the output stream */
void SS_SCNIMFT::ReadRestart(istream& in)
{
	ElementBaseT::ReadRestart(in);
}

/* contribution to the nodal residual forces */
//const dArray2DT& SS_SCNIMFT::InternalForce(int group)
//{
	/* check */
//	if (group != Group())
//		ExceptionT::GeneralFail("SS_SCNIMFT::InternalForce", 
//			"expecting solver group %d not %d", Group(), group);
//	return fForce;
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void SS_SCNIMFT::LHSDriver(GlobalT::SystemTypeT sys_type)
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
			ExceptionT::GeneralFail("FS_SCNIMFT::LHSDriver","Cannot get material\n");
		}
	
		AssembleParticleMass(fCurrMaterial->Density());
	}
	
	if (formK)
	{
		/* hold the smoothed strain */
		ArrayT<dSymMatrixT> strainList(1);
		strainList[0].Dimension(fSD);
		dSymMatrixT& strain = strainList[0];
		
		/* displacements */
		const dArray2DT& u = Field()(0,0);
	
		/* For now, just one material. Grab it */
		ContinuumMaterialT *mat = (*fMaterialList)[0];
		SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
		if (!fCurrMaterial)
		{
			ExceptionT::GeneralFail("SS_SCNIMFT::LHSDriver","Cannot get material\n");
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
		dMatrixT BJ(fSD == 2 ? 3 : 6, ndof), BK(fSD == 2 ? 3 : 6, ndof), K_JK;
		dMatrixT BJTCijkl(fSD, fSD == 2 ? 3 : 6);
		K_JK.Alias(fLHS);
		LinkedListT<dArrayT> bVectors_j;
		LinkedListT<int> nodeSupport_j;
		for (int i = 0; i < nNodes; i++)
		{	
			double w_i = fVoronoiCellVolumes[i]*constK; // integration weight
			
			int n_supp = nodalCellSupports.MinorDim(i);
			// Compute smoothed strain 
			strain = 0.0;
			
			dArrayT* bVec_i = bVectorArray(i);
			int* supp_i = nodalCellSupports(i);
			for (int j = 0; j < n_supp; j++) { 
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				bVec_i++;
				BJ.Multx(u(*supp_i++), strain.Pointer(), 1.0, dMatrixT::kAccumulate);
			}	
			fSSMatSupport->SetLinearStrain(&strainList);
		
			const dMatrixT& cijkl = fCurrMaterial->c_ijkl();
		
			// sum over pairs to get contribution to stiffness
			supp_i = nodalCellSupports(i);
			bVec_i = bVectorArray(i);
			for (int j = 0; j < n_supp; j++, supp_i++, bVec_i++)
			{
				bVectorToMatrix(bVec_i->Pointer(), BJ);
				BJTCijkl.MultATB(BJ, cijkl, 0);
				col_eqnos.Copy(field_eqnos(*supp_i));
				
				dArrayT* bVec_j = bVectorArray(i);
				int* supp_j = nodalCellSupports(i);
				for (int k = 0; k < n_supp; k++)
				{
					bVectorToMatrix(bVec_j->Pointer(), BK);
					bVec_j++;
					
					if (fSD == 2) 
					{
						BK(2,0) *= 2.; // I either have to do this here or on the RHS 
						BK(2,1) *= 2.; // It's the C_1212 e_12 + C_1221 e_21 factor of 2
					}
					else
						ExceptionT::GeneralFail("SS_SCNIMFT::LHSDriver","Not implemented for 3D yet\n");
						
					// K_JK = BT_J x Cijkl x B_K 
					K_JK.MultAB(BJTCijkl, BK, 0);
					
					K_JK *= w_i*constK;
					
					/* assemble */
					row_eqnos.Copy(field_eqnos(*supp_j++));
					support.AssembleLHS(group, fLHS, col_eqnos, row_eqnos);
				}
			}	
		}
	}
}

void SS_SCNIMFT::RHSDriver(void)
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
		ExceptionT::GeneralFail("SS_SCNIMFT::RHSDriver","Cannot get material\n");
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
	ArrayT<dSymMatrixT> strainList(1);
	strainList[0].Dimension(fSD);
	dSymMatrixT& strain = strainList[0];
	dMatrixT BJ(fSD == 2 ? 3 : 6, fSD);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	for (int i = 0; i < nNodes; i++)
	{
		double w_i = fVoronoiCellVolumes[i]; // integration weight

		int n_supp = nodalCellSupports.MinorDim(i);

		// Compute smoothed strain
		strain = 0.0;
			
		dArrayT* bVec_i = bVectorArray(i);
		int* supp_i = nodalCellSupports(i);
		for (int j = 0; j < n_supp; j++) { 
			bVectorToMatrix(bVec_i->Pointer(), BJ);
			bVec_i++;
			BJ.Multx(u(*supp_i++), strain.Pointer(), 1.0, dMatrixT::kAccumulate);
		}	
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
	
	fForce *= -constKd;

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

void SS_SCNIMFT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
#if __option(extended_errorcheck)
	if (BJ.Rows() != fSD*(fSD+1)/2) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.Cols() != fSD) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad minorDim");
#endif

	double* Bptr = BJ.Pointer();
	
	BJ = 0.;
	Bptr[0] = *bVector;
	if (fSD == 2)
	{
		Bptr[5] = .5 * *bVector++;
		Bptr[4] = *bVector;
		Bptr[2] = .5 * *bVector;
	}
	else // fSD == 3
	{
		Bptr[11] = Bptr[16] = 0.5 * *bVector++;
		Bptr[7] = *bVector;
		Bptr[5] = Bptr[15] = 0.5 * *bVector++;
		Bptr[14] = *bVector;
		Bptr[4] = Bptr[9] = 0.5 * *bVector;
	}
}


void SS_SCNIMFT::ReadMaterialData(ifstreamT& in)
{
	/* base class */
	SCNIMFT::ReadMaterialData(in);

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SSSolidMatT* mat = (SSSolidMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset] = mat->Need_Strain();
		needs[fNeedsOffset + 1] = mat->Need_Strain_last();
		
		/* consistency */
		needs[0] = needs[0] || needs[fNeedsOffset];
		needs[2] = needs[2] || needs[fNeedsOffset + 1];
	}

}

/* return a pointer to a new material list */
MaterialListT* SS_SCNIMFT::NewMaterialList(int nsd, int size)
{
	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = new SSMatSupportT(fSD, 1);
			
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SS_SCNIMFT::NewMaterialList","Could not instantiate material support\n");
		}

		if (nsd == 1)
			return new SSSolidMatList1DT(size, *fSSMatSupport);
		else if (nsd == 2)
			return new SSSolidMatList2DT(size, *fSSMatSupport);
		else if (nsd == 3)
			return new SSSolidMatList3DT(size, *fSSMatSupport);
		else
			return NULL;
	}
	else
	{
		if (nsd == 1)
			return new SSSolidMatList1DT;
		else if (nsd == 2)
			return new SSSolidMatList2DT;
		else if (nsd == 3)
			return new SSSolidMatList3DT;
		else
			return NULL;
	}	
}
	
// XML stuff below

/* describe the parameters needed by the interface */
void SS_SCNIMFT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

}

/* information about subordinate parameter lists */
void SS_SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

}

/* return the description of the given inline subordinate parameter list */
void SS_SCNIMFT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	/* inherited */
	ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SS_SCNIMFT::NewSub(const StringT& name) const
{
	/* inherited */
	return ElementBaseT::NewSub(name);
}
