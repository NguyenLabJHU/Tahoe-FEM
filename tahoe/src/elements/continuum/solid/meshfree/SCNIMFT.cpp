/* $Id: SCNIMFT.cpp,v 1.1 2004-01-27 01:26:19 cjkimme Exp $ */
#include "SCNIMFT.h"

#include "ArrayT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "fstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "iGridManagerT.h"
#include "iNodeT.h"
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "BasicFieldT.h"

#include "MeshFreeNodalShapeFunctionT.h"
#include "ContinuumMaterialT.h"
#include "SolidMaterialT.h"
#include "SolidMatSupportT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"


#ifdef __QHULL__
#include "CompGeomT.h"
#endif

using namespace Tahoe;

/* constructors */
SCNIMFT::SCNIMFT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	//MeshFreeFractureSupportT(ElementSupport().Input()),
	fSD(ElementSupport().NumSD()),
	fMaterialList(NULL),
	fSSMatSupport(NULL),
	fForce_man(0, fForce, field.NumDOF()),
	//fNodes(),
	fLocInitCoords(LocalArrayT::kInitCoords),
	//fLocDisp(LocalArrayT::kDisp),
	fFakeGeometry(NULL),
	fVoronoi(NULL)
	//fVoronoiVertices(),
	//fVoronoiCells(),
	//fVoronoiFacetIndices(),
	//fVoronoiFacetAreas(),
	//fVoronoiFacetNormals(),
	//fVoronoiCellVolumes()
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	
	/* check */
	//if (AutoBorder() && ElementSupport().Size() > 1)
	//	ExceptionT::BadInputValue("MeshFreeStrainSmoothedSST::MeshFreeStrainSmoothedSST", "auto-border not support in parallel");

}

SCNIMFT::SCNIMFT(const ElementSupportT& support):
	ElementBaseT(support),
	//MeshFreeFractureSupportT(ElementSupport().Input()),
	fNodes(),
	fFakeGeometry(NULL),
	fVoronoi(NULL)
{
	SetName("mfparticle");

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	
	/* check */
	//if (AutoBorder() && ElementSupport().Size() > 1)
	//	ExceptionT::BadInputValue("MeshFreeStrainSmoothedSST::MeshFreeStrainSmoothedSST", "auto-border not support in parallel");

}

/* destructor */
SCNIMFT::~SCNIMFT(void)
{
	/* free search grid */
	//delete fGrid;
	
	delete fConnectivities[0];
	
}

/* initialization */
void SCNIMFT::Initialize(void)
{
	const char caller[] = "SCNIMFT::Initialize";

	/* inherited */
	ElementBaseT::Initialize();
	
	/* re-dimension "element" force and stiffness contributions */
	fLHS.Dimension(fSD);
	
	/* allocate work space */
	fForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

	/* write parameters */
	ifstreamT& in = ElementSupport().Input();
	ostream& out = ElementSupport().Output();

	int qComputeVoronoiCell;
	in >> qComputeVoronoiCell;

	if (qComputeVoronoiCell)
	{
#ifndef __QHULL__
	ExceptionT::GeneralFail(caller,"Requires the QHull library\n");
#else 

		// Do the heavy lifting for the Voronoi Diagram now
		fVoronoi = new CompGeomT(fDeloneVertices);
		fVoronoi->ComputeVoronoiDiagram();
		dArray2DT& vorVerts = fVoronoi->VoronoiVertices();
		fVoronoiVertices.Alias(/*vorVerts.MajorDim(), vorVerts.MinorDim(),*/ vorVerts/*.Pointer()*/);
		fVoronoiCells.Alias(fVoronoi->VoronoiCells()); 		
		fVoronoiFacetIndices.Alias(fVoronoi->VoronoiFacetIndices());

		// Data for integration over boundary of each Voronoi region
		fVoronoiFacetAreas.Alias(fVoronoi->VoronoiFacetAreas());
		fVoronoiFacetNormals.Alias(fVoronoi->VoronoiFacetNormals());
		fVoronoiCellVolumes.Alias(fVoronoi->VoronoiCellVolumes());
		
		// Still need to worry about the boundary. Some cells might be clipped
        // I'm sidestepping that for now by reading in a good decomposition

		// Write output to file
		StringT vCellFile;
		in >> vCellFile;
		
		ofstreamT vout(vCellFile);

		if (vout.is_open())	
		{
			VoronoiDiagramToFile(vout);
			vout.close();
		}
		else 
  			cout  << " Unable to save data to file " << fileName << ". Ignoring error \n"; 
#endif
	} 
    else 
    {	// read in Voronoi information from a file
	  	StringT vCellFile;
	    in >> vCellFile;
	    
	    ifstreamT vin(vCellFile);

    	if (!vin.is_open())
			ExceptionT::GeneralFail(caller,"Unable to open file %s for reading", vCellFile);
	    
	    VoronoiDiagramFromFile(vin);  
	    
	    vin.close();
	}
	
	/* shape functions */
	/* only support single list of integration cells for now */
	if (fConnectivities.Length() > 1) {
		cout << "\n MeshFreeSSSolidT::SetShape: multiple element blocks within an\n"
		     <<   "     element group not supported. Number of blocks: " 
		     << fConnectivities.Length() << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* construct shape functions */
	fNodalShapes = new MeshFreeNodalShapeFunctionT(ElementSupport().NumSD(),
		fLocInitCoords, ElementSupport().InitialCoordinates(), *fElementConnectivities[0], 
		fVoronoiVertices, ElementSupport().Input());
	if (!fNodalShapes) throw ExceptionT::kOutOfMemory;
	
	/* echo parameters */
	fNodalShapes->WriteParameters(ElementSupport().Output());
	
	/* initialize */
//	fNodalShapes->Initialize();

	/* MLS stuff */
	fNodalShapes->SetSupportSize();

	/* exchange nodal parameters (only Dmax for now) */
	const ArrayT<int>* p_nodes_in = ElementSupport().ExternalNodes();
	if (p_nodes_in)
	{
		/* skip MLS fit at external nodes */
		iArrayT nodes_in;
		nodes_in.Alias(*p_nodes_in);
		fNodalShapes->SetSkipNodes(nodes_in);
		
		/* exchange */
		CommManagerT& comm = ElementSupport().CommManager();

		/* send all */
		dArray2DT& nodal_params = fNodalShapes->NodalParameters();

		/* initialize the exchange */
		int id = comm.Init_AllGather(nodal_params);
		
		/* do the exchange */
		comm.AllGather(id, nodal_params);
		
		/* clear the communication */
		comm.Clear_AllGather(id);
	}
	
	/* set nodal neighborhoods */
	fNodalShapes->SetNeighborData();

	/* initialize support data */
/*	iArrayT surface_nodes;
	if (fAutoBorder) {
		ArrayT<StringT> IDs;
		ElementBlockIDs(IDs);
		ElementSupport().Model().SurfaceNodes(IDs, surface_nodes,
			fFakeGeometry);
	}
	MeshFreeFractureSupportT::InitSupport(
		ElementSupport().Input(), 
		ElementSupport().Output(),
		fElementCards, 
		surface_nodes, 
		NumDOF(), 
		ElementSupport().NumNodes(),
		&ElementSupport().Model()
	);*/
	
	/* final MLS initializations */
	fNodalShapes->WriteStatistics(ElementSupport().Output());
	
	ComputeBMatrices();	
	
	/** Material Data */
	ReadMaterialData(in);
	WriteMaterialData(out);
	
	//TEMP - only works for one material right now, else would have to check
	//       for the material active within the integration cell (element)
	/*if (HasActiveCracks() && fMaterialList->Length() != 1)
	{	
		cout << "\n MeshFreeSSSolidT::Initialize: can only have 1 material in the group\n";
		cout <<   "     with active cracks" << endl;
		throw ExceptionT::kBadInputValue;
	}*/

//TEMP - write nodal parameters
#if 0
	ostream& out = FEManager().Output();
	out << "\n MeshFreeFDSolidT::Initialize: d_max:\n";
	const dArray2DT& nodal_params = fMFShapes->NodalParameters();
	const iArrayT* node_map = FEManager().NodeMap();
	int d_width = out.precision() + kDoubleExtra;
	out << setw(kIntWidth) << "node"
	<< setw(nodal_params.MinorDim()*d_width) << "d_max" << '\n';
	if (node_map)
		for (int i = 0; i < nodal_params.MajorDim(); i++)
		{
			out << setw(kIntWidth) << (*node_map)[i] + 1;
			nodal_params.PrintRow(i, out);
			out << '\n';
		}
	else
		for (int j = 0; j < nodal_params.MajorDim(); j++)
		{
			out << setw(kIntWidth) << j + 1;
			nodal_params.PrintRow(j, out);
			out << '\n';
		}
#endif

//TEMP - U connectivities
#if 0
	out << "\n MeshFreeFDSolidT::Initialize: connects U: 0, 1,...\n";
	fElemNodesEX->WriteNumbered(out);
	out.flush();
#endif

//TEMP - trace node
#if 0
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	int hit_node = -1;
	if (size == 1)
		hit_node = 13607 - 1;
	else if (size == 2 && rank == 0)
		hit_node = 5420 - 1;
	if (hit_node > 0) TraceNode(ElementSupport().Output(), hit_node, *this);
#endif
}

/* form of tangent matrix */
GlobalT::SystemTypeT SCNIMFT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void SCNIMFT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* writing output */
void SCNIMFT::RegisterOutput(void)
{
	/* "point connectivities" needed for output */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
	if (parition_nodes)
	{
		int num_nodes = parition_nodes->Length();
		fPointConnectivities.Alias(num_nodes, 1, parition_nodes->Pointer());
	}
	else /* ALL nodes */
	{
		fPointConnectivities.Dimension(ElementSupport().NumNodes(), 1);
		iArrayT tmp;
		tmp.Alias(fPointConnectivities);
		tmp.SetValueToPosition();				
	}

	/* block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* get output labels (per node) */
	ArrayT<StringT> n_labels, e_labels;
	GenerateOutputLabels(n_labels);

	/* set output specifier */
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
	OutputSetT output_set(GeometryT::kPoint, fPointConnectivities, n_labels, ChangingGeometry());
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

/* generate labels for output data */
void SCNIMFT::GenerateOutputLabels(ArrayT<StringT>& labels)
{
  	int ndof=NumDOF();
	if (ndof > 3) ExceptionT::GeneralFail("ParticlePairT::GenerateOutputLabels");

	/* displacement labels */
	const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	int num_labels = ndof; // displacements
	int num_stress=0;

	const char* stress[6];
	const char* strain[6];
	
	if (ndof==3)
	{
		num_stress=6;
		stress[2]="s33";
	  	stress[3]="s23";
	 	stress[4]="s13";
	  	stress[5]="s12";
	  	strain[0]="e11";
	  	strain[1]="e12";
	  	strain[2]="e13";
	  	strain[3]="e22";
	  	strain[4]="e23";
	  	strain[5]="e33";
	}
	else 
		if (ndof==2) 
		{
	   		num_stress=3;
	  		stress[0]="s11";
	  		stress[1]="s22";
	  		stress[2]="s12";
	  		strain[0]="e11";
	  		strain[1]="e12";
	  		strain[2]="e22";
	  	}
	  	else // ndof==1 
	   	{
	   		num_stress=1;
	  		stress[0] = "s11";
	  		strain[0] = "e11";
	  	}
		
	num_labels += 2*num_stress;
	labels.Dimension(num_labels);
	int dex = 0;
	for (dex = 0; dex < NumDOF(); dex++)
		labels[dex] = disp[dex];
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=stress[ns];
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=strain[ns];
}

void SCNIMFT::WriteOutput(void)
{
	/* max distance traveled since last reneighboring */
	ofstreamT& out = ElementSupport().Output();

}

/* compute specified output parameter and send for smoothing */
void SCNIMFT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/* trigger reconfiguration */
GlobalT::RelaxCodeT SCNIMFT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	return relax;
}

/* write restart data to the output stream */
void SCNIMFT::WriteRestart(ostream& out) const
{
	ElementBaseT::WriteRestart(out);
	
}

/* read restart data to the output stream */
void SCNIMFT::ReadRestart(istream& in)
{
	ElementBaseT::ReadRestart(in);
}

/* contribution to the nodal residual forces */
//const dArray2DT& SCNIMFT::InternalForce(int group)
//{
	/* check */
//	if (group != Group())
//		ExceptionT::GeneralFail("SCNIMFT::InternalForce", 
//			"expecting solver group %d not %d", Group(), group);
//	return fForce;
//}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return true if connectivities are changing */
bool SCNIMFT::ChangingGeometry(void) const
{
	return ElementSupport().CommManager().PartitionNodesChanging();
}

/* echo element connectivity data */
void SCNIMFT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
#pragma unused(out)
	const char caller[] = "SCNIMFT::EchoConnectivityData";
	
	/* access to the model database */
	ModelManagerT& model = ElementSupport().Model();

	/* read node set ids */
	ArrayT<StringT> ids;
	model.NodeSetList(in, ids);
	model.ManyNodeSets(ids, fNodes);

	// temporarily, read in connectivies of the element group. I need this
	// to piggyback on the search that finds the support of each node
	iArrayT matnums;
	model.ElementBlockList(in, ids, matnums);
	if (ids.Length() != 1)
		ExceptionT::GeneralFail(caller,"Not helping my kludge by having more than one element block\n");
	
	fElementConnectivities.Dimension(1);
	model.ReadConnectivity(ids[0]);

	/* set pointer to connectivity list */
	fElementConnectivities[0] = model.ElementGroupPointer(ids[0]);
	    
	/* store block data  ASSUMING bth block is material number b*/
	//fBlockData[0].Set(new_id, 0, elem_count, 0); 
	  
	/* set up ElementCards and equation arrays */
	//fElementCards.Dimension(elem_count);
	//fEqnos.Dimension(1);
	//fEqnos[0].Dimension(elem_count,fSD);

	/* set local array for coordinates */
	fLocInitCoords.Dimension(fNodes.Length(), fSD);
	fLocInitCoords.SetGlobal(ElementSupport().InitialCoordinates());
	fLocInitCoords.SetLocal(fNodes);

	// Get nodal coordinates to use in Initialize
	// really need to think about how many things I want to keep in memory
	fDeloneVertices.Dimension(fNodes.Length(), model.NumDimensions());
	fDeloneVertices.RowCollect(fNodes, model.Coordinates());
}

/* collecting element group equation numbers */
void SCNIMFT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

	/* dimension equations array */
	fEqnos.Configure(fNodalShapes->NodeNeighbors(), NumDOF());

	/* get local equations numbers */
	Field().SetLocalEqnos(fNodalShapes->NodeNeighbors(), fEqnos);

	/* add to list of equation numbers */
	eq_2.Append(&fEqnos);
}

/* assemble particle mass matrix into LHS of global equation system */
void SCNIMFT::AssembleParticleMass(const dArrayT& mass)
{
#pragma unused(mass)
	
	/* assemble all */
	ElementSupport().AssembleLHS(Group(), fLHS, Field().Equations());
}

/* form group contribution to the stiffness matrix */
void SCNIMFT::LHSDriver(GlobalT::SystemTypeT sys_type)
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
		//AssembleParticleMass(mass);
	}
	
	if (formK)
	{
	
		/* matrix format */
		dMatrixT::SymmetryFlagT format =
			(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
			dMatrixT::kWhole :
			dMatrixT::kUpperOnly;
	
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
			ExceptionT::GeneralFail("SCNIMFT::LHSDriver","Cannot get material\n");
		}

		const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
		int nNodes = fNodes.Length();

		/* assembly information */
		const ElementSupportT& support = ElementSupport();
		int group = Group();
		int ndof = NumDOF();
		fLHS.Dimension(ndof);
		
		iArrayT pair(2);
		const iArray2DT& field_eqnos = Field().Equations();
		iArray2DT pair_eqnos(2, ndof); 
		dMatrixT BJ, BK;
		BJ.Dimension(fSD == 2 ? 3 : 6, ndof);
		BK.Dimension(fSD == 2 ? 3 : 6, ndof);
		for (int i = 0; i < nNodes; i++)
		{	
			double w_i = fVoronoiCellVolumes[i]*constK; // integration weights
			
			// Compute smoothed strain 
			strain = 0.0;
			for (int j = 0; j < nodeSupport.MinorDim(i); j++)
			{
				bVectorToMatrix(bVectors[i](j), BJ);
				BJ.Multx(u(nodeSupport(i,j)), strain.Pointer(), 1.0, 1);
			}	
			fSSMatSupport->SetLinearStrain(&strainList);
		
			const dMatrixT& cijkl = fCurrMaterial->c_ijkl();
		
			// sum over pairs to get contribution to stiffness
			for (int j = 0; j < nodeSupport.MinorDim(i); j++)
			{
				bVectorToMatrix(bVectors[i](j), BJ);
				pair[0] = nodeSupport(i,j);
				for (int k = 0; k < nodeSupport.MinorDim(i); k++)
				{
					pair[1] = nodeSupport(i,k);
					bVectorToMatrix(bVectors[i](k), BK);
					
					fLHS = 0.;
					// K_JK = BT_J x Cijkl B_K 
					fLHS.MultATBC(BJ, cijkl, BK, format);
					
					fLHS *= w_i;
					
					/* assemble */
					pair_eqnos.RowCollect(pair, field_eqnos);
					support.AssembleLHS(group, fLHS, pair_eqnos);
				}
			}	
		}
	}
}

void SCNIMFT::RHSDriver(void)
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

	//TEMP - interial force not implemented
	if (formMa) ExceptionT::GeneralFail(caller, "inertial force not implemented");

	/* assembly information */
	const ElementSupportT& support = ElementSupport();
	int group = Group();
	int ndof = NumDOF();
	
	/* global coordinates */
	const dArray2DT& coords = support.CurrentCoordinates();
	
	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	int nNodes = fNodes.Length();
	
	fForce = 0.0;
	ArrayT<dSymMatrixT> strainList(1);
	strainList[0].Dimension(fSD);
	dSymMatrixT& strain = strainList[0];
	dMatrixT BJ(fSD == 2 ? 3 : 6, fSD);
	
	/* displacements */
	const dArray2DT& u = Field()(0,0);
	
	/* For now, just one material. Grab it */
	ContinuumMaterialT *mat = (*fMaterialList)[0];
	SolidMaterialT* fCurrMaterial = TB_DYNAMIC_CAST(SolidMaterialT*,mat);
	if (!fCurrMaterial)
	{
		ExceptionT::GeneralFail("SCNIMFT::RHSDriver","Cannot get material\n");
	}
	
	for (int i = 0; i < nNodes; i++)
	{
		double w_i = fVoronoiCellVolumes[i]; // integration weight

		if (!i) cout << "Node " << i;
	
		// Compute smoothed strain
		strain = 0.0;
		for (int j = 0; j < nodeSupport.MinorDim(i); j++)
		{
			if (!i) cout << " Node " << nodeSupport(i,j) << " in support " << bVectors[i](j)[0] << " " << bVectors[i](j)[1];
			bVectorToMatrix(bVectors[i](j), BJ);
			if (!i) cout << " u " << u(nodeSupport(i,j))[0] <<" "<<u(nodeSupport(i,j))[1] << " ";
			BJ.Multx(u(nodeSupport(i,j)), strain.Pointer(), 1.0, 1);
			if (!i) cout << strain[0] << " "<< strain[1] << " " << strain[2] << "\n";
		}	
		fSSMatSupport->SetLinearStrain(&strainList);
		
		const double* stress = fCurrMaterial->s_ij().Pointer();
		
		for (int j = 0; j < nodeSupport.MinorDim(i); j++)
		{
			bVectorToMatrix(bVectors[i](j), BJ);
			double* fint = fForce(nodeSupport(i,j));
			BJ.MultTx(stress, fint, w_i, 1);
		}	
	}

	/* assemble */
	ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());

}

void SCNIMFT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SCNIMFT::ReadMaterialData";

	/* construct material list */
	int size;
	in >> size;
	fMaterialList = NewMaterialList(NumSD(), size);
	if (!fMaterialList) ExceptionT::OutOfMemory(caller);

	/* read */
	fMaterialList->ReadMaterialData(in);
	
	/* check range */
	/*for (int i = 0; i < fBlockData.Length(); i++)
		if (fBlockData[i].MaterialID() < 0 || fBlockData[i].MaterialID() >= size)
			ExceptionT::BadInputValue(caller, "material number %d for element block %d is out of range",
				fBlockData[i].MaterialID()+1, i+1);*/
}

/* use in conjunction with ReadMaterialData */
void SCNIMFT::WriteMaterialData(ostream& out) const
{
	fMaterialList->WriteMaterialData(out);

	/* flush buffer */
	out.flush();
}

/* return a pointer to a new material list */
MaterialListT* SCNIMFT::NewMaterialList(int nsd, int size)
{
	/* full list */
	if (size > 0)
	{
		/* material support */
		if (!fSSMatSupport) {
			fSSMatSupport = new SSMatSupportT(fSD, fSD, 1);
			
			if (!fSSMatSupport)
				ExceptionT::GeneralFail("SCNIMFT::NewMaterialList","Could not instantiate material support\n");
				
			/* ElementSupportT sources */
			const ElementSupportT& e_support = ElementSupport();
			fSSMatSupport->SetRunState(e_support.RunState());
			fSSMatSupport->SetStepNumber(e_support.StepNumber());
//	p->SetIterationNumber(e_support.IterationNumber(Group()));
//TEMP - solvers not set up yet. For now, the source for the iteration number will
//       be set in the InitialCondition call for the subclass.
			fSSMatSupport->SetTime(e_support.Time());                              
			fSSMatSupport->SetTimeStep(e_support.TimeStep());
			fSSMatSupport->SetNumberOfSteps(e_support.NumberOfSteps());

			/* set pointer to local array */
			//p->SetLocalArray(fLocDisp);

			//fSSMatSupport = TB_DYNAMIC_CAST(SSMatSupportT*, p);
			//if (!fSSMatSupport)
			//	ExceptionT::GeneralFail("SCNIMFT::NewMaterialList","Cannot cast to SSMatSupport\n");
		}

		if (nsd == 1)
			return new SolidMatList1DT(size, *fSSMatSupport);
		else if (nsd == 2)
			return new SolidMatList2DT(size, *fSSMatSupport);
		else if (nsd == 3)
			return new SolidMatList3DT(size, *fSSMatSupport);
		else
			return NULL;
	}
	else
	{
		if (nsd == 1)
			return new SolidMatList1DT;
		else if (nsd == 2)
			return new SolidMatList2DT;
		else if (nsd == 3)
			return new SolidMatList3DT;
		else
			return NULL;
	}	
}
	
void SCNIMFT::bVectorToMatrix(double *bVector, dMatrixT& BJ)
{
#if __option(extended_errorcheck)
	if (BJ.MajorDim() != fSD*(fSD+1)/2) 
		ExceptionT::SizeMismatch("SCNIMFT::bVectorToMatrix","Matrix has bad majorDim");
	if (BJ.MinorDim() != fSD) 
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

void SCNIMFT::ComputeBMatrices(void)
{
	/* possible best implementation is to loop over all Delone edge
	 * centroids and compute all the necessary values only once per Voronoi
	 * facet. This approach minimizes number of times that the support of
	 * an arbitrary point in space (the Voronoi facet centroid) has to be
	 * found.
	 */

	const RaggedArray2DT<int>& nodeSupport = fNodalShapes->NodeNeighbors();
	int nNodes = fNodes.Length();
	
	/* allocate space for strain smoothing workspace */
	bVectors.Dimension(nNodes);
	for (int i = 0; i < nNodes; i++)
	{		
		/* allocate storage for a vector for each node in the support of node i */
		bVectors[i].Dimension(nodeSupport.MinorDim(i), fSD);
		bVectors[i] = 0.;
	}	
	
	dArrayT facetCentroid(fSD), facetNormal(fSD), facetIntegral(fSD);
	double* currentB, *currentI;
	int n_0, n_1, l_0, l_1;
	bool integralIsCurrent;
	for (int i = 0; i < fDeloneEdges.MajorDim(); i++)
	{
		facetCentroid = 0.; 
		n_0 = fDeloneEdges(i,0);
		n_1 = fDeloneEdges(i,1);
		if (i < nInteriorDeloneEdges)
		{
			facetCentroid.SumOf(fDeloneVertices(n_0),fDeloneVertices(n_1));
			facetCentroid /= 2.;
		}
		else
			facetCentroid.Set(fSD, fBoundaryDeloneCentroids(i-nInteriorDeloneEdges)); 
		facetNormal.DiffOf(fDeloneVertices(n_1), fDeloneVertices(n_0));
		facetNormal.UnitVector();
		
		if (!fNodalShapes->SetFieldAt(facetCentroid, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
		
		iArrayT supp_centroid(fNodalShapes->Neighbors());
		iArrayT supp_centroid_key(supp_centroid.Length());
		supp_centroid_key.SetValueToPosition();
		supp_centroid_key.SortAscending(supp_centroid);
		int n_supp_centroid = supp_centroid.Length();
		
		iArrayT supp_0(nodeSupport.MinorDim(n_0));
		supp_0.Copy(nodeSupport(n_0));
		iArrayT supp_0_key(supp_0.Length());
		supp_0_key.SetValueToPosition();
		supp_0_key.SortAscending(supp_0);
		l_0 = supp_0.Length();
		
		iArrayT supp_1(nodeSupport.MinorDim(n_1));
		supp_1.Copy(nodeSupport(n_1));
		iArrayT supp_1_key(supp_1.Length());
		supp_1_key.SetValueToPosition();
		supp_1_key.SortAscending(supp_1);
		l_1 = supp_1.Length();
		
		/* Simultaneously loop over support of the two nodes that are endpoints of the
		 * current Delone edge and the nodes in the support of the midpoint of this
		 * edge. They're sorted to make this loop march forward until there are no
		 * more integrals over facets to compute. (i.e. until the support the edge
		 * midpoint is exhausted---then the rest of the contributions are zero.
		 */
		int* c = supp_centroid.Pointer();
		int* c_i = supp_centroid_key.Pointer();
		int* s_0 = supp_0.Pointer();
		int* s_1 = supp_1.Pointer();
		int ctr_i, ctr_0, ctr_1;
		ctr_i = ctr_0 = ctr_1 = 0;
		while (ctr_i < n_supp_centroid && (ctr_0 < l_0 || ctr_1 < l_1))
		{ 
			bool okToAdvance = true;	
			while (okToAdvance && ctr_i < n_supp_centroid)
			{
				if (ctr_0 < l_0)  // still comparing against node 0
					if (*c >= *s_0) 
						okToAdvance = false; 
				if (okToAdvance && ctr_1 < l_1)  // still comparing against node 1
					if (*c >= *s_1) 
						okToAdvance = false;
				if (okToAdvance)
				{
					c++;
					c_i++;
					ctr_i++;
				}
			}
			if (ctr_i != n_supp_centroid)	
			{
				integralIsCurrent = false;
				if (ctr_0 < l_0 && *c == *s_0) // *c is in support of centroid and support of Delone vertex 0
				{
					facetIntegral = facetNormal;
					facetIntegral *= fDualAreas[i]*phiValues[*c_i];		
					integralIsCurrent = true;
					currentI = facetIntegral.Pointer();
					currentB = bVectors[n_0](supp_0_key[ctr_0]);
					for (int j = 0; j < fSD; j++)
						*currentB++ += *currentI++;
				}
				if (ctr_0 != l_0 && *s_0 <= *c)
				{
					ctr_0++;
					if (ctr_0 < l_0)
						s_0++;
				}
				if (ctr_1 < l_1 && *c == *s_1)
				{
					if (!integralIsCurrent) // facetIntegral hasn't been evaluated this iteration
					{
						facetIntegral = facetNormal;
						facetIntegral *= fDualAreas[i]*phiValues[*c_i];
					}
					currentI = facetIntegral.Pointer();
					currentB = bVectors[n_1](supp_1_key[ctr_1]);
					for (int j = 0; j < fSD; j++)
						*currentB++ -= *currentI++; //NB change in sign; facet normal is inverted!
				}
				if (ctr_1 != l_1 && *s_1 <= *c)
				{
					ctr_1++;
					if (ctr_1 < l_1)
						s_1++;
				}
			}
		}
	}
	
	/** Loop over remaining edges */
	for (int i = 0; i < fNonDeloneEdges.Length(); i++)
	{
		facetCentroid = 0.; 
		n_0 = fNonDeloneEdges[i];
		facetCentroid.Set(fSD, fNonDeloneCentroids(i)); 
		facetNormal.Set(fSD, fNonDeloneNormals(i));
		facetNormal.UnitVector();
		
		if (!fNodalShapes->SetFieldAt(facetCentroid, NULL)) // shift = 0 or not ?
			ExceptionT::GeneralFail("SCNIMFT::ComputeBMatrices","Shape Function evaluation"
				"failed at Delone edge %d\n",i);
				
		const dArrayT& phiValues = fNodalShapes->FieldAt();			
				
		iArrayT supp_centroid(fNodalShapes->Neighbors());
		iArrayT supp_centroid_key(supp_centroid.Length());
		supp_centroid_key.SetValueToPosition();
		supp_centroid_key.SortAscending(supp_centroid);
		int n_supp_centroid = supp_centroid.Length();
		
		iArrayT supp_0(nodeSupport.MinorDim(n_0));
		supp_0.Copy(nodeSupport(n_0));
		iArrayT supp_0_key(supp_0.Length());
		supp_0_key.SetValueToPosition();
		supp_0_key.SortAscending(supp_0);
		l_0 = supp_0.Length();
		
		/* Simultaneously loop over support of the node and support of the integration 
		 * point. They may be the same. Does that matter?
		 */
		int* c = supp_centroid.Pointer();
		int* c_i = supp_centroid_key.Pointer();
		int* s_0 = supp_0.Pointer();
		int ctr_i, ctr_0;
		ctr_i = ctr_0 = 0;
		while (ctr_i < n_supp_centroid && ctr_0 < l_0)
		{ 
			bool okToAdvance = true;	
			while (okToAdvance && ctr_i < n_supp_centroid)
			{
				if (ctr_0 < l_0)  // still comparing against node 0
					if (*c >= *s_0) 
						okToAdvance = false; 
				if (okToAdvance)
				{
					c++;
					c_i++;
					ctr_i++;
				}
			}
			if (ctr_i != n_supp_centroid)	
			{
				if (ctr_0 < l_0 && *c == *s_0) // *c is in support of centroid and support of Delone vertex 0
				{
					facetIntegral = facetNormal;
					facetIntegral *= fIntegrationWeights[i]*phiValues[*c_i];		
					integralIsCurrent = true;
					currentI = facetIntegral.Pointer();
					currentB = bVectors[n_0](supp_0_key[ctr_0]);
					for (int j = 0; j < fSD; j++)
						*currentB++ += *currentI++;
				}
				if (ctr_0 != l_0 && *s_0 <= *c)
				{
					ctr_0++;
					if (ctr_0 < l_0)
						s_0++;
				}
			}
		}
	}
	
	// scale integrals by volumes of Voronoi cells
	for (int i = 0; i < fNodes.Length(); i++)
		bVectors[i] *= 1./fVoronoiCellVolumes[i];

}

void SCNIMFT::VoronoiDiagramToFile(ofstreamT& vout)
{
	
    int nVertices = fVoronoiVertices.MajorDim();
    int nSD = fVoronoiVertices.MinorDim();

    vout << fNodes.Length() << "\n";
    vout << nSD << "\n";
    vout << nVertices << "\n";

    // write out vertices of the VoronoiDiagram
    for (int i = 0; i < nVertices; i++)
    {
      	for (int j = 0; j < nSD; j++)
			vout << fVoronoiVertices(i,j) << " "; 
      	vout << "\n";
    }
    
    // write out Voronoi cells for each node
    for (int i = 0; i < fNodes.Length(); i++)
    {
		vout << i <<"\n";

		// number of vertices on its boundary
		vout << fVoronoiCells[i].Length() << " ";
		for (int j = 0; j < fVoronoiCells[i].Length(); j++)
	    	vout << fVoronoiCells[i][j] << " ";
		vout << "\n";

		// number of facets
		vout << fVoronoiFacetIndices[i].Length() << "\n";
		
		// data for each facet
		for (int j = 0; j < fVoronoiFacetIndices[i].Length(); j++)
	  	{
	   		vout << j << " " << fVoronoiFacetIndices[i][j].Length() << " ";
	    	for (int k = 0; k < fVoronoiFacetIndices[i][j].Length(); k++)
	      		vout << fVoronoiFacetIndices[i][j][k] << " ";
	    
	    	// area of facet
	    	vout << fVoronoiFacetAreas[i][j] << " ";
	    
	    	// facet normal vector
	    	for (int k = 0; k < nSD; k++)
	      		vout << fVoronoiFacetNormals[i](j,k) << " ";
	    	vout << "\n";
	  	}
	
		// cell volume
		vout << fVoronoiCellVolumes[i] << "\n";
    }
}	
	
void SCNIMFT::VoronoiDiagramFromFile(ifstreamT& vin)
{
	const char caller[] = "SCNIMFT::VoronoiDiagramFromFile";	

    int nNodes, nSD, nVertices;
    vin >> nNodes;
    
    /* minor consistency checks */
    if (nNodes != fNodes.Length())
    {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match node number\n");
    }
    
    vin >> nSD;
    if (nSD != 2 )
    {
		vin.close();
		ExceptionT::GeneralFail(caller,"Input Voronoi file does not match SD\n");
    }
    
    vin >> nVertices;
    
    // allocate memory for Voronoi diagram data structures
    fVoronoiVertices.Dimension(nVertices, nSD);
    fVoronoiCells.Dimension(nNodes);
    fVoronoiFacetIndices.Dimension(nNodes);
    fVoronoiCellVolumes.Dimension(nNodes);
    fVoronoiFacetAreas.Dimension(nNodes);
    fVoronoiFacetNormals.Dimension(nNodes);

    for (int i = 0 ; i < nVertices; i++)
    	for (int j = 0; j < nSD; j++)
	  		vin >> fVoronoiVertices(i,j);

    for (int i = 0; i < nNodes; i++)
    {
		int itmp;
		vin >> itmp;
		if (itmp != i)
	  		ExceptionT::GeneralFail(caller,"Bad Input Voronoi file\n");

		// read in number of vertices
		vin >> itmp;
		fVoronoiCells[i].Dimension(itmp);
		for (int j = 0; j < itmp; j++)
	  		vin >> fVoronoiCells[i][j];
	  		
		// number of facets
		int nFacets;
		vin >> nFacets;
		fVoronoiFacetIndices[i].Dimension(nFacets);
		fVoronoiFacetAreas[i].Dimension(nFacets);
		fVoronoiFacetNormals[i].Dimension(nFacets, nSD);
		for (int j = 0; j < nFacets; j++)
	  	{
	    	vin >> itmp;
	    	if (itmp != j)
	      		ExceptionT::GeneralFail(caller,"Bad Input Voronoi File\n");
	  
	    	// number of vertices in facet
	    	vin >> itmp;
	    	fVoronoiFacetIndices[i][j].Dimension(itmp);
	    
	    	for (int k = 0; k < itmp; k++)
				vin >> fVoronoiFacetIndices[i][j][k];
	    	
	    	vin >> fVoronoiFacetAreas[i][j];
	    	
	    	for (int k = 0; k < nSD; k++)
	      		vin >> fVoronoiFacetNormals[i](j,k);
        }
	
		vin >> fVoronoiCellVolumes[i];

	}
	
	/* Read in all DeloneEdges, but store the number temporarily in nInteriorDeloneEdges */
	vin >> nInteriorDeloneEdges; 
	fDeloneEdges.Dimension(nInteriorDeloneEdges, 2);
	fDualAreas.Dimension(nInteriorDeloneEdges);
	
	for (int i = 0; i < nInteriorDeloneEdges; i++)
		vin >> fDeloneEdges(i,0) >> fDeloneEdges(i,1);
		
	for (int i = 0; i < nInteriorDeloneEdges; i++)
		vin >> fDualAreas[i];
		
	int nCentroids;
	vin >> nCentroids;
	
	fBoundaryDeloneCentroids.Dimension(nCentroids, fSD);
	nInteriorDeloneEdges -= nCentroids; // Compute number of facets whose centroids are known
	
	for (int i = 0; i < nCentroids; i++)
		for (int j = 0; j < fSD; j++)
			vin >> fBoundaryDeloneCentroids(i,j);
			
	vin >> nCentroids; // number of boundary facets
	fNonDeloneEdges.Dimension(nCentroids);
	fNonDeloneNormals.Dimension(nCentroids,fSD);
	fNonDeloneCentroids.Dimension(nCentroids,fSD);
	fIntegrationWeights.Dimension(nCentroids);
	
	for (int i = 0; i < nCentroids; i++)
		vin >> fNonDeloneEdges[i];	

	for (int i = 0; i < nCentroids; i++)
		for (int j = 0; j < fSD; j++)
			vin >> fNonDeloneCentroids(i,j);
			
	for (int i = 0; i < nCentroids; i++)
		for (int j = 0; j < fSD; j++)
			vin >> fNonDeloneNormals(i,j);
			
	for (int i = 0; i < nCentroids; i++)
		vin >> fIntegrationWeights[i];
}

// XML stuff below

/* describe the parameters needed by the interface */
void SCNIMFT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

}

/* information about subordinate parameter lists */
void SCNIMFT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

}

/* return the description of the given inline subordinate parameter list */
void SCNIMFT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	/* inherited */
	ElementBaseT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SCNIMFT::NewSub(const StringT& list_name) const
{
	/* inherited */
	return ElementBaseT::NewSub(list_name);
}
