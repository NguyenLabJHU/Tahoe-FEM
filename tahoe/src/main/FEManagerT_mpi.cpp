/* $Id: FEManagerT_mpi.cpp,v 1.35.2.5 2004-07-13 16:42:41 paklein Exp $ */
/* created: paklein (01/12/2000) */
#include "FEManagerT_mpi.h"
#include <time.h>

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "NodeManagerT.h"
#include "ElementBaseT.h"
#include "IOManager_mpi.h"
#include "GraphT.h"
#include "IOBaseT.h"
#include "PartitionT.h"
#include "CommunicatorT.h"
#include "CommManagerT.h"

#include "ModelFileT.h"
#include "ExodusT.h"

using namespace Tahoe;

/* constructor */
FEManagerT_mpi::FEManagerT_mpi(const StringT& input, ofstreamT& output, 
	CommunicatorT& comm, const ArrayT<StringT>& argv, PartitionT* partition, TaskT task):
	FEManagerT(input, output, comm, argv),
	fPartition(partition),
	fTask(task),
	fExternIOManager(NULL)
{	
	if (fTask == kRun)
	{
		const char caller[] = "FEManagerT_mpi::FEManagerT_mpi";

		/* revise input file name */
		StringT suffix;
		suffix.Suffix(fInputFile);
		fInputFile.Root();
		fInputFile.Append(".p", Rank());
		fInputFile.Append(suffix);

		/* checks */
		if (!fPartition)
			ExceptionT::BadInputValue(caller, "partition information required for task %d", kRun);
		else if (fPartition->ID() != Rank())
			ExceptionT::MPIFail(caller, "partition ID %d does not match process rank %d",
				fPartition->ID(), Rank());
		
		StringT log_file;
		log_file.Root(input);
		//log_file.Append(".p", Rank());
		log_file.Append(".log");
		flog.open(log_file);
		
		/* redirect log messages */
		fComm.SetLog(flog);

		/* log */
		TimeStamp(caller);
	}
	else if (fTask == kDecompose)
		fInitCode = kAllButSolver;
}

/* destructor */
FEManagerT_mpi::~FEManagerT_mpi(void)
{
	/* log */
	TimeStamp("FEManagerT_mpi::~FEManagerT_mpi");

	/* restore log messages */
	if (fTask == kRun) fComm.SetLog(cout);
}

ExceptionT::CodeT FEManagerT_mpi::InitStep(void)
{
	/* set default output time stamp */
	if (fExternIOManager) fExternIOManager->SetOutputTime(Time());

	/* inherited */
	ExceptionT::CodeT error = FEManagerT::InitStep();
	if (error != ExceptionT::kNoError) {
		cout << "\n FEManagerT_mpi::InitStep: error: " << error << endl;
	}
	return error;
}

ExceptionT::CodeT FEManagerT_mpi::SolveStep(void)
{
	/* inherited */
	ExceptionT::CodeT error = FEManagerT::SolveStep();
	if (error != ExceptionT::kNoError) {
		cout << "\n FEManagerT_mpi::SolveStep: return: " << error << endl;
	}
	return error;
}

/* system relaxation */
GlobalT::RelaxCodeT FEManagerT_mpi::RelaxSystem(int group) const
{
	/* inherited */
	GlobalT::RelaxCodeT relax = FEManagerT::RelaxSystem(group);

	/* gather codes */
	iArrayT all_relax(Size());
	fComm.AllGather(relax, all_relax);

	/* code precedence */
	for (int i = 0; i < all_relax.Length(); i++)
		relax = GlobalT::MaxPrecedence(relax, GlobalT::RelaxCodeT(all_relax[i]));
	
	/* report */
	if (relax != GlobalT::kNoRelax)
	{
		cout << "\n Relaxation code at time = " << Time() << '\n';
		cout << setw(kIntWidth) << "proc";	
		cout << setw(kIntWidth) << "code" << '\n';	
		for (int i = 0; i < all_relax.Length(); i++)
		{
			cout << setw(kIntWidth) << i;	
			cout << setw(kIntWidth) << all_relax[i];
			cout << '\n';	
		}
	}

	return relax;
}

/* update solution */
void FEManagerT_mpi::Update(int group, const dArrayT& update)
{
	/* give a heart beat */
	const char caller[] = "FEManagerT_mpi::Update";

	/* give heartbeat */
	TimeStamp(caller);
	
	/* check sum */
	if (fComm.Sum(ExceptionT::kNoError) != 0) 
		ExceptionT::BadHeartBeat(caller); /* must trigger try block in FEManagerT::SolveStep */

	/* inherited */
	FEManagerT::Update(group, update);
}

/* initiate the process of writing output from all output sets */
void FEManagerT_mpi::WriteOutput(double time)
{
	/* set output time for the external IO manager */
	if (fExternIOManager) fExternIOManager->SetOutputTime(time);

	/* inherited */
	FEManagerT::WriteOutput(time);
}
	
void FEManagerT_mpi::WriteOutput(int ID, const dArray2DT& n_values,
	const dArray2DT& e_values) const
{
	/* output assembly mode */
	if (!fExternIOManager)
		/* do local IO */
		FEManagerT::WriteOutput(ID, n_values, e_values);
	else
		/* distribute/assemble/write */
		fExternIOManager->WriteOutput(ID, n_values, e_values);
}

/* (temporarily) direct output away from main out */
void FEManagerT_mpi::DivertOutput(const StringT& outfile)
{
	/* do local IO */
	if (!fExternIOManager)
	{
		/* need processor designation for split output */
		StringT outfile_p = outfile;
		if (Size() > 1) outfile_p.Append(".p", Rank());

		FEManagerT::DivertOutput(outfile_p);
	}
	else 
	{
		/* external I/O */
		fExternIOManager->DivertOutput(outfile);	

		/* system output */
		if (fCurrentGroup != -1) /* resolved group */
			fSO_DivertOutput[fCurrentGroup]	= true;
		else /* all groups */
			fSO_DivertOutput = true;
	}
}

void FEManagerT_mpi::RestoreOutput(void)
{
	/* do local IO */
	if (!fExternIOManager)
		FEManagerT::RestoreOutput();
	else
	{
		/* external I/O */
		fExternIOManager->RestoreOutput();

		/* system output */
		if (fCurrentGroup != -1) /* resolved group */
			fSO_DivertOutput[fCurrentGroup]	= false;
		else /* all groups */
			fSO_DivertOutput = false;
	}
}

/* domain decomposition */
void FEManagerT_mpi::Decompose(ArrayT<PartitionT>& partition, GraphT& graphU,
	bool verbose, int method)
{
	const char caller[] = "FEManagerT_mpi::Decompose";

	//TEMP
	//TimeStamp("FEManagerT_mpi::Decompose");

	/* check */
	if (partition.Length() == 1)
		ExceptionT::GeneralFail(caller, "expecting more than 1 partition");

	/* geometry file must be ascii external */
	if (fInputFormat != IOBaseT::kTahoeII && fInputFormat != IOBaseT::kExodusII)
		ExceptionT::BadInputValue(caller, "expecting file format %d or %d, not %d",
			IOBaseT::kTahoeII, IOBaseT::kExodusII, fInputFormat);

	/* decomposition method */
	bool use_new_methods = false; //TEMP
	bool dual_graph = (InterpolantDOFs() == 0);
	if (dual_graph && use_new_methods)
		DoDecompose_2(partition, graphU, verbose, method);
	else
		DoDecompose_1(partition, graphU, verbose, method);

	if (fInputFormat == IOBaseT::kTahoeII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ModelFileT model_ALL;
			model_ALL.OpenRead(fModelFile);
			
			/* set number of element sets */
			iArrayT elementID;
			if (model_ALL.GetElementSetID(elementID) != ModelFileT::kOK) ExceptionT::GeneralFail(caller);
			ArrayT<StringT> IDlist(elementID.Length());
			for (int i = 0; i < IDlist.Length(); i++)
				IDlist[i].Append(elementID[i]);
			partition[j].InitElementBlocks(IDlist);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get element set */
				iArray2DT set;
				if (model_ALL.GetElementSet(elementID[i], set) != ModelFileT::kOK)
					ExceptionT::GeneralFail(caller);
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(IDlist[i], set);
			}
		}
	}
	else if (fInputFormat == IOBaseT::kExodusII)
	{
		/* label element sets in partition data */
		for (int j = 0; j < partition.Length(); j++)
		{
			/* original model file */
			ExodusT model_ALL(cout);
			model_ALL.OpenRead(fModelFile);
			
			/* set number of element sets */
			iArrayT elementID(model_ALL.NumElementBlocks());
			model_ALL.ElementBlockID(elementID);
			ArrayT<StringT> IDlist(elementID.Length());
			for (int i = 0; i < IDlist.Length(); i++)
				IDlist[i].Append(elementID[i]);
			partition[j].InitElementBlocks(IDlist);	
			for (int i = 0; i < elementID.Length(); i++)
			{
				/* get set dimensions */
				int num_elems;
				int num_elem_nodes;
				model_ALL.ReadElementBlockDims(elementID[i], num_elems, num_elem_nodes);

				/* get element set */
				iArray2DT set(num_elems, num_elem_nodes);
				GeometryT::CodeT geometry_code;
				model_ALL.ReadConnectivities(elementID[i], geometry_code, set);
					
				/* correct node numbering offset */
				set--;	
				
				/* set partition */
				partition[j].SetElements(IDlist[i], set);
			}
		}
	}
	else ExceptionT::GeneralFail(caller);
}

/* accept parameter list */
void FEManagerT_mpi::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FEManagerT::TakeParameterList(list);

	/* collect model file and input format from ModelManager */
	fInputFormat = fModelManager->DatabaseFormat();
	fModelFile = fModelManager->DatabaseName();

	/* correct restart file name */
	if (fReadRestart) {
		StringT suffix;
		suffix.Suffix(fRestartFile);
		fRestartFile.Root();
		fRestartFile.Append(".p", Rank());
		fRestartFile.Append(suffix);
	}
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* (re-)set system to initial conditions */
ExceptionT::CodeT FEManagerT_mpi::InitialCondition(void)
{
	/* inherited */
	ExceptionT::CodeT error = FEManagerT::InitialCondition();
	
	/* set I/O */
	if (error == ExceptionT::kNoError && fExternIOManager) 
		fExternIOManager->NextTimeSequence(SequenceNumber());
		
	return error;
}

/* global number of first local equation */
int FEManagerT_mpi::GetGlobalEquationStart(int group, int start_eq_shift) const
{
	if (Size() == 1)
		return FEManagerT::GetGlobalEquationStart(group, start_eq_shift);
	else
	{
		/* number of local equations */
		int num_eq = fNodeManager->NumEquations(group);

		/* collect from all */
		int size = Size();
		iArrayT all_num_eq(size);
		fComm.AllGather(num_eq, all_num_eq);

		/* compute offset to local equations */
		int offset = 0;
		for (int i = 0; i < Rank(); i++)
			offset += all_num_eq[i];

		/* equation start */
		return offset + 1 + start_eq_shift;
	}
}

int FEManagerT_mpi::GetGlobalNumEquations(int group) const
{
	if (Size() == 1)
		return FEManagerT::GetGlobalNumEquations(group);
	else
	{
		int loc_num_eq = fNodeManager->NumEquations(group);
		iArrayT all_num_eq(Size());
		fComm.AllGather(loc_num_eq, all_num_eq);
		return all_num_eq.Sum();
	}
}

/* construct a new CommManagerT */
CommManagerT* FEManagerT_mpi::New_CommManager(void) const
{
	/* inherited */
	CommManagerT* comm_man = FEManagerT::New_CommManager();
	if (!comm_man) ExceptionT::GeneralFail();

	/* set the partition data */
	comm_man->SetPartition(fPartition);
	return comm_man;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* collect computation effort for each node */
void FEManagerT_mpi::WeightNodalCost(iArrayT& weight) const
{
	weight.Dimension(fNodeManager->NumNodes());
	weight = 1;
	fNodeManager->WeightNodalCost(weight);
	for (int i = 0 ; i < fElementGroups->Length(); i++)
		(*fElementGroups)[i]->WeightNodalCost(weight);
}

/* write time stamp to log file */
void FEManagerT_mpi::TimeStamp(const char* message) const
{
	/* log */
	fComm.Log(CommunicatorT::kUrgent, message);
}

/* decomposition methods */
void FEManagerT_mpi::DoDecompose_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int
	method)
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;

	/* collect element groups */
	for (int s = 0 ; s < fElementGroups->Length(); s++)
		(*fElementGroups)[s]->ConnectsU(connects_1, connects_2);		
	
	/* dual graph partitioning graph */
	AutoArrayT<const iArray2DT*> connectsX_1;
	
	/* collect minimal connects */
	for (int s = 0 ; s < fElementGroups->Length(); s++)
		(*fElementGroups)[s]->ConnectsX(connectsX_1);

	/* initialize graph */
	GraphT& graphX = graph;
	for (int r = 0; r < connectsX_1.Length(); r++)
		graphX.AddGroup(*(connectsX_1[r]));
		
	/* make graph */
	if (verbose) cout << " FEManagerT_mpi::DoDecompose_2: constructing dual graph" << endl;
	clock_t t0 = clock();		
	graphX.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
	     << " sec: FEManagerT_mpi::DoDecompose_2: construct graph" << endl;
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	WeightNodalCost(weight);
	graphX.Partition(config, weight, connects_1, connects_2, partition, true, method);
}

void FEManagerT_mpi::DoDecompose_1(ArrayT<PartitionT>& partition, GraphT& graph, 
	bool verbose, int method)
{
	/* connectivities for partititioning */
	AutoArrayT<const iArray2DT*> connects_1;
	AutoArrayT<const RaggedArray2DT<int>*> connects_2;
	AutoArrayT<const iArray2DT*> equivalent_nodes;

	/* collect connectivies from all solver groups */
	for (int i = 0; i < NumGroups(); i++)
		fNodeManager->ConnectsU(i,connects_1,connects_2, equivalent_nodes);
	
	/* collect element groups */
	for (int s = 0 ; s < fElementGroups->Length(); s++)
		(*fElementGroups)[s]->ConnectsU(connects_1, connects_2);
		

	/* initialize graph */
	GraphT& graphU = graph;
	for (int r = 0; r < connects_1.Length(); r++)
		graphU.AddGroup(*(connects_1[r])); 

	for (int k = 0; k < connects_2.Length(); k++)
		graphU.AddGroup(*(connects_2[k]));

	for (int k = 0; k < equivalent_nodes.Length(); k++)
		graphU.AddEquivalentNodes(*(equivalent_nodes[k]));
		
	/* make graph */
	clock_t t0 = clock();
	if (verbose) cout << " FEManagerT_mpi::DoDecompose_1: constructing graph" << endl;
	graphU.MakeGraph();
	clock_t t1 = clock();
	if (verbose)
		cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: FEManagerT_mpi::DoDecompose_1: construct graph" << endl;
	
	/* dual graph partitioning graph */
	int dual_graph = (InterpolantDOFs() == 0) ? 1 : 0;
	AutoArrayT<const iArray2DT*> connectsX_1;
	GraphT graphX;
	if (dual_graph == 1)
	{
		if (verbose) cout << " FEManagerT_mpi::DoDecompose_1: constructing dual graph" << endl;
		
		/* collect element groups */
		for (int s = 0 ; s < fElementGroups->Length(); s++)
			(*fElementGroups)[s]->ConnectsX(connectsX_1);

		/* initialize graph */
		for (int r = 0; r < connectsX_1.Length(); r++)
			graphX.AddGroup(*(connectsX_1[r]));
		
		/* make graph */
		clock_t t0 = clock();		
		graphX.MakeGraph();
		clock_t t1 = clock();
		if (verbose)
			cout << setw(kDoubleWidth) << double(t1 - t0)/CLOCKS_PER_SEC
		     << " sec: FEManagerT_mpi::DoDecompose_1: construct X graph" << endl;
	}
	
	/* generate partition */
	iArrayT config(1); //TEMP - will be scalar soon?
	config[0] = partition.Length();	
	iArrayT weight;
	WeightNodalCost(weight);
	if (dual_graph == 1)
		graphX.Partition(config, weight, graphU, partition, true, method);
	else
		graphU.Partition(config, weight, partition, true, method);
}
