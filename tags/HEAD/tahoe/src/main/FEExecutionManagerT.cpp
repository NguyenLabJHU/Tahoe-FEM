/* $Id: FEExecutionManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (09/21/1997)                                          */

#include "FEExecutionManagerT.h"

#include <iostream.h>
#include <time.h>

#include "fstreamT.h"

#include "Environment.h"
#include "Constants.h"
#include "ExceptionCodes.h"

#if defined(__MWERKS__) && __option(profile)
#include <profiler.h>
#endif

#ifdef __MPI__
#include "mpi.h"
#endif

#include "fstreamT.h"
#include "FEManagerT.h"
#include "FEManagerT_mpi.h"
#include "IOManager_mpi.h"
#include "StringT.h"
#include "GraphT.h"
#include "PartitionT.h"
#include "OutputSetT.h"

#include "ModelFileT.h"
#include "ExodusT.h"
#include "JoinOutputT.h"

/* Constructor */
FEExecutionManagerT::FEExecutionManagerT(int argc, char* argv[], char job_char,
	char batch_char):
	ExecutionManagerT(argc, argv, job_char, batch_char, 0)
{
#ifdef __CPLANT__ // not grouped output due to limited memory
	AddCommandLineOption("-split_io");
#endif

//TEMP - decomp testing
#ifdef __MACOS__
	bool do_decomp = false;
	if (do_decomp) AddCommandLineOption("-decomp");

	bool do_split_io = false;
	if (do_split_io) AddCommandLineOption("-split_io");
#endif
}

/* Prompt input files until "quit" */
void FEExecutionManagerT::Run(void)
{
	int size = 1;
	int rank = 0;
#ifdef __MPI__
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) throw eMPIFail;
#endif

	/* check command line arguments */
	if (size > 1)
	{
		for (int i = 0; i < fCommandLineOptions.Length(); i++)
			if (fCommandLineOptions[i] == "-decomp")
			{
				if (rank == 0)
					cout << "\n ::Run: run \"-decomp\" on one processor only" << endl;
				return;
			}
			else if (fCommandLineOptions[i] == "-join")
			{
				if (rank == 0)
					cout << "\n ::Run: run \"-join\" on one processor only" << endl;
				return;
			}
	}

	/* inherited on fall-through */
	 ExecutionManagerT::Run();
}

/**********************************************************************
* Protected
**********************************************************************/

/* MUST be overloaded */
void FEExecutionManagerT::RunJob(ifstreamT& in, ostream& status)
{
	/* run serial by default */
	int run_option = 0;

	/* check command line options */
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
	{
		if (fCommandLineOptions[i] == "-decomp")
			run_option = 3;
		else if (fCommandLineOptions[i] == "-join")
			run_option = 4;
	}

#ifdef __MPI__
	int size;
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;
	if (size > 1) run_option = 1;
#endif

	switch (run_option)
	{
		case 0:
			cout << "\n RunJob_serial: " << in.filename() << endl;
			RunJob_serial(in, status);
			break;
		case 1:
			cout << "\n RunJob_parallel: " << in.filename() << endl;
			RunJob_parallel(in, status);
			break;
		case 3:
			cout << "\n RunDecomp_serial: " << in.filename() << endl;
			RunDecomp_serial(in, status);
			break;	
		case 4:
			cout << "\n RunJoin_serial: " << in.filename() << endl;
			RunJoin_serial(in, status);
			break;	
		default:
			status << "\n FEExecutionManagerT::RunJob: unknown option:"
			     << run_option << endl;
			throw eGeneralFail;
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* standard serial driver */
void FEExecutionManagerT::RunJob_serial(ifstreamT& in,
	ostream& status) const
{
	/* open output file stream */
	StringT outfilename;
	ofstreamT out;	
	out.open(outfilename.DefaultName(in.filename()));

#ifdef __MWERKS__
	if (!out.is_open())
	{
		cout << "\n FEExecutionManagerT::RunJob_serial: could not open file: "
		     << outfilename << endl;
		return;
	}
#endif

	clock_t t0 = 0, t1 = 0, t2 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);

	int phase; // job phase
	try
	{
		t0 = clock();

		/* construction */
		phase = 0;
		in.set_marker('#');
		FEManagerT analysis1(in, out);
		analysis1.Initialize();

		t1 = clock();
		
		/* solution */
//ProfilerSetStatus(1);
		phase = 1;
		analysis1.Solve();
//ProfilerSetStatus(0);

		t2 = clock();
	}

	/* job failure */
	catch (int code)
	{
		status << "\n \"" << in.filename() << "\" exit on exception during the\n";
		if (phase == 0)
		{
			status << " construction phase. Check the input file for errors." << endl;
		
			/* echo some lines from input */
			if (code == eBadInputValue) Rewind(in, status);
		}
		else
		{
			status << " solution phase. See \"" << outfilename << "\" for a list";
			status << " of the codes.\n";
		}
		
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		

		out << endl;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: "   << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status << "    Stop time: " << ctime(&stoptime);
out   << "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}	

/* generate decomposition files */
void FEExecutionManagerT::RunDecomp_serial(ifstreamT& in, ostream& status) const
{
	/* prompt for decomp size */
	int size = 0;
	int count = 0;
	while (count == 0 || (count < 5 && size < 2))
	{
		count++;
		cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __MPI__)
		cout << '\n';
#endif					
		cin >> size;
		if (size == 0) break;
	}
	if (size < 2) return;

	/* set stream comment marker */
	in.set_marker('#');

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();
	
		/* get the model file name */
		StringT model_file, suffix;
		IOBaseT::FileTypeT format;
		GetModelFile(in, model_file, format);
	
		/* global output model file */
		StringT global_model_file;
		suffix.Suffix(model_file);
		global_model_file.Root(model_file);
		global_model_file.Append(".io", suffix);

		/* output map file */
		StringT map_file;
		map_file.Root(model_file);
		map_file.Append(".n", size);
		map_file.Append(".io.map");

		/* set output map and and generate decomposition */
		Decompose(in, size, model_file, global_model_file, format, map_file);
		t1 = clock();
	}

	/* job failure */
	catch (int code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
	status << "\n End Execution\n" << endl;
}

/* join parallel results files */
void FEExecutionManagerT::RunJoin_serial(ifstreamT& in, ostream& status) const
{
	/* set stream comment marker */
	in.set_marker('#');

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();
	
		/* get the model file name */
		StringT model_file, suffix;
		IOBaseT::FileTypeT format;
		GetModelFile(in, model_file, format);
	
		/* global output model file */
		StringT global_model_file;
		suffix.Suffix(model_file);
		global_model_file.Root(model_file);
		global_model_file.Append(".io", suffix);

		/* prompt for decomp size */
		int size;
		cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __MPI__)
		cout << '\n';
#endif					
		cin >> size;
		if (size < 2) return;

		/* construct joiner */
		JoinOutputT output_joiner(in, model_file, global_model_file, format, size);
		
		/* join files */
		output_joiner.Join();

		t1 = clock();
	}

	/* job failure */
	catch (int code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "         Join: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
	status << "\n End Execution\n" << endl;
}

/* testing for distributed execution */
#ifndef __MPI__
void FEExecutionManagerT::RunJob_parallel(ifstreamT& in, ostream& status) const
{
#pragma unused (in)
	status << "\n FEExecutionManagerT::RunJob_parallel: no mpi" << endl;
	throw;
}
#else
void FEExecutionManagerT::RunJob_parallel(ifstreamT& in, ostream& status) const
{
	/* set stream comment marker */
	in.set_marker('#');

	/* get rank and size */
	int rank, size;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS ||
	    MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) throw eMPIFail;

	/* time markers */
	clock_t t0 = 0, t1 = 0, t2 = 0, t3 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);

	/* output stream */
	StringT out_file;
	out_file.Root(in.filename());
	out_file.Append(".p", rank);
	out_file.Append(".out");
	ofstreamT out;
	out.open(out_file);

	int token; // for run time check sums
	int check_sum;

	try {
	t0 = clock();
	
	/* get the model file name */
	StringT model_file, suffix;
	IOBaseT::FileTypeT format;
	GetModelFile(in, model_file, format);
	
	/* global output model file */
	StringT global_model_file;
	suffix.Suffix(model_file);
	global_model_file.Root(model_file);
	global_model_file.Append(".io", suffix);

	/* output map file */
	StringT map_file;
	map_file.Root(model_file);
	map_file.Append(".n", size);
	map_file.Append(".io.map");

	/* set output map and and generate decomposition */
	token = 1;
	if (rank == 0)
	{
		try { Decompose(in, size, model_file, global_model_file, format, map_file); }
		catch (int code)
		{
			cout << " ::RunJob_parallel: exception on decomposition: " << code << endl;
			token = 0;
		}
	}

	/* synch and check status */
	check_sum = 0;
	if (MPI_Allreduce(&token, &check_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) !=
		MPI_SUCCESS) throw eMPIFail;
	if (check_sum != size) throw eGeneralFail;

	/* read partition information */
	PartitionT partition;
	StringT part_file;
	part_file.Root(model_file);
	part_file.Append(".n", size);
	part_file.Append(".part", rank);
	ifstreamT part_in(in.comment_marker(), part_file);
	token = 1;
	if (!part_in.is_open())
	{
		cout << "\n ::RunJob_parallel: could not open file: " << part_file << endl;
		token = 0;
	}
	else
	{
		part_in >> partition;
		if (partition.ID() != rank)
		{
			cout << "\n ::RunJob_parallel: rank and partition ID mismatch" << endl;
			token = 0;
		}
		
		/* set correct numbering scope */
		partition.SetScope(PartitionT::kLocal);
	}
	
	/* synch and check status */
	check_sum = 0;
	if (MPI_Allreduce(&token, &check_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) !=
		MPI_SUCCESS) throw eMPIFail;
	if (check_sum != size) throw eGeneralFail;
		
	/* write partial geometry files (if needed) */
	StringT partial_file;
	suffix.Suffix(model_file);
	partial_file.Root(model_file);
	partial_file.Append(".n", size);
	partial_file.Append(".p", rank);
	partial_file.Append(suffix);
	if (NeedModelFile(partial_file, format))
	{			
		cout << "\n ::RunJob_parallel: writing partial geometry file: " << partial_file << endl;
		EchoPartialGeometry(partition, model_file, partial_file, format);
		cout << " ::RunJob_parallel: writing partial geometry file: partial_file: "
		     << partial_file << ": DONE" << endl;
	}
	
	/* construct local problem (Initialize() changes the file name) */
	t1 = clock();
	token = 1;
	ifstreamT in_loc(in.comment_marker(), in.filename());
	if (!fJobCharPutBack)
	{
		char filetypechar;
		in_loc >> filetypechar;
	}
	FEManagerT_mpi FEman(in_loc, out, &partition, FEManagerT_mpi::kRun);
	try { FEman.Initialize(); }
	catch (int code)
	{
		status << "\n \"" << in_loc.filename() << "\" exit on exception during the\n";
		status << " construction phase. Check the input file for errors." << endl;
		
		/* echo some lines from input */
		if (code == eBadInputValue) Rewind(in_loc, status);
		token = 0;

		/* write exception codes to out file */
		FEman.WriteExceptionCodes(out);
	}
	
	check_sum = 0;
	if (MPI_Allreduce(&token, &check_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;
	if (check_sum != size)
	{
		/* gather tokens to rank 0 */
		iArrayT tokens;
		if (rank == 0) tokens.Allocate(size);
		if (MPI_Gather(&token, 1, MPI_INT, tokens.Pointer(), 1, MPI_INT, 0, MPI_COMM_WORLD)
			!= MPI_SUCCESS) throw eMPIFail;
		if (rank == 0)
		{
			status << "\n The following processes exit on exception during construction:\n";
			for (int i = 0; i < tokens.Length(); i++)
				if (tokens[i] != 1)
					status << setw(kIntWidth) << i << '\n';
			status.flush();
		}

		throw eGeneralFail;
	}
	
	/* external IO */
	IOManager_mpi* IOMan = NULL;
	if (!CommandLineOption("-split_io"))
	{
		/* read output map */
		iArrayT output_map;
		ReadOutputMap(in, map_file, output_map);

		/* set-up local IO */
		IOMan = new IOManager_mpi(in, output_map, *(FEman.OutputManager()), FEman.Partition(),
			global_model_file, format);
		if (!IOMan) throw eOutOfMemory;
		
		/* set external IO */
		FEman.SetExternalIO(IOMan);
	}

	/* solve */
	t2 = clock();
	token = 1;
	try { FEman.Solve(); }
	catch (int code)
	{
		status << "\n \"" << in.filename() << "\" exit on exception during the\n";
		status << " solution phase. See \"" << out_file << "\" for a list";
		status << " of the codes.\n";
		token = 0;

		/* write exception codes to out file */
		FEman.WriteExceptionCodes(out);		
	}

	/* free external IO */
	delete IOMan;

	check_sum = 0;
	if (MPI_Allreduce(&token, &check_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)
		!= MPI_SUCCESS) throw eMPIFail;
	if (check_sum != size)
	{
		/* gather tokens to rank 0 */
		iArrayT tokens;
		if (rank == 0) tokens.Allocate(size);
		if (MPI_Gather(&token, 1, MPI_INT, tokens.Pointer(), 1, MPI_INT, 0, MPI_COMM_WORLD)
			!= MPI_SUCCESS) throw eMPIFail;
		if (rank == 0)
		{
			status << "\n The following processes exit on exception during solution:\n";
			for (int i = 0; i < tokens.Length(); i++)
				if (tokens[i] != 1)
					status << setw(kIntWidth) << i << '\n';
			status.flush();
		}
		throw eGeneralFail;
	}
	t3 = clock(); } // end try

	/* job failure */
	catch (int code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		
		if (t3 == 0) t3 = clock();		
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   " Construction: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t3 - t2)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   " Construction: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: " << double(t3 - t2)/CLOCKS_PER_SEC << " sec.\n";
out <<   "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}	
#endif // __MPI__

/* print message on exception */
void FEExecutionManagerT::Rewind(ifstreamT& in, ostream& status) const
{
	/* reset stream */
	ifstream& istr = in;
	if (!istr.good()) istr.clear();
			
	/* rewind a couple of lines */
	in.rewind(4);
			
	status << " Error occurred while reading input near:\n\n";
	int i = 8;
	char line[255];
	istr.getline(line, 254);
	while(istr.good() && i-- > 0)
	{
		status << line << '\n';
		istr.getline(line, 254);
	}
	status << endl;
}

/* extract the model file name from the stream */
void FEExecutionManagerT::GetModelFile(ifstreamT& in, StringT& model_file,
	IOBaseT::FileTypeT& format) const
{
	/* partially construct FE manager */
	ifstreamT in_temp(in.comment_marker(), in.filename());
	if (!fJobCharPutBack)
	{
		char filetypechar;
		in_temp >> filetypechar;
	}

	ofstreamT out;
	FEManagerT fe_temp(in_temp, out);
	fe_temp.Initialize(FEManagerT::kParametersOnly);

	/* model file name*/
	model_file = fe_temp.ModelFile();
	
	/* file format */
	format = fe_temp.InputFormat();
}

/* initializations for rank 0 */
void FEExecutionManagerT::Decompose(ifstreamT& in, int size,
	const StringT& model_file, const StringT& global_model_file,
	IOBaseT::FileTypeT format, const StringT& output_map_file) const
{
	bool split_io = CommandLineOption("-split_io");
	bool need_output_map = NeedOutputMap(in, output_map_file, size) && !split_io;
	bool need_model_file = NeedModelFile(global_model_file, format) && !split_io;
	bool need_decomp = NeedDecomposition(in, model_file, size);
	if (need_output_map ||
	    need_decomp ||
	    need_model_file)
	{
		/* echo stream */
		StringT decomp_file;
		decomp_file.Root(in.filename());
		decomp_file.Append(".out");
		ofstreamT decomp_out;
		decomp_out.open(decomp_file);
	
		ifstreamT in_decomp(in.comment_marker(), in.filename());
		if (!fJobCharPutBack)
		{
			char filetypechar;
			in_decomp >> filetypechar;
		}

		/* construct global problem */
		FEManagerT_mpi global_FEman(in_decomp, decomp_out, NULL, FEManagerT_mpi::kDecompose);
		try { global_FEman.Initialize(FEManagerT::kAllButSolver); }
		catch (int code)
		{
			cout << "\n FEExecutionManagerT::Decompose: exception during construction: "
		         << global_FEman.Exception(code) << endl;
			throw code;
		}

		/* set output map */
		if (need_output_map)
		{
			cout << "\n Generating output map: " << output_map_file << endl;
			IOManager* global_IOman = global_FEman.OutputManager();
			iArrayT output_map;
			SetOutputMap(global_IOman->ElementSets(), output_map, size);
	
			/* write map file */
			ofstreamT map_out(output_map_file);
			map_out << "# number of processors\n";
			map_out << size << '\n';
			map_out << "# number of output sets\n";
			map_out << output_map.Length() << '\n';
			map_out << "# part to processor output map\n";
			map_out << output_map.wrap(8) << '\n';
			cout << " Generating output map: " << output_map_file << ": DONE"<< endl;
		}

		/* output model file */
		if (need_model_file)
		{
			cout << "\n Writing output model file: " << global_model_file << endl;
			global_FEman.WriteGeometryFile(global_model_file, format);
			cout << " Writing output model file: " << global_model_file << ": DONE" << endl;
		}
	
		/* decompose */
		ArrayT<PartitionT> partition(size);
		if (need_decomp)
		{
			/* graph object */
			GraphT graph;
			
			try
			{
				cout << "\n Decomposing: " << model_file << endl;
				global_FEman.Decompose(partition, graph, true);
				cout << " Decomposing: " << model_file << ": DONE"<< endl;
			}
			catch (int code)
			{
				cout << "\n FEExecutionManagerT::Decompose: exception during decomposition: "
			         << code << endl;
				throw code;
			}
			
			/* write decomposition map */
			if (format == IOBaseT::kExodusII)
			{
				/* file name */
				StringT map_file;
				map_file.Root(model_file);
				map_file.Append(".n", size);
				map_file.Append(".decomp.exo");
				cout << " Node map file: " << map_file << endl;
			
				/* database */
				ExodusT exo(cout);
				StringT str;
				ArrayT<StringT> str_list;
				const dArray2DT& coords = global_FEman.Coordinates();
				int nnd = coords.MajorDim();
				int nsd = coords.MinorDim();
				exo.Create(map_file, str, str_list, str_list, nsd, nnd,
					nnd, size, 0, 0);
			
				/* coordinates */
				exo.WriteCoordinates(coords);
				
				/* data in decomp file */
				dArrayT part(nnd); part = -1;
				dArrayT inex(nnd); inex = 0;
				for (int i = 0; i < size; i++)
				{
					/* "owned" nodes */
					const iArrayT& nd_i = partition[i].Nodes_Internal();
					const iArrayT& nd_b = partition[i].Nodes_Border();
					iArray2DT connects(nd_i.Length() + nd_b.Length(), 1);
					connects.CopyPart(0, nd_i, 0, nd_i.Length());
					connects.CopyPart(nd_i.Length(), nd_b, 0, nd_b.Length());
					
					/* convert to global numbering */
					if (partition[i].NumberScope() != PartitionT::kGlobal)
						partition[i].SetNodeScope(PartitionT::kGlobal, connects);
				
					/* label part */
					for (int j = 0; j < connects.Length(); j++)
						part[connects[j]] = i;
				
					/* label border nodes */
					int* pnd_b = connects.Pointer(nd_i.Length());
					for (int k = 0; k < nd_b.Length(); k++)
						inex[*pnd_b++] = 1;
				
					/* write to file */
					connects++;
					exo.WriteConnectivities(i+1, GeometryT::kPoint, connects);
				}
				
				/* degree of each node */
				iArrayT i_degree(nnd);
				int shift;
				graph.Degrees(i_degree, shift);
				if (shift != 0)
				{
					cout << "\n FEExecutionManagerT::Decompose: unexpected node number shift: "
					     << shift << endl;				
					throw eGeneralFail;
				}
				dArrayT degree(nnd);
				for (int j = 0; j < nnd; j++)
					degree[j] = i_degree[j];

				/* part labels */
				ArrayT<StringT> labels(3);
				labels[0] = "part";
				labels[1] = "inex";
				labels[2] = "degree";
				exo.WriteNodeLabels(labels);
				exo.WriteTime(1, 0.0);
				exo.WriteNodalVariable(1, 1, part);
				exo.WriteNodalVariable(1, 2, inex);
				exo.WriteNodalVariable(1, 3, degree);

				cout << " Node map file: " << map_file << ": DONE"<< endl;
			}
		}
		else
		{
			/* read partition information from stream */
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT part_file;
				part_file.Root(model_file);
				part_file.Append(".n", size);
				part_file.Append(".part", i);
				ifstreamT part_in(in.comment_marker(), part_file);
				
				part_in >> partition[i];
				partition[i].SetScope(PartitionT::kLocal);
			}
		}
		
		/* write partial geometry files */
		if (CommandLineOption("-decomp"))
		{
			//NOTE - not the most efficient way since the global data is
			//       read multiple times, but OK for now
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT partial_file, suffix;
				suffix.Suffix(model_file);
				partial_file.Root(model_file);
				partial_file.Append(".n", size);
				partial_file.Append(".p", i);
				partial_file.Append(suffix);
				
				if (NeedModelFile(partial_file, format))
				{			
					cout << "     Writing partial model file: " << partial_file << endl;
					try { EchoPartialGeometry(partition[i], model_file, partial_file, format); }
					catch (int error)
					{
						cout << "\n ::Decompose: exception writing file: " << partial_file << endl;
						throw error;
					}
				}
			}
		}
	}
	else
		cout << "\n ::Decompose: decomposition files exist" << endl;
}

/* returns 1 if a new decomposition is needed */
bool FEExecutionManagerT::NeedDecomposition(ifstreamT& in, const StringT& model_file,
	int size) const
{
	/* model file root */
	StringT root;
	root.Root(model_file);

	for (int i = 0; i < size; i++)
	{
		/* partition data file name */
		StringT part_file = root;
		part_file.Append(".n", size);
		part_file.Append(".part", i);
	
		/* open partition file */
		ifstreamT part_in(in.comment_marker(), part_file);
		if (part_in.is_open())
		{
			StringT version;
			part_in >> version;
		
			int num_parts, part_ID;
			part_in >> num_parts;
			part_in >> part_ID;

			if (!PartitionT::CheckVersion(version) ||
			     num_parts != size ||
			     part_ID != i) return true;
		}
		else
			return true;
	}
	return false;
}

/* returns true if the global output model file is not found */
bool FEExecutionManagerT::NeedModelFile(const StringT& model_file,
	IOBaseT::FileTypeT format) const
{
	switch (format)
	{
		case IOBaseT::kTahoeII:
		{
			ModelFileT file;
			if (file.OpenRead(model_file) == ModelFileT::kOK)
				return false;
			else
				return true;
		}	
		case IOBaseT::kExodusII:
		{
			ExodusT file(cout);
			if (file.OpenRead(model_file))
				return false;			
			else
				return true;
		}
		default:
		{
			ifstreamT test(model_file);
			if (test.is_open())
				return false;
			else
				return true;
		}
	}
}

/* returns true if a new output map is needed */
bool FEExecutionManagerT::NeedOutputMap(ifstreamT& in, const StringT& map_file,
	int size) const
{
	/* open map file */
	ifstreamT map_in(in.comment_marker(), map_file);
	if (map_in.is_open())
	{
		int num_parts;
		map_in >> num_parts;
		if (num_parts != size) return true;
	
		int num_sets;
		map_in >> num_sets;
		iArrayT map(num_sets);
		map_in >> map;
		
		int min, max;
		map.MinMax(min, max);
		if (min < 0 || max >= size)
			return true;
		else
			return false;
	}
	else
		return true;
}

void FEExecutionManagerT::ReadOutputMap(ifstreamT& in, const StringT& map_file,
	iArrayT& map) const
{
	/* map file */
	ifstreamT map_in(in.comment_marker(), map_file);
	if (!map_in.is_open())
	{
		cout << "\n FEExecutionManagerT::ReadOutputMap: could not open io map file: "
		     << map_file << endl;
		throw eGeneralFail;
	}

	/* read map */
	int size, num_sets;
	map_in >> size >> num_sets;
	map.Allocate(num_sets);
	map_in >> map;

	/* check */
	int min, max;
	map.MinMax(min, max);
	if (min < 0 || max >= size)
	{
		cout << "\n FEExecutionManagerT::ReadOutputMap: map error\n";
		cout << map.wrap(5) << '\n';
		cout.flush();
		throw eGeneralFail;
	}
}

/* set output map based on length of map */
void FEExecutionManagerT::SetOutputMap(const ArrayT<OutputSetT*>& output_sets,
	iArrayT& output_map, int size) const
{
	/* initialize map */
	output_map.Allocate(output_sets.Length());
	output_map = 0;
	
	iArrayT output_counts(size);
	output_counts = 0;
	
	/* output set size */
	iArrayT set_size(output_sets.Length());
	for (int i = 0; i < output_sets.Length(); i++)
	{
		const OutputSetT& set = *(output_sets[i]);
	
		int size = 0;
		size += set.NumNodes()*set.NodeOutputLabels().Length();
		size += set.NumElements()*set.ElementOutputLabels().Length();
		
		set_size[i] = size;
	}

	for (int j = 0; j < output_sets.Length(); j++)
	{
		int max_at;
		set_size.Max(max_at);
	
		int min_at;
		output_counts.Min(min_at);
	
		output_map[max_at] = min_at;
		output_counts[min_at] += set_size[max_at];
		set_size[max_at] = 0;
	}
}

/* write partial geometry files */
void FEExecutionManagerT::EchoPartialGeometry(const PartitionT& partition,
	const StringT& model_file, const StringT& partial_file,
	IOBaseT::FileTypeT format) const
{

	switch (format)
	{
		case IOBaseT::kExodusII:
			EchoPartialGeometry_ExodusII(partition, model_file, partial_file);
			break;
	
		case IOBaseT::kTahoeII:
			EchoPartialGeometry_TahoeII(partition, model_file, partial_file);
			break;

		default:
			cout << "\n ::WritePartialGeometry: unsupported file format: "
			     << format << endl;	
			throw eGeneralFail;
	}
}

void FEExecutionManagerT::EchoPartialGeometry_ExodusII(const PartitionT& partition,
	const StringT& model_file, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();

	/* original model file */
	ExodusT model_ALL(cout);
	if (!model_ALL.OpenRead(model_file))
	{
		cout << "\n FEExecutionManagerT::EchoPartialGeometry_ExodusII: error opening file: "
		     << model_file << endl;
		throw eGeneralFail;
	}		

	/* collect file creation information */
	StringT title = "partition file: ";
	title.Append(part);
	ArrayT<StringT> nothing;
	const iArrayT& node_map = partition.NodeMap();
	int num_node = node_map.Length();
	int num_dim  = model_ALL.NumDimensions();
	int num_blks = model_ALL.NumElementBlocks();
	int num_ns   = model_ALL.NumNodeSets();
	int num_ss   = model_ALL.NumSideSets();
	int num_elem = 0;

	iArrayT elementID(num_blks);
	model_ALL.ElementBlockID(elementID);
	for (int ii = 0; ii < num_blks; ii++)
		num_elem += (partition.ElementMap(elementID[ii])).Length();

	/* partial model file */
	ExodusT model(cout);
	model.Create(partial_file, title, nothing, nothing, num_dim, num_node,
		num_elem, num_blks, num_ns, num_ss);
	
	/* coordinates */
	dArray2DT coords_ALL(model_ALL.NumNodes(), num_dim);
	model_ALL.ReadCoordinates(coords_ALL);
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	model.WriteCoordinates(coords);
	coords.Free();
		
	/* element sets */
	for (int j = 0; j < elementID.Length(); j++)
	{
		/* read global block */
		int num_elems;
		int num_elem_nodes;
		model_ALL.ReadElementBlockDims(elementID[j], num_elems, num_elem_nodes);

		iArray2DT set_ALL(num_elems, num_elem_nodes);
		GeometryT::CodeT geometry_code;
		model_ALL.ReadConnectivities(elementID[j], geometry_code, set_ALL);

		/* collect connectivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elementID[j]);
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);

		/* map to local scope */
		set--;
		partition.SetNodeScope(PartitionT::kLocal, set);
		set++;

		/* write to file */
		model.WriteConnectivities(elementID[j], geometry_code, set);
	}
		
	/* node sets */
	iArrayT nodeID(num_ns);
	model_ALL.NodeSetID(nodeID);
	for (int k = 0; k < nodeID.Length(); k++)
	{
		/* whole node set */
		iArrayT nodeset_ALL(model_ALL.NumNodesInSet(nodeID[k]));
		model_ALL.ReadNodeSet(nodeID[k], nodeset_ALL);
		nodeset_ALL--;
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		
		/* non-empty set */
		if (1 || local_indices.Length() > 0)
		{
			iArrayT nodeset(local_indices.Length());
			nodeset.Collect(local_indices, nodeset_ALL);
			
			/* map to local numbering */
			partition.SetNodeScope(PartitionT::kLocal, nodeset);
				
			/* add */
			nodeset++;
			model.WriteNodeSet(nodeID[k], nodeset);
		}
	}

	/* side sets */		
	iArrayT sideID(num_ss);
	model_ALL.SideSetID(sideID);
	for (int l = 0; l < sideID.Length(); l++)
	{
		/* whole side set */
		int element_set_ID;
		iArray2DT sideset_ALL(model_ALL.NumSidesInSet(sideID[l]), 2);		
		model_ALL.ReadSideSet(sideID[l], element_set_ID, sideset_ALL);

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
			elements_ALL--;
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
						
			/* non-empty set */
			if (local_indices.Length() > 0)
			{
				sideset.Allocate(local_indices.Length(), sideset_ALL.MinorDim());
				sideset.RowCollect(local_indices, sideset_ALL);

				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
				elements++;
				sideset.SetColumn(0, elements);
			}

			/* add */
			model.WriteSideSet(sideID[l], element_set_ID, sideset);
		}			
	}
}

void FEExecutionManagerT::EchoPartialGeometry_TahoeII(const PartitionT& partition,
		const StringT& model_file, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();

	/* original model file */
	ModelFileT model_ALL;
	if (model_ALL.OpenRead(model_file) != ModelFileT::kOK)
	{
		cout << "\n FEExecutionManagerT::EchoPartialGeometry_TahoeII: error opening file: "
		     << model_file << endl;
		throw eGeneralFail;
	}		
	
	/* open model file */
	bool extern_file = true;
	ModelFileT model;
	model.OpenWrite(partial_file, extern_file);

	/* title */
	StringT title;
	if (model_ALL.GetTitle(title) != ModelFileT::kOK) throw eGeneralFail;
	title.Append(": partition ", part);
	if (model.PutTitle(title) != ModelFileT::kOK) throw eGeneralFail;

	/* nodal coordinates */
	const iArrayT& node_map = partition.NodeMap();
	dArray2DT coords_ALL;
	if (model_ALL.GetCoordinates(coords_ALL) != ModelFileT::kOK) throw eGeneralFail;
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	if (model.PutCoordinates(coords) != ModelFileT::kOK) throw eGeneralFail;
	coords.Free();	
		
	/* element sets */
	iArrayT elementID;
	if (model_ALL.GetElementSetID(elementID) != ModelFileT::kOK) throw eGeneralFail;
	for (int j = 0; j < elementID.Length(); j++)
	{
		/* collect connecitivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elementID[j]);
		iArray2DT set_ALL;
		if (model_ALL.GetElementSet(elementID[j], set_ALL) != ModelFileT::kOK) throw eGeneralFail;				
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);
			
		/* map to local scope */
		set--;
		partition.SetNodeScope(PartitionT::kLocal, set);
		set++;

		/* add */
		if (model.PutElementSet(elementID[j], set) != ModelFileT::kOK) throw eGeneralFail;
	}
		
	/* node sets */
	iArrayT nodeID;
	if (model_ALL.GetNodeSetID(nodeID) != ModelFileT::kOK) throw eGeneralFail;
	for (int k = 0; k < nodeID.Length(); k++)
	{
		/* whole node set */
		iArrayT nodeset_ALL;		
		if (model_ALL.GetNodeSet(nodeID[k], nodeset_ALL) != ModelFileT::kOK) throw eGeneralFail;
		nodeset_ALL--;
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		iArrayT nodeset(local_indices.Length());
		nodeset.Collect(local_indices, nodeset_ALL);
			
		/* map to local numbering */
		partition.SetNodeScope(PartitionT::kLocal, nodeset);
			
		/* add */
		nodeset++;
		if (model.PutNodeSet(nodeID[k], nodeset) != ModelFileT::kOK) throw eGeneralFail;
	}

	/* side sets */		
	iArrayT sideID;
	if (model_ALL.GetSideSetID(sideID) != ModelFileT::kOK) throw eGeneralFail;
	for (int l = 0; l < sideID.Length(); l++)
	{
		/* whole side set */
		int element_set_ID;
		iArray2DT sideset_ALL;		
		if (model_ALL.GetSideSet(sideID[l], element_set_ID, sideset_ALL) != ModelFileT::kOK) throw eGeneralFail;

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
			elements_ALL--;
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
			sideset.Allocate(local_indices.Length(), sideset_ALL.MinorDim());
			sideset.RowCollect(local_indices, sideset_ALL);
				
			/* map to (block) local numbering */
			if (sideset.MajorDim() > 0)
			{
				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
				elements++;
				sideset.SetColumn(0, elements);
			}
		}
			
		/* add */
		if (model.PutSideSet(sideID[l], element_set_ID, sideset) != ModelFileT::kOK) throw eGeneralFail;
	}

	/* close database */
	model.Close();
}
