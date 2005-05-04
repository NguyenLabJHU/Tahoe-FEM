/* $Id: ComparatorT.cpp,v 1.26 2005-05-04 06:40:38 paklein Exp $ */
#include "ComparatorT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#if defined (__GCC_3__) || defined (__GCC_4__)
#include <strstream>
#else
#include <strstream.h>
#endif

#include "ExceptionT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "dArrayT.h"

/* XML input files */
#include "ParameterListT.h"
#include "expat_ParseT.h"

const char kBenchmarkDirectory[] = "benchmark";
using namespace Tahoe;

/* default tolerances */
double abs_tol = 1.0e-10;
double rel_tol = 1.0e-08;

/* constructor */
ComparatorT::ComparatorT(int argc, char* argv[], char job_char, char batch_char):
	FileCrawlerT(argc, argv, job_char, batch_char),
	fAbsTol(abs_tol),
	fRelTol(rel_tol),
	fIsRoot(true),
	fAbsTol_batch(-1.0),
	fRelTol_batch(-1.0)
{
	/* tolerances reset on command line */
	int index;
	if (CommandLineOption("-abs_tol", index))
	{
		istrstream in((const char *) fCommandLineOptions[index + 1]);
		in >> abs_tol;
	}
	if (CommandLineOption("-rel_tol", index))
	{
		istrstream in((const char *) fCommandLineOptions[index + 1]);
		in >> rel_tol;
	}
}

/* prompt input files until "quit" */
void ComparatorT::Run(void)
{
	/* check for command line specification of the file names */
	int index;
	if (CommandLineOption("-f1", index))
	{
		/* first file */
		StringT file_1(fCommandLineOptions[index + 1]);
		if (CommandLineOption("-f2", index))
		{
			/* second file */
			StringT file_2(fCommandLineOptions[index + 1]);

			/* (re-)set default tolerances */
			fAbsTol = abs_tol;	
			fRelTol = rel_tol;

			/* compare */
			PassOrFail(file_1, file_2, false, true);
		}
		else
			ExceptionT::GeneralFail("ComparatorT::Run", "expecting command line option \"-f2\"");
	}
	else /* inherited */
		FileCrawlerT::Run();
}

/**********************************************************************
* Protected
**********************************************************************/

/* handle batch file command */
void ComparatorT::BatchFileCommand(const StringT& command, ifstreamT& batch)
{
	/* check */
	if (command[0] != '-') {
		cout << "\n ComparatorT::BatchFileCommand: commands must start with \"-\": \""
		     << command << '\"' << endl;
		throw ExceptionT::kGeneralFail;
	}

	if (command == "-abs_tol")
	{
		double tol = -99;
		batch >> tol;
		if (tol > 0) 
			fAbsTol_batch = tol;
		else
			cout << "\n ComparatorT::BatchFileCommand: error reading abs_tol" << endl;
	}
	else if (command == "-rel_tol")
	{
		double tol = -99;
		batch >> tol;
		if (tol > 0) 
			fRelTol_batch = tol;
		else
			cout << "\n ComparatorT::BatchFileCommand: error reading rel_tol" << endl;
	}
	else
		cout << "\n ComparatorT::BatchFileCommand: ignoring unrecognized command \"" 
		     << command << '\"' << endl;
}

void ComparatorT::RunJob(ifstreamT& in, ostream& status)
{
	const char caller[] = "ComparatorT::RunJob";

	/* check for multi-Tahoe files */
	char job_char = in.next_char(); /* peek at next char */

	if (job_char != fJobChar) /* standard job file */
	{
		/* append to list of files */
		fFiles.Append(in.filename());

		/* append to results */
		cout << "\nSTART: " << in.filename() << '\n';
		bool result = false;
		try 
		{ 
			int index; 
			// to handle batches from MakeCSE and Translate
			if (CommandLineOption ("-ext", index))
				result = PassOrFail_Extension (in, fCommandLineOptions[index + 1]);
			else
				result = PassOrFail(in); 
		}
		catch (ExceptionT::CodeT error) 
		{
	    	cout << in.filename() << ": " << "EXCEPTION: " << error << endl;
	    	result = false;
	  	}
		fPassFail.Append(result);
		cout << in.filename() << ": " << ((result) ? "PASS" : "FAIL") << '\n';
		cout << "\nEND: " << in.filename() << '\n';

		/* clear values passed in batch file */
		fAbsTol_batch = -1.0;
		fRelTol_batch = -1.0;
	}
	else /* multi-Tahoe jobs have embedded input files */
	{
		/* path */
		StringT path;
		path.FilePath(in.filename());
	
		/* store values passed in batch file */
		double abs_tol_batch = fAbsTol_batch;
		double rel_tol_batch = fRelTol_batch;

		/* clear job char */
		in >> job_char;
		
		/* files listed in the multi-Tahoe file */
		StringT atom_file, atom_bridge_file, cont_file, cont_bridge_file;
		in >> atom_file 
		   >> atom_bridge_file
		   >> cont_file
		   >> cont_bridge_file;
		
		/* atom file */
		atom_file.ToNativePathName();
		atom_file.Prepend(path);
		ifstreamT atom_in(in.comment_marker(), atom_file);
		atom_in >> job_char;
		if (job_char)
			RunJob(atom_in, status); /* check file */
		else
			ExceptionT::GeneralFail(caller, "expecting job char %c not %c from file \"%s\"",
				fJobChar, job_char, atom_in.filename());

		/* continuum file */
		cont_file.ToNativePathName();
		cont_file.Prepend(path);
		ifstreamT cont_in(in.comment_marker(), cont_file);
		cont_in >> job_char;
		if (job_char)
			RunJob(cont_in, status); /* check file */
		else
			ExceptionT::GeneralFail(caller, "expecting job char %c not %c from file \"%s\"",
				fJobChar, job_char, cont_in.filename());

		/* restore values passed in batch file */
		fAbsTol_batch = abs_tol_batch;
		fRelTol_batch = rel_tol_batch;
	}
}

/* batch file processing */
void ComparatorT::RunBatch(ifstreamT& in, ostream& status)
{
	/* keep path to the root batch file */
	bool is_root = false;
	if (fIsRoot) 
	{
		fIsRoot = false;
		is_root = true;
	}

	/* inherited */
	FileCrawlerT::RunBatch(in, status);
	
	/* write summary */
	if (is_root)
	{
		/* check */
		if (fFiles.Length() != fPassFail.Length()) throw ExceptionT::kGeneralFail;

		/* relative path */
		StringT path;
		path.FilePath(in.filename());

		/* open output stream */
		StringT file_name = "compare.summary";
		file_name.ToNativePathName();
		file_name.Prepend(path);
		ofstreamT summary;
		summary.open(file_name);
		if (!summary.is_open())
		{
			cout << "\n ComparatorT::RunBatch: ERROR: could not open file: " 
			     << file_name << endl;
			throw ExceptionT::kGeneralFail;
		}

		/* count pass/fail */
		int pass_count = 0;
		for (int j = 0; j < fPassFail.Length(); j++)
			if (fPassFail[j]) pass_count++;

		/* write summary results */
		const char* result[]= {"FAIL", "PASS"};
		summary << "number of tests: " << fFiles.Length() << '\n';
		summary << "         passed: " << pass_count << '\n';
		summary << "         failed: " << fFiles.Length() - pass_count << '\n';
		for (int i = 0; i < fFiles.Length(); i++)
			summary << result[fPassFail[i]] << ": " << fFiles[i] << endl;
			
		/* empty lists */
		fPassFail.Resize(0);
		fFiles.Resize(0);
		
		/* reset root */
		fIsRoot = true;		
	}
}

/* recursive dispatch */
void ComparatorT::JobOrBatch(ifstreamT& in, ostream& status)
{
	/* catch XML job files */
	StringT ext;
	ext.Suffix(in.filename());
	if (ext == ".xml") {

		/* read the input file */
		ParameterListT params;
		expat_ParseT parser;
		parser.Parse(in.filename(), params);

		/* plain tahoe simulation */
		if (params.Name() == "tahoe" || params.Name() == "tahoe_bridging" || params.Name() == "tahoe_THK")
			RunJob(in, status);
		else if (params.Name() == "tahoe_multi" || params.Name() == "tahoe_bridging_scale") /* multi-tahoe simulation */
		{
			/* add results to the list */
			cout << "\nSTART: " << in.filename() << '\n';

			StringT path, file;
			path.FilePath(in.filename());
			ifstreamT in_xml;
		
			/* continuum results */
			file = params.GetParameter("continuum_input");
			file.ToNativePathName();
			file.Prepend(path);
			in_xml.open(file);
			JobOrBatch(in_xml, status);
			in_xml.close();
			bool continuum_pass = fPassFail.Last();

			/* atomistic results */
			file = params.GetParameter("atom_input");
			file.ToNativePathName();
			file.Prepend(path);
			in_xml.open(file);
			JobOrBatch(in_xml, status);
			in_xml.close();
			bool atomistic_pass = fPassFail.Last();
			
			/* report results to the list */
			bool result = continuum_pass && atomistic_pass;
			fFiles.Append(in.filename());
			fPassFail.Append(result);
			cout << in.filename() << ": " << ((result) ? "PASS" : "FAIL") << '\n';
			cout << "\nEND: " << in.filename() << '\n';
		}
		else
			ExceptionT::GeneralFail("ComparatorT::JobOrBatch", "unrecognized name \"%s\"",
				params.Name().Pointer());
	}
	else /* inherited */
		FileCrawlerT::JobOrBatch(in, status);
}

/**********************************************************************
* Private
**********************************************************************/

/* compare results against benchmarks */
bool ComparatorT::PassOrFail(ifstreamT& in) //const
{
	/* path to the current directory */
	StringT path;
	path.FilePath(in.filename());

	/* reset default tolerances */
	bool do_relative, do_absolute;
	Tolerance (in, do_relative, do_absolute);
	
	/* job file name */
	StringT file_root = in.filename();
	file_root.Drop(path.StringLength());
	file_root.Root();
	
	/* look for results as [infile].run or [infile].io[N].run */
	bool pass_or_fail= true;
	bool OK = true;
	int count = -1;
	int found_count = 0;
	while (OK) {
	
		/* look for results file in current directory */
		StringT current = path;
		current.Append(file_root);
		if (count != -1) current.Append(".io", count);
		current.Append(".run");

		/* look for results file in benchmark directory */
		StringT benchmark(kBenchmarkDirectory);
		benchmark.ToNativePathName();
		benchmark.Prepend(path);
		if (file_root[0] != StringT::DirectorySeparator())
		{
			const char sep[2] = {StringT::DirectorySeparator(), '\0'};
			benchmark.Append(sep);
		}
		benchmark.Append(file_root);
		if (count != -1) benchmark.Append(".io", count);
		benchmark.Append(".run");
		
		/* do compare */
		if (fstreamT::Exists(benchmark))
		{
			found_count++;
			pass_or_fail = pass_or_fail && PassOrFail(current, benchmark, do_relative, do_absolute);
		}
		else if (count != -1)
			OK = false;

	 	/* next */
	 	count++;
	 }
	 
	 /* no results found */
	 if (found_count == 0) {
	 	cout << "no results found for input file: " << in.filename() << endl;
	 	pass_or_fail = false;
	 }
	 
	 /* return */
	 return pass_or_fail;
}

bool ComparatorT::PassOrFail_Extension (ifstreamT& in, const StringT& extension)
{
  /* set tolerance values */
  bool do_relative, do_absolute;
  Tolerance (in, do_relative, do_absolute);

  /* path to the current directory */
  StringT path;
  path.FilePath (in.filename());
  
  /* root */
  StringT file_root = in.filename();
  file_root.Drop (path.StringLength());
  file_root.Root();

  /* look for results as [infile][extension].geom */
  bool pass_or_fail = true;

  /* look for results file in current directory */
  StringT current = path;
  current.Append (file_root);
  current.Append (extension);
  current.Append (".geom");
  bool current_exists = false;
  if (fstreamT::Exists (current))
    current_exists = true;
  else
    {
      cout << "  no results found:" << current;
      cout << "\n    for input file: " << in.filename() << endl;
      pass_or_fail = false;
    }

  /* look for results file in benchmark directory */
  StringT benchmark(kBenchmarkDirectory);
  benchmark.ToNativePathName();
  benchmark.Prepend(path);
  if (file_root[0] != StringT::DirectorySeparator())
    {
      const char sep[2] = {StringT::DirectorySeparator(), '\0'};
      benchmark.Append(sep);
    }
  benchmark.Append (file_root);
  benchmark.Append (extension);
  benchmark.Append (".geom");
  bool benchmark_exists = false;
  if (fstreamT::Exists (benchmark))
    benchmark_exists = true;
  else
    {
      cout << "no benchmark found:" << benchmark;
      cout << "\n    for input file: " << in.filename() << endl;
      pass_or_fail = false;
    }

  if (current_exists && benchmark_exists)
    pass_or_fail = PassOrFail (current, benchmark, do_relative, do_absolute);
  return pass_or_fail;
}
void ComparatorT::Tolerance (ifstreamT& in, bool& do_relative, bool& do_absolute)
{
        /* (re-)set default tolerances */
        fAbsTol = abs_tol;
	fRelTol = rel_tol;

	/* default comparison method */
	do_relative = true;
	do_absolute = false;

	/* override local tolerance file */
	if (fAbsTol_batch > 0) {
		fAbsTol = fAbsTol_batch;
		do_absolute = true;
		do_relative = false;
	}
	else if (fRelTol_batch > 0) {
		fRelTol = fRelTol_batch;
		do_absolute = false;
		do_relative = true;
	}
	else
	{
		/* look for local tolerance file */
	        StringT path;
		path.FilePath(in.filename());
		StringT tol_file = "compare.tol";
		tol_file.ToNativePathName();
		tol_file.Prepend(path);
		ifstreamT tol_in(in.comment_marker(), tol_file);

		/* clear skip list */
		fSkipLabels.Resize(0);

		/* read */
		if (tol_in.is_open())
		{
			cout << "found local tolerance file: " << tol_file << '\n';
			StringT key;
			tol_in >> key;
			while (tol_in.good())
			{			
				if (key == "abs_tol")
				{
					tol_in >> fAbsTol;
			
					/* for now, disable relative tolerance if found */	
					do_relative = false;
					do_absolute = true;
				}
				else if (key == "rel_tol")
				{
					tol_in >> fRelTol;

					/* for now, disable relative tolerance if found */	
					do_relative = true;
					do_absolute = false;
				}
				else if (key == "skip")
				{
					StringT skip_label;
					tol_in >> skip_label;
					fSkipLabels.Append(skip_label);
				}
				tol_in >> key;
			}
		}
		else
			cout << "did not find local tolerance file: " << tol_file << '\n';
	}	
		
	cout << "using tolerances:\n";
	if (do_absolute) cout << "    abs_tol = " << fAbsTol << '\n';
	if (do_relative) cout << "    rel_tol = " << fRelTol << '\n';
	cout.flush();
}

/* compare results */
bool ComparatorT::PassOrFail(const StringT& file_1, const StringT& file_2,
	bool do_rel, bool do_abs)
{
	/* header */
	cout << "\n file 1: " << file_1 << '\n'
	     <<   " file 2: " << file_2 << " (reference)" << endl;

	/* open file_1 -> "current" */
	ModelManagerT res_1(cout);
	try {
		/* file format */
		IOBaseT::FileTypeT format = IOBaseT::name_to_FileTypeT(file_1);
		if (!res_1.Initialize(format, file_1, true)) throw ExceptionT::kDatabaseFail;			
	}
	catch (ExceptionT::CodeT error) {
		cout << "\n ComparatorT::PassOrFail: could not initialize file \"" << file_1 << '\"' << endl;	
		return false;
	}

	/* open file_2 -> "benchmark" */
	ModelManagerT res_2(cout);
	try {
		/* file format */
		IOBaseT::FileTypeT format = IOBaseT::name_to_FileTypeT(file_2);
		if (!res_2.Initialize(format, file_2, true)) throw ExceptionT::kDatabaseFail;			
	}
	catch (ExceptionT::CodeT error) {
		cout << "\n ComparatorT::PassOrFail: could not initialize file \"" << file_2 << '\"' << endl;	
		return false;
	}
	
	/* verify node dimensions */
	if (res_1.NumNodes() != res_2.NumNodes() ||
	    res_1.NumDimensions() != res_2.NumDimensions()) {
		cout << "\n ComparatorT::PassOrFail: node dimension mismatch" << endl;	    
	    return false;
	}
	if (res_1.NumNodeVariables() != res_2.NumNodeVariables()) {
		cout << "\n ComparatorT::PassOrFail: node variable dimension mismatch" << endl;	    
	    return false;
	}

	/* verify block IDs and dimensions */
	if (res_1.NumElementVariables() != res_2.NumElementVariables()) {
		cout << "\n ComparatorT::PassOrFail: element variable dimension mismatch" << endl;	    
	    return false;
	}
	const ArrayT<StringT>& block_1 = res_1.ElementGroupIDs();
	const ArrayT<StringT>& block_2 = res_2.ElementGroupIDs();
	if (block_1.Length() != block_2.Length()) {
		cout << "\n ComparatorT::PassOrFail: number of element block ID's don't match" << endl;	    
	    return false;		
	}
	for (int i = 0; i < block_1.Length(); i++) {
		if (block_1[i] != block_2[i]) {
			cout << "\n ComparatorT::PassOrFail: block ID mismatch" << endl;	    
	    	return false;		
		}
		else {
			int nel_1, nen_1;
			res_1.ElementGroupDimensions(block_1[i], nel_1, nen_1);
			int nel_2, nen_2;
			res_2.ElementGroupDimensions(block_2[i], nel_2, nen_2);
			if (nel_1 != nel_2 || nen_1 != nen_2) {
				cout << "\n ComparatorT::PassOrFail: block dimension mismatch" << endl;
	    		return false;		
			}
		}
	}
	
	/* verify number of results */
	dArrayT time_1, time_2;
	res_1.TimeSteps(time_1);
	res_2.TimeSteps(time_2);
	if (time_1.Length() != time_2.Length()) {
		cout << "\n ComparatorT::PassOrFail: number of output steps don't match" << endl;
		return false;
	}
	for (int i = 0; i < time_1.Length(); i++)
		if (fabs(time_1[i] - time_2[i]) > 1.0e-12) {
			cout << "\n ComparatorT::PassOrFail: output step time mismatch" << endl;
			return false;
		}
		
	/* verify node ids */
	iArrayT ids_1, ids_2;
	res_1.AllNodeIDs(ids_1);
	res_2.AllNodeIDs(ids_2);
	for (int i = 0; i < ids_1.Length(); i++)
		if (ids_1[i] != ids_2[i]) {
			cout << "\n ComparatorT::PassOrFail: node id mismatch" << endl;
			return false;
		}
	
	/* output labels */
	ArrayT<StringT> n_labels_1, n_labels_2;
	res_1.NodeLabels(n_labels_1);
	res_2.NodeLabels(n_labels_2);
	ArrayT<StringT> e_labels_1, e_labels_2;
	res_1.ElementLabels(e_labels_1);
	res_2.ElementLabels(e_labels_2);

	/* run through time steps */
	for (int i = 0; i < time_1.Length(); i++) {

		/* nodal values */
		dArray2DT n_values_1(res_1.NumNodes(), res_1.NumNodeVariables());
		dArray2DT n_values_2(res_2.NumNodes(), res_2.NumNodeVariables());
		res_1.AllNodeVariables(i, n_values_1);
		res_2.AllNodeVariables(i, n_values_2);

		/* compare nodal blocks */
		if (!CompareDataBlocks(n_labels_1, n_values_1, n_labels_2, n_values_2, do_rel, do_abs)) {
			cout << "nodal data check: FAIL" << '\n';
			return false;
		}
		else cout << "nodal data check: PASS" << '\n';

		/* element values */
		dArray2DT e_values_1(res_1.NumElements(), res_1.NumElementVariables());
		dArray2DT e_values_2(res_2.NumElements(), res_2.NumElementVariables());
		res_1.AllElementVariables(i, e_values_1);
		res_2.AllElementVariables(i, e_values_2);
		
		/* compare element blocks */
		if (!CompareDataBlocks(e_labels_1, e_values_1, e_labels_2, e_values_2, do_rel, do_abs)) {
			cout << "element data check: FAIL" << '\n';
			return false;
		}
		else cout << "element data check: PASS" << '\n';
	}

	/* pass if fails through */
	return true;
}

bool ComparatorT::PassOrFail_old(const StringT& file_1, const StringT& file_2,
	bool do_rel, bool do_abs)
{
	/* open file_1 -> "current" */
	ifstreamT current_in(file_1);
	cout << "looking for current results: " << file_1 << ": "
	     << ((current_in.is_open()) ? "found" : "not found") << '\n';
	if (!current_in.is_open()) return false;

	/* open file_2 -> "benchmark" */
	ifstreamT bench_in(file_2);
	cout << "looking for benchmark results: " << file_2 << ": "
	     << ((bench_in.is_open()) ? "found" : "not found") << '\n';
	if (!bench_in.is_open()) return false;

	/* init */
	StringT b_str, c_str;
	if (!bench_in.FindString("O U T P U T", b_str)) return false;
	if (!current_in.FindString("O U T P U T", c_str)) return false;
	
	/* block info */
	int b_num_block, c_num_block;
	double b_time, c_time;
	bool b_OK = ReadDataBlockInfo(bench_in, b_time, b_num_block);
	bool c_OK = ReadDataBlockInfo(current_in, c_time, c_num_block);

	/* compare blocks */
	while (b_OK && c_OK)
	{
		/* verify time */
		if (fabs(b_time - c_time) > kSmall) {
			cout << "time mismatch: " << b_time << " != " << c_time << '\n';
			return false;
		}
		else
			cout << "time: " << b_time << '\n';

		/* verify number of element blocks */
		if (b_num_block != c_num_block) {
			cout << "block count mismatch: " << b_num_block << " != " << c_num_block << '\n';
			return false;
		}
		else
			cout << "number of blocks: " << b_num_block << '\n';		
	
		/* read nodal data */
		ArrayT<StringT> b_node_labels;
		dArray2DT b_node_data;
		if (!ReadNodalData(bench_in, b_node_labels, b_node_data)) {
			cout << "error reading node data from: " << bench_in.filename() << '\n';
			return false;
		}

		ArrayT<StringT> c_node_labels;
		dArray2DT c_node_data;
		if (!ReadNodalData(current_in, c_node_labels, c_node_data)) {
			cout << "error reading node data from: " << current_in.filename() << '\n';
			return false;
		}

		/* compare nodal blocks */
		if (!CompareDataBlocks(b_node_labels, b_node_data, c_node_labels, c_node_data, 
			do_rel, do_abs)) {
			cout << "nodal data check: FAIL" << '\n';
			return false;
		}
		else cout << "nodal data check: PASS" << '\n';

		/* loop over element blocks */
		for (int i = 0; i < b_num_block; i++)
		{
			/* read element data */
			ArrayT<StringT> b_element_labels;
			dArray2DT b_element_data;
			StringT b_block_ID;
			if (!ReadElementData(bench_in, b_element_labels, b_element_data, b_block_ID)) {
				cout << "error reading element data from: " << bench_in.filename() << '\n';
				return false;
			}

			ArrayT<StringT> c_element_labels;
			dArray2DT c_element_data;
			StringT c_block_ID;
			if (!ReadElementData(current_in, c_element_labels, c_element_data, c_block_ID)) {
				cout << "error reading element data from: " << current_in.filename() << '\n';
				return false;
			}

			/* verify block ID's */
			if (b_block_ID != c_block_ID) {
				cout << "block ID mismatch: " << b_block_ID << " != " << c_block_ID << '\n';
				return false;
			}
			else
				cout << "block ID: " << b_block_ID << '\n';		

			/* compare element blocks */
			if (!CompareDataBlocks(b_element_labels, b_element_data, c_element_labels, c_element_data, 
				do_rel, do_abs)) {
				cout << "element data check: FAIL" << '\n';
				return false;
			}
			else cout << "element data check: PASS" << '\n';
		}

		/* next block */
		b_OK = ReadDataBlockInfo(bench_in, b_time, b_num_block);
		c_OK = ReadDataBlockInfo(current_in, c_time, c_num_block);
	}

	/* termination */
	if (b_OK == c_OK) 
		return true;
	else { 
		if (!b_OK) 	cout << "error reading block data from: " << bench_in.filename() << '\n';
		if (!c_OK) 	cout << "error reading block data from: " << current_in.filename() << '\n';
		return false;	
	}
}

/* read data block header */
bool ComparatorT::ReadDataBlockInfo(ifstreamT& in, double& time, int& num_blocks) const
{
	StringT str;

	/* time */
	if (!in.FindString("Time", str)) return false;
	if (!str.Tail('=', time)) return false;

	/* number of element blocks */
	if (!in.FindString("Number of blocks", str)) return false;
	if (!str.Tail('=', num_blocks)) return false;
	
	return true;
}

/* read block of nodal data */
bool ComparatorT::ReadNodalData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data) const
{
	/* advance nodal data */
	StringT str;
	if (!in.FindString("Nodal data:", str)) return false;

	/* get dimensions */
	int num_nodes, num_values;
	if (!in.FindString("Number of nodal points", str)) return false;
	str.Tail('=', num_nodes);
	if (!in.FindString("Number of values", str)) return false;
	str.Tail('=', num_values);
	
	if (num_values > 0)
	{
		/* read labels */
		labels.Allocate(num_values + 2); // +2 for index and node number
		for (int i = 0; i < labels.Length(); i++)
			in >> labels[i];
		if (!in.good()) return false;
	
		/* read values */
		data.Allocate(num_nodes, num_values + 2); // +2 for index and node number
		in >> data;
	}
	return true;
}

/* read block of element data */
bool ComparatorT::ReadElementData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data, StringT& block_ID) const
{
	/* advance nodal data */
	StringT str;
	if (!in.FindString("Element data:", str)) return false;

	/* get dimensions */
	int num_nodes, num_values;
	if (!in.FindString("Block ID", str)) return false;
	str.Tail('=', block_ID);

	if (!in.FindString("Number of elements", str)) return false;
	str.Tail('=', num_nodes);
	if (!in.FindString("Number of values", str)) return false;
	str.Tail('=', num_values);
	
	if (num_values > 0)
	{
		/* read labels */
		labels.Allocate(num_values + 2); // +2 for index and element number
		for (int i = 0; i < labels.Length(); i++)
			in >> labels[i];
		if (!in.good()) return false;
	
		/* read values */
		data.Allocate(num_nodes, num_values + 2); // +2 for index and element number
		in >> data;
	}
	return true;
}

/* compare blocks - normalized by data set 1 */
bool ComparatorT::CompareDataBlocks(const ArrayT<StringT>& labels_1, const dArray2DT& data_1,
	const ArrayT<StringT>& labels_2, const dArray2DT& data_2, bool do_rel, bool do_abs) const
{
	/* compare labels */
	if (labels_1.Length() != labels_2.Length()) {
		cout << "label length mismatch: " << labels_1.Length() << " != " << labels_2.Length() << '\n';
		return false;
	}
	for (int i = 0; i < labels_1.Length(); i++)
		if (labels_1[i] != labels_2[i]) {
			cout << "mismatch with label " << i+1 << ": " << labels_1[i] << " != " << labels_2[i] << '\n';
			return false;
		}

	/* compare data */
	if (data_1.MajorDim() != data_2.MajorDim() ||
	    data_1.MinorDim() != data_2.MinorDim()) {
	    cout << "data dimension mismatch" << '\n';
	    return false;
	}
	
	/* echo dimensions */
	cout << "  number of rows: " <<  data_1.MajorDim() << '\n';
	cout << "number of values: " <<  data_1.MinorDim() << '\n';

	double max_abs_error = 0.0, max_rel_error = 0.0;
	for (int j = 0; j < data_1.MinorDim(); j++)
	{
		cout << "comparing: \"" << labels_1[j] << "\"\n";

		if (fSkipLabels.HasValue(labels_1[j]))
			cout << "comparing: \"" << labels_1[j] << "\": SKIPPED\n";
		else
		{
			double error_norm = 0.0, bench_norm = 0.0;
			int max_rel_error_index = -1, max_abs_error_index = -1;
			double max_rel_error_var = 0.0, max_abs_error_var = 0.0;
			for (int i = 0; i < data_1.MajorDim(); i++)
			{
				/* error */
				double abs_error = data_2(i,j) - data_1(i,j);
				double rel_error = (fabs(data_2(i,j)) > kSmall) ? abs_error/data_2(i,j) : 0.0;

				/* norms */
				bench_norm += data_2(i,j)*data_2(i,j);
				error_norm += abs_error*abs_error;
			
				/* track maximums */
				if (fabs(abs_error) > fabs(max_abs_error_var)) {
					max_abs_error_var = abs_error;
					max_abs_error_index = i;
				}
				if (fabs(rel_error) > fabs(max_rel_error_var)) {
					max_rel_error_var = rel_error;
					max_rel_error_index = i;
				}		
			}
		
			/* compute norms */
			bench_norm = sqrt(bench_norm);
			error_norm = sqrt(error_norm);
		
			cout << "abs error norm: " << error_norm << '\n';
			cout << "rel error norm: ";
			if (bench_norm > kSmall)
				cout << error_norm/bench_norm << '\n';
			else
				cout << "0.0" << '\n';
			
			/* maximum errors */
			cout << "         max abs error: " << max_abs_error_var << '\n';
			cout << "index of max abs error: " << max_abs_error_index + 1 << '\n';
			cout << "         max rel error: " << max_rel_error_var << '\n';
			cout << "index of max rel error: " << max_rel_error_index + 1 << '\n';
		
			/* running max */
			max_abs_error = (fabs(max_abs_error_var) > fabs(max_abs_error)) ? max_abs_error_var : max_abs_error;
			max_rel_error = (fabs(max_rel_error_var) > fabs(max_rel_error)) ? max_rel_error_var : max_rel_error;
		}
	}
	
	/* results strings */
	const char* judge[] = {"PASS", "FAIL", "IGNORE", "OVERRIDE"};
	
	/* absolute tolerance test */
	const char* abs_judge = judge[0];
	if (do_abs)
	{
		if (fabs(max_abs_error) > fAbsTol)
		{
			if (fabs(max_rel_error) > fRelTol)
				abs_judge = judge[1];
			else
				abs_judge = judge[3];
		}
	}
	else
		abs_judge = judge[2];

	/* relative tolerance test */
	const char* rel_judge = judge[0];
	if (do_rel)
	{
		if (fabs(max_rel_error) > fRelTol)
		{
			if (fabs(max_abs_error) > fAbsTol)
				rel_judge = judge[1];
			else
				rel_judge = judge[3];
		}
	}
	else
		rel_judge = judge[2];
	
	/* report */
	cout << "absolute tolerance result: " << fAbsTol << " : " << max_abs_error << " : " << abs_judge << '\n';
	cout << "relative tolerance result: " << fRelTol << " : " << max_rel_error << " : " << rel_judge << '\n';

	if (abs_judge == judge[1] || rel_judge == judge[1])
		return false;
	else
		return true;
}
