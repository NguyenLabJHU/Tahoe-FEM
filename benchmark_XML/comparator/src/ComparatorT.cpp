/* $Id: ComparatorT.cpp,v 1.1 2001-06-11 02:10:12 paklein Exp $ */

#include "ComparatorT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#include <strstream.h>

#include "fstreamT.h"
#include "Constants.h"
#include "ExceptionCodes.h"
#include "StringT.h"
#include "dArray2DT.h"

const char kBenchmarkDirectory[] = "benchmark";

const double abs_tol = 1.0e-10;
const double rel_tol = 1.0e-08;

/* constructor */
ComparatorT::ComparatorT(int argc, char* argv[], char job_char, char batch_char,
	int jobcharputback):
	FileCrawlerT(argc, argv, job_char, batch_char, jobcharputback),
	fAbsTol(abs_tol),
	fRelTol(rel_tol)
{
	/* tolerances reset on command line */
	int index;
	if (CommandLineOption("-abs_tol", index))
	{
		istrstream in((const char *) fCommandLineOptions[index + 1]);
		in >> fAbsTol;
	}
	if (CommandLineOption("-rel_tol", index))
	{
		istrstream in((const char *) fCommandLineOptions[index + 1]);
		in >> fRelTol;
	}
}

/**********************************************************************
* Protected
**********************************************************************/

void ComparatorT::RunJob(ifstreamT& in, ostream& status)
{
#pragma unused(status)

	/* append to list of files */
	fFiles.Append(in.filename());

	/* append to results */
	cout << "\nSTART: " << in.filename() << '\n';
	fPassFail.Append(PassOrFail(in));
	cout << "\nEND: " << in.filename() << '\n';
}

/* batch file processing */
void ComparatorT::RunBatch(ifstreamT& in, ostream& status)
{
	/* keep path to the root batch file */
	bool is_root = false;
	if (fRoot.StringLength() == 0) 
	{
		fRoot.FilePath(in.filename());
		is_root= true;
	}

	/* inherited */
	FileCrawlerT::RunBatch(in, status);
	
	/* write summary */
	if (is_root)
	{
		/* check */
		if (fFiles.Length() != fPassFail.Length()) throw eGeneralFail;

		/* open output stream */
		StringT file_name = "compare.out";
		file_name.ToNativePathName();
		file_name.Prepend(fRoot);
		ofstreamT summary;
		summary.open(file_name);
		if (!summary.is_open())
		{
			cout << "\n ComparatorT::RunBatch: ERROR: could not open file: " 
			     << file_name << endl;
			throw eGeneralFail;
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
		
		/* clear root name */
		fRoot.Clear();		
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* compare results against benchmarks */
bool ComparatorT::PassOrFail(ifstreamT& in) const
{
	/* path to the current directory */
	StringT path;
	path.FilePath(in.filename());
	
	/* job file name */
	StringT file_root = in.filename();
	file_root.Drop(path.StringLength());
	file_root.Root();
	
	/* look for results file in current directory */
	StringT current = path;
	current.Append(file_root, ".run");
	ifstreamT current_in(current);
	cout << "looking for current results: " << current << ": "
	     << ((current_in.is_open()) ? "found" : "not found") << '\n';
	if (!current_in.is_open()) return false;

	/* look for results file in benchmark directory */
	StringT benchmark(kBenchmarkDirectory);
	benchmark.ToNativePathName();
	benchmark.Prepend(path);
	benchmark.Append(file_root, ".run");
	ifstreamT bench_in(benchmark);
	cout << "looking for benchmark results: " << benchmark << ": "
	     << ((bench_in.is_open()) ? "found" : "not found") << '\n';
	if (!bench_in.is_open()) return false;

	/* init */
	StringT b_str, c_str;
	if (!bench_in.FindString("O U T P U T", b_str)) return false;
	if (!current_in.FindString("O U T P U T", c_str)) return false;
	
	/* block info */
	int b_group, c_group;
	double b_time, c_time;
	bool b_OK = ReadDataBlockInfo(bench_in, b_group, b_time);
	bool c_OK = ReadDataBlockInfo(current_in, c_group, c_time);

	/* compare blocks */
	while (b_OK && c_OK)
	{
		/* verify block info */
		if (b_group != c_group) {
			cout << "group number mismatch: " << b_group << " != " << c_group << '\n';
			return false;
		}
		else
			cout << "group: " << b_group << '\n';		
		if (fabs(b_time - c_time) > kSmall) {
			cout << "time mismatch: " << b_time << " != " << c_time << '\n';
			return false;
		}
		else
			cout << "time: " << b_time << '\n';
	
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
		if (!CompareDataBlocks(b_node_labels, b_node_data, c_node_labels, c_node_data)) {
			cout << "nodal data fails check" << '\n';
			return false;
		}
		else cout << "nodal data passes check" << '\n';

		/* read element data */
		ArrayT<StringT> b_element_labels;
		dArray2DT b_element_data;
		if (!ReadElementData(bench_in, b_element_labels, b_element_data)) {
			cout << "error reading element data from: " << bench_in.filename() << '\n';
			return false;
		}

		ArrayT<StringT> c_element_labels;
		dArray2DT c_element_data;
		if (!ReadElementData(current_in, c_element_labels, c_element_data)) {
			cout << "error reading element data from: " << current_in.filename() << '\n';
			return false;
		}

		/* compare element blocks */
		if (!CompareDataBlocks(b_element_labels, b_element_data, c_element_labels, c_element_data)) {
			cout << "element data fails check" << '\n';
			return false;
		}
		else cout << "element data passes check" << '\n';

		/* next block */
		b_OK = ReadDataBlockInfo(bench_in, b_group, b_time);
		c_OK = ReadDataBlockInfo(current_in, c_group, c_time);
	}

	/* termination */
	if (b_OK == c_OK) 
		return true;
	else 
		return false;
}

/* read data block header */
bool ComparatorT::ReadDataBlockInfo(ifstreamT& in, int& group, double& time) const
{
	StringT str;

	/* group number */
	if (!in.FindString("Group number", str)) return false;
	if (!str.Tail('=', group)) return false;
	
	/* time */
	if (!in.FindString("Time", str)) return false;
	if (!str.Tail('=', time)) return false;
	
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
		labels.Allocate(num_values + 1); // +1 for node number
		for (int i = 0; i < labels.Length(); i++)
			in >> labels[i];
		if (!in.good()) return false;
	
		/* read values */
		data.Allocate(num_nodes, num_values + 1); // +1 for node number
		in >> data;
	}
	return true;
}

/* read block of element data */
bool ComparatorT::ReadElementData(ifstreamT& in, ArrayT<StringT>& labels, dArray2DT& data) const
{
	/* advance nodal data */
	StringT str;
	if (!in.FindString("Element data:", str)) return false;

	/* get dimensions */
	int num_nodes, num_values;
	if (!in.FindString("Number of elements", str)) return false;
	str.Tail('=', num_nodes);
	if (!in.FindString("Number of values", str)) return false;
	str.Tail('=', num_values);
	
	if (num_values > 0)
	{
		/* read labels */
		labels.Allocate(num_values + 1); // +1 for node number
		for (int i = 0; i < labels.Length(); i++)
			in >> labels[i];
		if (!in.good()) return false;
	
		/* read values */
		data.Allocate(num_nodes, num_values + 1); // +1 for node number
		in >> data;
	}
	return true;
}

/* compare blocks - normalized by data set 1 */
bool ComparatorT::CompareDataBlocks(const ArrayT<StringT>& labels_1, const dArray2DT& data_1,
	const ArrayT<StringT>& labels_2, const dArray2DT& data_2) const
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

		double error_norm = 0.0, bench_norm = 0.0;
		int max_rel_error_index = -1, max_abs_error_index = -1;
		double max_rel_error_var = 0.0, max_abs_error_var = 0.0;
		for (int i = 0; i < data_1.MajorDim(); i++)
		{
			/* error */
			double abs_error = data_2(i,j) - data_1(i,j);
			double rel_error = (fabs(data_1(i,j)) > kSmall) ? abs_error/data_1(i,j) : 0.0;

			/* norms */
			bench_norm += data_1(i,j)*data_1(i,j);
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
			cout << "-" << '\n';
			
		/* maximum errors */
		cout << "         max abs error: " << max_abs_error_var << '\n';
		cout << "index of max abs error: " << max_abs_error_index + 1 << '\n';
		cout << "         max rel error: " << max_rel_error_var << '\n';
		cout << "index of max rel error: " << max_rel_error_index + 1 << '\n';
		
		/* running max */
		max_abs_error = (fabs(max_abs_error_var) > fabs(max_abs_error)) ? max_abs_error_var : max_abs_error;
		max_rel_error = (fabs(max_rel_error_var) > fabs(max_rel_error)) ? max_rel_error_var : max_rel_error;
	}

	/* assess results */
	if (fabs(max_abs_error) > fAbsTol) cout << "absolute error exceeds tolerance " << fAbsTol << ": " << max_abs_error << '\n';
	if (fabs(max_rel_error) > fRelTol) cout << "relative error exceeds tolerance " << fRelTol << ": " << max_rel_error << '\n';
//	if (fabs(max_abs_error) > fAbsTol || fabs(max_rel_error) > fRelTol)
	if (fabs(max_abs_error) > fAbsTol)
		return false;
	else
	{
		cout << "absolute error within tolerance " << fAbsTol << ": " << max_abs_error << '\n';
		cout << "relative error within tolerance " << fRelTol << ": " << max_rel_error << '\n';
		return true;
	}
}
