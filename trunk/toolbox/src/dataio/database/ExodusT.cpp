/* $Id: ExodusT.cpp,v 1.20 2002-10-20 22:36:53 paklein Exp $ */
/* created: sawimme (12/04/1998)                                          */

#include "ExodusT.h"

/* ANSI headers */
#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "StringT.h"
#include "pArrayT.h"
#include "fstreamT.h"

#ifdef __ACCESS__ // with SEACAS support
#include "exodusII.h"

/* constructor for opening input or output file */

using namespace Tahoe;

ExodusT::ExodusT(ostream& message_out, int float_size):
	fOut(message_out),
	exoid(-1),
	comp_ws(float_size),
	io_ws(0)
{
	/* check */
	if (comp_ws != sizeof(float) && comp_ws != sizeof(double))
	{
		cout << "\n ExodusT::ExodusT: unexpected float size " << comp_ws << endl;
		throw ExceptionT::kBadInputValue;
	}

#if __option(extended_errorcheck)
	/* debugging messages */
	ex_opts(EX_VERBOSE);
#endif
}

/* destructor */
ExodusT::~ExodusT(void) { Close(); }

/* opening/closing files */
bool ExodusT::OpenRead(const StringT& filename)
{
	/* safe */
	Close();

	/* ex_open does not like to fail */
	ifstreamT test_in(filename);
	if (!test_in.is_open()) return false;

	/* open Exodus file */
	exoid = ex_open(filename, EX_READ, &comp_ws, &io_ws, &version);

	/* opened successfully */
	if (exoid < 0)
		return false;
	else
	{
		file_name = filename;

		/* check floating point number format */
		if (io_ws != sizeof(float) && io_ws != sizeof(double))
		{
			fOut << "\n ExodusT::OpenRead: known float size: " << io_ws;
			fOut << endl;
			throw ExceptionT::kBadInputValue;
		}
	
		/* get initialization data */
	 	ArrayT<char> title(MAX_LINE_LENGTH);
		Try("ExodusT::OpenRead",
			ex_get_init(exoid, title.Pointer(), &num_dim, &num_nodes, &num_elem,
			            &num_elem_blk, &num_node_sets, &num_side_sets),
			true);

		/* see if there are any quality assurance data strings */
//		ReadQA();

		/* may want to also read information records */
//		ReadInfo();

		return true;
	}
}

bool ExodusT::OpenWrite(const StringT& filename)
{
	/* safe */
	Close();

	/* ex_open does not like to fail */
	ifstreamT test_in(filename);
	if (!test_in.is_open()) return false;

	/* re-open Exodus file */
	// io_ws will be set to the size in the database, which should be double
	comp_ws = sizeof(double);
	exoid = ex_open(filename, EX_WRITE, &comp_ws, &io_ws, &version);
	if (exoid < 0)
		return false;
	else
	{
		file_name = filename;
		return true;
	}
}

bool ExodusT::Create(const StringT& filename, const StringT& title,
		ArrayT<StringT>& info, ArrayT<StringT>& QA, int dim, int nodes,
		int elem, int num_blks, int node_sets, int side_sets)
{
	/* safe */
	Close();

	/* create Exodus file */
	io_ws = comp_ws;
	if (strlen(filename) >= MAX_LINE_LENGTH)
		file_name.Take(filename, MAX_LINE_LENGTH - 1);
	else
		file_name = filename;
	exoid = ex_create(file_name, EX_CLOBBER, &comp_ws, &io_ws);
	if (exoid < 0)
		return false;
	{
		/* set dimensions */
		num_dim = dim;
		num_nodes = nodes;
		num_elem = elem;
		num_elem_blk = num_blks;
		num_node_sets = node_sets;
		num_side_sets = side_sets;
	
		/* write parameters to the database */
		Try("ExodusT::WriteParameters",
			ex_put_init(exoid, title, num_dim, num_nodes,
			       num_elem, num_elem_blk, num_node_sets, num_side_sets),
			true);

		/* write QA */
		WriteQA(QA);
		
		/* write info */
		WriteInfo(info);

		return true;
	}
}

void ExodusT::Close(void)
{
	if (exoid > 0)
	{
		ex_close(exoid);
		Clear();
	}
}

/* accessors */
void ExodusT::ElementBlockID(nArrayT<int>& ID) const
{
	/* check */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (ID.Length() != num_elem_blk) throw ExceptionT::kSizeMismatch;

	/* non-empty */
	if (num_elem_blk > 0)
	{
		/* access */	
		Try("ExodusT::ElementBlockID",
			ex_get_elem_blk_ids(exoid, ID.Pointer()),
			true);
	}
}

void ExodusT::NodeSetID(nArrayT<int>& ID) const
{
	/* check */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (ID.Length() != num_node_sets) throw ExceptionT::kSizeMismatch;

	/* non-empty */
	if (num_node_sets > 0)
	{
		/* access */	
		Try("ExodusT::NodeSetID",
			ex_get_node_set_ids(exoid, ID.Pointer()),
			true);
	}
}

void ExodusT::SideSetID(nArrayT<int>& ID) const
{
	/* check */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (ID.Length() != num_side_sets) throw ExceptionT::kSizeMismatch;

	/* non-empty */
	if (num_side_sets > 0)
	{
		/* access */	
		Try("ExodusT::SideSetID",
			ex_get_side_set_ids(exoid, ID.Pointer()),
			true);
	}
}

/* coordinates */
void ExodusT::ReadCoordinates(dArray2DT& coords) const
{
	/* checks */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (coords.MajorDim() != num_nodes ||
	    coords.MinorDim() != num_dim) throw ExceptionT::kSizeMismatch;

	dArray2DT xyz(num_dim, num_nodes);
	double* px = xyz(0);
	double* py = (num_dim > 1) ? xyz(1) : NULL;
	double* pz = (num_dim > 2) ? xyz(2) : NULL;
	
	/* read coordinate data */
	Try("ExodusT::ReadCoordinates",
		ex_get_coord(exoid, px, py, pz),
		true);

	/* copy to transpose */
	coords.Transpose(xyz);
}

void ExodusT::WriteCoordinates(const dArray2DT& coords,
	const nArrayT<int>* node_map) const
{
	/* checks */
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	dArrayT x(num_nodes), y(num_nodes), z(num_nodes);
	double *px = x.Pointer(), *py = y.Pointer(), *pz = z.Pointer();

	for (int i = 0; i < num_nodes; i++)
	  {
	    *px++ = coords(i,0);
	    if (num_dim > 1)
	      *py++ = coords(i,1);
	    else
	      *py++ = 0.0;
	    if (num_dim==3)
	      *pz++ = coords(i,2);
	    else
	      *pz++ = 0.0;
	  }

	Try("ExodusT::WriteCoordinates: ex_put_coord",
		ex_put_coord(exoid, x.Pointer(), y.Pointer(), z.Pointer()),
		true);

	/* optional node map */
	if (node_map && node_map->Length() > 0)
	{
		/* check */
		if (node_map->Length() != coords.MajorDim()) throw ExceptionT::kSizeMismatch;
	
		Try("ExodusT::WriteCoordinates: ex_put_node_map",
			ex_put_node_num_map(exoid, node_map->Pointer()),
			true);
	}
}

void ExodusT::ReadNodeMap(nArrayT<int>& node_map) const
{
	/* checks */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (node_map.Length() != num_nodes) throw ExceptionT::kSizeMismatch;

	Try("ExodusT::ReadNodeMap: ",
		ex_get_node_num_map(exoid, node_map.Pointer()),
		true);
}

/* element block */
void ExodusT::ReadElementBlockDims(int block_ID, int& num_elems, int& num_elem_nodes) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* read block parameters */
	ArrayT<char> type(MAX_STR_LENGTH);
	int num_attr;
	Try("ExodusT::ElementBlockDims",
		ex_get_elem_block(exoid, block_ID, type.Pointer(),
			&num_elems, &num_elem_nodes, &num_attr),
		true);
}

void ExodusT::ReadConnectivities(int block_ID, GeometryT::CodeT& code,
	iArray2DT& connects) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* read attributues */		
	char elem_type[MAX_STR_LENGTH];
	int num_elems;
	int num_elem_nodes;
	int num_attr;
	Try("ExodusT::ReadConnectivities: ex_get_elem_block",
		ex_get_elem_block(exoid, block_ID, elem_type, &num_elems, &num_elem_nodes, &num_attr),
		true);

	/* check dims */
	ReadElementBlockDims(block_ID, num_elems, num_elem_nodes);
	if (connects.MajorDim() != num_elems ||
	    connects.MinorDim() != num_elem_nodes)
	{
		cout << "\n ExodusT::ReadConnectivities: database dimensions of block ID "
		     << block_ID << " {" << num_elems << "," << num_elem_nodes
		     << "}\n"
		     <<   "     do not match the destination array {" << connects.MajorDim()
		     << "," << connects.MinorDim() << "}: " << file_name << endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* non-empty set */
	if (connects.MajorDim() > 0)
	{
		/* resolve geometry code */
		code = ToGeometryCode(elem_type);
		//NOTE: for some reason ExodusII does not store the element type,
		//      or any information, declared for empty groups

		/* read connectivity data */
		Try("ExodusT::ReadConnectivities: ex_get_elem_conn",
			ex_get_elem_conn(exoid, block_ID, connects.Pointer()),
			true);

		ConvertElementNumbering (connects, code);	
	}
	else
		code = GeometryT::kNone;
}

void ExodusT::WriteConnectivities(int block_ID, GeometryT::CodeT code,
	const iArray2DT& connects, const nArrayT<int>* elem_map)
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* check dimensions */
	if (connects.MajorDim() > num_elem)
	{
		cout << "\n ExodusT::WriteConnectivities: number of elements in block "
		     << block_ID << " exceeds\n"
		     <<   "     the total number of elements " << num_elem << endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* get element type name and output nodes */
	int num_output_nodes;
	StringT elem_type;
	GetElementName(connects.MinorDim(), code, elem_type, num_output_nodes);

	/* set connectivity size based on num_output_nodes */
	iArray2DT tempconn;
	if (num_output_nodes != connects.MinorDim())
	{
		tempconn.Dimension(connects.MajorDim(), num_output_nodes);
		for (int i = 0; i < tempconn.MajorDim(); i++)
	  		for (int j = 0; j < tempconn.MinorDim(); j++)
	    		tempconn (i,j) = connects (i,j);
	}
	else
		tempconn.Alias(connects);
	
	/* write element block attributes */
	int num_attr = 0; // set to zero for now
	Try ("ExodusT::WriteConnectivities: ex_put_elem_block",
		ex_put_elem_block(exoid, block_ID, elem_type.Pointer(), tempconn.MajorDim(),
			tempconn.MinorDim(), num_attr),
		true);

	/* non-empty */
	if (tempconn.MajorDim() > 0)
	{
		ConvertElementNumbering (tempconn, code);
	
		/* write connectivities */
		Try ("ExodusT::WriteConnectivities: ex_put_elem_conn",
			ex_put_elem_conn(exoid, block_ID, tempconn.Pointer()),
			true);
	}
		
	/* optional element map */
	if (elem_map && elem_map->Length() > 0)
	{
		/* check */
		if (elem_map->Length() != connects.MajorDim()) throw ExceptionT::kSizeMismatch;
	
		Try("ExodusT::WriteConnectivities: ex_put_elem_num_map",
			ex_put_elem_num_map(exoid, elem_map->Pointer()),
			true);
	}
}

/* node sets */
int ExodusT::NumNodesInSet(int set_ID) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* read set parameters */
	int num_set_nodes, num_dist;
	Try("ExodusT::NumNodesInSet",
		ex_get_node_set_param(exoid, set_ID, &num_set_nodes, &num_dist),
		true);
	
	return num_set_nodes;
}

void ExodusT::ReadNodeSet(int set_ID, nArrayT<int>& nodes) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* check dims */
	if (NumNodesInSet(set_ID) != nodes.Length()) throw ExceptionT::kSizeMismatch;

	/* non-empty */
	if (nodes.Length() > 0)
	{
		/* read node set */
		Try("ExodusT::ReadNodeSet",
			ex_get_node_set(exoid, set_ID, nodes.Pointer()),
			true);
	}
}

void ExodusT::ReadNodeSets(const nArrayT<int>& set_ID, nArrayT<int>& nodes) const
{
	/* get total number of nodes */
	int num_set_nodes = 0;
	for (int i = 0; i < set_ID.Length(); i++)
		num_set_nodes += NumNodesInSet(set_ID[i]);

	/* allocate */
	nodes.Dimension(num_set_nodes);

	/* read */
	int count = 0;
	nArrayT<int> tmp;
	for (int j = 0; j < set_ID.Length(); j++)
	{
		tmp.Set(NumNodesInSet(set_ID[j]), nodes.Pointer(count));
		ReadNodeSet(set_ID[j], tmp);
		count += tmp.Length();
	}
}

void ExodusT::WriteNodeSet(int set_ID, const nArrayT<int>& nodes) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	Try("ExoduT::WriteNodeSet: ex_put_node_set_param",
		ex_put_node_set_param(exoid, set_ID, nodes.Length(), nodes.Length()),
		true);

	/* non-empty */
	if (nodes.Length() > 0)
	{
		Try("ExodusT::WriteNodeSet: ex_put_node_set",
			ex_put_node_set(exoid, set_ID, nodes.Pointer()),
			true);
	}
}

/* side sets */
int ExodusT::NumSidesInSet(int set_ID) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* read set parameters */
	int num_sides, num_dist;
	Try("ExodusT::NumSidesInSet",
		ex_get_side_set_param(exoid, set_ID, &num_sides, &num_dist),
		true);
	
	return num_sides;
}

void ExodusT::ReadSideSet(int set_ID, int& block_ID, iArray2DT& sides) const
{
	/* checks */
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	if (NumSidesInSet(set_ID) != sides.MajorDim()) throw ExceptionT::kSizeMismatch;

	/* non-empty set */
	if (sides.MajorDim() > 0)
	{
		/* read data */
		iArray2DT temp(2, sides.MajorDim());
		Try("ExodusT::ReadSideSet: ex_get_side_set",
			ex_get_side_set(exoid, set_ID, temp(0), temp(1)),
			true);

		/* convert to element block numbering */
		nArrayT<int> elements;
		temp.RowAlias(0, elements);
		GlobalToBlockElementNumbers(block_ID, elements);

		/* found valid element block */
		if (block_ID > 0)
		{
			/* read block parameters */
			char type[MAX_STR_LENGTH];
			int nel, nen, num_attr;
			Try("ExodusT::ReadSideSet: ex_get_elem_block",
				ex_get_elem_block(exoid, block_ID, type,
					&nel, &nen, &num_attr),
				true);	

			/* convert facet numbering */
			nArrayT<int> facets;
			temp.RowAlias(1, facets);
			ConvertSideSetOut(type, facets);

			/* transpose data */
			sides.Transpose(temp);
		}
		else
		{
			cout << "\n ExodusT::ReadSideSet: could not read side set " << set_ID << endl;
			sides = 0;
		}
	}
	else
		block_ID = 0;
}

void ExodusT::WriteSideSet(int set_ID, int block_ID, const iArray2DT& sides) const
{
	/* check */
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* write parameters */
	Try("ExodusT::WriteSideSet: ex_put_side_set",
		ex_put_side_set_param(exoid, set_ID, sides.MajorDim(), 0),
		true);	

	/* non-empty set */
	if (sides.MajorDim() > 0)
	{
		/* transposed data */
		iArray2DT temp(2, sides.MajorDim());
		temp.Transpose(sides);

		/* read block parameters */
		char type[MAX_STR_LENGTH];
		int nel, nen, num_attr;
		Try("ExodusT::WriteSideSet: ex_get_elem_block",
			ex_get_elem_block(exoid, block_ID, type,
				&nel, &nen, &num_attr),
			true);

		/* convert facet numbering */
		nArrayT<int> facets;
		temp.RowAlias(1, facets);
		ConvertSideSetIn(type, facets);
	
		/* convert to global numbering */
		nArrayT<int> elements;
		temp.RowAlias(0, elements);
		BlockToGlobalElementNumbers(block_ID, elements);
	
		/* write */
		Try("ExodusT::WriteSideSet: ex_put_side_set",
			ex_put_side_set(exoid, set_ID, temp(0), temp(1)),
			true);
	}
}

/* variable results */
void ExodusT::WriteLabels(const ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	const char vartypes [3] = { 'n', 'e', 'g' };
	char type = vartypes[t];

	/* type = "n" or "e" or "g" for nodal, element or global */
	if (labels.Length() > 0)
	{
		/* write number of variables */
		Try("ExodusT::WriteLabels: ex_put_var_param",
			ex_put_var_param(exoid, &type, labels.Length()),
			true);

		/* change array of strings to array of char* */
		ArrayT<char*> var_names(labels.Length());
		for (int i = 0; i < labels.Length(); i++)
			var_names[i] = labels[i].Pointer();

		/* write variable names */
		Try("ExodusT::WriteLabels: ex_put_var_names",
			ex_put_var_names(exoid, &type, labels.Length(), var_names.Pointer()),
			true);
	}
}	

void ExodusT::WriteTime(int step, double time) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* write time value, each time value corresponds to a time_step
	 * the first time_step must be one and be incremented by one */
	Try("ExodusT::WriteTime",
		ex_put_time(exoid, step, &time),
		true);
}

void ExodusT::WriteNodalVariable(int step, int index, const dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* don't write empty lists */
	if (fValues.Length() == 0) return;

	/* the time_step must correspond to the time_value of the printed increment
	 * index corresponds to the variable name list */
	Try("ExodusT::WriteNodalVariable",
		ex_put_nodal_var(exoid, step, index, fValues.Length(), fValues.Pointer()),
		true);
}

void ExodusT::WriteElementVariable(int step, int block_ID, int index,
		const dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* don't write empty lists */
	if (fValues.Length() == 0) return;

	/* the time_step must correspond to the time_value of the printed increment
	 * index corresponds to the variable name list */
	Try("ExodusT::WriteElementVariable",
		ex_put_elem_var (exoid, step, index, block_ID, fValues.Length(), fValues.Pointer()),
		true);
}

void ExodusT::WriteGlobalVariable(int step, const dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* don't write empty lists */
	if (fValues.Length() == 0) return;

	/* the time_step must correspond to the time_value of the printed increment */
	Try("ExodusT::WriteGlobalVariable",
		ex_put_glob_vars (exoid, step, fValues.Length(), fValues.Pointer()),
		true);
}

/* read results data */
void ExodusT::ReadLabels(ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	const char vartypes [3] = { 'n', 'e', 'g' };
	char type = vartypes[t];

	/* read number of variables */
	int num_labels;
	Try("ExodusT::ReadLabels: ex_get_var_param",
		ex_get_var_param(exoid, &type, &num_labels),
		true);

	/* type = "n" or "e" or "g" for nodal, element or global */
	if (num_labels > 0)
	{
		/* allocate array of char* */
		pArrayT<char*> var_names(num_labels);
		for (int i = 0; i < num_labels; i++)
		{
			char* str = new char[MAX_STR_LENGTH];
			if (!str) throw ExceptionT::kOutOfMemory;
			var_names[i] = str;
		}

		/* read variable names */
		Try("ExodusT::ReadLabels: ex_get_var_names",
			ex_get_var_names(exoid, &type, var_names.Length(), var_names.Pointer()),
			true);
			
		/* copy in */
		labels.Dimension(num_labels);
		for (int j = 0; j < num_labels; j++)
		{	
			char* str = var_names[j];
			labels[j] = str;
		}
	}
}	

int ExodusT::NumTimeSteps(void) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* return values */
	int i_ret;
	float f_ret;
	char c_ret;

	/* inquire */
	Try("ExodusT::NumTimeSteps",
		ex_inquire(exoid, EX_INQ_TIME, &i_ret, &f_ret, &c_ret),
		true);

	return i_ret;
}

void ExodusT::ReadTime(int step, double& time) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* inquire */
	Try("ExodusT::ReadTime",
		ex_get_time(exoid, step, &time),
		true);
}

int ExodusT::NumVariables (ExodusT::VariableTypeT t) const
{
	const char vartypes [3] = { 'n', 'e', 'g' };
	char type = vartypes[t];
	int num_labels;
	Try("ExodusT::NumNodeVariables: ex_get_var_param",
		ex_get_var_param(exoid, &type, &num_labels),
		true);
	return num_labels;
}

void ExodusT::ReadNodalVariable(int step, int index, dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* don't try to read empty lists */
	if (fValues.Length() == 0) return;

	/* the time_step must correspond to the time_value of the printed increment
	 * index corresponds to the variable name list */
	Try("ExodusT::ReadNodalVariable",
		ex_get_nodal_var(exoid, step, index, fValues.Length(), fValues.Pointer()),
		true);
}

void ExodusT::ReadElementVariable(int step, int block_ID, int index,
	dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;
	
	/* check dimensions */
	int num_elems, num_elem_nodes;
	ReadElementBlockDims(block_ID, num_elems, num_elem_nodes);
	if (num_elems != fValues.Length()) throw ExceptionT::kSizeMismatch;

	/* don't read from empty block */
	if (num_elems == 0) return;

	/* the time_step must correspond to the time_value of the printed increment
	 * index corresponds to the variable name list */
	Try("ExodusT::WriteElementVariable",
		ex_get_elem_var (exoid, step, index, block_ID, fValues.Length(), fValues.Pointer()),
		true);
}

void ExodusT::ReadGlobalVariable(int step, dArrayT& fValues) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	/* don't read empty lists */
	if (fValues.Length() == 0) return;

	/* the time_step must correspond to the time_value of the printed increment */
	Try("ExodusT::WriteGlobalVariable",
		ex_get_glob_vars (exoid, step, fValues.Length(), fValues.Pointer()),
		true);
}

void ExodusT::ReadQA(ArrayT<StringT>& qa_records) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	int num_qa_rec;
	Try("ExodusBaseT::ReadInit ex_inq_qa",
	    ex_inquire(exoid, EX_INQ_QA, &num_qa_rec, NULL, NULL),
	    1);

	/* set up QA strings */
/* DEC will not allow allocation based on passed constant */
if (num_qa_rec > MAX_QA_REC)
	  cout << "\nExodusBaseT::ReadQA, num records > MAX_QA_REC\n";
char *recs[MAX_QA_REC][4];
	for (int i = 0; i < num_qa_rec; i++)
	  for (int j = 0; j < 4; j++)
	    recs[i][j] = new char[MAX_STR_LENGTH];

	/* read records */
	Try("ExodusBaseT::ReadQA", ex_get_qa(exoid, recs), 1);
	
	int m=0;
	qa_records.Dimension(4*num_qa_rec);
	for (int ii=0; ii < num_qa_rec; ii++)
	  for (int j=0; j < 4; j++)
	    qa_records[m++] = recs[ii][j];
}

void ExodusT::ReadInfo(ArrayT<StringT>& info_records) const
{
	if (exoid < 0) throw ExceptionT::kGeneralFail;

	int num_info;
	Try("ExodusT::ReadInfo ex_inq_info",
		ex_inquire(exoid, EX_INQ_INFO, &num_info, NULL, NULL),
		1);

/* DEC will not allow allocation based on passed constant */
if (num_info > MAX_INFO)
	  cout << "\nExodusT::ReadInfo, num records > MAX_INFO\n";
	char *info[MAX_INFO];
	for (int i = 0; i < num_info; i++)
	  info[i] = new char[MAX_LINE_LENGTH];

	/* read records */
	Try("ExodusT::ReadInfo", ex_get_info(exoid, info), 1);

	info_records.Dimension(num_info);
	for (int j=0; j < num_info; j++)
	  info_records[j] = info[j];
}

/* convert global element numbers to block local numbers. assumes
* global element numbers are continuous within and between block
* and that the order of blocks is set by the global number. side
* sets may not include elements from more than one block */
void ExodusT::GlobalToBlockElementNumbers(int& block_ID, nArrayT<int>& elements) const
{
	int start_num = elements[0];
	
	/* get element block ID's */
	nArrayT<int> all_block_ID(num_elem_blk);
	ex_get_elem_blk_ids(exoid, all_block_ID.Pointer());
	
	/* find element number offset */
	int curr_blk = 0;
	int shift = 0;
	int next_shift = 0;
	for (int i = 0; i < num_elem_blk && next_shift < start_num; i++)
	{
		/* fetch block dimensions */
		int num_elems, num_elem_nodes;
		ReadElementBlockDims(all_block_ID[i], num_elems, num_elem_nodes);
	
		/* next numbering shift */
		curr_blk = i;
		shift = next_shift;
		next_shift += num_elems;
	}
	
	/* shift to block local numbers */
	elements -= shift;
	
	/* check */
	block_ID = all_block_ID[curr_blk];
	int num_elems, num_elem_nodes;
	ReadElementBlockDims(block_ID, num_elems, num_elem_nodes);
	int min, max;
	elements.MinMax(min, max);
	if (min < 1 || max > num_elems)
	{
		cout << "\n ExodusT::BlockElementNumbers: side set specification error:\n";
		cout <<   "     element number {min, max} = {" << min << "," << max;
		cout << "} exceeds the bounds\n";
		cout <<   "     of element group " << block_ID + 1 << '\n';
		cout <<   "     Associated element group could not be determined" << endl;
		block_ID = -1;
	}
}

void ExodusT::BlockToGlobalElementNumbers(int block_ID, nArrayT<int>& elements) const
{
	/* get element block ID's */
	nArrayT<int> all_block_ID(num_elem_blk);
	ex_get_elem_blk_ids(exoid, all_block_ID.Pointer());

	/* find element number shift */
	int shift = 0;
	for (int i = 0; i < num_elem_blk && all_block_ID[i] != block_ID; i++)
	{
		/* fetch block dimensions */
		int num_elems, num_elem_nodes;
		ReadElementBlockDims(all_block_ID[i], num_elems, num_elem_nodes);
	
		/* accumulate shift */
		shift += num_elems;
	}
	
	/* apply shift */
	elements += shift;
}


/*************************************************************************
* Protected
*************************************************************************/

void ExodusT::WriteQA(const ArrayT<StringT>& qa_records_) const
{
	/* truncate QA records */
	ArrayT<StringT> qa_records(qa_records_.Length());
		for (int k = 0; k < qa_records.Length(); k++)
			if (strlen(qa_records_[k]) >= MAX_STR_LENGTH)
				qa_records[k].Take(qa_records_[k], MAX_STR_LENGTH - 1);
			else
				qa_records[k] = qa_records_[k];
			

	/* DEC will not allow allocation based on passed constant */
	int num_recs = qa_records.Length()/4;
	if (num_recs > MAX_QA_REC)
		cout << "\nExodusT::WriteQA, num records > MAX_QA_REC\n";
	char *recs[MAX_QA_REC][4];
	int m=0;
	for (int i=0; i < num_recs; i++)
		for (int j=0; j < 4; j++)
			recs[i][j] = qa_records[m++].Pointer();
	/* you must send char *(*)[] */
	Try("ExodusT::WriteQA", ex_put_qa(exoid, num_recs, recs), 1);
}

void ExodusT::WriteInfo(const ArrayT<StringT>& info_records_) const
{
	/* truncate info records */
	ArrayT<StringT> info_records(info_records_.Length());
	for (int k = 0; k < info_records.Length(); k++)
		if (strlen(info_records_[k]) >= MAX_LINE_LENGTH)
			info_records[k].Take(info_records_[k], MAX_LINE_LENGTH - 1);
		else
			info_records[k] = info_records_[k];
	
	/* DEC will not allow allocation based on passed constant */
	int num_recs = info_records.Length();
	if (num_recs > MAX_INFO)
		cout << "\nExodusT::WriteInfo, num records > MAX_INFO\n";
	char *recs[MAX_INFO];
	for (int i=0; i < num_recs; i++)
	  recs[i] = info_records[i++].Pointer();

	/* write info_records */
	if (info_records.Length() > 0)
	  Try("ExodusT::WriteInfo",
	      ex_put_info(exoid, info_records.Length(), recs),
	      1);
}

/* return the element name and number of output nodes for the given
* geometry and number of element nodes */
void ExodusT::GetElementName(int elemnodes, GeometryT::CodeT code,
	StringT& elem_name, int& num_output_nodes) const
{
	switch (code)
	{
		case GeometryT::kPoint:
	    	    if (num_dim == 1)
	      		elem_name = "CIRCLE";
			else
				elem_name = "SPHERE";	
			num_output_nodes = 1;
		    break;

	        case GeometryT::kLine:
		        elem_name = "BEAM";
		        num_output_nodes = (elemnodes < 3) ? 2 : 3;
		        break;

		case GeometryT::kTriangle:
			elem_name =  "TRIANGLE";
			num_output_nodes = (elemnodes < 6) ? 3 : 6;
	    	        break;

		case GeometryT::kQuadrilateral:
			elem_name =  "QUAD";
			num_output_nodes = (elemnodes < 8) ? 4 : 8;
			break;

		case GeometryT::kHexahedron:
	    	elem_name = "HEX";
			num_output_nodes = (elemnodes < 20) ? 8 : 20;
			break;

		case GeometryT::kTetrahedron:
			elem_name =  "TETRA";
			num_output_nodes = (elemnodes < 10) ? 4 : 10;
			break;

		case GeometryT::kPentahedron:
			elem_name = "WEDGE";
	    	num_output_nodes = (elemnodes < 15) ? 6 : 15;
			break;

	    //case GeometryT::kSHELL:
	    //elem_name = "SHELL";

		default:
			cout << "\nExodusT::GetElementName cannot find name\n\n";
			throw ExceptionT::kGeneralFail;
	  }
}

/* return the geometry code for the given element name */
GeometryT::CodeT ExodusT::ToGeometryCode(const StringT& elem_name) const
{
	const char *elem_names[8] = {
		"CIRCLE",
		"SPHERE",
	        "BEAM",
		"TRIANGLE",
		"QUAD",
		"HEX",
		"TETRA",
		"WEDGE"};
	GeometryT::CodeT geom_codes[8] = {
		GeometryT::kPoint,
		GeometryT::kPoint,
		GeometryT::kLine,
		GeometryT::kTriangle,
		GeometryT::kQuadrilateral,
		GeometryT::kHexahedron,
		GeometryT::kTetrahedron,
		GeometryT::kPentahedron};

	StringT elem_name_upper(elem_name);
	elem_name_upper.ToUpper();
	
	GeometryT::CodeT code = GeometryT::kNone;
	for (int i = 0; i < 7 && code == GeometryT::kNone; i++)
		if (strncmp(elem_name_upper, elem_names[i], 3) == 0)
			code = geom_codes[i];

	/* could not resolve element name */
	if (code == GeometryT::kNone)
	{
		cout << "\n ExodusT::ToGeometryCode: code not resolve element name "
		     << elem_name << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	return code;
}

/*************************************************************************
* Private
*************************************************************************/

/* convert side set information to fe++ numbering convention */
void ExodusT::ConvertSideSetOut(const char* elem_type, nArrayT<int>& sides) const
{
	if (strncmp(elem_type, "HEX", 3) == 0)
	{
		/* facet numbering map: tahoe[exodusII] */
		int hex_map[7] = {-1,3,4,5,6,1,2};
		
		int*  p = sides.Pointer();
		int dim = sides.Length();
		for (int i = 0; i < dim; i++)
		{
			*p = hex_map[*p];
			p++;
		}
	}
}

void ExodusT::ConvertSideSetIn(const char* elem_type, nArrayT<int>& sides) const
{
	if (strncmp(elem_type, "HEX", 3) == 0)
	{
		/* facet numbering map: exodusII[tahoe] */
		int hex_map[7] = {-1,5,6,1,2,3,4};
		
		int*  p = sides.Pointer();
		int dim = sides.Length();
		for (int i = 0; i < dim; i++)
		{
			*p = hex_map[*p];
			p++;
		}
	}
}

/* convert element numbering to/from tahoe/ensight/abaqus convention */
void ExodusT::ConvertElementNumbering (iArray2DT& conn, int fcode) const
{
int convert = 0;
int dimension = 4;
int start = 0;

/* hex20 is divided into 5 sets of 4 numbers */
if (fcode == GeometryT::kHexahedron && conn.MinorDim() >= 20)
{
convert = 1;
dimension = 4;
start = 12;
}
/* wedge15 is divided into 5 sets of 3 numbers */
if (fcode == GeometryT::kPentahedron && conn.MinorDim() >= 15)
{
convert = 1;
dimension = 3;
start = 9;
}

/* basically swapping last two sets of numbers, set length = dimension */
if (convert)
{
nArrayT<int> temp (2*dimension);
for (int i=0; i < conn.MajorDim(); i++)
	{
	  /* transfer to temp space */
	  int *ptemp = temp.Pointer();
	  int *pconn = conn (i) + start;
	  for (int j=0; j < 2*dimension; j++)
	    *ptemp++ = *(pconn + j);

	  /* rewrite connectivity */
	  for (int k=0; k < dimension; k++)
	    *pconn++ = temp [k + dimension];
	  for (int k2=0; k2 < dimension; k2++)
	    *pconn++ = temp [k2];
	}
}
}

/* clear all parameter data */
void ExodusT::Clear(void)
{
	file_name = "\0";

	/* opening parameters */
	exoid   =-1;
	io_ws   = 0;
	version = 0.0;

	/* initialization parameters */
	num_dim       = 0;
	num_nodes     = 0;
	num_elem      = 0;
	num_elem_blk  = 0;
	num_node_sets = 0;
	num_side_sets = 0;
}

/* process return values - (do_warning == 1) implies print
* warning message for (code > 0) */
void ExodusT::Try(const char* caller, int code, bool do_warning) const
{
	if (code < 0)
	{
		fOut << "\n " << caller << ": returned error: " << code << '\n'
		     <<   "     file: " << file_name << endl;
		throw ExceptionT::kGeneralFail;
	}

	if (code > 0 && do_warning)
		fOut << "\n " << caller << ": returned warning: " << code << '\n'
		     <<   "     file: " << file_name << endl;
}
#else
#pragma warn_unusedarg off
/* constructor for opening input or output file */
ExodusT::ExodusT(ostream& out, int float_size):
	fOut(out),
	exoid(-1),
	comp_ws(float_size),
	io_ws(0)
{ }
ExodusT::~ExodusT(void) { }
bool ExodusT::OpenRead(const StringT& filename) 
{
  cout << "\n ExodusT::OpenRead: file format not available" << endl;
  throw ExceptionT::kGeneralFail;
  return false; 
}
bool ExodusT::OpenWrite(const StringT& filename)
{
  cout << "\n ExodusT::OpenWrite: file format not available" << endl;
  throw ExceptionT::kGeneralFail;
  return false; 
}
bool ExodusT::Create(const StringT& filename, const StringT& title,
		ArrayT<StringT>& info, ArrayT<StringT>& QA, int dim, int nodes,
		int elem, int num_blks, int node_sets, int side_sets)
{
  cout << "\n ExodusT::OpenCreate: file format not available" << endl;
  throw ExceptionT::kGeneralFail;
  return false; 
}
void ExodusT::Close(void) { throw ExceptionT::kGeneralFail; }
void ExodusT::ElementBlockID(nArrayT<int>& ID) const { throw ExceptionT::kGeneralFail; }
void ExodusT::NodeSetID(nArrayT<int>& ID) const { throw ExceptionT::kGeneralFail; }
void ExodusT::SideSetID(nArrayT<int>& ID) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadCoordinates(dArray2DT& coords) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteCoordinates(const dArray2DT& coords, const nArrayT<int>* node_map) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadNodeMap(nArrayT<int>& node_map) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadElementBlockDims(int block_ID, int& num_elems, int& num_elem_nodes) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadConnectivities(int block_ID, GeometryT::CodeT& code, iArray2DT& connects) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteConnectivities(int block_ID, GeometryT::CodeT code, const iArray2DT& connects, const nArrayT<int>* elem_map) { throw ExceptionT::kGeneralFail; }
int ExodusT::NumNodesInSet(int set_ID) const { return 0; }
void ExodusT::ReadNodeSet(int set_ID, nArrayT<int>& nodes) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadNodeSets(const nArrayT<int>& set_ID, nArrayT<int>& nodes) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteNodeSet(int set_ID, const nArrayT<int>& nodes) const { throw ExceptionT::kGeneralFail; }
int ExodusT::NumSidesInSet(int set_ID) const { return 0; }
void ExodusT::ReadSideSet(int set_ID, int& block_ID, iArray2DT& sides) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteSideSet(int set_ID, int block_ID, const iArray2DT& sides) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteTime(int step, double time) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteNodalVariable(int step, int index, const dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteElementVariable(int step, int block_ID, int index, const dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteGlobalVariable(int step, const dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadLabels(ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteLabels(const ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const { throw ExceptionT::kGeneralFail; }
int ExodusT::NumTimeSteps(void) const { return 0; }
void ExodusT::ReadTime(int step, double& time) const { throw ExceptionT::kGeneralFail; }
int ExodusT::NumVariables (ExodusT::VariableTypeT t) const { return 0; }
void ExodusT::ReadNodalVariable(int step, int index, dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadElementVariable(int step, int block_ID, int index, dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadGlobalVariable(int step, dArrayT& fValues) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadQA(ArrayT<StringT>& qa_records) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ReadInfo(ArrayT<StringT>& info_records) const { throw ExceptionT::kGeneralFail; }
void ExodusT::GlobalToBlockElementNumbers(int& block_ID, nArrayT<int>& elements) const { throw ExceptionT::kGeneralFail; }
void ExodusT::BlockToGlobalElementNumbers(int block_ID, nArrayT<int>& elements) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteQA(const ArrayT<StringT>& qa_records) const { throw ExceptionT::kGeneralFail; }
void ExodusT::WriteInfo(const ArrayT<StringT>& info_records) const { throw ExceptionT::kGeneralFail; }
void ExodusT::GetElementName(int elemnodes, GeometryT::CodeT code, StringT& elem_name, int& num_output_nodes) const { throw ExceptionT::kGeneralFail; }
GeometryT::CodeT ExodusT::ToGeometryCode(const StringT& elem_name) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ConvertSideSetOut(const char* elem_type, nArrayT<int>& sides) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ConvertSideSetIn(const char* elem_type, nArrayT<int>& sides) const { throw ExceptionT::kGeneralFail; }
void ExodusT::ConvertElementNumbering (iArray2DT& conn, int fcode) const { throw ExceptionT::kGeneralFail; }
void ExodusT::Clear(void) { throw ExceptionT::kGeneralFail; }
void ExodusT::Try(const char* caller, int code, bool do_warning) const { throw ExceptionT::kGeneralFail; }
#pragma warn_unusedarg reset
#endif /* __ACCESS__ */
