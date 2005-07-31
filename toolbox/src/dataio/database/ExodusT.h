/* $Id: ExodusT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: sawimme (12/04/1998)                                          */

#ifndef _EXODUS_T_H_
#define _EXODUS_T_H_

#include "Environment.h"

/* direct members */
#include "iArrayT.h"
#include "StringT.h"
#include "GeometryT.h"

/* forward declarations */
class iArray2DT;
class dArrayT;
class dArray2DT;

class ExodusT
{
public:

	/* constructors */
	ExodusT(ostream& out, int float_size = sizeof(double));

	/* destructor */
	~ExodusT(void);	

	/* opening/creating/closing files - returns the file ID */
	bool OpenRead(const StringT& filename);
	bool OpenWrite (const StringT& filename);
	bool Create(const StringT& filename, const StringT& title,
		ArrayT<StringT>& info, ArrayT<StringT>& QA, int dim, int nodes,
		int elem, int num_blks, int node_sets, int side_sets);
	void Close(void);
	int FileID(void) const;

	/* accessors */
	const StringT& Filename(void) const;
	int NumNodes(void) const;
	int NumDimensions(void) const;
	int NumNodeSets(void) const;
	int NumElementBlocks(void) const;
	int NumSideSets(void) const;
	void ElementBlockID(iArrayT& ID) const;
	void NodeSetID(iArrayT& ID) const;
	void SideSetID(iArrayT& ID) const;
	
	/* coordinates */
	void ReadCoordinates(dArray2DT& coords) const;
	void WriteCoordinates(const dArray2DT& coords, const iArrayT* node_map = NULL) const;
	void ReadNodeMap(iArrayT& node_map) const;

	/* element blocks */
	void ReadElementBlockDims(int block_ID, int& num_elems, int& num_elem_nodes) const;
	void ReadConnectivities(int block_ID, GeometryT::CodeT& code,
		iArray2DT& connects) const;
	void WriteConnectivities(int block_ID, GeometryT::CodeT code,
		const iArray2DT& connects, const iArrayT* elem_map = NULL);

	/* node sets */
	int  NumNodesInSet(int set_ID) const;
	void ReadNodeSet(int set_ID, iArrayT& nodes) const;
	void ReadNodeSets(const iArrayT& set_ID, iArrayT& nodes) const;
	void WriteNodeSet(int set_ID, const iArrayT& nodes) const;
	
	/* side sets */
	int  NumSidesInSet(int set_ID) const;
	void ReadSideSet(int set_ID, int& block_ID, iArray2DT& sides) const;
	void WriteSideSet(int set_ID, int block_ID, const iArray2DT& sides) const;

	/* write results */
	void WriteNodeLabels(const ArrayT<StringT>& labels) const;
	void WriteElementLabels(const ArrayT<StringT>& labels) const;
	void WriteGlobalLabels(const ArrayT<StringT>& labels) const;
	void WriteTime(int step, double time) const;
	void WriteNodalVariable(int step, int index, const dArrayT& fValues) const;
	void WriteElementVariable(int step, int block_ID, int index,
		const dArrayT& fValues) const;
	void WriteGlobalVariable(int step, const dArrayT& fValues) const;

	/* read results data */
	void ReadNodeLabels(ArrayT<StringT>& labels) const;
	void ReadElementLabels(ArrayT<StringT>& labels) const;
	void ReadGlobalLabels(ArrayT<StringT>& labels) const;
	int NumTimeSteps(void) const;
	void ReadTime(int step, double& time) const;
	void ReadNodalVariable(int step, int index, dArrayT& fValues) const;
	void ReadElementVariable(int step, int block_ID, int index,
		dArrayT& fValues) const;
	void ReadGlobalVariable(int step, dArrayT& fValues) const;

	/* convert global element numbers to block local numbers. assumes
	 * global element numbers are continuous within and between block
	 * and that the order of blocks is set by the global number. side
	 * sets may not include elements from more than one block */
	void GlobalToBlockElementNumbers(int& block_ID, iArrayT& elements) const;
	void BlockToGlobalElementNumbers(int  block_ID, iArrayT& elements) const;

	/* Read and echo Quality Assurance and Information data strings */
	void ReadQA(ArrayT<StringT>& records) const;
	void ReadInfo(ArrayT<StringT>& info_records) const;

protected:

	/* return the element name and number of output nodes for the given
	 * geometry and number of element nodes */
	void GetElementName(int elemnodes, GeometryT::CodeT code,
		StringT& elem_name, int& num_output_nodes) const;

	/* return the geometry code for the given element name */
	GeometryT::CodeT ToGeometryCode(const StringT& elem_name) const;

private:

	enum DimensionsT { MAX_QA_REC = 5, MAX_INFO = 5 };

	/* Write Quality Assurance data strings */
	void WriteQA(const ArrayT<StringT>& records) const;
	void WriteInfo(const ArrayT<StringT>& info_records) const;

	/* convert side set information to fe++ numbering convention */
	void ConvertSideSetOut(const char* elem_type, iArrayT& sides) const;
	void ConvertSideSetIn(const char* elem_type, iArrayT& sides) const;

	void ConvertElementNumbering (iArray2DT& conn, int code) const;

	/* labels */
	void WriteLabels(const ArrayT<StringT>& labels, const char *type) const;
	void ReadLabels(ArrayT<StringT>& labels, char *type) const;
			
	/* clear all parameter data */
	void Clear(void);
	
	/* process return values - (do_warning == 1) prints
	 * warning message for (code > 0) */
	void Try(const char* caller, int code, bool do_warning) const;

private:

	/* message stream */
	ostream& fOut;

	/* opening parameters */
	StringT file_name;
	int   exoid;
	int   comp_ws; // floating point size in program
	int   io_ws;   // floating point size in file data
	float version;

	/* initialization parameters */
	int   num_dim;
	int   num_nodes;
	int   num_elem;
	int   num_elem_blk;
	int   num_node_sets;
	int   num_side_sets;
};

/* inlines */

/* accessors */
inline int ExodusT::FileID(void) const { return exoid; }
inline const StringT& ExodusT::Filename(void) const { return file_name; }
inline int ExodusT::NumNodes(void) const { return num_nodes; }
inline int ExodusT::NumDimensions(void) const { return num_dim; }
inline int ExodusT::NumNodeSets(void) const { return num_node_sets; }
inline int ExodusT::NumElementBlocks(void) const { return num_elem_blk; }
inline int ExodusT::NumSideSets(void) const { return num_side_sets; }

#endif /* _EXODUS_T_H_ */
