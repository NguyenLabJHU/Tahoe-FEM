/* $Id: ExodusT.h,v 1.7 2002-07-05 22:26:25 paklein Exp $ */
/* created: sawimme (12/04/1998) */

#ifndef _EXODUS_T_H_
#define _EXODUS_T_H_

#include "Environment.h"

/* direct members */
#include "nArrayT.h"
#include "StringT.h"
#include "GeometryT.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;
class dArrayT;
class dArray2DT;

/** interface for reading and writing ExodusII files */
class ExodusT
{
public:

  /** variable type for reading variables and labels */
  enum VariableTypeT { kNode=0, kElement, kGlobal };

  /** constructors */
  ExodusT(ostream& message_out, int float_size = sizeof(double));

  /** destructor */
  ~ExodusT(void);	

  /** opens an existing file for reading */
  bool OpenRead(const StringT& filename);

  /** opens an existing file for appending/updating */
  bool OpenWrite (const StringT& filename);

  /** create a new output database 
   * \param filename name of file created
   * \param title title of simulation
   * \param info optional destriptors
   * \param QA quality assurance strings
   * \param dim spatial degrees of freedom
   * \param nodes total number of nodes in file
   * \param elem total number of elements in file
   * \param num_blks number of element blocks
   * \param node_sets number of node sets
   * \param side_sets number of side sets */
  bool Create(const StringT& filename, const StringT& title,
	      ArrayT<StringT>& info, ArrayT<StringT>& QA, int dim, int nodes,
	      int elem, int num_blks, int node_sets, int side_sets);

  /** closes the database and clears parameters */
  void Close(void);

  int FileID(void) const; /**< returns the File ID */
  const StringT& Filename(void) const;/**< returns the filename */
  int NumNodes(void) const;/**< returns the total number of files */
  int NumDimensions(void) const;/**< returns the spatial degree of freedom */
  int NumNodeSets(void) const;/**< returns the number of node sets */
  int NumElementBlocks(void) const;/**< returns the number of element blocks */
  int NumSideSets(void) const;/**< returns the number of side sets */
  void ElementBlockID(nArrayT<int>& ID) const;/**< returns a list of element block ID numbers */
  void NodeSetID(nArrayT<int>& ID) const;/**< returns a list of node set ID numbers */
  void SideSetID(nArrayT<int>& ID) const;/**< returns a list of side set ID numbers */
	
  void ReadCoordinates(dArray2DT& coords) const; /**< returns the coordinates */
  void WriteCoordinates(const dArray2DT& coords, const nArrayT<int>* node_map = NULL) const; /**< write coordinates and node map (list of node IDs) */
  void ReadNodeMap(nArrayT<int>& node_map) const; /**< returns list of node IDs */

  /** returns dimensions of element block */
  void ReadElementBlockDims(int block_ID, int& num_elems, int& num_elem_nodes) const; 
  /** read an element block and its geometry code */
  void ReadConnectivities(int block_ID, GeometryT::CodeT& code, iArray2DT& connects) const;
  /** write an element block and its element map (list of element IDs) */
  void WriteConnectivities(int block_ID, GeometryT::CodeT code, const iArray2DT& connects, const nArrayT<int>* elem_map = NULL); 

  /** return the number of nodes in a node set */
  int  NumNodesInSet(int set_ID) const;
  /** return the list of nodes in the set, consecutive numbering, starting at 1 */
  void ReadNodeSet(int set_ID, nArrayT<int>& nodes) const;
  /** return the list of nodes in multiple sets */
  void ReadNodeSets(const nArrayT<int>& set_ID, nArrayT<int>& nodes) const;
  /** write a node set, use consecutive numbering starting at 1 */
  void WriteNodeSet(int set_ID, const nArrayT<int>& nodes) const;
	
  /** return number of element facets in a side set */
  int  NumSidesInSet(int set_ID) const;
  /** read the side set and the element block the facets are contained within */
  void ReadSideSet(int set_ID, int& block_ID, iArray2DT& sides) const;
  /** write a side set that is contained within a specified element block */
  void WriteSideSet(int set_ID, int block_ID, const iArray2DT& sides) const;

  /** write variable labels */
  void WriteLabels(const ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const;
  /** write time at the print step */
  void WriteTime(int step, double time) const;
  /** write data for one node variable (index) at a certain time step for all nodes */
  void WriteNodalVariable(int step, int index, const dArrayT& fValues) const;
  /** write data for one element variable (index) at a certain time step for a certain element block */
  void WriteElementVariable(int step, int block_ID, int index, const dArrayT& fValues) const;
  /** write data for one global variable at a certain time step */
  void WriteGlobalVariable(int step, const dArrayT& fValues) const;

  /** read variable labels */
  void ReadLabels(ArrayT<StringT>& labels, ExodusT::VariableTypeT t) const;
  /** return the number of time steps */
  int NumTimeSteps(void) const;
  /** return time step value */
  void ReadTime(int step, double& time) const;
  /** return the number of variables for the type specified */
  int NumVariables (ExodusT::VariableTypeT t) const;
  /** read data for one node variable (index) at a certain time step for all nodes */
  void ReadNodalVariable(int step, int index, dArrayT& fValues) const;
  /** read data for one element variable (index) at a certain time step for a certain element block */
  void ReadElementVariable(int step, int block_ID, int index, dArrayT& fValues) const;
  /** read data for one global variable at a certain time step */
  void ReadGlobalVariable(int step, dArrayT& fValues) const;

  /** convert global element numbers to block local numbers. assumes
   * global element numbers are continuous within and between block
   * and that the order of blocks is set by the global number. side
   * sets may not include elements from more than one block */
  void GlobalToBlockElementNumbers(int& block_ID, nArrayT<int>& elements) const;
  /** convert block local numbers to global element numbers. assumes
   * global element numbers are continuous within and between block
   * and that the order of blocks is set by the global number. side
   * sets may not include elements from more than one block */
  void BlockToGlobalElementNumbers(int  block_ID, nArrayT<int>& elements) const;

  /** Read and echo Quality Assurance and Information data strings */
  void ReadQA(ArrayT<StringT>& records) const;
  /** read optional descriptions */
  void ReadInfo(ArrayT<StringT>& info_records) const;

protected:

  /** return the element name and number of output nodes for the given
   * geometry and number of element nodes */
  void GetElementName(int elemnodes, GeometryT::CodeT code, StringT& elem_name, int& num_output_nodes) const;

  /** return the geometry code for the given element name */
  GeometryT::CodeT ToGeometryCode(const StringT& elem_name) const;

private:

  /** maximum limits */
  enum DimensionsT { MAX_QA_REC = 5, MAX_INFO = 5 };

  /** Write Quality Assurance data strings */
  void WriteQA(const ArrayT<StringT>& records) const;
  /** write optional descriptions */
  void WriteInfo(const ArrayT<StringT>& info_records) const;

  /** convert side set information to fe++ numbering convention */
  void ConvertSideSetOut(const char* elem_type, nArrayT<int>& sides) const;
  /** convert side set information from fe++ numbering convention */
  void ConvertSideSetIn(const char* elem_type, nArrayT<int>& sides) const;

  /* convert element numbering from/to fe++ numbering convention */
  void ConvertElementNumbering (iArray2DT& conn, int code) const;
			
  /** clear all parameter data */
  void Clear(void);
	
  /** process return values - (do_warning == 1) prints
   * warning message for (code > 0) */
  void Try(const char* caller, int code, bool do_warning) const;

private:

  /** message stream */
  ostream& fOut;

  /* opening parameters */
  StringT file_name; /**< file opened */
  int   exoid; /**< exodus file ID */
  int   comp_ws; /**< floating point size in program */
  int   io_ws;   /**< floating point size in file data */
  float version; /**< database library version */

  /* initialization parameters */
  int   num_dim; /**< spatial DOF */
  int   num_nodes; /**< number of nodes */
  int   num_elem; /**< total number of elements */
  int   num_elem_blk; /**< number of element blocks */
  int   num_node_sets; /**< number of node sets */
  int   num_side_sets; /**< number of side sets */
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

} // namespace Tahoe 
#endif /* _EXODUS_T_H_ */
