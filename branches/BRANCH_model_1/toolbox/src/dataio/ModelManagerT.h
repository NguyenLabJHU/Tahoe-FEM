/* $Id: ModelManagerT.h,v 1.4.2.11 2001-10-30 16:43:14 sawimme Exp $ */
/* created: sawimme July 2001 */

#ifndef _MODELMANAGER_T_H_
#define _MODELMANAGER_T_H_

/* direct */
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "StringT.h"
#include "GeometryT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "InputBaseT.h"

/* forward */
#include "ios_fwd_decl.h"
class ifstreamT;

/** Interface for any code to input data. Data types: Coordinates, Element 
 * Connectivity, NodeSets and SideSets are stored. Node and Element Maps 
 * along with variable data can be retrieved but is not stored in this class. 
 * All array data is referred to by its StringT ID or its int index with 
 * the data type. Special functions are added to read card data and arrays of
 * ID values from Tahoe parameter files. */
class ModelManagerT
{
 public:
  /** constructor.
   * \param message error messages ostream */
  ModelManagerT (ostream& message);
  ~ModelManagerT (void);

  /** initialize
   * Use only 1 of the 3 initialize functions. This function is meant to be
   * used in conjuction with a Tahoe parameter file. The readonly option controls
   * whether or not the InputBaseT* is initialized. If readonly is true, the pointer
   * is not cast and no data from the model file is registered. This option is
   * used when decomposing or joining parallel jobs.
   * \param in stream containing input file format and model file name 
   * \param readonly flag to only read data, if false the InputBaseT class is initialized */
  void Initialize (ifstreamT& in, bool readonly);
  /** initialize
   * Use only 1 of the 3 initialize functions. This function initializes the
   * InputBaseT* and register array data found in the model file 
   * \param format IO database format.
   * \param database Name of model file (null for inline text files */
  void Initialize (const IOBaseT::FileTypeT format, const StringT& database);
  /** initialize
   * Use only 1 of the 3 initialize functions. This function queriues the user
   * interactively for database format and model file name. It is meant to be used
   * with translator programs. It initializes the InputBaseT* and registers array 
   * data found in the model file. */
  void Initialize (void);

  /** echo database format and model file name to the ostream in the print format
   * used by Tahoe */
  void EchoData (ostream& o) const;
  /** access function
   * \return format database format
   * \return name model file name */
  void Format (IOBaseT::FileTypeT& format, StringT& name) const;

  /** InputBaseT node registration
   * called by InputBaseT* to register data
   * array data will be read later from model file, null array is allocated
   * \param length number of node points
   * \param dof spatial degrees of freedom  */
  bool RegisterNodes (int length, int dof);
  /** InputBaseT element registration
   * called by InputBaseT* to register data
   * array data will be read later from model file, null array is allocated 
   * \param name ID value 
   * \param numelems number of elements in connectivity set
   * \param numelemnodes number of nodes per element in set
   * \param code geometry code of elements in set */
  bool RegisterElementGroup (const StringT& name, int numelems, int numelemnodes, GeometryT::CodeT code);
  /** InputBaseT nodeset registration
   * called by InputBaseT* to register data
   * array data will be read later from model file, null array is allocated 
   * \param name ID value
   * \param length number of nodes in the set */
  bool RegisterNodeSet (const StringT& name, int length); 
  /** InputBaseT sideset registration
   * called by InputBaseT* to register data
   * array data will be read later from model file, null array is allocated 
   * \param name ID value
   * \param length number of facets in the set
   * \param local true if the side set elements are locally numbered within an element group
   * \param groupindex element group index which contains this side set, or -1 for unknown globally numbered sets */
  bool RegisterSideSet (const StringT& name, int length, bool local, int groupindex);

  /** external node registration 
   * register data not found in model file, data is copied into a storage array
   * \param coords array of coordinates */
  bool RegisterNodes (dArray2DT& coords);
  /** external element registration 
   * register data not found in model file, data is copied into a storage array 
   * \param name ID value 
   * \param conn array of connectivities
   * \param code geometry code of elements in set */
  bool RegisterElementGroup (const StringT& name, iArray2DT& conn, GeometryT::CodeT code);
  /** external node set registration 
   * register data not found in model file, data is copied into a storage array
   * \param name ID value
   * \param set array of nodes */
  bool RegisterNodeSet (const StringT& name, iArrayT& set);
  /** external side set registration 
   * register data not found in model file, data is copied into a storage array
   * \param name ID value
   * \param set array of facets
   * \param local true if the side set elements are locally numbered within an element group
   * \param groupindex element group index which contains this side set, or -1 for unknown globally numbered sets */
  bool RegisterSideSet (const StringT& name, iArray2DT& set, bool local, int groupindex);

  /** parameter file node registration
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name */
  bool RegisterNodes (ifstreamT& in);
  /** parameter file element registration
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param name ID value
   * \param code geometry code of elements in set */
  bool RegisterElementGroup (ifstreamT& in, const StringT& name, GeometryT::CodeT code);
  /** parameter file node set registration
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param name ID value */
  bool RegisterNodeSet (ifstreamT& in, const StringT& name);
  /** parameter file side set registration
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param name ID value
   * \param local true if the side set elements are locally numbered within an element group
   * \param groupindex element group index which contains this side set, or -1 for unknown globally numbered sets */
  bool RegisterSideSet (ifstreamT& in, const StringT& name, bool local, int groupindex);

  /* reads from input file the coordinate dimensions and coords if inline database format
   * \param in stream containing data or external file name */
  void ReadInlineCoordinates (ifstreamT& in);
  /** reads element block data from Tahoe parameter file 
   * Element block data: the number of blocks, list of names, and matnums
   * if the database format is inline, it will read connectivity and register it
   * \param in stream containing element block data or external file name
   * \return indexes array of element group indexes
   * \return matnums array of corresponding material IDs */
  void ElementBlockList (ifstreamT& in, iArrayT& indexes, iArrayT& matnums);
  /** reads node set block data from Tahoe paramter file
   * Node Set block data: number of sets and set IDs
   * if the datbase format is inline, it will read number of nodes and register them
   * \param in stream containing node set block data or external file name
   * \return indexes array of node set indexes */
  void NodeSetList (ifstreamT& in, iArrayT& indexes);
  /** reads side set block data from Tahoe paramter file
   * Side Set block data: number of sets and set IDs
   * if the datbase format is inline, it will read number of sides and register them
   * set multidatabasesets = false for places where the number of sets is not in the
   * parameter file and it is assumed that there is only one set to read
   * \param in stream containing side set block data or external file name
   * \return indexes array of node set indexes
   * \param multidatabasesets flag for slight parameter file inconsistency, see more info above */
  void SideSetList (ifstreamT& in, iArrayT& indexes, bool multidatabasesets);

  /** read IC/KBC/FBC card type data from Tahoe parameter file, including number of cards
   * \param in stream for parameter file
   * \param out error messaging stream
   * \return nodes array of node set arrays 
   * \return data array of integer card data, each row corresponds to a node set
   * \return value array of double values, each member corresponds to a node set */
  int ReadCards (ifstreamT& in, ostream& out, ArrayT<iArrayT>& nodes, iArray2DT& data, dArrayT& value);
  /** read traction card overall dimensions from Tahoe parameter file
   * Call this function in conjuction with ReadTractionSetData and ReadTractionSideSet
   * if inline database format, number tractions and number of sets are read
   * else number of sets are read
   * \param in stream for parameter file
   * \return numlines number of cards in parameter file
   * \return numsets number of sets */
  void ReadNumTractionLines (ifstreamT& in, int& numlines, int& numsets);
  /** read a set of data
   * if inline database format, the element block index and dimensions for a set of data is read
   * else the set of data is dimensioned as 1, element block index is set later
   * \param in stream for parameter file
   * \return blockindex element block index the set is contained within for inline text
   * \return setsize number of cards to read from parameter file for this set */
  void ReadTractionSetData (ifstreamT& in, int& blockindex, int& setsize);
  /** reads a set of traction cards
   * this read a partial traction card, the rest must be read by an element class, 
   * if inline database format, the element and facet is read
   * else the set name is read and the element block index is determined
   * \param in stream for parameter file
   * \return blockindex element block index the set is contained within for model file data
   * \return localsides array of facets from model file or just one facet from inline text, locally numbered */
  void ReadTractionSideSet (ifstreamT& in, int& blockindex, iArray2DT& localsides);

  /** access coordinate dimensions
   * \return length number of nodes
   * \return dof spatial degree of freedom */
  void CoordinateDimensions (int& length, int& dof) const;
  /** return a reference to the coordinate array, whether it is filled or empty */
  const dArray2DT& CoordinateReference (void) const;
  /** read the coordinate array if not yet read from the model file and returns a reference to the array */
  const dArray2DT& Coordinates (void);
  /** reads the coordinate array if not yet read from the model file, no return accessor */
  void ReadCoordinates (void);

  /** determine if coordinates are written wth 3 DOF for 2D elements
   * Patran, Abaqus, EnSight, etc. always store coordinates in 3D */
  bool AreElements2D (void) const;

  /** returns the number of element groups/blocks/sets */
  int NumElementGroups (void) const;
  /** returns an array of element groups names */
  void ElementGroupNames (ArrayT<StringT>& names) const;
  /** returns the index for the element group name */
  int ElementGroupIndex (const StringT& name) const;
  /** returns the dimensions for the element group */
  void ElementGroupDimensions (int index, int& numelems, int& numelemnodes) const;
  /** returns the geometry code for the element group */
  GeometryT::CodeT ElementGroupGeometry (int index) const;
  /** reads the elements if not yet read from the model file and returns a reference to the array
   * \note element node numbering is global, continuous, and offset to zero  */
  const iArray2DT& ElementGroup (int index);
  /** reads the elements if not yet read from the model, no return accessor */
  void ReadConnectivity (int index);
  /** returns the pointer to the element group array, whether it is filled or empty
   * \note element node numbering is global, continuous, and offset to zero */
  const iArray2DT* ElementGroupPointer (int index) const;

  /** access node map data
   * \note data is not offset and may not be continuous */
  void AllNodeMap (iArrayT& map);
  /** access element map data
   * \note data is not offset and may not be continuous */
  void AllElementMap (iArrayT& map);
  /** access element map data for a given element group name
   * \note data is not offset and may not be continuous */
  void ElementMap (StringT& name, iArrayT& map);

  /** return number of node sets */
  int NumNodeSets (void) const;
  /** return array of node set IDs */
  void NodeSetNames (ArrayT<StringT>& names) const;
  /** return index for the node set name */
  int NodeSetIndex (const StringT& name) const;
  /** return node set length */
  int NodeSetLength (int index) const;
  /** return reference to node set array
      \note node numbering is global, continuous, and offset to zero */
  const iArrayT& NodeSet (int index);
  /** compile the set of node sets indicated by indexes into one sorted array called nodes
   * \note node numbering is global, continuous, and offset to zero */
  void ManyNodeSets (const iArrayT& indexes, iArrayT& nodes);

  /** returns the number of side sets */
  int NumSideSets (void) const;
  /** return the side set IDs */
  void SideSetNames (ArrayT<StringT>& names) const;
  /** returns index for the side set name */
  int SideSetIndex (const StringT& name) const;
  /** returns side set length */
  int SideSetLength (int index) const;
  /** returns reference to side set array
   *\note elements in the array may be numbered locally or globally, but are offset to zero and continuous
   *\note facets numbers are offset to zero */
  const iArray2DT& SideSet (int index) const;
  /** determines if model file storage is local or global */
  bool IsSideSetLocal (int index) const;
  /** determines element group that contains the side set, -1 is returned for globally numbered sets */
  int SideSetGroupIndex (int sidesetindex) const;
  /** convert locally numbered set to globally numbered */
  void SideSetLocalToGlobal (const int localelemindex, const iArray2DT& local, iArray2DT& global);
  /** convert globally numbered set to locally numbered and determine element group containing set */
  void SideSetGlobalToLocal (int& localelemindex, iArray2DT& local, const iArray2DT& global);

  /** add nodes to the coordinate array
   * \param newcoords array of coordinates to add
   * \return new_node_tags node tags, globally numbered, continuous, offset to zero
   * \return newtotalnumnodes number of nodes after adding */
  void AddNodes (const dArray2DT& newcoords, iArrayT& new_node_tags, int& newtotalnumnodes);
  /** dupicate nodes to expand the coordinate array
   * \param nodes array of node tags that will be duplicated
   * \return new_node_tags node tags, globally numbered, continuous, offset to zero
   * \return newtotalnumnodes number of nodes after adding */
  void DuplicateNodes (const iArrayT& nodes, iArrayT& new_node_tags, int& newtotalnumnodes);
  /** adjust the DOF of the coordinate array from 3D to 2D by dropping the 3rd coordiante value */
  void AdjustCoordinatesto2D (void);

  /** call this function if the connectivity group/block/set is altered through the element group pointer
   * the number of elements and element nodes is updated */
  void UpdateConnectivityDimensions (int index);
  /** add elements to an element group array
   * \param index element group index
   * \param connects connectivity of elements to add
   * \return new_elem_tags element tags, locally numbered, continuous, offset to zero
   * \return newtotalnumelems number of elements in the group after adding */
  void AddElement (int index, const iArray2DT& connects, iArrayT& new_elem_tags, int& newtotalnumelems);

  /** This closes the link to InputBaseT, it does not clear any stored data */
  void CloseModel (void);

  /** checks to see if external file or actual data */
  ifstreamT& OpenExternal (ifstreamT& in, ifstreamT& in2, ostream& out, bool verbose, const char* fail) const;

  /** returns the number of time steps */
  int NumTimeSteps (void);
  /** return time values */
  void TimeSteps (dArrayT& steps);
  
  /** return number of node variables found */
  int NumNodeVariables (void);
  /** returns node variable labels */
  void NodeLabels (ArrayT<StringT>& labels);
  /** return node variable values for the stepindex for all node points */
  void AllNodeVariables (int stepindex, dArray2DT& values);
  /** return node variable values for the stepindex for all nodes in the element set */
  void NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values);
  /** return node variable values for the setpindex for all nodes in the node set */
  void NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values);

  /** return number of element variables found */
  int NumElementVariables (void);
  /** returns element variable labels */
  void ElementLabels (ArrayT<StringT>& labels);
  /** returns element variable values for the stepindex for all elements */
  void AllElementVariables (int stepindex, dArray2DT& values);
  /** returns element variable values for the stepindex for all elements in the element set */
  void ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values);

  /** return the number of quadrature points per element for a specified element set */
  int NumElementQuadPoints (StringT& elsetname);
  /** return number of quadrature variables found */
  int NumQuadratureVariables (void);
  /** returns quadrature variable labels */
  void QuadratureLabels (ArrayT<StringT>& labels);
  /** returns quadrature variable values for the stepindex for all elements */
  void AllQuadratureVariables (int stepindex, dArray2DT& values);
  /** returns quadrature variable values for the stepindex for all elements in the element set */
  void QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values); 

 private:
  /** sets the InputBaseT pointer, scans the model file, registers array data found */
  void ScanModel (const StringT& database);
  /** scans the model file for element groups and registers them */
  bool ScanElements (void);
  /** scans the model file for node sets and registers them */
  bool ScanNodeSets (void);
  /** scans the model file for side sets and registers them */
  bool ScanSideSets (void);
  
  /** checks the name of the set being registered against names already in the registry */
  bool CheckName (const ArrayT<StringT>& list, const StringT& name, const char *type) const;

 protected:
  ostream& fMessage; /**< where to write error messages */
  IOBaseT::FileTypeT fFormat; /**< database format */
  InputBaseT *fInput; /**< database class */
  StringT fInputName; /**< model file name or basis */

 private:
  /* dimensional information */
  iArrayT fCoordinateDimensions; /**< num nodes and dof */
  iAutoArrayT fElementLengths; /**< number of elements */
  iAutoArrayT fElementNodes; /**< number of element nodes */
  iAutoArrayT fNodeSetDimensions; /**< number of nodes in set */
  iAutoArrayT fSideSetDimensions; /**< number of sides in set */
  AutoArrayT<bool> fSideSetIsLocal; /**< flag for globally or locally numbered */
  iAutoArrayT fSideSetGroupIndex; /**< -1 for globally numbered or element group that contains set */

  /* set parameters */
  AutoArrayT<StringT> fElementNames; /**< element group IDs */
  AutoArrayT<StringT> fNodeSetNames; /**< node set IDs */
  AutoArrayT<StringT> fSideSetNames; /**< side set IDs */
  AutoArrayT<GeometryT::CodeT> fElementCodes; /** element group geometry codes */

  /* data */
  int fNumElementSets; /**< number element groups registered */
  int fNumNodeSets; /**< number of node sets registered */
  int fNumSideSets; /**< number of side sets registered */
  dArray2DT fCoordinates; /**< coordinates */
  AutoArrayT<iArray2DT> fElementSets; /**< connectivities */ 
  AutoArrayT<iArrayT> fNodeSets; /**< node sets */
  AutoArrayT<iArray2DT> fSideSets; /**< side sets */
};

inline void ModelManagerT::Format (IOBaseT::FileTypeT& format, StringT& name) const
{ 
  format = fFormat;
  name = fInputName;
}

inline int ModelManagerT::NumElementGroups (void) const { return fNumElementSets; }
inline int ModelManagerT::NumNodeSets (void) const { return fNumNodeSets; }
inline int ModelManagerT::NumSideSets (void) const { return fNumSideSets; }

inline int ModelManagerT::NumTimeSteps (void) { return fInput->NumTimeSteps(); }
inline void ModelManagerT::TimeSteps (dArrayT& steps) { fInput->ReadTimeSteps (steps); }

inline int ModelManagerT::NumNodeVariables (void) { return fInput->NumNodeVariables (); }
inline void ModelManagerT::NodeLabels (ArrayT<StringT>& labels) { fInput->ReadNodeLabels (labels); }
inline void ModelManagerT::AllNodeVariables (int stepindex, dArray2DT& values) { fInput->ReadAllNodeVariables (stepindex, values); }
inline void ModelManagerT::NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadNodeVariables (stepindex, elsetname, values); }
inline void ModelManagerT::NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values) { fInput->ReadNodeSetVariables (stepindex, nsetname, values); }

inline int ModelManagerT::NumElementVariables (void) { return fInput->NumElementVariables (); }
inline void ModelManagerT::ElementLabels (ArrayT<StringT>& labels) { fInput->ReadElementLabels (labels); }
inline void ModelManagerT::AllElementVariables (int stepindex, dArray2DT& values) { fInput->ReadAllElementVariables (stepindex, values); }
inline void ModelManagerT::ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadElementVariables (stepindex, elsetname, values); }

inline int ModelManagerT::NumElementQuadPoints (StringT& name) { return fInput->NumElementQuadPoints(name); }
inline int ModelManagerT::NumQuadratureVariables (void) { return fInput->NumQuadratureVariables (); }
inline void ModelManagerT::QuadratureLabels (ArrayT<StringT>& labels) { fInput->ReadQuadratureLabels (labels); }
inline void ModelManagerT::AllQuadratureVariables (int stepindex, dArray2DT& values) { fInput->ReadAllQuadratureVariables (stepindex, values); }
inline void ModelManagerT::QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadQuadratureVariables (stepindex, elsetname, values); }

#endif
