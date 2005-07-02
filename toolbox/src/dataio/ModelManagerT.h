/* $Id: ModelManagerT.h,v 1.37 2005-07-02 22:44:41 paklein Exp $ */
/* created: sawimme July 2001 */
#ifndef _MODELMANAGER_T_H_
#define _MODELMANAGER_T_H_

#include "ios_fwd_decl.h"

/* direct members */
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "StringT.h"
#include "GeometryT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "InputBaseT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
template <class TYPE> class nVariArray2DT;
class InverseMapT;

/** Interface for any code to input data. Data types: Coordinates, Element 
 * Connectivity, NodeSets and SideSets are stored. Node and Element Maps 
 * along with variable data can be retrieved but is not stored in this class. 
 * All array data is referred to by its StringT ID or its int index with 
 * the data type. Special functions are added to read card data and arrays of
 * ID values from Tahoe parameter files. */
class ModelManagerT
{
 public:
 
 	/** class constants */
 	enum ConstantsT {
 		kNotFound =-1 /**< returned by element/node/set queries if ID not found */
 		};

	/** side set scopes */
	enum SideSetScopeT {
		 kLocal = 0, /**< elements numbered within the element block */
		kGlobal = 1  /**< elements numbered with respect to all elements */
		};

  /** constructor.
   * \param message error messages ostream */
  ModelManagerT (ostream& message);
  
  /** destructor */
  ~ModelManagerT (void);
  
  /** \name initialization */
  /** This function initializes the
   * InputBaseT* and register array data found in the model file 
   * \param format IO database format.
   * \param database Name of model file (null for inline text files)
   * \param scan_model if true model dimensions are scanned
   * \return true if model database is open, false otherwise */
  bool Initialize(const IOBaseT::FileTypeT format, const StringT& database, bool scan_model);

  /** echo database format and model file name to the ostream in the print format
   * used by Tahoe */
  void EchoData (ostream& o) const;

	/** This closes the link to InputBaseT, it does not clear any stored data */
	void CloseModel(void);

	/** checks to see if external file or actual data */
	ifstreamT& OpenExternal (ifstreamT& in, ifstreamT& in2, ostream& out, bool verbose, const char* fail) const;
  
	/** \name basic accessors */
	/*@{*/
	/** returns QA records */
	void QARecords (ArrayT<StringT>& records);

	/** return the database file format */
	IOBaseT::FileTypeT DatabaseFormat(void) const { return fFormat; };

	/** return the database file name */
	const StringT& DatabaseName(void) const { return fInputName; };

	/** returns the number of output time steps */
	int NumTimeSteps (void);

	/** return times for the output data */
	void TimeSteps (dArrayT& steps);
	/*@}*/
  
	/** \name nodal coordinate information */
	/*@{*/
	/** number of nodes */
	int NumNodes(void) const;

	/** number of spatial dimensions */
	int NumDimensions (void) const;

	/** reads from input file the coordinate dimensions and coords if inline database format
	 * \param in stream containing data or external file name */
	void ReadInlineCoordinates (ifstreamT& in);

	/** return the coordinates.
	 * reads the coordinate array if not yet read from the model file and returns a 
	 * reference to the array */
	const dArray2DT& Coordinates(void);
  
	/** reads the coordinate array if not yet read from the model file, no return accessor */
	void ReadCoordinates (void);

	/** access node map data
	 * \param ids not offset and may not be continuous */
	void AllNodeIDs(iArrayT& ids);

	/* access node data 
	 * \param map offset, globally numbered, and continuous */
	void AllNodeMap (iArrayT& map);
	/*@}*/

	/** \name element information */
	/*@{*/
	/** returns the number of element groups */
	int NumElementGroups(void) const;

	/** returns the total number of elements in groups */
	int NumElements(void) const;

	/** determine if coordinates are written wth 3 DOF for 2D elements;
	 * Patran, Abaqus, EnSight, etc. always store coordinates in 3D */
	bool AreElements2D (void) const;
	
	/** reads element block data from Tahoe parameter file.
	 * Element block data: the number of blocks, list of names, and matnums
	 * if the database format is inline, it will read connectivity and register it
	 * \param in stream containing element block data or external file name
	 * \param ID returned array of element group ID's
	 * \param matnums returned array of corresponding material IDs */
	void ElementBlockList (ifstreamT& in, ArrayT<StringT>& ID, iArrayT& matnums);

	/** return the ID of the element group at the given index. The name of the
	 * element group is the string form of the database-specific element block
	 * identifier. */
	const StringT& ElementGroupID(int index) const;

	/** returns an array of element groups names. The names of the
	 * element group are the string form of the database-specific element block
	 * identifiers. */
	 const ArrayT<StringT>& ElementGroupIDs(void) const { return fElementNames; };

	/** returns the index for the element group name.
	 * \return the element group index or ModelManagerT::kNotFound if ID is not found. */
	int ElementGroupIndex (const StringT& ID) const;

	/** return an unused element ID 
	 * \param prefix prefix for ID. The prefix can be empty. */
	StringT FreeElementID(const StringT& prefix) const;

	/** returns the dimensions for the element group */
	void ElementGroupDimensions (const StringT& ID, int& numelems, int& numelemnodes) const;

	/** returns the geometry code for the element group */
	GeometryT::CodeT ElementGroupGeometry (const StringT& ID) const;

	/** reads the elements if not yet read from the model file and returns a 
	 * reference to the array
	 * \note element node numbering is global, continuous, and offset to zero  */
	const iArray2DT& ElementGroup (const StringT& ID);

	/** reads the elements if not yet read from the model, no return accessor */
	void ReadConnectivity (const StringT& ID);

	/** returns the pointer to the element group array, whether it is filled or empty
	 * \note element node numbering is global, continuous, and offset to zero */
	const iArray2DT* ElementGroupPointer (const StringT& ID) const;

	/** collect pointers to the element groups with the ID's given
	 * \param IDs list of element block ID's
	 * \param blocks returns with pointers to the element blocks with the
	 *        given ID's. Array is re-dimensioned to the length of IDs */
	void ElementGroupPointers(const ArrayT<StringT>& IDs, 
		ArrayT<const iArray2DT*>& blocks) const;

	/** access element map data
	 * \param ids not offset and may not be continuous */
	void AllElementIDs(iArrayT& ids);

	/** access element map data for a given element group name. The names of the
	 * element group are the string form of the database-specific element block
	 * identifiers.
	 * \param ids not offset and may not be continuous */
	void ElementIDs (const StringT& ID, iArrayT& ids);

	/** access element set data for a given element group name. The names of the
	 * element group are the strig form of the database-specific element block
	 * identifiers.
	 * \param map offset, globally numbered, and continuous */
	void ElementMap (const StringT& ID, iArrayT& map);

	/** return the "bounding" elements. For the body comprised of the given
	 * list of element blocks, determine the elements with faces on the body
	 * boundary and the corresponding neighbors. All returned arrays are
	 * dimensioned during the call. Function is non-const because element
	 * blocks that have not been read yet will be read.
	 * \param IDs list of element blocks comprising the body
	 * \param elements elements within each block that have faces on the body 
	 *        boundary. Element are numbered sequentially through all blocks
	 * \param neighbors neighoring elements for bounding element 
	 * \param geometry use the supplied GeometryBaseT is non-NULL; otherwise
	 *        construct a temporary */
	void BoundingElements(const ArrayT<StringT>& IDs, iArrayT& elements, 
		iArray2DT& neighbors, const GeometryBaseT* geometry = NULL);

	/** return the neighoring elements for each element. For the body comprised of the given
	 * list of element blocks, determine the neighboring elements for each element . 
	 * The returned array is dimensioned during the call. Function is non-const because element
	 * blocks that have not been read yet will be read.
	 * \param IDs list of element blocks comprising the body
	 * \param neighbors neighoring elements for each element in the blocks passed in. The major
	 *        dimension of this array is the total number of elements in the given list of blocks.
	 *        The minor dimension is the number of neighboring faces for the given element geometry.
	 *        The element associated with each row in the array and the neighbor numbers in each row
	 *        are elements numbered with respect to their position in the collection of elements defined
	 *        by the list of ID's passed in. A neighbor of -1 indicates the element has no neighbor
	 *        on the given face. The faces are number based on the cannonical face number defined for
	 *        the given element geometry.
	 * \param geometry use the supplied GeometryBaseT is non-NULL; otherwise
	 *        construct a temporary */
	void ElementNeighbors(const ArrayT<StringT>& IDs, iArray2DT& neighbors, 
		const GeometryBaseT* geometry = NULL);

	/** return the list of element block ID's containing the given list of nodes 
	 * \param nodes list of node numbers
	 * \param returns with a list of element block ID's that contain the give nodes */
	void ElementGroupIDsWithNodes(const ArrayT<int>& nodes, ArrayT<StringT>& element_ids);

	/** element faces on the group "surface" grouped into contiguous patches */
	void SurfaceFacets(const ArrayT<StringT>& IDs,
		GeometryT::CodeT& geometry_code,
		ArrayT<iArray2DT>& surface_facet_sets,
		iArrayT& surface_nodes,
		const GeometryBaseT* geometry = NULL);
		
	/** finds element facets, facet nodes, and facet numbers */
	void SurfaceFacets(const ArrayT<StringT>& IDs,
		GeometryT::CodeT& geometry_code,
		iArray2DT& surface_facets,
		iArrayT& surface_nodes,
		iArrayT& facet_numbers,
		iArrayT& elem_numbers,
		const GeometryBaseT* geometry = NULL);
	
	/** generate a list of nodes on the "surface" of the element group
	 * based in the group connectivities */
	void SurfaceNodes(const ArrayT<StringT>& IDs, 
		iArrayT& surface_nodes,
		const GeometryBaseT* geometry = NULL);
	/*@}*/

	 /** compute the nodal area associated with each striker node */
	 void ComputeNodalArea(const iArrayT& node_tags, 
		const ArrayT<StringT>& surface_blocks, 
		dArrayT& nodal_area, InverseMapT& inverse_map);


	/** \name node set information */
	/*@{*/
	/** return number of node sets */
	int NumNodeSets (void) const;

	/** return an unused node set ID 
	 * \param prefix prefix for ID. The prefix can be empty. */
	StringT FreeNodeSetID(const StringT& prefix) const;

	/** reads node set block data from Tahoe paramter file.
	 * Node Set block data: number of sets and set IDs
	 * if the datbase format is inline, it will read number of nodes and register them
	 * \param in stream containing node set block data or external file name
	 * \param ID returned array of node set ID's */
	void NodeSetList (ifstreamT& in, ArrayT<StringT>& ID);

	/** return array of node set IDs. The names of the
	 * node set are the string form of the database-specific element block
	 * identifiers. */
	const ArrayT<StringT>& NodeSetIDs(void) const { return fNodeSetNames; };

	/** return index for the node set name
	 * \return the node set index or ModelManagerT::kNotFound if ID is not found. */  
	int NodeSetIndex (const StringT& ID) const;

	/** return node set length */
	int NodeSetLength (const StringT& ID) const;

	/** return reference to node set array
	 * \note node numbering is global, continuous, and offset to zero */
	const iArrayT& NodeSet (const StringT& ID);

	/** return mapped node set array.
	 * compile the set of node sets indicated by indexes into one sorted array called nodes
	 * \note node numbering is global, continuous, and offset to zero */
	void ManyNodeSets (const ArrayT<StringT>& ID, iArrayT& nodes);

	/** return mapped node set array.
	 * compile the set of node sets indicated by indexes into one sorted array called nodes
	 * \note node numbering is global, continuous, and offset to zero */
	void ManyNodeSets (const ArrayT<StringT>& ID, AutoArrayT<int>& nodes);
	/*@}*/

	/** \name side set information */
	/*@{*/
	/** returns the number of side sets */
	int NumSideSets (void) const;

	/** reads side set block data from Tahoe paramter file.
	 * Side Set block data: number of sets and set IDs
	 * if the datbase format is inline, it will read number of sides and register them
	 * set multidatabasesets = false for places where the number of sets is not in the
	 * parameter file and it is assumed that there is only one set to read
	 * \param in stream containing side set block data or external file name
	 * \param ID returned array of side set ID's
	 * \param multidatabasesets flag for slight parameter file inconsistency, see more info above */
	void SideSetList (ifstreamT& in, ArrayT<StringT>& ID, bool multidatabasesets);

	/** return the side set names. The names of the
	 * side sets are the string form of the database-specific element block
	 * identifiers. */
	const ArrayT<StringT>& SideSetIDs(void) const { return fSideSetNames; };

	/** returns index for the side set name. The names of the
	 * side sets are the string form of the database-specific element block
	 * identifiers.
	 * \return the side set index or ModelManagerT::kNotFound if ID is not found. */  
	int SideSetIndex (const StringT& ID) const;

	/** returns side set length */
	int SideSetLength (const StringT& ID) const;

	/** returns reference to side set array. Element numbering in the array may be numbered locally 
	 * or globally, but are offset to zero and continuous facets numbers are offset to zero. The
	 * scope of the element numbering can be determined from ModelManagerT::IsSideSetLocal, and
	 * can be transformed using ModelManagerT::SideSetLocalToGlobal and ModelManagerT::SideSetGlobalToLocal. */
	const iArray2DT& SideSet (const StringT& ID);

	/** return side set as nodes on faces. Limited to element geometries for
	 * which all facets have the same number of nodes.
	 * \param ID side set ID
	 * \param facet_geom geometry of each face
	 * \param facet_nodes number of nodes in each face
	 * \param faces returns with the global node numbers on each face: [nface] x [nfn] */
	void SideSet(const StringT& ID, ArrayT<GeometryT::CodeT>& facet_geom,
		iArrayT& facet_nodes, iArray2DT& faces);

	/** determines if model file storage is local or global */
	bool IsSideSetLocal (const StringT& ID) const;

	/** determines element group that contains the side set, -1 is returned for 
	 * globally numbered sets */
	const StringT& SideSetGroupID (const StringT& ss_ID) const;

	/** convert locally numbered set to globally numbered
	 * \param element_ID ID of the associated element group */
	void SideSetLocalToGlobal (const StringT& element_ID, const iArray2DT& local, iArray2DT& global);

	/** convert globally numbered set to locally numbered and determine element group containing set */
	void SideSetGlobalToLocal(const iArray2DT& global, iArray2DT& local, StringT& element_ID);
	/*@}*/

	/** \name model transformations */
	/*@{*/
	/** add nodes to the coordinate array
	 * \param newcoords array of coordinates to add
	 * \param new_node_tags returned node tags, globally numbered, continuous, offset to zero
	 * \param newtotalnumnodes returned number of nodes after adding */
	void AddNodes (const dArray2DT& newcoords, iArrayT& new_node_tags, int& newtotalnumnodes);

	/** dupicate nodes to expand the coordinate array
	 * \param nodes array of node tags that will be duplicated
	 * \param new_node_tags returned node tags, globally numbered, continuous, offset to zero
	 * \param newtotalnumnodes returned number of nodes after adding */
	void DuplicateNodes (const iArrayT& nodes, iArrayT& new_node_tags, int& newtotalnumnodes);

	/** resize the coordinate array. Excess values at the tail of the coordinate list are
	 * discarded. Additional space added to the array is not initialized */
	void ResizeNodes(int num_nodes);

	/** adjust the DOF of the coordinate array from 3D to 2D by dropping the 3rd coordiante value */
	void AdjustCoordinatesto2D (void);

	/** replace the coordinate list */
	void UpdateNodes(dArray2DT& coordinates, bool keep) { RegisterNodes(coordinates, keep); };

	/** overwrite the given the coordinates */
	void UpdateNodes(const dArray2DT& coordinates, const ArrayT<int>& nodes);

	/** overwrite the coordinates of the given node */
	void UpdateNode(const dArrayT& coordinates, int node);

	/** call this function if the connectivity group/block/set is altered and replacement is needed
	 * the number of elements and element nodes is updated
	 * \param conn updated connectivities
	 * \param keep flag to control memory allocation. If keep is true and conn
	 *        is not an alias, the model manager will take ownership of the
	 *        memory in conn. conn will be passed back as an alias. Otherwise
	 *        the data in conn will be copied. */
	void UpdateElementGroup(const StringT& ID, iArray2DT& conn, bool keep);

	/** change the number of elements in the element group */
	void ResizeElementGroup(const StringT& ID, int num_elements);

	/** update the nodes in an existing node set */
	void UpdateNodeSet(const StringT& ID, iArrayT& node_set, bool keep);

	/** update the nodes in an existing side set */
	void UpdateSideSet(const StringT& ID, iArray2DT& side_set, bool keep);

	/** add elements to an element group array
	 * \param index element group index
	 * \param connects connectivity of elements to add
	 * \param new_elem_tags returned element tags, locally numbered, continuous, offset to zero
	 * \param newtotalnumelems returned number of elements in the group after adding */
	void AddElement (const StringT& ID, const iArray2DT& connects, iArrayT& new_elem_tags, int& newtotalnumelems);
	/*@}*/

	/** \name registration functions
	 * Registering model elements not read by the model manager. */
	 /*@{*/
  /** external node registration.
   * register data not found in model file, data is copied into a storage array
   * \param coords array of coordinates 
   * \param keep flag to control memory allocation. If keep is true and coords
   *        is not an alias, the model manager will take ownership of the
   *        memory in coords. coords will be passed back as an alias. Otherwise
   *        the data in coords will be copied. */
  bool RegisterNodes(dArray2DT& coords, bool keep);

  /** parameter file node registration.
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name */
  bool RegisterNodes (ifstreamT& in);

  /** external element registration.
   * register data not found in model file, data is copied into a storage array 
   * \param ID ID value. Names should be the string form of the database-specific
   *        element block identifiers.
   * \param conn array of connectivities
   * \param code geometry code of elements in set
   * \param keep flag to control memory allocation. If keep is true and conn
   *        is not an alias, the model manager will take ownership of the
   *        memory in conn. conn will be passed back as an alias. Otherwise
   *        the data in conn will be copied.
   *
   * \note RegisterElementGroup puts an element block in the geometry data
   *       that was not present in the original geometry file. Currently, this
   *       means values at the nodes may be written but output of element values is 
   *       not supported.
   */
  bool RegisterElementGroup (const StringT& ID, iArray2DT& conn, GeometryT::CodeT code, 
  	bool keep);

  /** parameter file element registration.
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param ID ID value. Names should be the string form of the database-specific
   *        element block identifiers.
   * \param code geometry code of elements in set */
  bool RegisterElementGroup (ifstreamT& in, const StringT& ID, GeometryT::CodeT code);

  /** create an element group. Create an element group with the given ID, geometry code,
   * and number of element nodes. Return true if the group was registered, or false if
   * registration fails for any reason. The group is initialized with zero elements. The
   * group can be dimensioned with a call to ResizeElementGroup. */
  bool RegisterElementGroup(const StringT& ID, GeometryT::CodeT code, int numelemnodes);

  /** external node set registration. 
   * register data not found in model file, data is copied into a storage array
   * \param ID ID value. Names should be the string form of the database-specific
   *        node set identifiers.
   * \param set array of nodes
   * \param keep flag to control memory allocation. If keep is true and set
   *        is not an alias, the model manager will take ownership of the
   *        memory in set. set will be passed back as an alias. Otherwise
   *        the data in set will be copied. */
  bool RegisterNodeSet (const StringT& ID, iArrayT& set, bool keep);

  /** parameter file node set registration.
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param ID ID value. Names should be the string form of the database-specific
   *        node set identifiers. */
  bool RegisterNodeSet (ifstreamT& in, const StringT& ID);
  
  /** external side set registration.
   * register data not found in model file, data is copied into a storage array
   * \param ss_ID ID value. Names should be the string form of the database-specific
   *        side set identifiers.
   * \param set array of facets
   * \param local true if the side set elements are locally numbered within an element group
   * \param element_ID element group ID which contains this side set, or -1 for unknown globally numbered sets
   * \param keep flag to control memory allocation. If keep is true and set
   *        is not an alias, the model manager will take ownership of the
   *        memory in set. set will be passed back as an alias. Otherwise
   *        the data in set will be copied. */
  bool RegisterSideSet (const StringT& ss_ID, iArray2DT& set, SideSetScopeT scope, 
  	const StringT& element_ID, bool keep);

  /** parameter file side set registration.
   * called by functions that are reading the Tahoe parameter file, if the database format
   * is inline or external file data. Data is read into a storage array.
   * \param in stream containing data or external file name
   * \param ss_ID ID value. Names should be the string form of the database-specific
   *        side set identifiers.
   * \param local true if the side set elements are locally numbered within an element group
   * \param element_ID element group ID which contains this side set, or -1 for unknown globally numbered sets */
  bool RegisterSideSet (ifstreamT& in, const StringT& ss_ID, SideSetScopeT local, 
  	const StringT& element_ID);
	 /*@}*/

	/** \name Tahoe-specific boundary condition specifications */
	/*@{*/
  /** read IC/KBC/FBC card type data from Tahoe parameter file, including number of cards.
   * \param in stream for parameter file
   * \param out error messaging stream
   * \param nodes returned array of node set arrays 
   * \param return data returned array of integer card data, each row corresponds to a node set
   * \param return value returned array of double values, each member corresponds to a node set */
  int ReadCards (ifstreamT& in, ostream& out, ArrayT<iArrayT>& nodes, iArray2DT& data, dArrayT& value);

  /** read traction card overall dimensions from Tahoe parameter file.
   * Call this function in conjuction with ReadTractionSetData and ReadTractionSideSet
   * if inline database format, number tractions and number of sets are read
   * else number of sets are read
   * \param in stream for parameter file
   * \param numlines returned number of cards in parameter file
   * \param numsets returned number of sets */
  void ReadNumTractionLines (ifstreamT& in, int& numlines, int& numsets);

  /** read a set of data.
   * if inline database format, the element block ID and dimensions for a set of data is read
   * else the set of data is dimensioned as 1, element block ID is set later
   * \param in stream for parameter file
   * \param element_ID returned element block ID the set is contained within for inline text
   * \param setsize returned number of cards to read from parameter file for this set */
  void ReadTractionSetData (ifstreamT& in, StringT& element_ID, int& setsize);

  /** reads a set of traction cards.
   * this read a partial traction card, the rest must be read by an element class, 
   * if inline database format, the element and facet is read
   * else the set ID is read and the element block ID is determined
   * \param in stream for parameter file
   * \param element_ID returned element block ID the set is contained within for model file data
   * \param localsides returned array of facets from model file or just one facet from inline text, locally numbered */
  void ReadTractionSideSet (ifstreamT& in, StringT& element_ID, iArray2DT& localsides);
	/*@}*/
  
	/** \name nodal results data */
	/*@{*/
  /** return number of node variables found */
  int NumNodeVariables (void);

  /** returns node variable labels.
   * \param labels returns with the labels for the nodal output data. Arrays is
   *        (re-)dimensioned by the call. */
  void NodeLabels (ArrayT<StringT>& labels);

  /** returns an array of 0 or >0 values to indicate if the element set had data for 
   * each node variable */
  void NodeVariablesUsed (const StringT& ID, iArrayT& used);

  /** returns one variable for the stepindex for all nodes */
  void AllNodeVariable (int step, int varindex, dArrayT& values);

  /** returns one variable for nodes in an element set */
  void NodeVariable (int step, const StringT& ID, int varindex, dArrayT& values);

  /** return node variable values for the stepindex for all node points */
  void AllNodeVariables (int stepindex, dArray2DT& values);

  /** return node variable values for the stepindex for all nodes in the element set. 
   * \param elsetname name of the element group. The names are the string form of the
   *        database-specific element block identifiers. */
  void NodeVariables (int stepindex, const StringT& ID, dArray2DT& values);

  /** return node variable values for the setpindex for all nodes in the node set.
   * \param ID ID of the node set. The names are the string form of the
   *        database-specific element block identifiers. */
  void NodeSetVariables (int stepindex, const StringT& ID, dArray2DT& values);
	/*@}*/

	/** \name element results data */
	/*@{*/
  /** return number of element variables found */
  int NumElementVariables (void);

  /** returns element variable labels 
   * \param labels returns with the labels for the element output variables. Array is
   *        (re-)dimensioned by the call. */
  void ElementLabels (ArrayT<StringT>& labels);

  /** returns an array of 0 or >0 values to indicate if the element set had data for each element variable */
  void ElementVariablesUsed (const StringT& ID, iArrayT& used);

  /** returns one variable for the stepindex for all nodes */
  void AllElementVariable (int step, int varindex, dArrayT& values);

  /** returns one variable for nodes in an element set */
  void ElementVariable (int step, const StringT& ID, int varindex, dArrayT& values);

  /** returns element variable values for the stepindex for all elements */
  void AllElementVariables (int stepindex, dArray2DT& values);

  /** returns element variable values for the stepindex for all elements in the element set */
  void ElementVariables (int stepindex, const StringT& ID, dArray2DT& values);

  /** return the number of quadrature points per element for a specified element set */
  int NumElementQuadPoints (const StringT& ID);

  /** return number of quadrature variables found */
  int NumQuadratureVariables (void);

  /** returns quadrature variable labels */
  void QuadratureLabels (ArrayT<StringT>& labels);

  /** returns an array of 0 or >0 values to indicate if the element set had data for each quadrature variable */
  void QuadratureVariablesUsed (const StringT& ID, iArrayT& used);

  /** returns one variable for the stepindex for all nodes */
  void AllQuadratureVariable (int step, int varindex, dArrayT& values);

  /** returns one variable for nodes in an element set */
  void QuadratureVariable (int step, const StringT& ID, int varindex, dArrayT& values);

  /** returns quadrature variable values for the stepindex for all elements */
  void AllQuadratureVariables (int stepindex, dArray2DT& values);

  /** returns quadrature variable values for the stepindex for all elements in the element set */
  void QuadratureVariables (int stepindex, const StringT& ID, dArray2DT& values); 
	/*@}*/

 private:
 
 	/** return a reference to the input class. Throws eDataBaseFail if
 	 * the ModelMahagerT::fInput is NULL */
 	InputBaseT& Input(const char* caller = NULL) const;
 
 	/** return true of the ID's match */
 	bool ID_Match(const StringT& a, const StringT& b) const;
 
	/** return the length of the ID string not including any trailing white-space padding */
 	int ID_Length(const StringT& ID) const;

  /** sets the InputBaseT pointer, scans the model file, registers array data found
   * \return true if successful, false otherwise */
  bool ScanModel (const StringT& database);

  /** scans the model file for element groups and registers them */
  bool ScanElements (void);
  /** scans the model file for node sets and registers them */
  bool ScanNodeSets (void);
  /** scans the model file for side sets and registers them */
  bool ScanSideSets (void);
  
  /** checks the name of the set being registered against names already in the registry */
  bool CheckID (const ArrayT<StringT>& list, const StringT& ID, const char *type) const;

	/** clear database parameters */
	void Clear(void);

	/** return a pointer to a new GeometryBaseT. User is responsible for deleting class. */
	GeometryBaseT* NewGeometry(GeometryT::CodeT geometry, int nen) const;

 protected:
 
  ostream& fMessage; /**< where to write error messages */
  IOBaseT::FileTypeT fFormat; /**< database format */
  InputBaseT *fInput; /**< database class */
  StringT fInputName; /**< model file name or basis */

 private:
 
  /** \name dimensions */
  /*@{*/
  iArrayT fCoordinateDimensions; /**< num nodes and dof */
  iAutoArrayT fElementLengths; /**< number of elements */
  iAutoArrayT fElementNodes; /**< number of element nodes */
  iAutoArrayT fNodeSetDimensions; /**< number of nodes in set */
  iAutoArrayT fSideSetDimensions; /**< number of sides in set */
  AutoArrayT<SideSetScopeT> fSideSetScope; /**< flag for globally or locally numbered */
  iAutoArrayT fSideSetGroupIndex; /**< -1 for globally numbered or element group that contains set */
  /*@}*/	

  /** \name set parameters */
  /*@{*/	
  AutoArrayT<StringT> fElementNames; /**< element group IDs */
  AutoArrayT<StringT> fNodeSetNames; /**< node set IDs */
  AutoArrayT<StringT> fSideSetNames; /**< side set IDs */
  AutoArrayT<GeometryT::CodeT> fElementCodes; /** element group geometry codes */
  /*@}*/	

  /** \name model data */
  /*@{*/	
  dArray2DT fCoordinates; /**< coordinates */
  AutoArrayT<iArray2DT*> fElementSets; /**< connectivities */ 
  AutoArrayT<iArrayT*> fNodeSets; /**< node sets */
  AutoArrayT<iArray2DT*> fSideSets; /**< side sets */
  /*@}*/
  
  /** memory manager for the coordinate array */
  nVariArray2DT<double> fCoordinates_man;
  
  /** memory managers for the connectivities */
  AutoArrayT<nVariArray2DT<int>* > fElementSets_man;
};

/* return a reference to the input class */
inline InputBaseT& ModelManagerT::Input(const char* caller) const {
	const char this_caller[] = "ModelManagerT::Input"; 
	if (!fInput) 
		ExceptionT::DatabaseFail((caller) ? caller : this_caller, 
			"input source is not initialized");
	return *fInput;
}

inline int ModelManagerT::NumNodes(void) const { return fCoordinateDimensions[0]; };
inline int ModelManagerT::NumDimensions (void) const { return fCoordinateDimensions[1]; };

inline int ModelManagerT::NumElementGroups (void) const { return fElementSets.Length(); }
inline int ModelManagerT::NumNodeSets (void) const { return fNodeSets.Length(); }
inline int ModelManagerT::NumSideSets (void) const { return fSideSets.Length(); }

inline int ModelManagerT::NumElements(void) const
{
	return Input("NumTimeSteps").NumGlobalElements();
}

inline int ModelManagerT::NumTimeSteps(void) { 
	return Input("NumTimeSteps").NumTimeSteps(); 
}

inline void ModelManagerT::TimeSteps(dArrayT& steps) {
	Input("TimeSteps").ReadTimeSteps(steps); 
}

inline int ModelManagerT::NumNodeVariables(void) { 
	return Input("NumNodeVariables").NumNodeVariables (); 
}

inline void ModelManagerT::NodeLabels(ArrayT<StringT>& labels) { 
	labels.Dimension(NumNodeVariables());
	Input("NodeLabels").ReadNodeLabels(labels); 
}

inline void ModelManagerT::NodeVariablesUsed(const StringT& ID, iArrayT& used) { 
	Input("NodeVariablesUsed").NodeVariablesUsed(ID, used); 
}

inline void ModelManagerT::AllNodeVariable(int step, int varindex, dArrayT& values) { 
	Input("AllNodeVariable").ReadAllNodeVariable (step, varindex, values); 
}

inline void ModelManagerT::NodeVariable(int step, const StringT& ID, int varindex, dArrayT& values) { 
	Input("NodeVariable").ReadNodeVariable(step, ID, varindex, values); 
}

inline void ModelManagerT::AllNodeVariables(int stepindex, dArray2DT& values) { 
	Input("AllNodeVariables").ReadAllNodeVariables (stepindex, values); 
}

inline void ModelManagerT::NodeVariables(int stepindex, const StringT& ID, dArray2DT& values) {
	Input("NodeVariables").ReadNodeVariables (stepindex, ID, values); 
}

inline void ModelManagerT::NodeSetVariables(int stepindex, const StringT& ID, dArray2DT& values) {
	Input("NodeSetVariables").ReadNodeSetVariables (stepindex, ID, values); 
}

inline int ModelManagerT::NumElementVariables(void) { 
	return Input("NumElementVariables").NumElementVariables (); 
}

inline void ModelManagerT::ElementLabels(ArrayT<StringT>& labels) { 
	labels.Dimension(NumElementVariables());
	Input("ElementLabels").ReadElementLabels(labels); 
}

inline void ModelManagerT::ElementVariablesUsed(const StringT& ID, iArrayT& used) {
	Input("ElementVariablesUsed").ElementVariablesUsed(ID, used); 
}	

inline void ModelManagerT::AllElementVariable(int step, int varindex, dArrayT& values) { 
	Input("AllElementVariable").ReadAllElementVariable (step, varindex, values); 
}

inline void ModelManagerT::ElementVariable(int step, const StringT& ID, int varindex, dArrayT& values) {
	Input("ElementVariable").ReadElementVariable (step, ID, varindex, values); 
}

inline void ModelManagerT::AllElementVariables (int stepindex, dArray2DT& values) { 
	Input("AllElementVariables").ReadAllElementVariables (stepindex, values); 
}

inline void ModelManagerT::ElementVariables(int stepindex, const StringT& ID, dArray2DT& values) {
	Input("ElementVariables").ReadElementVariables (stepindex, ID, values); 
}

inline int ModelManagerT::NumElementQuadPoints(const StringT& ID) { 
	return Input("NumElementQuadPoints").NumElementQuadPoints(ID); 
}

inline int ModelManagerT::NumQuadratureVariables(void) { 
	return Input("NumQuadratureVariables").NumQuadratureVariables (); 
}

inline void ModelManagerT::QuadratureLabels(ArrayT<StringT>& labels) { 
	Input("QuadratureLabels").ReadQuadratureLabels (labels); 
}

inline void ModelManagerT::QuadratureVariablesUsed(const StringT& ID, iArrayT& used) {
	Input("QuadratureVariablesUsed").QuadratureVariablesUsed (ID, used); 
}

inline void ModelManagerT::AllQuadratureVariable(int step, int varindex, dArrayT& values) { 
	Input("AllQuadratureVariable").ReadAllQuadratureVariable (step, varindex, values); 
}

inline void ModelManagerT::QuadratureVariable(int step, const StringT& ID, int varindex, dArrayT& values) {
	Input("QuadratureVariable").ReadQuadratureVariable (step, ID, varindex, values); 
}

inline void ModelManagerT::AllQuadratureVariables(int stepindex, dArray2DT& values) { 
	Input("AllQuadratureVariables").ReadAllQuadratureVariables (stepindex, values); 
}

inline void ModelManagerT::QuadratureVariables(int stepindex, const StringT& ID, dArray2DT& values) {
	Input("QuadratureVariables").ReadQuadratureVariables (stepindex, ID, values); 
}

inline void ModelManagerT::QARecords (ArrayT<StringT>& records) { 
	Input("QARecords").QARecords(records); 
}

inline const StringT& ModelManagerT::ElementGroupID(int index) const {
	if (index < 0 || index >= fElementNames.Length()) {
		cout << "\n ModelManagerT::ElementGroupID: index " << index << " out of range {0," 
		     << fElementNames.Length() - 1 << "}" << endl;
		throw ExceptionT::kOutOfRange;
	}
	return fElementNames[index]; 
};

inline void ModelManagerT::UpdateNode(const dArrayT& coordinates, int node)
{
	/* write in */
	fCoordinates.SetRow(node, coordinates);
}

} // namespace Tahoe 
#endif
