/* $Id: InputBaseT.h,v 1.13 2002-07-05 22:26:26 paklein Exp $ */
/* created: sawimme (08/12/1999) */

#ifndef _INPUTBASE_T_H_
#define _INPUTBASE_T_H_

#include "IOBaseT.h"

#include "GeometryT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* foward declaration */
class iArrayT;
class iArray2DT;
class dArrayT;
class dArray2DT;
class StringT;
template <class TYPE> class ArrayT;
class ModelManagerT;

/** derived classes must:\n
 * 1. offset connectivity to start node numbering at zero\n
 * 2. consecutively number elements\n
 * 3. consecutively number coordinates\n
 * 4. offset node sets to start node numbering at zero\n
 * 5. offset side sets to start element and facet numbering at zero\n
 * 6. side sets use offset global element numbers and offset local facets\n
 * 7. Node and Element Maps are not offset to start at zero and are global.\n
 * 8. Connectivity uses offset global consecutive node numbering.\n
 * 9. The derived class must check dimension arrays before setting arrays.\n
 *10. ReadVariable functions return node values for all nodes
 *        or nodes in node set or just nodes used by an element group.\n
 *11. ReadVariable functions return only elements in group for element
 *        and quadrature data.\n
 *12. For ReadVariable functions, the step is an integer starting at zero,
 *        increasing consecutively.\n
 *13. ReadVariable functions return values for all labels.\n
 *
 * Most functions are not const due to the fact that databases like 
 * ABAQUS need to keep track of the size of the buffer read. */
class InputBaseT : public IOBaseT
{
public:

  /** constructor */
  InputBaseT(ostream& out);

  /** destructor */
  virtual ~InputBaseT (void);

  /** return the database string version of the element block ID */

  /** open the input source. \return true if successful, false otherwise */
  virtual bool Open (const StringT& filename) = 0;

  /** close the input source */
  virtual void Close (void) = 0;

  /** return names, Array must be preallocated */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const = 0;
  /** return names, Array must be preallocated */
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const = 0;
  /** return names, Array must be preallocated */
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const = 0;

  /** return dimenesions */
  virtual int  NumElementGroups (void) const = 0;
  /** return dimenesions */
  virtual int  NumSideSets (void) const = 0;
  /** return dimenesions */
  virtual int  NumNodeSets (void) const = 0;

  /* NODES */
  virtual int  NumNodes (void) const = 0;
  virtual int  NumDimensions (void) const = 0; /**< should return num dims to be used */

	/** ids for all nodes. ids for all nodes in the coordinate list. ids may not be
	 * compact or ordered */
	virtual void ReadNodeID(iArrayT& node_id) = 0; 

	/** read coordinates */
	virtual void ReadCoordinates(dArray2DT& coords) = 0;
	
	/** read coordinates and ids. ids for all nodes in the coordinate list. 
	 * ids may not be compact or ordered */
	virtual void ReadCoordinates(dArray2DT& coords, iArrayT& node_id) = 0;

  /* ELEMENTS */
  virtual int  NumGlobalElements (void) const = 0; /**< for all element sets */
  virtual int  NumElements (const StringT& name) = 0; /**< for the set specified */
  virtual int  NumElementNodes (const StringT& name) = 0; /**< typically for the first element in the set */
  virtual int  NumElementQuadPoints (const StringT& name) = 0; /**< typically for the first element in the set */
  virtual void ReadAllElementMap (iArrayT& elemmap) = 0; /**< all elements, not offset, can be discontinuous */
  virtual void ReadGlobalElementMap (const StringT& name, iArrayT& elemmap) = 0; /**< set elements, not offset, can be discontinuous */
  virtual void ReadGlobalElementSet (const StringT& name, iArrayT& set) = 0; /**< offset, continuous */
  virtual void ReadConnectivity (const StringT& name, iArray2DT& connects) = 0; /**< offset nodes, continuous */
  virtual void ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode) = 0;

  virtual int  NumNodesInSet (const StringT& name) = 0;
  virtual void ReadNodeSet (const StringT& name, iArrayT& nodes) = 0; /**< offset nodes, continuous */

  virtual bool AreSideSetsLocal (void) const = 0;
  virtual int  NumSidesInSet (const StringT& setname) const = 0;
  virtual StringT SideSetGroupName (const StringT& setname) const = 0;
  virtual void ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const = 0; /**< offset elements & facets, continuous */
  virtual void ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const = 0; /**< offset elements & facets, continuous */
  
  /** record[0] = progname, record[1] = version, record[2] = date, record[3] = time */
  virtual void QARecords (ArrayT<StringT>& records) = 0;

  virtual int  NumTimeSteps (void) const = 0;
  virtual void ReadTimeSteps (dArrayT& steps) = 0;

  /** for all nodes or elements */
  virtual int  NumNodeVariables (void) const = 0;
  virtual int  NumElementVariables (void) const = 0;
  virtual int  NumQuadratureVariables (void) const = 0;

  /** for only the nodes or elements in the block */
  virtual void NodeVariablesUsed (const StringT& name, iArrayT& used) = 0;
  virtual void ElementVariablesUsed (const StringT& name, iArrayT& used) = 0;
  virtual void QuadratureVariablesUsed (const StringT& name, iArrayT& used) = 0;  

  /** for all nodes or elements */
  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const = 0;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const = 0;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const = 0;  

  /** step starts at zero and increases by one,
   * varindex refers to the variables index position in the Label list */
  virtual void ReadAllNodeVariable (int step, int varindex, dArrayT& values) = 0; /**< one variables for all nodes */
  virtual void ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values) = 0; /**< one variable for nodes in an element set */
  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues) = 0; /**< all variables for all nodes */
  virtual void ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues) = 0; /**< all variables for nodes in an element set */
  virtual void ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues) = 0; /**< all variables for nodes in a node set */

  virtual void ReadAllElementVariable (int step, int varindex, dArrayT& values) = 0; /** < one variable for all elements */
  virtual void ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values) = 0; /** < one variable for an element set */
  virtual void ReadAllElementVariables (int step, dArray2DT& evalues) = 0; /**< all variables for all elements */
  virtual void ReadElementVariables (int step, const StringT& name, dArray2DT& evalues) = 0; /**< all variables for elements in set */

  virtual void ReadAllQuadratureVariable (int step, int varindex, dArrayT& values) = 0; /** < one variable for all elements */
  virtual void ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values) = 0; /** < one variable for an element set */
  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues) = 0; /**< all variables for all quad points */
  virtual void ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues) = 0; /**< all variables for elements in set */

};

inline InputBaseT::InputBaseT (ostream& out) : IOBaseT (out) { }
inline InputBaseT::~InputBaseT (void) { }

} // namespace Tahoe 
#endif
