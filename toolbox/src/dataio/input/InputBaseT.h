/* $Id: InputBaseT.h,v 1.5.2.1 2001-10-15 19:01:27 sawimme Exp $ */
/* created: sawimme (08/12/1999) */

#ifndef _INPUTBASE_T_H_
#define _INPUTBASE_T_H_

#include "IOBaseT.h"

#include "GeometryT.h"

/* foward declaration */
#include "ios_fwd_decl.h"
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
 *
 * Most functions are not const due to the fact that databases like 
 * ABAQUS need to keep track of the size of the buffer read. */
class InputBaseT : public IOBaseT
{
public:
  InputBaseT (ostream& out);
  virtual ~InputBaseT (void);

  virtual void Open (const StringT& filename) = 0;
  virtual void Close (void) = 0;

  /* return names, Array must be preallocated */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const = 0;
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const = 0;
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const = 0;

  /* return dimenesions */
  virtual int  NumElementGroups (void) const = 0;
  virtual int  NumSideSets (void) const = 0;
  virtual int  NumNodeSets (void) const = 0;

  /* NODES */
  virtual int  NumNodes (void) const = 0;
  virtual int  NumDimensions (void) const = 0; /* should return num dims to be used */
  virtual void ReadNodeMap (iArrayT& nodemap) = 0; /* all nodes, not offset, can be discontinuous */
  virtual void ReadCoordinates (dArray2DT& coords) = 0;
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap) = 0;

  /* ELEMENTS */
  virtual int  NumGlobalElements (void) const = 0; /* for all element sets */
  virtual int  NumElements (StringT& name) = 0; /* for the set specified */
  virtual int  NumElementNodes (StringT& name) = 0; /* typically for the first element in the set */
  virtual int  NumElementQuadPoints (StringT& name) = 0; /* typically for the first element in the set */
  virtual void ReadAllElementMap (iArrayT& elemmap) = 0; /* all elements, not offset, can be discontinuous */
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap) = 0; /* set elements, not offset, can be discontinuous */
  virtual void ReadGlobalElementSet (StringT& name, iArrayT& set) = 0; /* offset, continuous */
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects) = 0; /* offset nodes, continuous */
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode) = 0;

  virtual int  NumNodesInSet (StringT& name) = 0;
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes) = 0; /* offset nodes, continuous */

  virtual bool AreSideSetsLocal (void) const = 0;
  virtual int  NumSidesInSet (StringT& setname) const = 0;
  virtual StringT SideSetGroupName (StringT& setname) const = 0;
  virtual void ReadSideSetLocal (StringT& setname, iArray2DT& sides) const = 0; /* offset elements & facets, continuous */
  virtual void ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const = 0; /* offset elements & facets, continuous */
  
  /* record[0] = progname, record[1] = version, record[2] = date, record[3] = time */
  virtual void QARecords (ArrayT<StringT>& records) = 0;

  virtual int  NumTimeSteps (void) const = 0;
  virtual void ReadTimeSteps (dArrayT& steps) = 0;

  /* for all nodes or elements */
  virtual int  NumNodeVariables (void) const = 0;
  virtual int  NumElementVariables (void) const = 0;
  virtual int  NumQuadratureVariables (void) const = 0;

  /* for all nodes or elements */
  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const = 0;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const = 0;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const = 0;  

  /* step starts at zero and increases by one */
  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues) = 0; /* all nodes */
  virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues) = 0; /* nodes in an element set */
  virtual void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues) = 0; /* nodes in a node set */

  virtual void ReadAllElementVariables (int step, dArray2DT& evalues) = 0; /* all elements */
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues) = 0; /* elements in set */

  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues) = 0; /* all quad points */
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues) = 0; /* elements in set */

};

inline InputBaseT::InputBaseT (ostream& out) : IOBaseT (out) { }
inline InputBaseT::~InputBaseT (void) { }

#endif
