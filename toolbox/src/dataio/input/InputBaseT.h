/* $Id: InputBaseT.h,v 1.2 2001-08-03 19:16:43 sawimme Exp $ */
/* created: sawimme (08/12/1999)                                          */

/* A derived classes must
   1. offset connectivity to start node numbering at zero 
   2. consecutively number elements within a group
   3. consecutively number coordinates
   4. offset node sets to start node numbering at zero
   5. offset side sets to start element and facet numbering at zero
   6. side sets use offset global element numbers and offset local facets
   7. Node and Element Maps are not offset to start at zero and are global.
   8. Connectivity uses offset global consecutive node numbering.
   9. The derived class must check dimension arrays before setting arrays.
  10. ReadVariable functions return node values for all nodes
          not just nodes used by a particular element group.
  11. ReadVariable functions return only elements in group for element
          and quadrature data.
  12. For ReadVariable functions, the step is an integer starting at zero,
          increasing consecutively.

   Most functions are not const due to the fact that databases like 
   ABAQUS need to keep track of the size of the buffer read. 
*/

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

class InputBaseT : public IOBaseT
{
public:
  InputBaseT (ostream& out);
  virtual ~InputBaseT (void);

  virtual void Open (const StringT& filename) = 0;
  virtual void Close (void) = 0;

  /* virtual with derived classes */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames);
  virtual void SideSetNames (ArrayT<StringT>& sidenames);
  virtual void NodeSetNames (ArrayT<StringT>& nodenames);

  virtual int  NumElementGroups (void) = 0;
  virtual int  NumSideSets (void) = 0;
  virtual int  NumNodeSets (void) = 0;

  virtual int  NumNodes (void) = 0;
  virtual int  NumDimensions (void) = 0;
  virtual void ReadNodeMap (iArrayT& nodemap) = 0;
  virtual void ReadCoordinates (dArray2DT& coords) = 0;
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap) = 0;

  virtual int  NumGlobalElements (void) = 0;
  virtual int  NumElements (StringT& name);
  virtual int  NumElementNodes (StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap) = 0;
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (StringT& name);
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void) = 0;
  virtual int  NumSidesInSet (StringT& name);
  virtual int  SideSetGroupIndex (StringT& name);
  virtual void ReadSideSetLocal (StringT& name, iArray2DT& sides);
  virtual void ReadSideSetGlobal (StringT& name, iArray2DT& sides);
  
  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void);
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void);
  virtual int  NumElementVariables (void);
  virtual int  NumQuadratureVariables (void);

  //virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) = 0;
  virtual void ReadElementLabels (StringT& name, ArrayT<StringT>& elabels);
  virtual void ReadQuadratureLabels (StringT& name, ArrayT<StringT>& qlabels);  

  //virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues) = 0;
  //virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues);
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);

 private:
  virtual void ElementGroupIDs (iArrayT& groupnums);
  virtual void SideSetIDs (iArrayT& sidenums);
  virtual void NodeSetIDs (iArrayT& nodenums);

  virtual int  NumElements_ID (int ID);
  virtual int  NumElementNodes_ID (int ID);
  virtual void ReadGlobalElementMap_ID (int ID, iArrayT& elemmap);
  virtual void ReadConnectivity_ID (int ID, iArray2DT& connects);
  virtual void ReadGeometryCode_ID (int ID, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet_ID (int ID);
  virtual void ReadNodeSet_ID (int ID, iArrayT& nodes);

  virtual int  NumSidesInSet_ID (int ID);
  virtual int  SideSetGroupIndex_ID (int ID);
  virtual void ReadSideSetLocal_ID (int ID, iArray2DT& sides);
  virtual void ReadSideSetGlobal_ID (int ID, iArray2DT& sides);

  virtual void ReadElementLabels_ID (int ID, ArrayT<StringT>& elabels);
  virtual void ReadQuadratureLabels_ID (int ID, ArrayT<StringT>& qlabels);  

  virtual void ReadElementVariables_ID (int step, int ID, dArray2DT& evalues);
  virtual void ReadQuadratureVariables_ID (int step, int ID, dArray2DT& qvalues);

};

#endif
