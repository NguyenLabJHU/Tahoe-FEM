/* $Id: EnSightInputT.h,v 1.2 2001-08-03 19:16:43 sawimme Exp $ */
/* created: sawimme (05/18/1998)                                          */

#ifndef _ENSIGHTINPUT_T_H_
#define _ENSIGHTINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "EnSightT.h"
#include "StringT.h"
#include "iArray2DT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
template <class TYPE> class ArrayT;
class dArray2DT;
class iArray2DT;

class EnSightInputT : public InputBaseT
{
public:
  EnSightInputT (ostream& out, bool binary);

  virtual void Open (const StringT& file);
  virtual void Close (void);

  /* virtual with InputManager base class */
  virtual int  NumElementGroups (void);
  virtual int  NumSideSets (void);
  virtual int  NumNodeSets (void);

  virtual int  NumNodes (void);
  virtual int  NumDimensions (void);
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual bool AreSideSetsLocal (void);

  virtual int  NumGlobalElements (void);
  virtual void ReadAllElementMap (iArrayT& elemmap);

  virtual int  NumTimeSteps (void);
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void);
  virtual int  NumElementVariables (void);
  
 protected:
  virtual void ElementGroupIDs (iArrayT& groupnums);
  
  virtual int NumElements_ID (int ID);
  virtual int NumElementNodes_ID (int ID);
  virtual void ReadGlobalElementMap_ID (int ID, iArrayT& elemmap);
  virtual void ReadConnectivity_ID (int ID, iArray2DT& connects);
  virtual void ReadGeometryCode_ID (int ID, GeometryT::CodeT& code);

  virtual void ReadElementLabels_ID (int ID, ArrayT<StringT>& elabels);
  virtual void ReadElementVariables_ID (int step, int ID, dArray2DT& evalues);

 private:
  bool AdvanceStream (istream& in, const char* key) const;
  void ScanGeometryFile (void);
  
  StringT CreateVariableFile (const StringT& old, int inc) const;
  void ReadVariableData (ArrayT<bool>& vector, ArrayT<StringT>& labels, int group_id, dArray2DT& values, int currentinc, bool nodal) const;
  
 private:
  EnSightT fData;
  StringT fGeometryFile;
  StringT fCaseFile;
  iArray2DT fPartDimensions; // num_nodes, num_elems, partID
  int fStartIncrement;
  int fIncrement;
};

inline void EnSightInputT::Close (void) { }
inline int EnSightInputT::NumSideSets (void) { return 0; }
inline int EnSightInputT::NumNodeSets (void) { return 0; }
inline int EnSightInputT::NumDimensions (void) { return 3; }
inline bool EnSightInputT::AreSideSetsLocal (void) { return true; }

#endif
