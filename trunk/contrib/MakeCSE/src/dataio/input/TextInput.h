// file: TextInput.h

// created:       SAW (08/11/99)

#ifndef _TEXTINPUT_H_
#define _TEXTINPUT_H_

#include "InputBaseT.h"

/* direct members */
#include "ModelFileT.h"

class TextInput : public InputBaseT
{
 public:
  TextInput (ostream& out, const StringT& filename);

  /* virtual with derived classes */
  virtual int  NumElementGroups (void) const;
  virtual int  NumSideSets (void) const;
  virtual int  NumNodeSets (void) const;
  virtual void GroupNumbers (iArrayT& groupnums) const;
  virtual void SideSetNumbers (iArrayT& sidenums) const;
  virtual void NodeSetNumbers (iArrayT& nodenums) const;

  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);
  virtual void ReadConnectivity (int group, GeometryT::GeometryCode& geocode, iArray2DT& connects, iArrayT& elementmap);
  virtual void ReadNodeSet (int set_num, iArrayT& nodes) const;
  virtual void ReadSideSet (int set_num, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (int set_num, iArray2DT& sides) const;
  virtual void Close (void);
  virtual void QARecords (ArrayT<StringT>& records) const;
  virtual void ReadTimeSteps (dArrayT& steps);
  virtual void ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id);
  virtual void ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues);

 private:
  ModelFileT fData;
};

#endif
