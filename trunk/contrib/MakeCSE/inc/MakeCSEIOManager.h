// file: MakeCSEIOManager.h

// created: SAW 12/21/99

#ifndef _MAKECSEIOMANAGER_H_
#define _MAKECSEIOMANAGER_H_

#include "IOManager.h"

class MakeCSEIOManager : public IOManager
{
 public:

  enum DataType { kFacet,
		  kZone,
		  kZoneEdge,
		  kZoneEdgeNSets,
		  kBoundary,
		  kElementSplit,
		  kSingleNodes,
		  kContactData,
		  kRenumber,
		  kExecution,
		  kMapNodes,
		  kCopySide,
		  kBlockToNode,
		  kNumDataTypes };

  MakeCSEIOManager (ostream& out);

  virtual void Interactive (void);

  virtual void InputData (int& data, int key) const;
  virtual void InputData (iArrayT& data, int key) const;

 private:
  virtual void Parse (ifstream_x& in, StringT& word1);
  void ReadIDColumnal (ifstream_x& in, int key, int numoptions = 1);
  void Read2IDColumnal (ifstream_x& in, int key, int numoptions = 1);
  bool ReadID (ifstream_x& in, int& start, int& stop) const;
  void ReadMultiID (ifstream_x& in, int key);

  void InteractiveCSE (void);
  void InteractiveSplitElement (void);
  void Read (char* first, int key, int num);
  void Read2D (char* first, char *second, int key, int num);
  void Read3D (char* first, char *second, char *third, int key, int num);
  

 private:
  int fVerbose;
  int fRenumber;
  int fZoneEdge;
  ArrayT<iArrayT> fData;
};

#endif
