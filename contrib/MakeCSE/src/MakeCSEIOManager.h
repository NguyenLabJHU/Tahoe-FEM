// file: MakeCSEIOManager.h

// created: SAW 12/21/99

#ifndef _MAKECSEIOMANAGER_H_
#define _MAKECSEIOMANAGER_H_

#include "ModelManagerT.h"
#include "OutputBaseT.h"

using namespace Tahoe;

class MakeCSEIOManager : public ModelManagerT
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

  void ReadParameters (ifstreamT& in, bool interactive, const StringT& program, const StringT& version);
  void Interactive (void);

  void InputData (int& data, int key) const;
  void InputData (iArrayT& data, int key) const;

  void WriteGeometry (void);

 private:
  void InteractiveIO (void);
  void ReadInputFile (ifstreamT& in);

  void Parse (ifstreamT& in, StringT& word1);
  bool ReadWord1 (ifstreamT& in, StringT& word1) const;
  void ReadIDColumnal (ifstreamT& in, int key, int numoptions = 1);
  void Read2IDColumnal (ifstreamT& in, int key, int numoptions = 1);
  bool ReadID (ifstreamT& in, int& start, int& stop) const;
  void ReadMultiID (ifstreamT& in, int key);

  void InteractiveCSE (void);
  void InteractiveSplitElement (void);
  void Read (const char* first, int key, int num);
  void Read2D (const char* first, const char *second, int key, int num);
  void Read3D (const char* first, const char *second, const char *third, int key, int num);
  
  void ReadInputFormat (ifstreamT& in);
  void ReadOutputFormat (ifstreamT& in);
  void SetInput (void);
  void SetOutput (const StringT& program_name, 
		  const StringT& version, const StringT& title, 
		  const StringT& input_file, 
		  IOBaseT::FileTypeT output_format);

 private:
  StringT fTitle;
  IOBaseT::FileTypeT fOutputFormat;
  OutputBaseT* fOutput;
	ofstream fEchoInput;
	bool fEcho;
	bool fExternTahoeII;

  int fVerbose;
  int fRenumber;
  int fZoneEdge;
  ArrayT<iArrayT> fData;
};

#endif
