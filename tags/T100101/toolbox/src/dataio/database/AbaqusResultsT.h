/* $Id: AbaqusResultsT.h,v 1.3 2001-09-21 13:49:59 sawimme Exp $ */
/*
   CREATED: S. Wimmer 9 Nov 2000

   To add varible update these:
   1. NVT
   2. VariableKeyT
   3. VariableName ()
   4. VariableKey ()
   5. IntToVariableKey ()
   6. VariableKeyIndex ()
   7. VariableWrittenWithNodeNumber ()

*/

#ifndef _ABAQUSRESULTS_T_H_
#define _ABAQUSRESULTS_T_H_

/* direct members */
#include <fstream.h>
#include "iArrayT.h"
#include "dArrayT.h"
#include "StringT.h"
#include "dArray2DT.h"
#include "iAutoArrayT.h"
#include "AutoArrayT.h"
#include "GeometryT.h"

class AbaqusResultsT
{
 public:
  enum ElementType { kUnknown = -11,
		     kTriangle = 12,
		     kQuad = 13,
		     kHex = 14,
		     kTet = 15,
		     kWedge = 16,
		     kShell = 17 };

  enum VariableType { kQuadVar = 0,
		      kElemVar = 1,
		      kNodeVar = 2 };

  enum NumVariables { NVT = 26 };
  enum VariableKeyT { kNone = -1,
		      kTEMP = 2, kLOADS = 3, kFLUXS = 4,
		      kSDV = 5, kCOORD = 8,
		      kS = 11, kSINV = 12,
		      kE = 21, kPE = 22, kCE = 23, kIE = 24, kEE = 25,
		      kU = 101, kV = 102, kA = 103, kNCOORD = 107,
		      kSP = 401,
		      kEP = 403, kNEP = 404, kLEP = 405, kERP = 406,
		      kEEP = 408, kIEP = 409, kTHEP = 410, kPEP = 411, kCEP = 412 };

  AbaqusResultsT (ostream& message);
  void Initialize (char *filename);
  void Close (void);
  void ScanFile (int &numelems, int &numnodes, int &numtimesteps, int &nummodes);

  void ElementSetNames (ArrayT<StringT>& names) const;
  void NodeSetNames (ArrayT<StringT>& names) const;
  void NodeMap (iArrayT& n) const;
  int  NumNodesInSet (StringT& name);
  void NodeSet (StringT& name, iArrayT& nset);

  void ElementMap (iArrayT& e) const;
  int  NumElements (StringT& name);
  int  NumElementNodes (StringT& name);
  int  NumElementQuadPoints (StringT& name);
  void ElementSet (StringT& name, iArrayT& elset);
  void GeometryCode (StringT& name, GeometryT::CodeT& code);

  int  NumNodeVariables (void) const;
  int  NumElementVariables (void) const;
  int  NumQuadratureVariables (void) const;

  void ModeData (int index, int &number, double &mode) const;
  void TimeData (int index, int &number, double &time) const;

  void NodeVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const;
  void ElementVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const;
  void QuadratureVariables (ArrayT<VariableKeyT>& keys, iArrayT& dims) const;

  void ReadVariables (VariableType vt, int step, dArray2DT& values, StringT& name);

  const char* VariableName (int index) const;
  VariableKeyT IntToVariableKey (int key) const;
  VariableKeyT VariableKey (int index) const;
  int VariableKeyIndex (VariableKeyT key) const;

  bool NextCoordinate (int &number, dArrayT& nodes);
  bool NextElement (int &number, GeometryT::CodeT &type, iArrayT &nodes);

  void VersionNotes (ArrayT<StringT>& records);
  void ResetFile (void);
 
  int NumElements (void) const;
  int ElementNumber (int index) const;
  int VariableDimension (int index) const;
  int NodeNumber (int index) const;
  int NumElementSets (void) const;
  int NumNodeSets (void) const;

 private:
  bool ReadVersion (void);
  void NextMode (int &number, double &mode);
  void NextTimeSteps (int &number, double &time);
  void ScanElement (void);
  void ReadOutputDefinitions (int &outputmode);
  void ReadElementHeader (int& objnum, int& intpt, int& secpt, int &location);
  void ScanVariable (AbaqusResultsT::VariableKeyT key, int outputmode, int location);

  bool VariableWrittenWithNodeNumber (AbaqusResultsT::VariableKeyT key) const;
  bool CorrectType (int outputmode, int objnum, int intpt, int location, VariableType vt, int& ID) const;

  int TranslateElementName (char *, GeometryT::CodeT &, int &);
  int TranslateContinuum (char *, GeometryT::CodeT &, int &);
  int Translate2D (char *, GeometryT::CodeT &, int &);
  int Translate3D (char *, GeometryT::CodeT &, int &);
  int TranslateShell (char *, GeometryT::CodeT &, int &);
  
  void AdvanceTo (int target);
  bool SkipAttributes (void);
  int  ReadNextRecord (int &key);

  bool Read (StringT& s, int n);
  bool Read (int& i);
  bool Read (double& d);
  bool CheckBufferSize (istream& in, int numchars);
  void CheckBufferSize (istream& in);

  enum OutputType { kElementOutput = 0,
		    kNodalOutput = 1,
		    kModalOutput = 2,
		    kElemSetOutput = 3 };

  enum GeneralKeys { ELEMENTHEADER = 1,
		     ELEMENT = 1900,
		     NODE = 1901,
		     ACTIVEDOF = 1902,
		     OUTPUTDEFINE = 1911,
		     VERSION = 1921,
		     HEADING = 1922,
		     NODESET = 1931,
		     NODESETCONT = 1932,
		     ELEMENTSET = 1933,
		     ELEMSETCONT = 1934,
		     LABELREF = 1940,
		     MODAL = 1980,
		     STARTINCREMENT = 2000,
		     ENDINCREMENT = 2001 };

  enum OutputParamsT { dprecision = 15 };

  enum StatusT { OKAY = -101, BAD = -102, END = -103 };

 private:
  ifstream fIn;
  ostream& fMessage;
  StringT fFileName;

  bool fBinary;
  int fBufferDone;
  int fBufferSize;
  StringT fBuffer;
  int fCurrentLength;

  int fNumNodes;
  int fNumElements;

  int fNumNodeSets;
  AutoArrayT<StringT> fNodeSetNames;
  int fNumElementSets;
  AutoArrayT<StringT> fElementSetNames;

  int fStartCount;
  int fEndCount;
  int fModalCount;

  iArrayT fVarDimension;
  iArrayT fVarType;
  iArrayT fVarArrayColumn;
  int fNumNodeVars;
  int fNumElemVars;
  int fNumQuadVars;

  iAutoArrayT fElementNumber;
  iAutoArrayT fNodeNumber;

  iAutoArrayT fTimeIncs;
  AutoArrayT<double> fTimeSteps;

  iAutoArrayT fModeIncs;
  AutoArrayT<double> fModeSteps; 
};

inline int AbaqusResultsT::NumElements (void) const { return fNumElements; }
inline int AbaqusResultsT::ElementNumber (int index) const { return fElementNumber[index]; }
inline int AbaqusResultsT::VariableDimension (int index) const { return fVarDimension[index]; }
inline int AbaqusResultsT::NodeNumber (int index) const { return fNodeNumber[index]; }
inline int AbaqusResultsT::NumElementSets (void) const { return fNumElementSets; }
inline int AbaqusResultsT::NumNodeSets (void) const { return fNumNodeSets; }
inline int AbaqusResultsT::NumNodeVariables (void) const { return fNumNodeVars; }
inline int AbaqusResultsT::NumElementVariables (void) const { return fNumElemVars; }
inline int AbaqusResultsT::NumQuadratureVariables (void) const { return fNumQuadVars; }

#endif
