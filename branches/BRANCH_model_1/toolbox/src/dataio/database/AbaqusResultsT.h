/* $Id: AbaqusResultsT.h,v 1.3.2.1 2001-11-06 14:20:39 sawimme Exp $ */
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
#include "iArray2DT.h"
#include "AbaqusVariablesT.h"

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

  enum AnalysisTypeT { kStatic = 1,
		       kDynamic = 12 };

  AbaqusResultsT (ostream& message);

  void Initialize (char *filename);

  void Create (char *filename, bool binary, int numelems, int numnodes, 
	       double elemsize);
  void OpenWrite (char *filename, bool binary, int bufferwritten);

  /** close file and return amount of buffer written */
  int Close (void);

  void ScanFile (int &numelems, int &numnodes, int &numtimesteps, 
		 int &nummodes);

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

  void NodeVariables (iArrayT& keys, iArrayT& dims) const;
  void ElementVariables (iArrayT& keys, iArrayT& dims) const;
  void QuadratureVariables (iArrayT& keys, iArrayT& dims) const;

  void ReadVariables (AbaqusVariablesT::TypeT vt, int step, dArray2DT& values,
		      StringT& name);

  const char* VariableName (int index) const;
  int VariableKey (const char *name) const;
  int VariableKey (int index) const;
  int VariableKeyIndex (int key) const;

  bool NextCoordinate (int &number, dArrayT& nodes);
  bool NextElement (int &number, GeometryT::CodeT &type, iArrayT &nodes);

  void WriteConnectivity (GeometryT::CodeT code, int startnumber, 
			  const iArray2DT& connects);
  void WriteCoordinates (const iArrayT& nodes_used, const dArray2DT& coords);
  void WriteElementSet (const StringT& name, const iArrayT& elms);
  void WriteNodeSet (const StringT& name, const iArrayT& nodes);
  void WriteActiveDOF (const iArrayT& active);
  void WriteHeading (const StringT& heading);

  void WriteStartIncrement (int step, int inc, double totaltime, 
     double time, double timeincrement, AbaqusResultsT::AnalysisTypeT atype);

  void WriteOutputDefinition (int key, const StringT& setname, GeometryT::CodeT code, 
			      int numelemnodes);
  void WriteNodeVariables (int &index, const iArrayT& keys, const dArray2DT& values, 
			   const iArrayT& nodes_used, int numdir, int numshear);
  void WriteEndIncrement (void);

  void VersionNotes (ArrayT<StringT>& records);
  void ResetFile (void);
 
  int NumElements (void) const;
  int ElementNumber (int index) const;
  int VariableDimension (int index) const;
  int NodeNumber (int index) const;
  int NumElementSets (void) const;
  int NumNodeSets (void) const;

 private:

  enum NumVariables { NVT = 23 };

  enum OutputType { kElementOutput = 0,
		    kNodalOutput = 1,
		    kModalOutput = 2,
		    kElemSetOutput = 3 };

  enum ElementVarType { kElementQuadrature = 0,
			kElementCentroidal = 1,
			kElementNodal = 2,
			kElementRebar = 3,
			kElementNodeAveraged = 4,
			kElementWhole = 5};
  
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

  bool ReadVersion (void);
  void NextMode (int &number, double &mode);
  void NextTimeSteps (int &number, double &time);
  void ScanElement (void);
  void ReadOutputDefinitions (int &outputmode);
  void ReadElementHeader (int& objnum, int& intpt, int& secpt, int &location);
  void ScanVariable (int key, int outputmode, int location);

  void WriteOutputDefinition (int key, const StringT& setname, 
			      GeometryT::CodeT code, int numnodes,
			      int count, int& numdir, int& numshear);
  void WriteElementHeader (int key, int number, int intpt, int secpt, 
			   AbaqusResultsT::ElementVarType flag, int numdirect, 
			   int numshear, int numdir, int numsecforc); 

  bool VariableWrittenWithNodeNumber (int key) const;
  bool CorrectType (int outputmode, int objnum, int intpt, int location, 
		    AbaqusVariablesT::TypeT vt, int& ID) const;

  int TranslateElementName (char *, GeometryT::CodeT &, int &);
  int TranslateContinuum (char *, GeometryT::CodeT &, int &);
  int Translate2D (char *, GeometryT::CodeT &, int &);
  int Translate3D (char *, GeometryT::CodeT &, int &);
  int TranslateShell (char *, GeometryT::CodeT &, int &);
  void GetElementName (GeometryT::CodeT geometry_code, int elemnodes, 
		       int& num_output_nodes, StringT& elem_name) const;
  
  void AdvanceTo (int target);
  bool SkipAttributes (void);
  int  ReadNextRecord (int &key);

  bool Read (StringT& s, int n);
  bool Read (int& i);
  bool Read (double& d);
  bool CheckBufferSize (istream& in, int numchars);
  void CheckBufferSize (istream& in);

  void Write (int i);
  void Write (double d);
  void Write (const StringT& s, int blocks = 1);
  void WriteASCII (const StringT& s);
  void CheckBufferSize (ostream& out);

  void SetVariableNames (void);

 private:
  ifstream fIn;
  ofstream fOut;
  ostream& fMessage;
  StringT fFileName;
  const StringT fMarker;
  const StringT fOutVersion;

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

  ArrayT<AbaqusVariablesT> fVariableTable;
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
inline int AbaqusResultsT::VariableDimension (int index) const { return fVariableTable[index].Dimension(); }
inline int AbaqusResultsT::NodeNumber (int index) const { return fNodeNumber[index]; }
inline int AbaqusResultsT::NumElementSets (void) const { return fNumElementSets; }
inline int AbaqusResultsT::NumNodeSets (void) const { return fNumNodeSets; }
inline int AbaqusResultsT::NumNodeVariables (void) const { return fNumNodeVars; }
inline int AbaqusResultsT::NumElementVariables (void) const { return fNumElemVars; }
inline int AbaqusResultsT::NumQuadratureVariables (void) const { return fNumQuadVars; }

#endif
