/* $Id: AbaqusResultsT.h,v 1.12 2002-07-02 19:57:00 cjkimme Exp $ */
/*
   CREATED: S. Wimmer 9 Nov 2000

   To add varible update these:
   1. enum NumVariables { NVT= }
   2. SetVariableNames

*/

#ifndef _ABAQUSRESULTS_T_H_
#define _ABAQUSRESULTS_T_H_

/* direct members */
#include <fstream.h>
#include "ifstreamT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "StringT.h"
#include "dArray2DT.h"
#include "iAutoArrayT.h"
#include "AutoArrayT.h"
#include "GeometryT.h"
#include "iArray2DT.h"
#include "AbaqusVariablesT.h"


namespace Tahoe {

class AbaqusResultsT
{
 public:

  /** geometry code, until Tahoe gets shell elements */
  enum ElementType { kUnknown = -11,
		     kTriangle = 12,
		     kQuad = 13,
		     kHex = 14,
		     kTet = 15,
		     kWedge = 16,
		     kShell = 17 };

  /** analysis types */
  enum AnalysisTypeT { kStatic = 1,
		       kDynamic = 12 };

  /** constructor */
  AbaqusResultsT (ostream& message);

  /** used to open a database, either binary or ASCII, it will determine.
   * \return true if successful, false otherwise. */
  bool Initialize (const char *filename);

  /** used to create a new database
   * \param filename name of file
   * \param binary binary or ASCII format
   * \param numelems number of elements
   * \param numnodes number of nodes 
   * \param elemsize characteristic element size */
  void Create (const char *filename, bool binary, int numelems, int numnodes, 
	       double elemsize);

  /** used to append variable data to preexisting file 
   * \param filename name of file
   * \param binary binary or ASCII format
   * \param bufferwritten amount of buffer already written to file */
  void OpenWrite (const char *filename, bool binary, int bufferwritten);

  /** close file and return amount of buffer written */
  int Close (void);

  /** scan input database after initializing. 
      \return true if scan successful, false otherwise */
  bool ScanFile (int &numelems, int &numnodes, int &numtimesteps, int &nummodes);

  /** access element set names */
  void ElementSetNames (ArrayT<StringT>& names) const;

  /** access node set names */
  void NodeSetNames (ArrayT<StringT>& names) const;

  /** acces node map */
  void NodeMap (iArrayT& n) const;

  /** acces num nodes in a set */
  int  NumNodesInSet (const StringT& name) const;

  /** return the node set */
  void NodeSet (const StringT& name, iArrayT& nset) const;

  void ElementMap (iArrayT& e) const; /**< access element map */
  int  NumElements (const StringT& name) const; /**< number of elements in a set */
  int  NumElementNodes (const StringT& name); /**< number of element nodes in first element in the set */
  int  NumElementQuadPoints (const StringT& name) const; /**< number of quad pts in first element in the set */
  void ElementSet (const StringT& name, iArrayT& elset) const; /**< access element set */
  void GeometryCode (const StringT& name, GeometryT::CodeT& code); /**< access geometry code */

  int  NumNodeVariables (void) const;
  int  NumElementVariables (void) const;
  int  NumQuadratureVariables (void) const;

  void ModeData (int index, int &number, double &mode) const; /**< mode data access */
  void TimeData (int index, int &number, double &time) const; /**< time step data access */

  void NodeVariables (iArrayT& keys, iArrayT& dims) const; /**< access node keys and dimensions  */
  void ElementVariables (iArrayT& keys, iArrayT& dims) const; /**< access elem keys and dimensions */
  void QuadratureVariables (iArrayT& keys, iArrayT& dims) const; /**< access quad keys and dimensions */

  /** returns dimensions of those variables used by this set name */
  void VariablesUsed (const StringT& name, AbaqusVariablesT::TypeT vt, iArrayT& used);

  /** read all variables of the type specified for the step specified for either the node or element set */
  void ReadVariables (AbaqusVariablesT::TypeT vt, int step, dArray2DT& values, const StringT& name);

  const char* VariableName (int index) const; /**< convert index to name */
  int VariableKey (const char *name) const; /**< convert name to key */
  int VariableKey (int index) const; /**< convert index to key */
  int VariableKeyIndex (int key) const; /**< convert key to index */

  bool NextCoordinate (int &number, dArrayT& nodes); /**< read next coordinate in file */
  bool NextElement (int &number, GeometryT::CodeT &type, iArrayT &nodes); /**< read next element in file */

  /** write connectivity set, number elements consecutively from startnumber */
  void WriteConnectivity (GeometryT::CodeT code, int startnumber, const iArray2DT& connects);

  /** write coordinates, number nodes using nodes_used */
  void WriteCoordinates (const iArrayT& nodes_used, const dArray2DT& coords);

  /** write element set */
  void WriteElementSet (const StringT& name, const iArrayT& elms);

  /** write node set */
  void WriteNodeSet (const StringT& name, const iArrayT& nodes);

  /** describe the active degrees of freedom */
  void WriteActiveDOF (const iArrayT& active);

  /** write a title */
  void WriteHeading (const StringT& heading);

  /** start a print increment */
  void WriteStartIncrement (int step, int inc, double totaltime, 
     double time, double timeincrement, AbaqusResultsT::AnalysisTypeT atype);

  /** for each variable for each time inc, define the setname, etc */
  void WriteOutputDefinition (int key, const StringT& setname, GeometryT::CodeT code, 
			      int numelemnodes);

  /** write variables that are nodal, whether they are element data averaged at the nodes
      or node point data */
  void WriteNodeVariables (int &index, const iArrayT& keys, const dArray2DT& values, 
			   const iArrayT& nodes_used, int numdir, int numshear);

  /** write element variables */
  void WriteElementVariables (int &index, const iArrayT& keys, const dArray2DT& values, 
			      const iArrayT& els_used, int numdir, int numshear);

  /** close the printed increment */
  void WriteEndIncrement (void);

  /** access version notes */
  void VersionNotes (ArrayT<StringT>& records);

  /** reset the input file to the start of the file */
  void ResetFile (void);
 
  int NumElements (void) const; /**< number of elements in file */
  int ElementNumber (int index) const; /**< element ID */
  int VariableDimension (int index) const; /**< number of variable components */
  int NodeNumber (int index) const; /**< node ID */
  int NumElementSets (void) const; /**< number of element sets */
  int NumNodeSets (void) const; /**< number of node sets */

 private:

  /** number of variables defined */
  enum NumVariables { NVT = 24 };

  /** Type of element output records */
  enum OutputType { kElementOutput = 0,
		    kNodalOutput = 1,
		    kModalOutput = 2,
		    kElemSetOutput = 3 };

  /** storage of element variable */
  enum ElementVarType { kElementQuadrature = 0,
			kElementCentroidal = 1,
			kElementNodal = 2,
			kElementRebar = 3,
			kElementNodeAveraged = 4,
			kElementWhole = 5};
  
  /** record keys */
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

  enum OutputParamsT { dprecision = 15, kDoubleSize = 8, kIntSize = 4};

  /** return status
      \note OKAY is item found, 
      \note BAD means cannot read or incorrect item found
      \note END means EOF found */
  enum StatusT { OKAY = -101, BAD = -102, END = -103 };

  /** read version record */
  bool ReadVersion (void);

  /** read next mode record found in file */
  bool NextMode (int &number, double &mode);

  /** read next time increment found in file */
  bool NextTimeSteps (int &number, double &time);

  /** scan the element for dimensions */
  void ScanElement (void);

  /** store element or node set */
  void StoreSet (iArrayT& set);

  /** read output definition record for output mode */
  void ReadOutputDefinitions (int &outputmode, StringT& setname);

  /** read element header for either element number or node number and location */
  void ReadElementHeader (int& objnum, int& intpt, int& secpt, int &location);

  /** scan variable record and save data about variables found */
  void ScanVariable (int key, int outputmode, int location);

  /** write element header record */
  void WriteElementHeader (int key, int number, int intpt, int secpt, 
			   AbaqusResultsT::ElementVarType flag, int numdirect, 
			   int numshear, int numdir, int numsecforc); 

  /** is the variable written with the node number preceeding the component values */
  bool VariableWrittenWithNodeNumber (int key) const;

  /** determines node or element list for ReadVariables
   * \return true if subset of all nodes or elements */
  bool DataPoints (AbaqusVariablesT::TypeT vt, const StringT& name, iArrayT& set) const;

  /** advances fIn to the time index given */
  void AdvanceToTimeIncrement (int step);

  /** determine variable type */
  bool CorrectType (int outputmode, int objnum, int intpt, int location, 
		    AbaqusVariablesT::TypeT vt, int& ID) const;

  int TranslateElementName (const char *, GeometryT::CodeT &, int &) const; /**< translate element name into geometry code */
  int TranslateContinuum (const char *, GeometryT::CodeT &, int &) const; /**< translate continuum elements */
  int TranslateSpring (const char *name, GeometryT::CodeT &type, int &numintpts) const; /**< translate spring elements */
  int Translate2D (const char *, GeometryT::CodeT &, int &) const; /** translate 2D continuum elements */
  int Translate3D (const char *, GeometryT::CodeT &, int &) const; /** translate 3D continuum elements */
  int TranslateShell (const char *, GeometryT::CodeT &, int &) const; /** translate shell elements to quad or tris */
  int TranslateRigid (const char *name, GeometryT::CodeT &type, int &numintpts) const; /**< translate rigid elements */
  /** generate element name and determine number of element nodes to write */
  void GetElementName (GeometryT::CodeT geometry_code, int elemnodes, 
		       int& num_output_nodes, StringT& elem_name) const;
  
  /** advance to the next record of type target in the file */
  bool AdvanceTo (int target);
  /** skip the rest of the record */
  bool SkipAttributes (void);
  
  /** read the start of the next record to figure out what it is
      \note END is returned if EOF is encountered */
  int ReadNextRecord (int &key);

  bool ReadSetName (StringT& s, int n); /**< read the string attribute and delete leading/trailing spaces */
  bool Read (StringT& s, int n); /**< read the string attribute, n blocks */
  bool Read (int& i); /**< read the integer attribute */
  bool Read (double& d); /**< read the double attribute */
  bool CheckBufferSize (istream& in, int numchars); /**< check the ASCII buffer size, read more if necessary */
  void CheckBufferSize (istream& in); /**< check the binary buffer size, read header/footer */

  void Write (int i); /**< write the integer attribute */
  void Write (double d); /**< write the double attribute */
  void Write (const StringT& s, int blocks = 1); /**< write the string attribute */
  void WriteASCII (const StringT& s); /**< write the string */
  void CheckBufferSize (ostream& out); /**< check binary buffer, write header/footer */

  void SetVariableNames (void); /**< set up the table of variable definitions */

 private:
  ifstreamT fIn; /**< input file stream */
//ifstream fIn; /**< input file stream */
  ofstream fOut; /**< output file stream */
  ostream& fMessage; /**< error message stream */
  StringT fFileName; /**< file name */
  const StringT fMarker; /**< ASCII record key marker */
  const StringT fOutVersion; /**< version of output file written */

  bool fBinary; /**< binary or ASCII file format */
  int fBufferDone; /**< amount of buffer read or written */
  int fBufferSize; /**< amount of buffer to read or write between header/footer or per line */
  StringT fBuffer; /**< the ASCII buffer */
  int fCurrentLength; /**< current length of record left to read */

  int fNumNodes; /**< num nodes in input file */
  int fNumElements; /**< num element in input file */

  int fNumNodeSets; /**< num node sets in input file */
  AutoArrayT<StringT> fNodeSetNames; /**< node set names in input file */
  AutoArrayT<iArrayT> fNodeSets; /**< node sets in input file */
  int fNumElementSets; /**< num element sets in input file */
  AutoArrayT<StringT> fElementSetNames; /**< element set names in input file */
  AutoArrayT<iArrayT> fElementSets; /**< element sets in input file */
  int fStartCount; /**< number of start increments in input file */
  int fEndCount; /**< number of end increments in input file */
  int fModalCount; /**< number of mode records in input file */

  ArrayT<AbaqusVariablesT> fVariableTable; /**< table of variable definitions */
  int fNumNodeVars; /**< number of node vars in input file */
  int fNumElemVars; /**< number of elem vars in input file */
  int fNumQuadVars; /**< number of quad vars in input file */

  iAutoArrayT fElementNumber; /**< element IDs */
  iAutoArrayT fNumElementQuadPoints; /**< num quad pts per element */
  iAutoArrayT fNodeNumber; /**< node IDs */

  iAutoArrayT fTimeIncs; /**< time incs */
  AutoArrayT<double> fTimeSteps; /**< time values */

  iAutoArrayT fModeIncs; /**< mode numbers */
  AutoArrayT<double> fModeSteps;  /**< mode values */
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

} // namespace Tahoe 
#endif
