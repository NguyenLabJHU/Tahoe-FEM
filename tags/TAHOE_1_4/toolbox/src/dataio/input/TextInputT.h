/* $Id: TextInputT.h,v 1.2 2003-11-10 22:14:22 cjkimme Exp $ */
#ifndef _TEXT_INPUT_T_H_
#define _TEXT_INPUT_T_H_

/* base classes */
#include "InputBaseT.h"

/* direct members */
#include "StringT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArrayT;

class TextInputT : public InputBaseT
{
public:
  TextInputT (ostream& out);

  virtual bool Open (const StringT& filename);
  virtual void Close (void);

  /* return names, Array must be preallocated */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const;
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const;
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const;

  /* return dimenesions */
  virtual int  NumElementGroups (void) const;
  virtual int  NumSideSets (void) const;
  virtual int  NumNodeSets (void) const;

  /* NODES */
  virtual int  NumNodes (void) const;
  virtual int  NumDimensions (void) const;

	/** ids for all nodes. ids for all nodes in the coordinate list. ids may not be
	 * compact or ordered */
	virtual void ReadNodeID(iArrayT& node_id); 

	/** read coordinates */
	virtual void ReadCoordinates(dArray2DT& coords);
	
	/** read coordinates and ids. ids for all nodes in the coordinate list. 
	 * ids may not be compact or ordered */
	virtual void ReadCoordinates(dArray2DT& coords, iArrayT& node_id);

  /* ELEMENTS */
  virtual int  NumGlobalElements (void) const;
  virtual int  NumElements (const StringT& name);
  virtual int  NumElementNodes (const StringT& name);
  virtual int  NumElementQuadPoints (const StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (const StringT& name, iArrayT& elemmap);
  virtual void ReadGlobalElementSet (const StringT& name, iArrayT& set);
  virtual void ReadConnectivity (const StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (const StringT& name);
  virtual void ReadNodeSet (const StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void) const;
  virtual int  NumSidesInSet (const StringT& setname) const;
  virtual StringT SideSetGroupName (const StringT& setname) const;
  virtual void ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const;
  
  /* record[0] = progname, record[1] = version, record[2] = date, record[3] = time */
  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  /* for all nodes or elements */
  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;

  virtual void NodeVariablesUsed (const StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (const StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (const StringT& name, iArrayT& used);  

  virtual void ReadNodeLabels (ArrayT<StringT>& labels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;

  /* step starts at zero and increases by one */
  virtual void ReadAllNodeVariable (int step, int varindex, dArrayT& values);
  virtual void ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues);

  virtual void ReadAllElementVariable (int step, int varindex, dArrayT& values);
  virtual void ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadAllElementVariables (int step, dArray2DT& evalues);
  virtual void ReadElementVariables (int step, const StringT& name, dArray2DT& evalues);

  virtual void ReadAllQuadratureVariable (int step, int varindex, dArrayT& values);
  virtual void ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  virtual void ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues);

 private:
  bool OpenFile (ifstreamT& in, const char *ext) const;
  bool ScanGeometryFile (ifstreamT& in);
  bool ScanResultsFile (ifstreamT& in);
  bool ScanResultsFile_old (ifstreamT& in);
  bool AdvanceToBlock (ifstreamT& in, const StringT& name, const char *tname) const;
  void DataBlock (ifstreamT& in, iArrayT& used, iArrayT& ids, dArray2DT& vals, bool nodal) const;

	/** \name backward compatibility */
	/*@{*/
	/** return true if results file uses the pre-TOC format */
	bool is_old_format(const StringT& file) const;
 	void ReadNodeVariables_old(int step, const StringT& name, dArray2DT& nvalues);
	void ReadAllNodeVariables_old(int step, dArray2DT& nvalues);
	void ReadAllElementVariables_old(int step, dArray2DT& evalues);
	void ReadElementVariables_old(int step, const StringT& name, dArray2DT& evalues);
	/*@}*/

	/** return the results file name for the given output step */
	void ResultsFile(const StringT& toc_file, int step, StringT& file) const;

 private:

	StringT fFileRoot;
	StringT fFilePath;

  AutoArrayT<StringT> fBlockID;
  AutoArrayT<int> fBlockNumElem;
  AutoArrayT<int> fBlockNumElemNode;
  AutoArrayT<GeometryT::CodeT> fBlockGeometry;
  int fNumNodes;
  int fNumElements;
  int fNumDOF;

  AutoArrayT<double> fTimeSteps;
  AutoArrayT<StringT> fNodeVariable;
  AutoArrayT<StringT> fElementVariable;
};

inline void TextInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{ sidenames.Free(); }
inline void TextInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{ nodenames.Free(); }
inline int TextInputT::NumElementGroups (void) const
{ return fBlockID.Length(); }
inline int TextInputT::NumSideSets (void) const
{ return 0; }
inline int TextInputT::NumNodeSets (void) const
{ return 0; }
inline int TextInputT::NumNodes (void) const
{ return fNumNodes; }
inline int TextInputT::NumDimensions (void) const
{ return fNumDOF; }
inline int TextInputT::NumElementQuadPoints(const StringT& name)
{ 
#ifdef __MWERKS__
#pragma unused(name)
#endif
	return 0; 
}
inline int TextInputT::NumGlobalElements (void) const
{ return fNumElements; }
inline int TextInputT::NumNodesInSet (const StringT& name)
{
#ifdef __MWERKS__
#pragma unused (name)
#endif
  return 0;
}
inline bool TextInputT::AreSideSetsLocal (void) const
{ return true; }
inline int TextInputT::NumSidesInSet (const StringT& setname) const
{ 
#ifdef __MWERKS__
#pragma unused (setname)
#endif
  return 0; 
}

inline StringT TextInputT::SideSetGroupName (const StringT& setname) const
{ 
#ifdef __MWERKS__
#pragma unused (setname)
#endif
  StringT s = "";
  return s;
}

inline void TextInputT::QARecords (ArrayT<StringT>& records)
{
#ifdef __MWERKS__
#pragma unused (records)
#endif
}
inline int TextInputT::NumTimeSteps (void) const { return fTimeSteps.Length(); }
inline int TextInputT::NumNodeVariables (void) const
{ return fNodeVariable.Length(); }
inline int TextInputT::NumElementVariables (void) const
{ return fElementVariable.Length(); }
inline int TextInputT::NumQuadratureVariables (void) const
{ return 0; }


} // namespace Tahoe 
#endif
