/* $Id: TranslateIOManager.h,v 1.17 2003-11-18 19:36:27 paklein Exp $ */
#ifndef _TRANSLATE_IOMANAGER_H_
#define _TRANSLATE_IOMANAGER_H_

#include "ModelManagerT.h"
#include "OutputBaseT.h"
#include "StringT.h"
#include "ArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "ofstreamT.h"

namespace Tahoe {

class TranslateIOManager
{
 public:

	/** constructor */
	TranslateIOManager(ostream& message, istream& in, bool write);

	/** destructor */
	virtual ~TranslateIOManager(void);

	/** set up to echo input parameters */
	void SetEcho (int selection, const StringT& file);

	/** run program */
	virtual void Translate (const StringT& program, const StringT& version, const StringT& title);

 protected:
  void SetInput (void);
  virtual void SetOutput (const StringT& program, const StringT& version, const StringT& title);

  void InitializeVariables (void);
  void InitializeNodeVariables (void);
  void InitializeElementVariables(void);
  void InitializeQuadVariables (void);

  void InitializeElements (int& elementgroup, StringT& name);
  void InitializeNodePoints (iArrayT& nodes, iArrayT& index);
  void SelectElements(StringT& ID, iArrayT& elements, iArrayT& index);

  void InitializeTime (void);

  virtual void TranslateVariables (void);

  virtual void WriteGeometry (void);
  void WriteNodes (void);
  void WriteNodeSets (void);
  void WriteElements (void);
  void WriteSideSets (void);

  void VariableQuery (const ArrayT<StringT>& names, iArrayT& list);

	void ReNameLabels(const StringT& data_type, ArrayT<StringT>& labels);
	void FilterLabels(const StringT& data_type, iArrayT& filter, ArrayT<StringT>& labels);

 protected:
  ostream& fMessage;
  istream& fIn;
  bool fWrite;
  ModelManagerT fModel;

  int fOutputFormat;
  StringT fOutputName;
  int fCoords;  // flag used to determine if coords are treated as a variable

  // variable information
  int fNumNV;
  int fNumEV;
  int fNumQV;
  ArrayT<StringT> fNodeLabels;
  ArrayT<StringT> fElementLabels;
  ArrayT<StringT> fQuadratureLabels;
  iArrayT fNVUsed;
  iArrayT fEVUsed;
  iArrayT fQVUsed;
  
  iArrayT fKeepNodeLabels;
  iArrayT fKeepElementLabels;

  // time step information
  int fNumTS;
  dArrayT fTimeSteps;
  iArrayT fTimeIncs;

  OutputBaseT* fOutput;
  bool fEcho;
  ofstreamT fEchoOut;

 private:
  iArrayT fOutputID;
  bool fOneOutputSet;

  dArray2DT fCoordinates;
  iArrayT fNodeMap; // not stored in ModelManager
  iArrayT fNodeID;
  ArrayT<iArray2DT> fGlobalSideSets; 
  // need to store here instead of in ModelManager so that they are all global

};

} // namespace Tahoe

#endif
