/* $Id: TranslateIOManager.h,v 1.6 2002-02-13 18:08:32 sawimme Exp $ */

#ifndef _TRANSLATE_IOMANAGER_H_
#define _TRANSLATE_IOMANAGER_H_

#include "ModelManagerT.h"
#include "OutputBaseT.h"
#include "StringT.h"
#include "ArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"

class TranslateIOManager
{
 public:

  TranslateIOManager (ostream& message, istream& in, bool write);
  virtual void Translate (const StringT& program, const StringT& version, const StringT& title);

 protected:
  void SetInput (void);
  virtual void SetOutput (const StringT& program, const StringT& version, const StringT& title);

  void InitializeVariables (void);
  void InitializeNodeVariables (void);
  void InitializeQuadVariables (void);

  void InitializeElements (int& elementgroup, StringT& name) const;
  void InitializeNodePoints (iArrayT& nodes, iArrayT& index);

  void InitializeTime (void);

  virtual void TranslateVariables (void);

  virtual void WriteGeometry (void);
  void WriteNodes (void);
  void WriteNodeSets (void);
  void WriteElements (void);
  void WriteSideSets (void);

  void VariableQuery (const ArrayT<StringT>& names, iArrayT& list);

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

  // time step information
  int fNumTS;
  dArrayT fTimeSteps;
  iArrayT fTimeIncs;

 private:
  OutputBaseT* fOutput;
  iArrayT fOutputID;

  iArrayT fNodeMap; // not stored in ModelManager
  ArrayT<iArray2DT> fGlobalSideSets; 
  // need to store here instead of in ModelManager so that they are all global

};

#endif
