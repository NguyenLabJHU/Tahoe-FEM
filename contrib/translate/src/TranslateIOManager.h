/* $Id: TranslateIOManager.h,v 1.3 2001-09-10 16:47:22 sawimme Exp $ */

#ifndef _TRANSLATE_IOMANAGER_H_
#define _TRANSLATE_IOMANAGER_H_

#include "DataManagerT.h"
#include "OutputBaseT.h"
#include "StringT.h"
#include "ArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"

class TranslateIOManager
{
 public:

  TranslateIOManager (ostream& message);
  virtual void Translate (const StringT& program, const StringT& version, const StringT& title);

 protected:
  virtual void SetOutput (const StringT& program, const StringT& version, const StringT& title);

  void InitializeVariables (void);
  void InitializeTime (void);

  virtual void TranslateVariables (void);

  virtual void WriteGeometry (void);
  void WriteNodes (void);
  void WriteNodeSets (void);
  void WriteElements (void);
  void WriteSideSets (void);

 protected:
  ostream& fMessage;

  DataManagerT fModel;

  OutputBaseT* fOutput;
  StringT fOutputName;
  iArrayT fOutputID;

  int fNumNV;
  int fNumEV;
  int fNumQV;
  ArrayT<StringT> fNodeLabels;
  ArrayT<StringT> fElementLabels;
  ArrayT<StringT> fQuadratureLabels;

  int fNumTS;
  dArrayT fTimeSteps;
  iArrayT fTimeIncs;

  iArrayT fNodeMap; // not stored in ModelManager
  ArrayT<iArray2DT> fGlobalSideSets; 
  // need to store here instead of in ModelManager so that they are all global
};

#endif
