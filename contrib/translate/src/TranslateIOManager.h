/* $Id: TranslateIOManager.h,v 1.2 2001-09-07 13:25:29 sawimme Exp $ */

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
  void Translate (const StringT& program, const StringT& version, const StringT& title);

 private:
  void SetOutput (const StringT& program, const StringT& version, const StringT& title);

  void InitializeVariables (void);
  void InitializeTime (void);

  void TranslateVariables (bool xy);

  void WriteGeometry (void);
  void WriteNodes (void);
  void WriteNodeSets (void);
  void WriteElements (void);
  void WriteSideSets (void);

 private:
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

  int fNumTS;
  dArrayT fTimeSteps;

  iArrayT fNodeMap;

  ArrayT<iArray2DT> fGlobalSideSets;
};

#endif
