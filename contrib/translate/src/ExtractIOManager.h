
#ifndef _EXTRACT_IOMANAGER_H_
#define _EXTRACT_IOMANAGER_H_

#include "TranslateIOManager.h"

class ExtractIOManager : public TranslateIOManager
{
 public:
  ExtractIOManager (ostream& message);
  virtual void Translate (const StringT& program, const StringT& version, const StringT& title);

 private:
  virtual void SetOutput (const StringT& program, const StringT& version, const StringT& title);
  virtual void InitializeVariables (void);
  virtual void InitializeNodePoints (void);
  virtual void TranslateVariables (void);

 private:
  int fOutputFormat;
  ofstream fOutFile;

  int fNumNP;
  iArrayT fNodePoints;
  iArrayT fNodePointIndex;
};

#endif
