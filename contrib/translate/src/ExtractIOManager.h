
#ifndef _EXTRACT_IOMANAGER_H_
#define _EXTRACT_IOMANAGER_H_

#include "TranslateIOManager.h"
#include "ofstreamT.h"

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

  void PrepFiles (StringT& ext, int digits) const;
  void OpenFile (ofstreamT& o, int index, int digits, StringT& ext, bool append) const;

 private:
  int fOutputFormat;

  int fCoords;
  int fNumNP;
  iArrayT fNodePoints;
  iArrayT fNodePointIndex;
};

#endif
