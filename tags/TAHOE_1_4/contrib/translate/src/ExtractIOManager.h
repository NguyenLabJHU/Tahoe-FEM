
#ifndef _EXTRACT_IOMANAGER_H_
#define _EXTRACT_IOMANAGER_H_

#include "TranslateIOManager.h"
#include "ofstreamT.h"

namespace Tahoe {

class ExtractIOManager : public TranslateIOManager
{
 public:
  ExtractIOManager (ostream& message, istream& in, bool write);
  void Translate (const StringT& program, const StringT& version, const StringT& title);

 protected:
  virtual void Initialize (void) = 0;
  virtual void TranslateVariables (void) = 0;

  void SetOutput (const StringT& program, const StringT& version, const StringT& title);
  void PrepFiles (iArrayT& varsused, ArrayT<StringT>& labels);
  void WriteVarData (iArrayT& varsused, int ts) const;

 private:
  void OpenFile (ofstreamT& o, const StringT& name, bool append) const;

 protected:
  int fNumItems;
  ArrayT<StringT> fItemNames;
  iArrayT fItemIndex;

  dArray2DT fVarData;
  dArray2DT fCoordinates;

 private:
  int fCheck;
  int fNumDigits;
  StringT fOutfileExtension;
};

} // namespace Tahoe

#endif
