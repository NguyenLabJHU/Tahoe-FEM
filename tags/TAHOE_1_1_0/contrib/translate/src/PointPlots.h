
#ifndef _POINT_PLOTS_H_
#define _POINT_PLOTS_H_

#include "TranslateIOManager.h"
#include "ofstreamT.h"

namespace Tahoe {

class PointPlots : public TranslateIOManager
{
 public:
  PointPlots (ostream& message, istream& in, bool write);
  virtual void Translate (const StringT& program, const StringT& version, const StringT& title);

 private:
  virtual void SetOutput (const StringT& program, const StringT& version, const StringT& title);
  virtual void TranslateVariables (void);
  void OpenFile (ofstreamT& o, int index, int digits, StringT& ext) const;

 private:
  int fElementGroup;
};

} // namespace Tahoe

#endif
