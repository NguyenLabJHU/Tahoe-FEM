
#ifndef _EXTRACT_QUAD_H_
#define _EXTRACT_QUAD_H_

#include "ExtractIOManager.h"
#include "ofstreamT.h"

namespace Tahoe {

class ExtractQuad : public ExtractIOManager
{
 public:
  ExtractQuad (ostream& message, istream& in, bool write);

 protected:
  void Initialize (void);
  void TranslateVariables (void);

 private:
  StringT fElementName;
  int fElementGroup;
};

} //namespace Tahoe

#endif
