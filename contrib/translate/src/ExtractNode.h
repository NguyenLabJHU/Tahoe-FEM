#ifndef _EXTRACT_NODE_H_
#define _EXTRACT_NODE_H_

#include "ExtractIOManager.h"

class ExtractNode : public ExtractIOManager
{
 public:
  ExtractNode (ostream& message, istream& in, bool write);
  
 protected:
  void Initialize (void);
  void TranslateVariables (void);
};

#endif
