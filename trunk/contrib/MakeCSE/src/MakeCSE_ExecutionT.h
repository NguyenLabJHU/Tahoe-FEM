#ifndef _MAKECSE_EXECUTION_T_H_
#define _MAKECSE_EXECUTION_T_H_

#include "StringT.h"

namespace Tahoe {

class ifstreamT;

class MakeCSE_ExecutionT
{
 public:
  MakeCSE_ExecutionT (void); /**< constructor */
  void Run (void); /**< run program */
 private:
  void RunBatchOrJob (ifstreamT& in); /**< batch or job */
  void RunJob (ifstreamT& in); /**< run individual job */
 private:
  bool fInteractive;
  const StringT fProgram;
  const StringT fVersion;
};
} // namespace

#endif
