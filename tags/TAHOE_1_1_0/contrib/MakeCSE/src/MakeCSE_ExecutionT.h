#ifndef _MAKECSE_EXECUTION_T_H_
#define _MAKECSE_EXECUTION_T_H_

#include "StringT.h"

namespace Tahoe {

class ifstreamT;
class sArrayT;

class MakeCSE_ExecutionT
{
 public:
  MakeCSE_ExecutionT (void); /**< constructor */
  void Run (const sArrayT& lineoptions); /**< run program */
 private:
  void RunInteractive (void); /**< prompt user for batch/job file */
  void RunBatchOrJob (ifstreamT& in); /**< batch or job */
  void RunJob (ifstreamT& in); /**< run individual job */
 private:
  bool fInteractive;
  const StringT fProgram;
  const StringT fVersion;
};
} // namespace

#endif
