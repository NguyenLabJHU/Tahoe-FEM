/* created: sawimme (05/11/2001) */

#ifndef _AVSOUTPUT_T_H_
#define _AVSOUTPUT_T_H_

/* direct members */
#include "OutputBaseT.h"

namespace Tahoe {

/* forward declarations */
class AVST;

class AVSOutputT : public OutputBaseT
{
 public:
  AVSOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary);
  void WriteGeometry (void);
  void WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

 private:
  StringT CreateFileName (int ID) const;
  void CountVariables (int &num_nvars, const ArrayT<StringT>& labels) const;
  void WriteCoordinates (ostream &avsout, AVST &avs, int index, iArrayT &nodesused) const;
  void WriteConnectivity (ostream &avsout, AVST &avs, int index, iArrayT &nodesused) const;
  void WriteVariable (ostream &avsout, AVST &avs, const ArrayT<StringT>& labels, const dArray2DT& values, int num_vars) const;
  bool IsVector (const ArrayT<StringT>& inlabels, int index, ArrayT<StringT>& extension, int dof) const;

 private:
  bool fBinary;
};
 
} // namespace Tahoe 
#endif
