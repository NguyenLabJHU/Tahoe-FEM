/* created: sawimme April 2002 */

#ifndef _PATRANOUTPUT_T_H_
#define _PATRANOUTPUT_T_H_

/* direct members */
#include "OutputBaseT.h"

/* forward declarations */


class PatranOutputT : public OutputBaseT
{
 public:
  PatranOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary);
  void WriteGeometry (void);
  void WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

 private:
  void FileName (int ID, StringT& filename, const char* ext) const;
  int GetPatranElementType (GeometryT::CodeT geom) const;

 private:
  bool fBinary;
};

#endif
