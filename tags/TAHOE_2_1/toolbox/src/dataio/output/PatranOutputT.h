/* created: sawimme April 2002 */

#ifndef _PATRANOUTPUT_T_H_
#define _PATRANOUTPUT_T_H_

/* direct members */
#include "OutputBaseT.h"
#include "PatranT.h"

/* forward declarations */

namespace Tahoe {

class PatranOutputT : public OutputBaseT
{
 public:
  PatranOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary);
  void WriteGeometry (void);
  void WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values);

 private:
  void FileName (int ID, StringT& filename, const char* ext) const;
  void WriteConnectivity (ostream& patout, int& firstID, int ID, iArrayT& nodes_used) const;
  void WriteNamedComponents (ostream& patout, int& firstID, int ID) const;
  PatranT::NamedTypes GetPatranNamedType (GeometryT::CodeT geom) const;
  PatranT::ElementTypes GetPatranElementType (GeometryT::CodeT geom) const;

 private:
  bool fBinary;
  PatranT fPatran;
};

} // namespace Tahoe 
#endif
