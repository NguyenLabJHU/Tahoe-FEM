/* created: saw (05.11.2001)                                              */

#ifndef _AVS_T_H_
#define _AVS_T_H_

/* direct members */
#include "GeometryT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class StringT;
class iArray2DT;
class dArray2DT;

class AVST
{
 public:
  AVST (ostream& log, bool binary);
  void WriteHeader (ostream& out, int numpoints, int numcells, int numnodedata, int numcelldata, int nummodeldata) const;
  void WriteCoordinates (ostream& out, dArray2DT& coords, int firstID) const;
  void WriteCells (ostream& out, GeometryT::CodeT code, iArray2DT& connects, int matid, int firstID) const;
  void WriteDataHeader (ostream& out, const ArrayT<StringT>& labels) const;
  void WriteDataHeader (ostream &out, const ArrayT<StringT>& labels, const ArrayT<int>& dimensions) const;
  void WriteData (ostream& out, const dArray2DT& data, int firstID) const;

  /* future 
     void WriteCoordinates (ostream& out, dArray2DT& coords, iArrayT& nodes_used) const;
     void WriteCells (ostream& out, GeometryT::CodeT code, iArray2DT& connects, int matid, iArrayT& elemIDs) const;
     void WriteData (ostream& out, dArray2DT& data, iArrayT& nodes_used) const;
  */

 private:
  enum WidthsT { iwidth = 10, dwidth = 14, dprecision = 5};
  void GetElementName (GeometryT::CodeT code, StringT& name, int &numnodes) const;
  void WriteForward (ostream& out, int *pc, int numnodes) const;
  void WriteBackward (ostream& out, int *pc, int numnodes, int numfacenodes) const;
  void WriteArray2DT (ostream& out, const dArray2DT& data, int firstID) const;

 private:
  ostream &fOut;
  const bool fBinary;
};

} // namespace Tahoe 
#endif
