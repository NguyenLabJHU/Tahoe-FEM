// file: TecPlotT.h

// created: saw (06.06.2000)

#ifndef _TECPLOT_T_H_
#define _TECPLOT_T_H_

/* direct members */
#include "iosfwd.h"
#include "GeometryT.h"
#include "ArrayT.h"

/* forward declarations */
class StringT;
class dArray2DT;
class iArray2DT;
class iArrayT;

class TecPlotT
{
 public:
  TecPlotT (ostream& out, bool point);

  void WriteHeader (ostream& out, const StringT& title, const ArrayT<StringT>& variablenames) const;

  void WriteIJKZone (ostream& out, const StringT& title, const iArrayT& ijk) const;
  
  // must use WriteConnecitivity with this
  void WriteFEZone (ostream& out, const StringT& title, int numnodes, int numelems, GeometryT::GeometryCode code, bool connectivity) const;

  // write data can only be call once if using point format
  // but may be called repeatly, in proper order, for block format
  void WriteData (ostream& out, const dArray2DT& data) const;

  // only used with WriteFEZone
  void WriteConnectivity (ostream& out, GeometryT::GeometryCode code, const iArray2DT& connects) const;

 private:
  ostream& fOut;
  bool fPoint; // data can be written in point or block format
};

#endif
