// file: CSEBaseT.h

// created: SAW 5/2/2000

#ifndef _MakeCSE_CSEBASET_H_
#define _MakeCSE_CSEBASET_H_

#include "MakeCSE_ElementBaseT.h"

namespace Tahoe {

class MakeCSE_CSEBaseT : public MakeCSE_ElementBaseT
{
 public:
  MakeCSE_CSEBaseT (ostream& fMainOut, const StringT& ID);

  virtual void Initialize (GeometryT::CodeT code, int numregfacenodes);

  virtual void SetNodes (int e1local, const iArrayT& regelemnodes);

  virtual void CSElemFaces (iArrayT& faces) const;

 private:
  void CSEType (GeometryT::CodeT code, int numfacenodes);
  void Set3DNodes (int linearfacenodes);

 private:
  iArrayT fSurfaceFacets;
  iArrayT fNodeNumberingMap;
};

inline void MakeCSE_CSEBaseT::CSElemFaces (iArrayT& faces) const { faces = fSurfaceFacets; }

} //namespace Tahoe
#endif
