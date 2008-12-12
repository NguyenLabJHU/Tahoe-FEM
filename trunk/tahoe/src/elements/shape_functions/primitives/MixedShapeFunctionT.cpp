//
// $Id: MixedShapeFunctionT.cpp,v 1.1 2008-12-12 00:40:26 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

#include "MixedShapeFunctionT.h"

namespace Tahoe {

  MixedShapeFunctionT::MixedShapeFunctionT(ShapeFunctionT& shapeFunction,
      GeometryT::CodeT geometry_code, int numIP, const LocalArrayT& coords,
      bool is_open_set) :
    ShapeFunctionT(geometry_code, numIP, coords, is_open_set),
    fShapeFunction(shapeFunction)
  {

  }

  MixedShapeFunctionT::~MixedShapeFunctionT()
  {
    //
    // TODO Auto-generated destructor stub
    //

  }

  //
  // compute local shape functions and derivatives
  //
  void MixedShapeFunctionT::SetDerivatives(void)
  {

    const int nmn = fCoords.NumberOfNodes();

    if (nmn != 1) {
      ShapeFunctionT::SetDerivatives();
    } else {

      // Compute primary interpolation function Jacobian determinants
      dArrayT dets;
      fShapeFunction.JacobianDets(dets);

      // Average
      const int nmi = NumIP();
      const int nsd = NumSD();
      const int nip = dets.Length();

      dArrayT mixedDets(nmn);
      mixedDets[0] = 0.0;

      for (int i = 0; i < nip; ++i) {
        mixedDets[0] += dets[i];
      }
      mixedDets[0] /= nip;

      // Store
      fDet = mixedDets;
      fDomain->SetStoredJacobianDets(mixedDets);

    }

  }

} // namespace Tahoe
