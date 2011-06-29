//
// $Id: MixedShapeFunctionT.h,v 1.1 2008-12-12 00:40:26 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

#ifndef _MIXED_SHAPE_FUNCTION_T_H_
#define _MIXED_SHAPE_FUNCTION_T_H_

#include <cassert>
#include "ShapeFunctionT.h"

namespace Tahoe {

  class MixedShapeFunctionT: public ShapeFunctionT {
  public:
    MixedShapeFunctionT(ShapeFunctionT& shapeFunction,
        GeometryT::CodeT geometry_code, int numIP, const LocalArrayT& coords,
        bool is_open_set = false);

    virtual ~MixedShapeFunctionT();

    // compute jacobians of the nodal values at integration points
    virtual void JacobianDets(dArrayT& dets);

    virtual void SetDerivatives(void);

  private:

    // Reference to primary interpolation functions.
    ShapeFunctionT& fShapeFunction;

  };

}

#include "MixedShapeFunctionT.i.h"

#endif // _MIXED_SHAPE_FUNCTION_T_H_
