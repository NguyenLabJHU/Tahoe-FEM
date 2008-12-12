//
// $Id: MixedShapeFunctionT.i.h,v 1.1 2008-12-12 00:40:26 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//
namespace Tahoe {

  // compute jacobians of the nodal values at integration points
  inline void MixedShapeFunctionT::JacobianDets(dArrayT& dets)
  {
    ShapeFunctionT::JacobianDets(dets);
  }

} // namespace Tahoe
