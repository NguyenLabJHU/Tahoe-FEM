/* $Id: VTKBodyT.cpp,v 1.13 2001-12-08 00:17:19 recampb Exp $ */

#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

#include "dArrayT.h"


/* array behavior */
const bool ArrayT<VTKBodyT>::fByteCopy = true;

/* constructor */
VTKBodyT::VTKBodyT(VTKBodyDataT* body_data)
{
  body = body_data;
}

#if 0
/* conversion */
VTKBodyT::operator VTKBodyDataT()
{
  if (!body) throw eGeneralFail;
  return *body;
}
#endif
