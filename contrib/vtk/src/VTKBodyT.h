/* $Id: VTKBodyT.h,v 1.10 2001-12-08 00:17:19 recampb Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"
#include "VTKBodyDataT.h"


/* forward declarations */

class VTKBodyDataT;

class VTKBodyT: public iConsoleObjectT
{
 public:

  /** default constuctor */
  VTKBodyT(void) { body = NULL; };

  /** constructor */
  VTKBodyT(VTKBodyDataT* body_data);
  
  /** return pointer to the body data */
  VTKBodyDataT* BodyData(void) { return body; };

  /** comparison operator */
  bool operator==(const VTKBodyT& rhs) { return body == rhs.body; };

  /** rvalue - smart pointer */
  VTKBodyDataT* operator->(); //CW wouldn't call functions with conversion

 private:

  VTKBodyDataT* body; 
};

inline VTKBodyDataT* VTKBodyT::operator->()
{
  return body;
}


#endif
