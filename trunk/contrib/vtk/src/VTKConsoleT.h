/* $Id: VTKConsoleT.h,v 1.17 2001-10-26 02:14:53 paklein Exp $ */

#ifndef _VTK_CONSOLE_T_H_
#define _VTK_CONSOLE_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "VTKFrameT.h"
#include "StringT.h"
#include "AutoArrayT.h"

/* VTK forward declarations */
class vtkTIFFWriter;

/* forward declarations */
class VTKBodyT;

class VTKConsoleT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKConsoleT(void);

  /** destructor */
  ~VTKConsoleT(void);

  /* execute given command - returns false on fail */
  virtual bool iDoCommand(const StringT& command, StringT& line);

 private:

  int test;
  double numRen;
  //StringT inFile;

  vtkRenderWindow *renWin;
  vtkRenderWindowInteractor *iren;
  //vtkRendererSource *renSrc;
 
  vtkTIFFWriter *writer;
 
  /** list of frame in the display window */
  ArrayT<VTKFrameT> fFrames;

  /** list of display objects */
  AutoArrayT<VTKBodyT*> fBodies;
};

#endif
