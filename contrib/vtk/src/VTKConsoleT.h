/* $Id: VTKConsoleT.h,v 1.16 2001-10-25 21:40:19 recampb Exp $ */

#ifndef _VTK_CONSOLE_T_H_
#define _VTK_CONSOLE_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "VTKFrameT.h"

/* forward declarations */
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkPoints;
class vtkUnstructuredGrid;
class vtkUnstructuredGridReader;
class vtkDataSetMapper;
class vtkActor;
class vtkScalarBarActor;
class vtkCubeAxesActor2D;
class vtkRendererSource;
class vtkTIFFWriter;
class vtkWindowToImageFilter;
class vtkLookupTable;
class vtkIdFilter;
class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;
class vtkScalars;
class vtkCamera;
class vtkCubeAxesActor2D;
class StringT;
class vtkWarpVector;
class vtkVectors;


class VTKConsoleT: public iConsoleObjectT
{
 public:

  /* constructor */
  VTKConsoleT(void);

  /* execute given command - returns false on fail */
  virtual bool iDoCommand(const StringT& command, StringT& line);

 private:


  int test;
  double numRen;
  // StringT source_file;
  StringT inFile;

  vtkRenderWindow *renWin;
  vtkRenderWindowInteractor *iren;
  vtkRendererSource *renSrc;
 
  vtkTIFFWriter *writer;
 

  ArrayT<VTKFrameT> fFrames;
};

#endif
