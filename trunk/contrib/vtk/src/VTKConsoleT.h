/* $Id: VTKConsoleT.h,v 1.15 2001-10-24 00:52:08 recampb Exp $ */

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

  double time;

  int test;
  double numRen;
  // StringT source_file;

  vtkRenderWindow *renWin;
  vtkRenderWindowInteractor *iren;
 
  vtkTIFFWriter *writer;
   
  /* vtkScalars *scalars[num_time_steps]; */
 

  ArrayT<VTKFrameT> fFrames;
};

#endif
