/* $Id: VTKConsoleT.h,v 1.3 2001-09-26 18:07:58 recampb Exp $ */

#ifndef _VTK_CONSOLE_T_H_
#define _VTK_CONSOLE_T_H_

/* base class */
#include "iConsoleObjectT.h"

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



class VTKConsoleT: public iConsoleObjectT
{
 public:

  /* constructor */
  VTKConsoleT(void);

  /* execute given command - returns false on fail */
  virtual bool iDoCommand(const StringT& command, StringT& line);

  /* void dataUpdate(void); */

 private:
  double valRange1, valRange2;
  double hueRange1,hueRange2;
  double satRange1, satRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1, scalarRange2;
  int numColors;
  StringT source_file;
  StringT output_file;
  vtkUnstructuredGridReader *ugr;
  vtkRenderer *renderer;
  vtkRenderWindow *renWin;
  vtkRenderWindowInteractor *iren;
  vtkLookupTable *lut;
  vtkDataSetMapper *ugridMapper;
  vtkActor *ugridActor;
  vtkActor *wireActor;
  vtkScalarBarActor *scalarBar;
  vtkRendererSource *renSrc;
  vtkTIFFWriter *writer;
  

};

#endif
