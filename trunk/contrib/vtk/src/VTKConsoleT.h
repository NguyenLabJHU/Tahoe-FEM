/* $Id: VTKConsoleT.h,v 1.7 2001-10-02 18:40:31 recampb Exp $ */

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
class vtkIdFilter;
class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;
class vtkScalars;


class VTKConsoleT: public iConsoleObjectT
{
 public:

  /* constructor */
  VTKConsoleT(void);

  /* execute given command - returns false on fail */
  virtual bool iDoCommand(const StringT& command, StringT& line);

 private:
  double valRange1, valRange2;
  double hueRange1,hueRange2;
  double satRange1, satRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1, scalarRange2;
  int numColors;
  StringT source_file;
  StringT output_file;
  StringT outFileName;
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
  vtkIdFilter *ids;
  vtkSelectVisiblePoints *visPts;
  vtkLabeledDataMapper *ldm;
  vtkActor2D *pointLabels;
  vtkUnstructuredGrid *ugrid;  
  /* vtkScalars *scalars[num_time_steps]; */
  vtkScalars *scalars[10000];
  int num_time_steps;

};

#endif
