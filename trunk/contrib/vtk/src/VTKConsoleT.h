/* $Id: VTKConsoleT.h,v 1.14 2001-10-23 00:23:08 recampb Exp $ */

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
  double valRange1, valRange2;
  double hueRange1,hueRange2;
  double satRange1, satRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1[100], scalarRange2[100];
  double time;
  double scale_factor;
  int numColors;
  int num_node_variables;
  int currentVarNum;
  int test;
  double numRen;
  StringT source_file;
  StringT output_file;
  StringT outFileName;
  StringT sbTitle;
  StringT varList;
  vtkUnstructuredGridReader *ugr;
  vtkRenderer *renderer;
  vtkRenderer *renderer2;
  vtkRenderer *renderer3;
  vtkRenderer *renderer4;
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
  vtkScalars *scalars [1000][100];
  vtkVectors *vectors [1000][100];
  int num_time_steps;
  vtkCamera *cam;
  vtkCubeAxesActor2D *axes;
  vtkPoints *points;
  ArrayT<StringT> node_labels;
  vtkUnstructuredGrid *warpGrid;
  vtkPoints *dPoints;
  vtkWarpVector *warp;
  int frameNum;

  ArrayT<VTKFrameT> fFrames;
};

#endif
