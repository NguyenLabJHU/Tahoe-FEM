/* $Id: VTKFrameT.h,v 1.2 2001-10-24 00:52:08 recampb Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"


/* forward declarations */
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkPoints;
class vtkUnstructuredGrid;
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


class VTKFrameT: public StringT, public iConsoleObjectT
{
 public:

  /* constructor */
  VTKFrameT(void);

  //private:

  /* dummy variable */
  //int fAge;
  
  int num_nodes;
  int num_dim;
  StringT inFile;
  double valRange1, valRange2;
  double hueRange1,hueRange2;
  double satRange1, satRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1[100], scalarRange2[100];
  double scale_factor;
  int numColors;
  int num_node_variables;
  int num_time_steps;
  ArrayT<StringT> node_labels;
  int currentVarNum;
  StringT output_file;
  StringT outFileName;
  StringT sbTitle;
  StringT varList;
  vtkRenderer *renderer;
  vtkLookupTable *lut;
  vtkDataSetMapper *ugridMapper;
  vtkActor *ugridActor;
  vtkActor *wireActor;
  vtkScalarBarActor *scalarBar;
  vtkRendererSource *renSrc;
  vtkIdFilter *ids;
  vtkSelectVisiblePoints *visPts;
  vtkLabeledDataMapper *ldm;
  vtkActor2D *pointLabels;
  vtkUnstructuredGrid *ugrid;
  vtkScalars *scalars [1000][100];
  vtkVectors *vectors [1000][100];
  vtkCamera *cam;
  vtkCubeAxesActor2D *axes;
  vtkPoints *points;
  vtkWarpVector *warp;
  int frameNum;

};

#endif
