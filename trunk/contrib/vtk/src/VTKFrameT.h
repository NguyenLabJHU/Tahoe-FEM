/* $Id: VTKFrameT.h,v 1.4 2001-10-25 21:40:19 recampb Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "VTKBodyT.h"

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
class VTKBodyT;

class VTKFrameT: public iConsoleObjectT
{
 public:

  /* constructor */
  VTKFrameT(void);
  //VTKBodyT Body(const int bodyNum) { return bodies(bodyNum);};
  
  //vtkRenderer* Renderer(void) {return renderer;};
  //private:
 
  vtkRenderer *renderer;
  /* vtkScalarBarActor *scalarBar; */
  vtkRendererSource *renSrc;
 /*  vtkIdFilter *ids; */
/*   vtkSelectVisiblePoints *visPts; */
/*   vtkLabeledDataMapper *ldm; */
/*   vtkActor2D *pointLabels; */
/*   vtkUnstructuredGrid *ugrid; */
/*   vtkScalars *scalars [1000][100]; */
/*   vtkVectors *vectors [1000][100]; */
/*   vtkCamera *cam; */
/*   vtkCubeAxesActor2D *axes; */
/*   vtkPoints *points; */
/*   vtkWarpVector *warp; */
/*   int frameNum; */
  ArrayT<VTKBodyT> bodies;
};

#endif
