/* $Id: VTKFrameT.h,v 1.5 2001-10-26 02:14:53 paklein Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "AutoArrayT.h"

/* VTK forward declarations */
class vtkRenderer;
class vtkRendererSource;
class vtkRenderWindow;
class vtkRenderWindowInteractor;

/* forward declarations */
class VTKBodyT;

class VTKFrameT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKFrameT(void);

  /** destructor */
  ~VTKFrameT(void);

  //VTKBodyT Body(const int bodyNum) { return bodies(bodyNum);};
  
  /** return a pointer to the frame's renderer */
  vtkRenderer* Renderer(void) { return renderer; };

private:

  vtkRenderer *renderer;
  vtkRendererSource *renSrc;

  /* vtkScalarBarActor *scalarBar; */
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

  /** pointers to bodies displayed in the frame */
  AutoArrayT<VTKBodyT*> bodies;
};

#endif
