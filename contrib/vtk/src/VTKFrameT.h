/* $Id: VTKFrameT.h,v 1.7 2001-11-01 19:16:44 recampb Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "AutoArrayT.h"
#include "VTKBodyT.h"

/* VTK forward declarations */
class vtkRenderer;
class vtkRendererSource;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkTIFFWriter;

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
  
  vtkRenderer* getRen(void) {return renderer;};
  //private:
  
  void ResetView(void);

  /** return a pointer to the frame's renderer */
  vtkRenderer* Renderer(void) { return renderer; };

   virtual bool iDoCommand(const StringT& command, StringT& line);
   vtkRenderWindow* getRenWin(void) {return fRenWin;};
   vtkRenderWindowInteractor* getIren(void) {return fIren;};
   StringT getName(void) {return bodies[0]->inFile;};
   
  vtkRenderWindow *fRenWin;
  vtkRenderWindowInteractor *fIren;
  AutoArrayT<VTKBodyT*> bodies;
 private:
  
  vtkRenderer *renderer;
  vtkRendererSource *renSrc;
  vtkTIFFWriter *writer;

  


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

  //AutoArrayT<VTKBodyT*> pBodies;

  /** pointers to bodies displayed in the frame */


};

#endif
