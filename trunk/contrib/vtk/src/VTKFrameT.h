/* $Id: VTKFrameT.h,v 1.8 2001-11-06 02:39:51 recampb Exp $ */

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
//class vtkRendererSource;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkTIFFWriter;

/* forward declarations */
class VTKBodyT;
class VTKConsoleT;

class VTKFrameT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKFrameT(void);

  /** destructor */
  ~VTKFrameT(void);

  /** return a pointer to the specified frame body */
  VTKBodyT* Body(int bodyNum) { return bodies[bodyNum]; };
  
  /** add body to the frame.
   * returns true if the body was added the frame, false
   * otherwise */
  bool AddBody(VTKBodyT* body);

  /** delete body from the frame. returns true if body was found and
   * and removed, false otherwise. */
  bool RemoveBody(VTKBodyT* body);

  vtkRenderer* getRen(void) {return renderer;};
  //private:
  
  void ResetView(void);

  /** return a pointer to the frame's renderer */
  vtkRenderer* Renderer(void) { return renderer; };

   virtual bool iDoCommand(const StringT& command, StringT& line);
   
   /** set controlling console object */
   void setConsole(VTKConsoleT* console) { fConsole = console; };

   /** set the renderer window */
   void setRenWin(vtkRenderWindow* renWin) { fRenWin = renWin; };
   
   /** set the window interactor */
   void setIren(vtkRenderWindowInteractor* iren) {fIren = iren; };

   StringT getName(void) {return bodies[0]->inFile;};
   
 private:

   /** controlling console object */
   VTKConsoleT* fConsole;
  
  vtkRenderer *renderer;
  //  vtkRendererSource *renSrc;
  vtkTIFFWriter *writer;

  vtkRenderWindow *fRenWin;
  vtkRenderWindowInteractor *fIren;
  AutoArrayT<VTKBodyT*> bodies;
  

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
