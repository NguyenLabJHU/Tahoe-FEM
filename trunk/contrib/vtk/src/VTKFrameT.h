/* $Id: VTKFrameT.h,v 1.14 2001-11-29 21:22:43 recampb Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "AutoArrayT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

/* VTK forward declarations */
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkActor2D;
class vtkSelectVisisblePoints;
class vtkLabeledDataMapper;
class vtkCubeAxesActor2D;

/* forward declarations */
class VTKBodyT;
class VTKBodyDataT;
class VTKConsoleT;

class VTKFrameT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKFrameT(void);

  /** destructor */
  ~VTKFrameT(void);

  /** return a pointer to the specified frame body */
  VTKBodyT* Body(int bodyNum) { return &bodies[bodyNum]; };
  
  /** add body to the frame.
   * returns true if the body was added the frame, false
   * otherwise */
  bool AddBody(VTKBodyDataT* body);

  /** delete body from the frame. returns true if body was found and
   * and removed, false otherwise. */
  bool RemoveBody(VTKBodyDataT* body);

  vtkRenderer* getRen(void) {return renderer;};
  //private:
  
  void ResetView(void);

  void ShowFrameNum(StringT);

  /** return a pointer to the frame's renderer */
  vtkRenderer* Renderer(void) { return renderer; };

   virtual bool iDoCommand(const CommandSpecT& command, StringT& line);
   
   /** set controlling console object */
   void setConsole(VTKConsoleT* console) { fConsole = console; };

   /** set the renderer window */
   void setRenWin(vtkRenderWindow* renWin) { fRenWin = renWin; };
   
   /** set the window interactor */
   void setIren(vtkRenderWindowInteractor* iren) {fIren = iren; };

   StringT getName(void) {return bodies[0]->inFile;};

 protected:

   /** write prompt for the specific argument of the command */
   virtual void ValuePrompt(const CommandSpecT& command, int index, ostream& out) const;  

 private:

   /** controlling console object */
   VTKConsoleT* fConsole;
  
  vtkRenderer *renderer;

  vtkRenderWindow *fRenWin;
  vtkRenderWindowInteractor *fIren;

  AutoArrayT<VTKBodyT> bodies;

  vtkActor2D* pointLabels;
  vtkSelectVisiblePoints* visPts;
  vtkLabeledDataMapper* ldm;
  vtkCubeAxesActor2D* axes;


};

#endif
