/* $Id: VTKFrameT.h,v 1.15 2001-12-10 12:44:08 paklein Exp $ */

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
class vtkCubeAxesActor2D;
class vtkScalarBarActor;
class vtkTextMapper;
class vtkActor2D;

/* forward declarations */
class VTKBodyT;
class VTKBodyDataT;
class VTKConsoleT;

class VTKFrameT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKFrameT(VTKConsoleT& console);

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
  
  void ResetView(void);

	/** label in the frame */
  	void ShowFrameLabel(const StringT& label);

	/** remove the frame label */
  	void HideFrameLabel(void);

  /** return a pointer to the frame's renderer */
  vtkRenderer* Renderer(void) { return fRenderer; };

	/** execute console command. \return true is executed normally */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);
   
   /** set controlling console object */
//   void setConsole(VTKConsoleT* console) { fConsole = console; };

   /** set the renderer window */
//   void setRenWin(vtkRenderWindow* renWin) { fRenWin = renWin; };
   
   /** set the window interactor */
//   void setIren(vtkRenderWindowInteractor* iren) {fIren = iren; };

//   const StringT& getName(void) const {return bodies[0]->SourceFile(); };
//what was this for????

 protected:

	/** call to re-render the window contents. Call goes through the console,
	 * so all window contents will be brought up to date and re-drawn. */
	void Render(void) const;

   /** write prompt for the specific argument of the command */
   virtual void ValuePrompt(const CommandSpecT& command, int index, ostream& out) const;  

 private:

	/** controlling console object */
	VTKConsoleT& fConsole;
  
  	/** renderer for this frame */
	vtkRenderer *fRenderer;
	
	/** window where frame appears */
//	vtkRenderWindow *fRenWin;

	/** interactor for the window in which the frame appears */
//	vtkRenderWindowInteractor *fIren;

	/** bodies appearing in the frame */
	AutoArrayT<VTKBodyT> bodies;

	/** scalar bar appearing in the frame */
	vtkScalarBarActor* scalarBar;
	
	/** frame label */
	vtkTextMapper* fLabelMapper;
	vtkActor2D*    fLabelActor;
};

#endif
