/* $Id: VTKFrameT.h,v 1.18 2002-09-03 07:10:36 paklein Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "AutoArrayT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

using namespace Tahoe; 

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
  VTKBodyT* Body(int bodyNum) { return bodies[bodyNum]; };
  
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

	/** update data. Only updates data managed directly by the frame. */
	void UpdateData(void);

 protected:

	/** call to re-render the window contents. Call goes through the console,
	 * so all window contents will be brought up to date and re-drawn. */
	void Render(void);

	/** write prompt for the specific argument of the command */
	virtual void ValuePrompt(const CommandSpecT& command, int index, ostream& out) const;  

	/** transform from display to world coordinates.
	 * \param worldPt has to be allocated as 4 vector */
	void ComputeDisplayToWorld(double x, double y, double z, double *worldPt);

	/** transform from world to display coordinates.
	 * \param displayPt has to be allocated as 3 vector */
	void ComputeWorldToDisplay(double x, double y, double z, double *displayPt);

	/** get the current time string */
	void GetCurrentTime(StringT& string) const;

 private:

	/** controlling console object */
	VTKConsoleT& fConsole;
  
  	/** renderer for this frame */
	vtkRenderer *fRenderer;

	/** bodies appearing in the frame */
	AutoArrayT<VTKBodyT*> bodies;

	/** scalar bar appearing in the frame */
	vtkScalarBarActor* scalarBar;
	
	/** frame label */
	StringT        fFrameLabel;
	vtkTextMapper* fLabelMapper;
	vtkActor2D*    fLabelActor;

	/** time label */
	StringT        fTimeLabel;
	vtkTextMapper* fTimeLabelMapper;
	vtkActor2D*    fTimeLabelActor;
};

#endif
