
/* $Id: VTKConsoleT.h,v 1.31 2002-09-22 19:54:59 paklein Exp $ */

#ifndef _VTK_CONSOLE_T_H_
#define _VTK_CONSOLE_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

/* direct members */
#include "VTKFrameT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"
#include "StringT.h"
#include "AutoArrayT.h"
#include "Array2DT.h"

using namespace Tahoe; 

/* forward declarations */
class VTKBodyT;
class VTKBodyDataT;
class vtkPointPicker;
class vtkCellPicker;
class vtkActor;
class vtkActorCollection;

class VTKConsoleT: public iConsoleObjectT
{
 public:

  /** constructor */
  VTKConsoleT(const ArrayT<StringT>& arguments);

  /** destructor */
  ~VTKConsoleT(void);

  /** execute given command. \return true if OK, false on fail */
  virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

  /** return the list of bodies */
  const ArrayT<VTKBodyDataT*> Bodies(void) const { return fBodies; };
  
  static AutoArrayT<vtkActor*> pickedPoints;

 protected:
 
    /** write prompt for the specific argument of the command */
   virtual void ValuePrompt(const CommandSpecT& command, int index, ostream& out) const;  

 private:
  
  /** construct body from the given file path. If successful,
   * new body is appended to the end of VTKConsoleT::fBodies.
   * \return true if successful, false otherwise */
  bool AddBody(const StringT& format, const StringT& file);

  /** reset the frame layout. Frames that fit will be transferred
   * to the new layout. Others will be discarded. 
   * \param num_x number of frames in the horizontal direction 
   * \param num_y number of frames in the vertical direction */
  void SetFrameLayout(int num_x, int num_y);

  /** prints picked point and scalar value to screen */
  static void PickPoints(void *);

  static void PickCells(void *);



  /** returns the index of the requested option */
  bool CommandLineOption(const char* str, int& index, int start = 0) const;

 private:

  /** argument list passed into constructor */
  const ArrayT<StringT> fArguments;

  /** the display window */
  vtkRenderWindow *renWin;
  bool fRenderHold;

  /** interactor to rotate, translate, zoom, etc */
  vtkRenderWindowInteractor *iren;

  /** list of frame in the display window */
  Array2DT<VTKFrameT*> fFrames;

  /** list of display objects */
  int fBodyCount;
  AutoArrayT<VTKBodyDataT*> fBodies;
  vtkPointPicker* pointPicker;
  vtkCellPicker* cellPicker;


  //static vtkActorCollection* pickedPoints;
 
};

#endif
