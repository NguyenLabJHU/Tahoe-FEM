/* $Id: VTKBodyT.h,v 1.13 2001-12-13 02:57:59 paklein Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"
#include "VTKBodyDataT.h"

/* forward declarations */
class VTKBodyDataT;
class VTKFrameT;
class vtkCubeAxesActor2D;

class vtkIdFilter;
class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;

/** interface to console graphics object. Each appearance of an
 * object on screen has its own VTKBodyT which may share the
 * underlying model data in the VTKBodyDataT. */
class VTKBodyT: public iConsoleObjectT
{
 public:

	/** default constuctor */
	VTKBodyT(void);

	/** constructor */
	VTKBodyT(VTKFrameT* frame, VTKBodyDataT* body_data);
  
	/** destructor */
	~VTKBodyT(void);
  
	/** return pointer to the body data */
	VTKBodyDataT* BodyData(void) { return fBodyData; };

	/** comparison operator */
	bool operator==(const VTKBodyT& rhs) { return fBodyData == rhs.fBodyData; };

	/** rvalue - smart pointer */
	VTKBodyDataT* operator->();

	/** execute console command. \return true is executed normally */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

 	/** add actors in self to the given renderer */
 	void AddToFrame(void);

 	/** add actors in self to the given renderer */
 	void RemoveFromFrame(void);

 	/** show node numbers */
 	//void ShowNodeNumbers(vtkRenderer* renderer);

 	/** hide node numbers */
 	//void HideNodeNumbers(vtkRenderer* renderer);
 	
 	/** change the plot variable */
	bool ChangeVars(const StringT& var);

 private:

	/** frame where body is displayed */
	VTKFrameT* fFrame;

	/** body data */
	VTKBodyDataT* fBodyData;
	
	/** coordinate axes */
	ArrayT<vtkCubeAxesActor2D*> fAxes;
	
	/* node numbers */
	ArrayT<vtkIdFilter*> fIDFilter;
	ArrayT<vtkSelectVisiblePoints*> fVisPoints;
	ArrayT<vtkLabeledDataMapper*> fNodeLabelMapper;
	ArrayT<vtkActor2D*> fNodeLabelActor;	
};

inline VTKBodyDataT* VTKBodyT::operator->()
{
  return fBodyData;
}

#endif
