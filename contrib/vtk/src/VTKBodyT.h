/* $Id: VTKBodyT.h,v 1.9 2001-11-29 21:22:43 recampb Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"
#include "ExodusT.h"
#include "VTKBodyDataT.h"


/* forward declarations */
class vtkPoints;
class vtkCellArray;
class vtkUnstructuredGrid;
class vtkDataSetMapper;
class vtkActor;
class vtkLookupTable;
class vtkIdFilter;
class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;
class vtkScalars;
class vtkWarpVector;
class vtkVectors;
class vtkScalarBarActor;
class ExodusT;
class VTKBodyDataT;

class VTKBodyT: public iConsoleObjectT
{
 public:

  /** default constuctor */
  VTKBodyT(void) { body = NULL; };

  /** constructor */
  VTKBodyT(VTKBodyDataT* body_data);
  
  /** return pointer to the body data */
  VTKBodyDataT* BodyData(void) { return body; };

  /** comparison operator */
  bool operator==(const VTKBodyT& rhs) { return body == rhs.body; };

  /** rvalue - smart pointer */
  VTKBodyDataT* operator->(); //CW wouldn't call functions with conversion

 private:

  VTKBodyDataT* body; 
};

inline VTKBodyDataT* VTKBodyT::operator->()
{
  return body;
}


#endif
