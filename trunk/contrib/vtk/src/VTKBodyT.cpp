/* $Id: VTKBodyT.cpp,v 1.12 2001-11-29 21:22:43 recampb Exp $ */

#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

#include "vtkPoints.h"

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkRendererSource.h"
#include "vtkWindowToImageFilter.h"
#include "vtkLookupTable.h"
#include "vtkIdFilter.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkLabeledDataMapper.h"
#include "vtkActor2D.h"
#include "vtkFieldData.h"
#include "vtkCamera.h"
#include "vtkWarpVector.h"
#include "vtkVectors.h"
#include "vtkScalarBarActor.h"
#include "StringT.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"

/* array behavior */
const bool ArrayT<VTKBodyT>::fByteCopy = true;

/* constructor */
VTKBodyT::VTKBodyT(VTKBodyDataT* body_data)
{
  body = body_data;
}

#if 0
/* conversion */
VTKBodyT::operator VTKBodyDataT()
{
  if (!body) throw eGeneralFail;
  return *body;
}
#endif
