/* $Id: VTKBodyT.cpp,v 1.11 2001-11-20 01:04:03 recampb Exp $ */

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
const bool ArrayT<VTKBodyT*>::fByteCopy = true;

/* constructor */
VTKBodyT::VTKBodyT(VTKBodyDataT* body_data)
{



}

/* destructor */
VTKBodyT::~VTKBodyT(void)
{
  /* clean up */
  points->Delete();
  vtk_cell_array->Delete();
  ugrid->Delete();
  warp->Delete();
  ugridMapper->Delete();
  ugridActor->Delete();
  lut->Delete();
  for (int i = 0; i<num_time_steps; i++)
    for (int j = 0; j<num_node_variables; j++)
      {
      scalars[i][j]->Delete();
      vectors[i][j]->Delete();
      }
}

void VTKBodyT::SetLookupTable(void)
{
  if (node_labels.Length() > 0) {
	sbTitle.Append(node_labels[currentVarNum]); 
	sbTitle.Append(" for time step 000");
  }
  scalarBar->SetLookupTable(ugridMapper->GetLookupTable());

  scalarBar->SetTitle(sbTitle);
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
  scalarBar->SetOrientationToHorizontal();
  scalarBar->SetWidth(0.8); 
  scalarBar->SetHeight(0.17);  
}

void VTKBodyT::UpdateData(void)
{
  if (num_node_variables >0){
    lut->SetHueRange(hueRange1, hueRange2);
    lut->SetSaturationRange(satRange1, satRange2);
    lut->SetValueRange(valRange1, valRange2);
    lut->SetAlphaRange(alphaRange1, alphaRange2);
    lut->SetNumberOfColors(numColors);
    warp->SetScaleFactor(scale_factor);
    ugridMapper->SetScalarRange(scalarRange1[currentVarNum], scalarRange2[currentVarNum]);
    ugridMapper->SetLookupTable(lut);
  }
}

void VTKBodyT::DefaultValues(void)
{
  numColors = 256;
  hueRange1 = 0.6667; hueRange2 = 0;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
}

bool VTKBodyT::ChangeVars(const StringT& var)
{
  /* find variable number */
  int varNum = -1;
  for (int i = 0; varNum == -1 && i < node_labels.Length(); i++)
	if (node_labels[i] == var)
	  varNum = i;

  /* change if found */
  if (varNum == -1)
	return false;
  else {
	ugrid->GetPointData()->SetScalars(scalars[currentStepNum][varNum]);
	ugridMapper->SetScalarRange(scalarRange1[varNum],scalarRange2[varNum]);
	if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	  ugrid->GetPointData()->SetVectors(vectors[currentStepNum][varNum]);

	sbTitle = "";
	sbTitle.Append(node_labels[varNum]); 
	sbTitle.Append(" for time step 000 ");
	//sbTitle.Append(currentStepNum,3);
	scalarBar->SetTitle(sbTitle);
	currentVarNum = varNum;
	return true;
  }
}

void VTKBodyT::SelectTimeStep(int stepNum)
{
  if (num_node_variables >0){
    sbTitle.Drop(-3);
    sbTitle.Append(stepNum,3);
    scalarBar->SetTitle(sbTitle);
    // ugrid->GetPointData()->SetScalars(scalars[currentStepNum]);
    ugrid->GetPointData()->SetScalars(scalars[stepNum][currentVarNum]);
    if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
      ugrid->GetPointData()->SetVectors(vectors[stepNum][currentVarNum]);
    currentStepNum= stepNum;
  }
}

void VTKBodyT::ChangeDataColor(int color)
{
  ugridMapper->ScalarVisibilityOff();
  if (color ==1)
    ugridActor->GetProperty()->SetColor(1,0,0);
  else if (color==2)
    ugridActor->GetProperty()->SetColor(0,1,0);
  else if (color==3)
    ugridActor->GetProperty()->SetColor(0,0,1);
  else
    cout << "invalid color";
}
