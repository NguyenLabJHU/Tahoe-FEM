/* $Id: VTKConsoleT.cpp,v 1.3 2001-09-26 18:07:57 recampb Exp $ */

#include "VTKConsoleT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkRendererSource.h"
#include "vtkTIFFWriter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkLookupTable.h"


//void VTKConsoleT::dataUpdate(void) {
  

VTKConsoleT::VTKConsoleT(void)
{
  /* set console name */
  iSetName("vtk");

  /* add variables to the console */


  iAddVariable("min_Hue_Range", hueRange1);
  iAddVariable("max_Hue_Range", hueRange2);
  iAddVariable("min_Value_Range", valRange1);
  iAddVariable("max_Value_Range", valRange2);
  iAddVariable("min_Saturation_Range", satRange1);
  iAddVariable("max_Saturation_Range", satRange2);
  iAddVariable("min_Alpha_Range", alphaRange1);
  iAddVariable("max_Alpha_Range", alphaRange2);
  iAddVariable("min_Scalar_Range", scalarRange1);
  iAddVariable("max_Scalar_Range", scalarRange2);
  iAddVariable("numColors", numColors);
  iAddVariable("source_file", source_file);
  iAddVariable("output_file", output_file);
//   iAddVariable("ugr", ugr);
//   iAddVariable("renderer", renderer);
//   iAddVariable("renWin", renWin);
//   iAddVariable("iren", iren);
//   iAddVariable("lut", lut);
//   iAddVariable("ugridMapper", ugridMapper);
//   iAddVariable("ugridActor", ugridActor);
//   iAddVariable("wireActor", wireActor);
//   iAddVariable("scalarBar", scalarBar);
//   iAddVariable("renSrc", renSrc);
  
  /* add console commands */
//   iAddCommand("Integer_Print");
//   iAddCommand("String_Print");
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Save");

  source_file = "../../example_files/heat/data_file1.vtk";
  output_file = "A.tif";
  numColors = 256;
  hueRange1 = 0; hueRange2 = 0.6667;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
  scalarRange1 = 6; scalarRange2 = 17;
  ugr = vtkUnstructuredGridReader::New();
  renderer = vtkRenderer::New();
  renWin = vtkRenderWindow::New();
  iren = vtkRenderWindowInteractor::New();
  lut = vtkLookupTable::New();
  ugridMapper = vtkDataSetMapper::New();
  ugridActor = vtkActor::New();
  wireActor = vtkActor::New();
  scalarBar = vtkScalarBarActor::New();
  writer = vtkTIFFWriter::New();
  renSrc = vtkRendererSource::New();

  ugr->SetFileName(source_file);

  renWin->AddRenderer(renderer);
 
  iren->SetRenderWindow(renWin);

  lut->SetHueRange(hueRange1, hueRange2);
  lut->SetSaturationRange(satRange1,satRange2);
  lut->SetValueRange(valRange1,valRange2);
  lut->SetAlphaRange(alphaRange1,alphaRange2);
  lut->SetNumberOfColors(numColors);
  lut->Build();

  
  ugridMapper->SetInput(ugr->GetOutput());
  ugridMapper->SetScalarRange(scalarRange1,scalarRange2);
  ugridMapper->SetLookupTable(lut);
   ugridMapper->ImmediateModeRenderingOn();


  ugridActor->SetMapper(ugridMapper);
  ugridActor->AddPosition(0,0.001,0);

  
  wireActor->SetMapper(ugridMapper);
  wireActor->GetProperty()->SetRepresentationToWireframe();
  wireActor->GetProperty()->SetColor(0,0,0);

  renderer->AddActor(ugridActor);
  renderer->AddActor(wireActor);
  renderer->SetBackground(0,0,0);

  renWin->SetPosition(668, 0);
  renWin->SetSize(600, 700);

  
  scalarBar->SetLookupTable(ugridMapper->GetLookupTable());
  scalarBar->SetTitle("Temperature for time .50000");
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
  scalarBar->SetOrientationToHorizontal();
  scalarBar->SetWidth(0.8);
  scalarBar->SetHeight(0.17);
  renderer->AddActor(scalarBar);



  // dataUpdate();
}



/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{
//   if (command == "Integer_Print")
//     {
//       cout << "int = " << fInteger << endl;
//       return true;
//     }
//   else if (command == "String_Print")
//     {
//       cout << "string = " << fString << endl;
//       return true;
//     }

  if (command == "Start_Rendering")
  {
    renWin->Render();
    iren->Start();
    return true;
  }

  else if (command == "Update_Rendering")
  {

    ugridMapper->SetScalarRange(scalarRange1,scalarRange2);

    lut->SetHueRange(hueRange1, hueRange2);
    lut->SetSaturationRange(satRange1,satRange2);
    lut->SetValueRange(valRange1,valRange2);
    lut->SetAlphaRange(alphaRange1,alphaRange2);
    lut->SetNumberOfColors(numColors);
    ugr->SetFileName(source_file);
    renWin->Render();
    iren->Start();
    // dataUpdate();
    return true;
  }
    
  else if (command == "Reset_to_Default_Values")
    {
      source_file = "../../example_files/heat/data_file1.vtk";
      numColors = 256;
      hueRange1 = 0; hueRange2 = 0.6667;
      valRange1 = 1; valRange2 = 1;
      satRange1 = 1; satRange2 = 1;
      alphaRange1 = 1; alphaRange2 = 1;
      scalarRange1 = 6; scalarRange2 = 17;
      return true;
    }

  else if (command == "Save")
    {
      //  cout << "Enter name for file to be saved to: ";
      //  cin >> output_file;
    
      renSrc->SetInput(renderer);
      renSrc->WholeWindowOn();
      writer->SetInput(renSrc->GetOutput());
      writer->SetFileName(output_file);
      writer->Write();
      return true;
    }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
