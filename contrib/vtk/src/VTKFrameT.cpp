/* $Id: VTKFrameT.cpp,v 1.8 2001-11-06 02:39:51 recampb Exp $ */

#include "VTKFrameT.h"
#include "VTKConsoleT.h"
#include "VTKBodyT.h"
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
//#include "vtkTIFFWriter.h"
//#include "vtkWindowToImageFilter.h"
#include "vtkLookupTable.h"
#include "vtkIdFilter.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkLabeledDataMapper.h"
#include "vtkActor2D.h"
#include "vtkFieldData.h"
#include "vtkCamera.h"
#include "vtkWarpVector.h"
#include "vtkVectors.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"
#include "VTKConsoleT.h"

/* constructor */
VTKFrameT::VTKFrameT(void):
  renderer(NULL)
{
  /* set up display classes */
  renderer = vtkRenderer::New();
 
  /* set up console modifiable variables */
//   iAddVariable("age", fAge);
//   iAddVariable("str_len", (const int) fLength);
//   iAddVariable("string", (const char*) fArray);
 /* add variables to the console */
//   iAddVariable("min_Hue_Range", hueRange1);
//   iAddVariable("max_Hue_Range", hueRange2);
//   iAddVariable("min_Value_Range", valRange1);
//   iAddVariable("max_Value_Range", valRange2);
//   iAddVariable("min_Saturation_Range", satRange1);
//   iAddVariable("max_Saturation_Range", satRange2);
//   iAddVariable("min_Alpha_Range", alphaRange1);
//   iAddVariable("max_Alpha_Range", alphaRange2);
//   iAddVariable("numColors", numColors);
//   //iAddVariable("source_file", source_file);
//   iAddVariable("min_Scalar_Range", scalarRange1);
//   iAddVariable("max_Scalar_Range", scalarRange2);
//   iAddVariable("scale_factor", scale_factor);

  /* add console commands */
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Reset_view");
  iAddCommand("Save");
  //iAddCommand("Save_flip_book_images");

  iAddCommand("Show_Node_Numbers");
  iAddCommand("Hide_Node_Numbers");
  iAddCommand("Color_bar_on");
  iAddCommand("Color_bar_off");
  iAddCommand("X_axis_rotation");
  iAddCommand("Y_axis_rotation");
  iAddCommand("Z_axis_rotation");
  //iAddCommand("Flip_book");
  iAddCommand("Change_background_color");
  iAddCommand("Select_time_step");
  iAddCommand("Show_axes");
  iAddCommand("Hide_axes");
  iAddCommand("Choose_variable");

  iAddCommand("Add_body");
  iAddCommand("Remove_body");
}

/* destructor */
VTKFrameT::~VTKFrameT(void)
{
  /* free display classes */
  renderer->Delete(); renderer = NULL;
}

void VTKFrameT::ResetView(void)
{
  renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
  renderer->GetActiveCamera()->SetPosition(0,0,1);
  renderer->GetActiveCamera()->ComputeViewPlaneNormal();
  renderer->GetActiveCamera()->SetViewUp(0,1,0);
  renderer->GetActiveCamera()->OrthogonalizeViewUp();
  renderer->ResetCamera();
}

bool VTKFrameT::AddBody(VTKBodyT* body)
{
  if (bodies.AppendUnique(body))
    {
      renderer->AddActor(body->Actor());
      renderer->AddActor(body->SBActor());
      // frame needs to keep track of the current visible scalar bar
      // so that at most one is visible and SBActors() can be exhanged
      // in renderer
      ResetView();
      return true;
    }
  else
    return false;
}

bool VTKFrameT::RemoveBody(VTKBodyT* body)
{
  int index = bodies.PositionOf(body);
  if (index == -1)
    return false;
  else
    {
      VTKBodyT* body = bodies[index];
      
      /* remove from renderer */
      renderer->RemoveActor(body->Actor());
      renderer->RemoveActor(body->SBActor()); // if added to the renderer earlier
      ResetView();

      /* remove from body list */
      bodies.DeleteAt(index);
      return true;
    }
}

/* execute given command - returns false on fail */
bool VTKFrameT::iDoCommand(const StringT& command, StringT& line)
{
  int sfbTest;
  int  varNum;
  double xRot, yRot, zRot;
  double timeStep;

  if (command == "Start_Rendering")
    {
      for (int i = 0; i < bodies.Length(); i++)
	bodies[i]->UpdateData();
      fRenWin->Render();
      cout << getName() << endl;
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      fIren->Start();
      return true;
    }

  else if (command == "Update_Rendering")
  {
    for (int i = 0; i < bodies.Length(); i++)
      bodies[i]->UpdateData();
    fRenWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;   
    fIren->Start();
    return true;
  }
  else if (command == "Add_body")
    {
      /* list of bodies */
      const ArrayT<VTKBodyT*>& bodies = fConsole->Bodies();

      cout << "body to add (0," << bodies.Length()-1 << "): ";
      int body = -99;
      cin >> body;
      char line[255];
      cin.getline(line, 254);
      if (body >= 0 && body < bodies.Length())
	{
	  if (AddBody(bodies[body]))
	    {
	      
	      fRenWin->Render();
	      return true;
	    }
	  else
	    {
	      cout << "body already added to frame" << endl;
	      return false;
	    }
	}
      else
	{
	  cout << "body number out of range" << endl;
	  return false;
	}
    }
  else if (command == "Remove_body")
    {
      cout << "body to remove (0," << bodies.Length()-1 << "): ";
      int body = -99;
      cin >> body;
      char line[255];
      cin.getline(line, 254);
      if (body >= 0 && body < bodies.Length())
	{
	  if (RemoveBody(bodies[body]))
	    {
	      fRenWin->Render();
	      return true;
	    }
	  else
	    {
	      cout << "frame did not contain the body" << endl;
	      return false;
	    }
	}
      else
	{
	  cout << "body number out of range" << endl;
	  return false;
	}
    }
    
//   else if (command == "Reset_to_Default_Values")
//     {
//       // source_file = "../../example_files/heat/data_file1.vtk";
//       numColors = 256;
//       hueRange1 = 0; hueRange2 = 0.6667;
//       valRange1 = 1; valRange2 = 1;
//       satRange1 = 1; satRange2 = 1;
//       alphaRange1 = 1; alphaRange2 = 1;
//       // scalarRange1 = 6; scalarRange2 = 17;
//       //  ugridMapper->SetScalarRange(scalarRange1,scalarRange2);
//       lut->SetHueRange(hueRange1, hueRange2);
//       lut->SetSaturationRange(satRange1,satRange2);
//       lut->SetValueRange(valRange1,valRange2);
//       lut->SetAlphaRange(alphaRange1,alphaRange2);
//       lut->SetNumberOfColors(numColors);
      
//       renWin->Render();
//       // iren->Start();
//       return true;
//     }

  else if (command == "Reset_view")
    {
      ResetView();
      renderer->ResetCamera();
      fRenWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      fIren->Start();
      return true;
    }

//   else if (command == "Save")
//     {
//       cout << "Enter name for file to be saved to: ";
//       cin >> output_file;
//       char line[255];
//       cin.getline(line, 254);
//       cout << "Save image at: \n 1: current view\n 2: default view: ";
//       cin >> sfbTest;
//       cin.getline(line, 254);
//       /* if default camera desired */
//       if (sfbTest == 2) {
// 	  cam->SetFocalPoint(0,0,0);
// 	  cam->SetPosition(0,0,1);
// 	  cam->ComputeViewPlaneNormal();
// 	  cam->SetViewUp(0,1,0);
// 	  cam->OrthogonalizeViewUp();
// 	  renderer->SetActiveCamera(cam);
// 	  renderer->ResetCamera();
// 	  renWin->Render();
// 	}
//       int saveOpt;
//       renSrc->SetInput(renderer);
//       if (numRen == 4){ 
// 	cout << "Save options:\n 1: single plot\n 2: all 4 plots: ";
// 	cin >> saveOpt;
// 	cin.getline(line, 254);
// 	if (saveOpt == 2)
// 	  renSrc->WholeWindowOn();
	
//       }
//       writer->SetInput(renSrc->GetOutput());
//       writer->SetFileName(output_file);
//       writer->Write();
//       renWin->Render();
//       cout << "File " << output_file << " has been saved." << endl;
//       //  iren->Start();
//       return true;
//     }

//   else if (command == "Show_Node_Numbers")
//     {
//       //ids->SetInput(ugrid);
//       //ids->PointIdsOn();
//       //ids->CellIdsOn();
//       //ids->FieldDataOn();
//       //visPts->SetInput(ids->GetOutput());
//       visPts->SetInput(warp->GetOutput());
//       visPts->SetRenderer(renderer);

//       //ldm->SetInput(warp->GetOutput());
//       ldm->SetInput(visPts->GetOutput());
//       //ldm->SetInput(ids->GetOutput());
//       ldm->SetLabelModeToLabelIds();
//       ldm->ShadowOff();
//       // ldm->SetLabelModeToLabelFieldData();
//       pointLabels->SetMapper(ldm);
//       pointLabels->VisibilityOn();
//       renderer->AddActor2D(pointLabels);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;      

//     }


//   else if (command == "Color_bar_off")
//     {
//       renderer->RemoveActor(scalarBar);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;
//     }

//   else if (command == "Color_bar_on")
//     {
//       renderer->AddActor(scalarBar);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;
//     }


//   else if (command == "X_axis_rotation")
//     {
 
//       cout << "Using the right-hand rule, rotate by how many degrees?: ";
//       cin >> xRot;
//       char line[255];
//       cin.getline(line, 254);
//       renderer->GetActiveCamera()->Elevation(xRot);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       //iren->Start();
//       return true;
//     }


//   else if (command == "Y_axis_rotation")
//     {

//       cout << "Using the right-hand rule, rotate by how many degrees?: ";
//       cin >> yRot;
//       char line[255];
//       cin.getline(line, 254);
//       renderer->GetActiveCamera()->Azimuth(yRot);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       //iren->Start();
//       return true;
//     }

//   else if (command == "Z_axis_rotation")
//     {

//       cout << "Using the right-hand rule, rotate by how many degrees?: ";
//       cin >> zRot;
//       char line[255];
//       cin.getline(line, 254);
//       renderer->GetActiveCamera()->Roll(zRot);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       //iren->Start();
//       return true;
//     }

  //  else if (command == "Flip_book")
  //{
  //}

  //  else if (command== "Save_flip_book_images")
  //{
  //}
	 
//   else if (command=="Change_background_color")
//     {
//       int bgColor;
//       do {
// 	cout << "choose background color:\n 1: black\n 2: white\n 3: red\n 4: green\n 5: blue\n";
// 	cin >> bgColor;
// 	char line[255];
// 	cin.getline(line, 254);
//       } while (bgColor != 1 && bgColor !=2 && bgColor != 3 && bgColor !=4 && bgColor !=5);
      
//       switch (bgColor) {
//       case 1: 
// 	renderer->SetBackground(0,0,0);
// 	break;
//       case 2:
// 	renderer->SetBackground(1,1,1);
// 	break;
//       case 3:
// 	renderer->SetBackground(1,0,0);
// 	break;
//       case 4:
// 	renderer->SetBackground(0,1,0);
// 	break;
//       default:
// 	renderer->SetBackground(0,0,1);
// 	break;
//       }
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;
//     }

  else if (command == "Select_time_step")
    {
      int step;
      cout << "choose frame number from 0 to " << bodies[0]->num_time_steps-1 <<" to be displayed: ";
      cin >> step;
      bodies[0]->SelectTimeStep(step);
      char line[255];
      cin.getline(line, 254);

     
      fRenWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      fIren->Start();
      return true;
      
    }

//    else if (command == "Hide_Node_Numbers")

//     {

// // //       ids->PointIdsOff();
// // //       ids->Update();
//       pointLabels->VisibilityOff();
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;   
 
//     }
  
//   else if (command == "Show_axes")
//   {

//   // x,y,z axes
      
//       axes->SetInput(ugrid);
//      axes->SetCamera(renderer->GetActiveCamera());
//     //  axes->SetLabelFormat("%6.4g");
//      // axes->ShadowOn();
//      // axes->SetFlyModeToOuterEdges();
//       axes->SetFlyModeToClosestTriad();
//       axes->SetFontFactor(3.8);
//       axes->GetProperty()->SetColor(0,1,1);
//       // axes->SetBounds(0,1,0,1,0,1);
//       //axes->ZAxisVisibilityOff();
//       axes->VisibilityOn();
//       renderer->AddActor2D(axes);
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;
//   }

//   else if (command == "Hide_axes")
//     {
//       axes->VisibilityOff();
//       renWin->Render();
//       cout << "type 'e' in the graphics window to exit interactive mode" << endl;
//       iren->Start();
//       return true;
//     }

  else if (command == "Choose_variable")
    {
      cout << "choose variable number from 0 to " << bodies[0]->num_node_variables-1 <<" to be displayed\n" << bodies[0]->varList;      
      cin >> varNum;
      char line[255];
      cin.getline(line, 254);
      bodies[0]->ChangeVars(varNum);
      fRenWin->Render();
      return true;

    }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
