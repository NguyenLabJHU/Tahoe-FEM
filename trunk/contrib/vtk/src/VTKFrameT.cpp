/* $Id: VTKFrameT.cpp,v 1.14 2001-11-20 01:04:04 recampb Exp $ */

#include "VTKFrameT.h"
#include "VTKConsoleT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
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
#include "vtkWarpVector.h"
#include "vtkVectors.h"
#include "vtkTextMapper.h"

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
 
  /* add console commands */
  iAddCommand("Interactive");
  iAddCommand("Update");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Reset_view");
  iAddCommand("Save");
  iAddCommand("Show_Node_Numbers");
  iAddCommand("Hide_Node_Numbers");
  iAddCommand("Color_bar_on");
  iAddCommand("Color_bar_off");
  iAddCommand("X_axis_rotation");
  iAddCommand("Y_axis_rotation");
  iAddCommand("Z_axis_rotation");
  iAddCommand("Zoom");
  iAddCommand("Change_background_color");
  iAddCommand("ChangeDataColor");
  iAddCommand("Select_time_step");
  iAddCommand("Show_axes");
  iAddCommand("Hide_axes");
  iAddCommand("Choose_variable");

  iAddCommand("AddBody");
  iAddCommand("RemoveBody");
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
  renderer->GetActiveCamera()->Zoom(0.85);
  renderer->ResetCamera();
}

bool VTKFrameT::AddBody(VTKBodyDataT* body)
{
  if (bodies.AppendUnique(body))
    {
      renderer->AddActor(body->Actor());
      // renderer->AddActor(body->SBActor());
      // frame needs to keep track of the current visible scalar bar
      // so that at most one is visible and SBActors() can be exhanged
      // in renderer
      ResetView();
      StringT name = "body";
      int index = bodies.PositionOf(body);
      name.Append(index);
      bodies[index]->iSetName(name);
      iAddSub(*bodies[index]);
      return true;
    }
  else
    return false;
}

bool VTKFrameT::RemoveBody(VTKBodyDataT* body)
{
  int index = bodies.PositionOf(body);
  if (index == -1)
    return false;
  else
    {
      VTKBodyDataT* body = bodies[index];
      iDeleteSub(*bodies[index]);
      /* remove from renderer */
      renderer->RemoveActor(body->Actor());
      renderer->RemoveActor(body->SBActor()); // if added to the renderer earlier
      ResetView();

      /* remove from body list */
      bodies.DeleteAt(index);
      return true;
    }
}

void VTKFrameT::ShowFrameNum(StringT Name)
{
  vtkTextMapper* textMapper = vtkTextMapper::New();
  textMapper->SetInput(Name);
  textMapper->SetFontSize(16);
  vtkActor2D* text = vtkActor2D::New();
  text->SetMapper(textMapper);
  text->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  text->GetPositionCoordinate()->SetValue(0.4,.95);
  text->GetProperty()->SetColor(0,1,0);
  renderer->AddActor(text);
  fRenWin->Render();
  
}

/* execute given command - returns false on fail */
bool VTKFrameT::iDoCommand(const StringT& command, StringT& line)
{
  int sfbTest;
  int  varNum;
  double xRot, yRot, zRot, zoom;
  double timeStep;

  if (command == "Update")
    {
      for (int i = 0; i < bodies.Length(); i++)
		bodies[i]->UpdateData();
      fRenWin->Render();
      return true;
    }
  else if (command == "Interactive")
  {
    for (int i = 0; i < bodies.Length(); i++)
      bodies[i]->UpdateData();
    fRenWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;   
    fIren->Start();
    return true;
  }
  else if (command == "AddBody")
    {
      /* list of bodies */
      const ArrayT<VTKBodyDataT*>& bodies = fConsole->Bodies();

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
  else if (command == "RemoveBody")
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

  else if (command == "Show_Node_Numbers")
    {
      pointLabels = vtkActor2D::New();
      visPts = vtkSelectVisiblePoints::New();
      ldm = vtkLabeledDataMapper::New();
      if (bodies[0]->num_node_variables > 0)
	visPts->SetInput(bodies[0]->warp->GetOutput());
      else
	visPts->SetInput(bodies[0]->ugrid);
      visPts->SetRenderer(renderer);
      ldm->SetInput(visPts->GetOutput());
      ldm->SetLabelModeToLabelIds();
      ldm->ShadowOff();
      pointLabels->SetMapper(ldm);
      pointLabels->VisibilityOn();
      renderer->AddActor2D(pointLabels);
      fRenWin->Render();
      return true;      

    }

  else if (command == "Hide_Node_Numbers")
    {
      pointLabels->VisibilityOff();
      pointLabels->Delete();
      visPts->Delete();
      ldm->Delete();
      fRenWin->Render();
      return true;

    }



  else if (command == "Color_bar_off")
    {
      renderer->RemoveActor(bodies[0]->SBActor());
      fRenWin->Render();
      return true;
    }

  else if (command == "Color_bar_on")
    {
      renderer->AddActor(bodies[0]->SBActor());
      fRenWin->Render();

      return true;
    }

  else if (command == "ChangeDataColor")
    {
      int color;
      char line[255];
      cout << "Choose color: \n 1: Red\n 2: Green\n 3: Blue: ";
      cin >> color;
      cin.getline(line, 254);
      bodies[0]->ChangeDataColor(color);
      fRenWin->Render();
      return true;
    }

  else if (command == "X_axis_rotation")
    {
 
      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> xRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Elevation(xRot);
      fRenWin->Render();
      return true;
    }


  else if (command == "Y_axis_rotation")
    {

      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> yRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Azimuth(-yRot);
      fRenWin->Render();
      return true;
    }

  else if (command == "Z_axis_rotation")
    {

      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> zRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Roll(-zRot);
      fRenWin->Render();
      return true;
    }

  else if (command == "Zoom")
    {

      cout << "Zoom by what magnitude?: ";
      cin >> zoom;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Zoom(zoom);
      fRenWin->Render();
      return true;
    }
	 
  else if (command=="Change_background_color")
    {
      int bgColor;
      do {
	cout << "choose background color:\n 1: black\n 2: white\n 3: red\n 4: green\n 5: blue\n";
	cin >> bgColor;
	char line[255];
	cin.getline(line, 254);
      } while (bgColor != 1 && bgColor !=2 && bgColor != 3 && bgColor !=4 && bgColor !=5);
      
      switch (bgColor) {
      case 1: 
	renderer->SetBackground(0,0,0);
	break;
      case 2:
	renderer->SetBackground(1,1,1);
	break;
      case 3:
	renderer->SetBackground(1,0,0);
	break;
      case 4:
	renderer->SetBackground(0,1,0);
	break;
      default:
	renderer->SetBackground(0,0,1);
	break;
      }
      fRenWin->Render();
      return true;
    }

  else if (command == "Select_time_step")
    {
      int step;
      cout << "choose frame number from 0 to " << bodies[0]->num_time_steps-1 <<" to be displayed: ";
      cin >> step;
      for (int i = 0; i < bodies.Length(); i++)
	bodies[i]->SelectTimeStep(step);
      char line[255];
      cin.getline(line, 254);
 
     
      fRenWin->Render();
      return true;
      
    }

  
  else if (command == "Show_axes")
    {
      
      // x,y,z axes
      axes = vtkCubeAxesActor2D::New();
      axes->SetInput(bodies[0]->ugrid);
      axes->SetCamera(renderer->GetActiveCamera());
      axes->SetLabelFormat("%6.4g");
      //axes->SetCornerOffset(.2);
      // axes->ShadowOn();
      // axes->SetFlyModeToOuterEdges();
      axes->SetFlyModeToClosestTriad();
      //axes->SetFontFactor(1.8);
      axes->GetProperty()->SetColor(0,1,1);
      // axes->SetBounds(0,1,0,1,0,1);
      //axes->ZAxisVisibilityOff();
      axes->VisibilityOn();
      renderer->AddActor2D(axes);
      fRenWin->Render();
      return true;
    }
  
  else if (command == "Hide_axes")
    {
      axes->VisibilityOff();
      axes->Delete();
      fRenWin->Render();
      return true;
    }

  else if (command == "Choose_variable")
    {
	  //NOTE: collect list of variables from all bodies???

      cout << "choose variable number from 0 to " << bodies[0]->num_node_variables-1 <<" to be displayed\n" << bodies[0]->varList;      
      cin >> varNum;
      char line[255];
      cin.getline(line, 254);
	  const StringT& var = (bodies[0]->NodeLabels())[varNum];
	  for (int i = 0; i < bodies.Length(); i++)
		bodies[i]->ChangeVars(var); // will not return true if body does not have the var
      fRenWin->Render();
      return true;

    }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
