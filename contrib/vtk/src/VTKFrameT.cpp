/* $Id: VTKFrameT.cpp,v 1.15 2001-11-29 21:22:43 recampb Exp $ */

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
#include "CommandSpecT.h"
#include "ArgSpecT.h"

/* constructor */
VTKFrameT::VTKFrameT(void):
  renderer(NULL)
{
  /* set up display classes */
  renderer = vtkRenderer::New();
 
  /* add console commands */
  iAddCommand(CommandSpecT("Interactive"));
  iAddCommand(CommandSpecT("Update"));
  iAddCommand(CommandSpecT("Reset_to_Default_Values"));
  iAddCommand(CommandSpecT("Reset_view"));
  iAddCommand(CommandSpecT("Save"));
  iAddCommand(CommandSpecT("Show_Node_Numbers"));
  iAddCommand(CommandSpecT("Hide_Node_Numbers"));
  iAddCommand(CommandSpecT("Color_bar_on"));
  iAddCommand(CommandSpecT("Color_bar_off"));

  CommandSpecT rotate("Rotate", false);
  ArgSpecT rot_x(ArgSpecT::double_, "x");
  rot_x.SetDefault(0.0);
  rot_x.SetPrompt("x-axis rotation");
  ArgSpecT rot_y(ArgSpecT::double_, "y");
  rot_y.SetDefault(0.0);
  rot_y.SetPrompt("y-axis rotation");
  ArgSpecT rot_z(ArgSpecT::double_, "z");
  rot_z.SetDefault(0.0);
  rot_z.SetPrompt("y-axis rotation");
  rotate.AddArgument(rot_x);
  rotate.AddArgument(rot_y);
  rotate.AddArgument(rot_z);
  iAddCommand(rotate);

  iAddCommand(CommandSpecT("Zoom"));
  iAddCommand(CommandSpecT("Change_background_color"));
  iAddCommand(CommandSpecT("ChangeDataColor"));
  iAddCommand(CommandSpecT("Select_time_step"));
  iAddCommand(CommandSpecT("Show_axes"));
  iAddCommand(CommandSpecT("Hide_axes"));
  iAddCommand(CommandSpecT("Choose_variable"));

  ArgSpecT body_num(ArgSpecT::int_);
  body_num.SetPrompt("body number");

  CommandSpecT add_body("AddBody");
  add_body.AddArgument(body_num);
  iAddCommand(add_body);

  CommandSpecT rem_body("RemoveBody");
  rem_body.AddArgument(body_num);
  iAddCommand(rem_body);
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

bool VTKFrameT::AddBody(VTKBodyDataT* body_data)
{
  /* add body to the list */
  if (bodies.AppendUnique(body_data))
    {
      renderer->AddActor(body_data->Actor());
      // renderer->AddActor(body->SBActor());
      // frame needs to keep track of the current visible scalar bar
      // so that at most one is visible and SBActors() can be exhanged
      // in renderer
      ResetView();
      StringT name = "body";
      int index = bodies.PositionOf(body_data);
      name.Append(index);
      bodies[index]->iSetName(name);
      iAddSub(bodies[index]);
      return true;
    }
  else
    return false;
}

bool VTKFrameT::RemoveBody(VTKBodyDataT* body_data)
{
  /* look for body data */
  int index = bodies.PositionOf(body_data);
  if (index == -1)
    return false;
  else
    {
      VTKBodyT& body = bodies[index];
      iDeleteSub(body);

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
bool VTKFrameT::iDoCommand(const CommandSpecT& command, StringT& line)
{
  int sfbTest;
  int  varNum;
  double xRot, yRot, zRot, zoom;
  double timeStep;

  if (command.Name() == "Update")
    {
      for (int i = 0; i < bodies.Length(); i++)
		bodies[i]->UpdateData();
      fRenWin->Render();
      return true;
    }
  else if (command.Name() == "Interactive")
  {
    for (int i = 0; i < bodies.Length(); i++)
      bodies[i]->UpdateData();
    fRenWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;   
    fIren->Start();
    return true;
  }
  else if (command.Name() == "AddBody")
    {
      int body;
      command.Argument(0).GetValue(body);

      const ArrayT<VTKBodyDataT*>& all_bodies = fConsole->Bodies();
      if (body >= 0 && body < all_bodies.Length())
	{
	  if (AddBody(all_bodies[body]))
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
  else if (command.Name() == "RemoveBody")
    {
      int body;
      command.Argument(0).GetValue(body);

      if (body >= 0 && body < bodies.Length())
	{
	  if (RemoveBody(bodies[body].BodyData()))
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

  else if (command.Name() == "Reset_view")
    {
      ResetView();
      renderer->ResetCamera();
      fRenWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      fIren->Start();
      return true;
    }

  else if (command.Name() == "Show_Node_Numbers")
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

  else if (command.Name() == "Hide_Node_Numbers")
    {
      pointLabels->VisibilityOff();
      pointLabels->Delete();
      visPts->Delete();
      ldm->Delete();
      fRenWin->Render();
      return true;

    }



  else if (command.Name() == "Color_bar_off")
    {
      renderer->RemoveActor(bodies[0]->SBActor());
      fRenWin->Render();
      return true;
    }

  else if (command.Name() == "Color_bar_on")
    {
      renderer->AddActor(bodies[0]->SBActor());
      fRenWin->Render();

      return true;
    }

  else if (command.Name() == "ChangeDataColor")
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

  else if (command.Name() == "Rotate")
    {
      double x, y, z;
      command.Argument(0).GetValue(x);
      command.Argument(1).GetValue(y);
      command.Argument(2).GetValue(z);

      if (fabs(x) > 1.0e-6)
	renderer->GetActiveCamera()->Elevation(x);

      if (fabs(y) > 1.0e-6)
	renderer->GetActiveCamera()->Azimuth(y);

      if (fabs(z) > 1.0e-6)
	renderer->GetActiveCamera()->Roll(z);

      fRenWin->Render();
      return true;
    }

  else if (command.Name() == "Zoom")
    {

      cout << "Zoom by what magnitude?: ";
      cin >> zoom;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Zoom(zoom);
      fRenWin->Render();
      return true;
    }
	 
  else if (command.Name() =="Change_background_color")
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

  else if (command.Name() == "Select_time_step")
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

  
  else if (command.Name() == "Show_axes")
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
  
  else if (command.Name() == "Hide_axes")
    {
      axes->VisibilityOff();
      axes->Delete();
      fRenWin->Render();
      return true;
    }

  else if (command.Name() == "Choose_variable")
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

/* write prompt for the specific argument of the command */
void VTKFrameT::ValuePrompt(const CommandSpecT& command, int index, ostream& out) const
{
  if (command.Name() == "AddBody")
    {
      /* list of bodies */
      const ArrayT<VTKBodyDataT*>& bodies = fConsole->Bodies();
      out << "body to add (0," << bodies.Length()-1 << ")\n";
    }
  else if (command.Name() == "RemoveBody")
    {
      out << "body to remove (0," << bodies.Length()-1 << ")\n";
    }
  else /* inherited */
    iConsoleObjectT::ValuePrompt(command, index, out);
}
