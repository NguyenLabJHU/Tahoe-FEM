/* $Id: VTKFrameT.cpp,v 1.17 2001-12-10 12:44:08 paklein Exp $ */

#include "VTKFrameT.h"
#include "VTKConsoleT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPoints.h"
//#include "vtkUnstructuredGrid.h"
//#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
//#include "vtkCubeAxesActor2D.h"
//#include "vtkRendererSource.h"
#include "vtkLookupTable.h"
//#include "vtkWarpVector.h"
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
VTKFrameT::VTKFrameT(VTKConsoleT& console):
	fConsole(console),
	fRenderer(NULL),
	fLabelMapper(NULL),
	fLabelActor(NULL)
{
	/* set up display classes */
	fRenderer = vtkRenderer::New();
	
	/* borrow some commands from the controlling console */
	CommandSpecT* com_spec;
	com_spec = fConsole.iCommand("Update");
	if (!com_spec) throw eGeneralFail;
	iAddCommand(*com_spec);

	com_spec = fConsole.iCommand("Interactive");
	if (!com_spec) throw eGeneralFail;
	iAddCommand(*com_spec);

	com_spec = fConsole.iCommand("SelectTimeStep");
	if (!com_spec) throw eGeneralFail;
	iAddCommand(*com_spec);

  /* own console commands */
//  iAddCommand(CommandSpecT("Reset_to_Default_Values"));
  iAddCommand(CommandSpecT("ResetView"));
  iAddCommand(CommandSpecT("Save"));
  iAddCommand(CommandSpecT("ShowNodeNumbers"));
  iAddCommand(CommandSpecT("HideNodeNumbers"));
//  iAddCommand(CommandSpecT("Color_bar_on"));
//  iAddCommand(CommandSpecT("Color_bar_off"));

  CommandSpecT rotate("Rotate", false);
  ArgSpecT rot_x(ArgSpecT::double_, "x");
  rot_x.SetDefault(0.0);
  rot_x.SetPrompt("x-axis rotation");
  ArgSpecT rot_y(ArgSpecT::double_, "y");
  rot_y.SetDefault(0.0);
  rot_y.SetPrompt("y-axis rotation");
  ArgSpecT rot_z(ArgSpecT::double_, "z");
  rot_z.SetDefault(0.0);
  rot_z.SetPrompt("z-axis rotation");
  rotate.AddArgument(rot_x);
  rotate.AddArgument(rot_y);
  rotate.AddArgument(rot_z);
  iAddCommand(rotate);

  iAddCommand(CommandSpecT("Zoom"));
  iAddCommand(CommandSpecT("Change_background_color"));
  iAddCommand(CommandSpecT("ChangeDataColor"));
//  iAddCommand(CommandSpecT("Select_time_step"));
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
  	fRenderer->Delete();
  	
  	/* frame label */
	if (fLabelMapper) fLabelMapper->Delete();
	if (fLabelActor) fLabelActor->Delete();

	/* free subs - assuming all added with VTKFrameT::AddBody */
	const ArrayT<iConsoleObjectT*>& isubs = iSubs();
	for (int i = 0; i < isubs.Length(); i++)
		delete isubs[i];
}

void VTKFrameT::ResetView(void)
{
  fRenderer->GetActiveCamera()->SetFocalPoint(0,0,0);
  fRenderer->GetActiveCamera()->SetPosition(0,0,1);
  fRenderer->GetActiveCamera()->ComputeViewPlaneNormal();
  fRenderer->GetActiveCamera()->SetViewUp(0,1,0);
  fRenderer->GetActiveCamera()->OrthogonalizeViewUp();
  fRenderer->GetActiveCamera()->Zoom(0.85);
  fRenderer->ResetCamera();
}

bool VTKFrameT::AddBody(VTKBodyDataT* body_data)
{
	/* add body to the list */
	VTKBodyT* new_body = new VTKBodyT(this, body_data);
	if (bodies.AppendUnique(*new_body))
	{
		new_body->AddToFrame();
		iAddSub(*new_body);
		
		/* light back side for 2D */
		if (body_data->NumSD() == 2 && !fRenderer->GetTwoSidedLighting())
			fRenderer->TwoSidedLightingOn();
		
		ResetView();
		return true;
    }
	else
		return false;
}

bool VTKFrameT::RemoveBody(VTKBodyDataT* body_data)
{
  	/* look for body data */
	int index = -1;
	for (int i = 0; index == -1 && bodies.Length(); i++)
		if (bodies[i].BodyData() == body_data)
			index = i;

  if (index == -1)
  {
  	cout << "body not found" << endl;
    return false;
  }
  else
    {
      VTKBodyT& body = bodies[index];
      iDeleteSub(body);

      /* remove from fRenderer */
      body.RemoveFromFrame();
      ResetView();

      /* remove from body list */
      bodies.DeleteAt(index);
      return true;
    }
}

void VTKFrameT::ShowFrameLabel(const StringT& label)
{
	if (!fLabelMapper)
	{
		fLabelMapper = vtkTextMapper::New();
#if __DARWIN__
  		fLabelMapper->SetFontSize(11);
#else
		fLabelMapper->SetFontSize(16);
#endif
	}
	fLabelMapper->SetInput(label);
	
	if (!fLabelActor) {
		fLabelActor = vtkActor2D::New();
  		fLabelActor->SetMapper(fLabelMapper);
		fLabelActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
		fLabelActor->GetPositionCoordinate()->SetValue(0.4,.95);
		fLabelActor->GetProperty()->SetColor(0,1,0);
  		fRenderer->AddActor(fLabelActor);
  	}
}

/* remove the frame label */
void VTKFrameT::HideFrameLabel(void)
{
	if (fLabelActor)
	{
  		fRenderer->RemoveActor(fLabelActor);
  		fLabelActor->Delete();
  		fLabelActor = NULL;
		fLabelMapper->Delete();
		fLabelMapper = NULL;
	}
}

/* execute given command - returns false on fail */
bool VTKFrameT::iDoCommand(const CommandSpecT& command, StringT& line)
{
  int sfbTest;
  int  varNum;
  double xRot, yRot, zRot, zoom;
  double timeStep;

  if (command.Name() == "Update")
      return fConsole.iDoCommand(command, line);
  else if (command.Name() == "Interactive")
      return fConsole.iDoCommand(command, line);
  else if (command.Name() == "SelectTimeStep")
      return fConsole.iDoCommand(command, line);
  else if (command.Name() == "AddBody")
    {
      int body;
      command.Argument(0).GetValue(body);

      const ArrayT<VTKBodyDataT*>& all_bodies = fConsole.Bodies();
      if (body >= 0 && body < all_bodies.Length())
	{
	  if (AddBody(all_bodies[body]))
	    {
	      
	      Render();
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
	      Render();
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

  else if (command.Name() == "ResetView")
    {
      ResetView();
      fRenderer->ResetCamera();
      Render();
      return true;
    }
  	else if (command.Name() == "ShowNodeNumbers")
    {
    	if (bodies.Length() == 0)
    		return false;
    	else
    	{
    		/* command spec */
    		CommandSpecT* show = bodies[0].iCommand("ShowNodeNumbers");
    		if (!show)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels ON */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i].iDoCommand(*show, tmp);

			Render();
			return true;
    	}    
    }
	else if (command.Name() == "HideNodeNumbers")
    {
    	if (bodies.Length() == 0)
    		return false;
    	else
    	{
    		/* command spec */
    		CommandSpecT* hide = bodies[0].iCommand("HideNodeNumbers");
    		if (!hide)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels OFF */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i].iDoCommand(*hide, tmp);

			Render();
			return true;
    	}    
	}
  else if (command.Name() == "Color_bar_off")
    {
 //     fRenderer->RemoveActor(bodies[0]->SBActor());
      Render();
      return true;
    }
  else if (command.Name() == "Color_bar_on")
    {
//      fRenderer->AddActor(bodies[0]->SBActor());
      Render();

      return true;
    }
  else if (command.Name() == "ChangeDataColor")
    {
    	cout << "not updated" << endl;
#if 0
      int color;
      char line[255];
      cout << "Choose color: \n 1: Red\n 2: Green\n 3: Blue: ";
      cin >> color;
      cin.getline(line, 254);
      bodies[0]->ChangeDataColor(color);
      Render();
#endif
      return true;
    }
	else if (command.Name() == "Rotate")
	{
		double x, y, z;
		command.Argument("x").GetValue(x);
		command.Argument("y").GetValue(y);
		command.Argument("z").GetValue(z);

		if (fabs(x) > 1.0e-6) fRenderer->GetActiveCamera()->Elevation(x);
		if (fabs(y) > 1.0e-6) fRenderer->GetActiveCamera()->Azimuth(y);
		if (fabs(z) > 1.0e-6) fRenderer->GetActiveCamera()->Roll(z);

		fRenderer->GetActiveCamera()->OrthogonalizeViewUp();
		Render();
		return true;
	}
  else if (command.Name() == "Zoom")
    {

      cout << "Zoom by what magnitude?: ";
      cin >> zoom;
      char line[255];
      cin.getline(line, 254);
      fRenderer->GetActiveCamera()->Zoom(zoom);
      Render();
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
	fRenderer->SetBackground(0,0,0);
	break;
      case 2:
	fRenderer->SetBackground(1,1,1);
	break;
      case 3:
	fRenderer->SetBackground(1,0,0);
	break;
      case 4:
	fRenderer->SetBackground(0,1,0);
	break;
      default:
	fRenderer->SetBackground(0,0,1);
	break;
      }
      Render();
      return true;
    }

  else if (command.Name() == "Select_time_step")
    {
      int step;
      cout << "choose frame number from 0 to " << bodies[0]->NumTimeSteps()-1 <<" to be displayed: ";
      cin >> step;
      for (int i = 0; i < bodies.Length(); i++)
	     bodies[i]->SelectTimeStep(step);
      char line[255];
      cin.getline(line, 254);
      Render();
      return true; 
    }
  else if (command.Name() == "Choose_variable")
	{
		cout << "choose variable number from 0 to " << bodies[0]->NumNodeVariables()-1 <<" to be displayed\n";
		const ArrayT<StringT>& labels = bodies[0]->NodeLabels();
		for (int i = 0; i < labels.Length(); i++)
			cout << setw(5) << i << ": " << labels[i] << '\n';
		cout << " variable: ";
		cin >> varNum;
		char line[255];
		cin.getline(line, 254);
		const StringT& var = (bodies[0]->NodeLabels())[varNum];
		for (int i = 0; i < bodies.Length(); i++)
			bodies[i]->ChangeVars(var); // will not return true if body does not have the var
		Render();
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
      const ArrayT<VTKBodyDataT*>& bodies = fConsole.Bodies();
      out << "body to add (0," << bodies.Length()-1 << ")\n";
    }
  else if (command.Name() == "RemoveBody")
    {
      out << "body to remove (0," << bodies.Length()-1 << ")\n";
    }
  else /* inherited */
    iConsoleObjectT::ValuePrompt(command, index, out);
}

/* call to re-render the window contents */
void VTKFrameT::Render(void) const
{
	CommandSpecT* spec = fConsole.iCommand("Update");
	if (!spec)
	{
		cout << "VTKFrameT::Render: \"Update\" command not found" << endl;
		throw eGeneralFail;
	}
	StringT line;
	fConsole.iDoCommand(*spec, line);
}
