/* $Id: VTKFrameT.cpp,v 1.20 2002-01-16 00:33:57 cjkimme Exp $ */

#include "VTKFrameT.h"
#include "VTKConsoleT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"
#include "VTKUGridT.h"

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
//#include "vtkVectors.h"
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
	fLabelActor(NULL),
	scalarBar(NULL)
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

	com_spec = fConsole.iCommand("FlipBook");
	if (!com_spec) throw eGeneralFail;
	iAddCommand(*com_spec);

	com_spec = fConsole.iCommand("Save");
	if (!com_spec) throw eGeneralFail;
	iAddCommand(*com_spec);

  /* own console commands */
//  iAddCommand(CommandSpecT("Reset_to_Default_Values"));
  iAddCommand(CommandSpecT("ResetView"));
  iAddCommand(CommandSpecT("ShowNodeNumbers"));
  iAddCommand(CommandSpecT("HideNodeNumbers"));
  iAddCommand(CommandSpecT("ShowAxes"));
  iAddCommand(CommandSpecT("HideAxes"));

	iAddCommand(CommandSpecT("ShowColorBar"));
	iAddCommand(CommandSpecT("HideColorBar"));

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

  CommandSpecT pan("Pan", false);
  ArgSpecT pan_x(ArgSpecT::double_, "x");
  pan_x.SetDefault(0.0);
  pan_x.SetPrompt("x-direction pan (dx/w)");
  ArgSpecT pan_y(ArgSpecT::double_, "y");
  pan_y.SetDefault(0.0);
  pan_y.SetPrompt("y-direction pan (dy/h)");
  pan.AddArgument(pan_x);
  pan.AddArgument(pan_y);
  iAddCommand(pan);

  CommandSpecT zoom("Zoom");
  ArgSpecT factor(ArgSpecT::double_);
  factor.SetPrompt("zoom factor");
  zoom.AddArgument(factor);
  iAddCommand(zoom);
  
  CommandSpecT change_bg("ChangeBackgroundColor");
  change_bg.SetPrompter(this);
  ArgSpecT bg_color(ArgSpecT::int_);
  bg_color.SetPrompt("background color");
  change_bg.AddArgument(bg_color);
  iAddCommand(change_bg);

//  iAddCommand(CommandSpecT("ChangeDataColor"));

  CommandSpecT choosevar("ChooseVariable");
  choosevar.SetPrompter(this);
  ArgSpecT varnum(ArgSpecT::int_);
  varnum.SetPrompt("variable number");
  choosevar.AddArgument(varnum);
  iAddCommand(choosevar);

  ArgSpecT body_num(ArgSpecT::int_);
  body_num.SetPrompt("body number");

  CommandSpecT add_body("AddBody");
  add_body.AddArgument(body_num);
  add_body.SetPrompter(this);
  iAddCommand(add_body);

  CommandSpecT rem_body("RemoveBody");
  add_body.SetPrompter(this);
  rem_body.AddArgument(body_num);
  iAddCommand(rem_body);
}

/* destructor */
VTKFrameT::~VTKFrameT(void)
{
	/* free display classes */
  	fRenderer->Delete();
  	
  	/* scalar bar */
  	if (scalarBar) scalarBar->Delete();
  	
  	/* frame label */
	if (fLabelMapper) fLabelMapper->Delete();
	if (fLabelActor) fLabelActor->Delete();
	
	/* free bodies */
	for (int i = 0; i < bodies.Length(); i++)
		delete bodies[i];
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
	if (bodies.AppendUnique(new_body))
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
		if (bodies[i]->BodyData() == body_data)
			index = i;

  if (index == -1)
  {
  	cout << "body not found" << endl;
    return false;
  }
  else
    {
      VTKBodyT* body = bodies[index];
      iDeleteSub(*body);

      /* remove from fRenderer */
      body->RemoveFromFrame();
      ResetView();

      /* remove from body list */
      bodies.DeleteAt(index);
      
      /* free body */
      delete body;
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
  else if (command.Name() == "FlipBook")
      return fConsole.iDoCommand(command, line);
  else if (command.Name() == "Save")
      return fConsole.iDoCommand(command, line);
  else if (command.Name() == "ShowColorBar")
  {
  	if (bodies.Length() == 0)
		return false;
	else if (scalarBar)
	{
		cout << "hide scalar bar first" << endl;
		return false;
	}
	else
  	{
  		/* use the first body */
  		VTKBodyDataT* body_data = bodies[0]->BodyData();
  		
  		/* get grids */
  		const ArrayT<VTKUGridT*>& ugrids = body_data->UGrids();
  		if (ugrids.Length() == 0) return false;
  	
		/* new scalar bar */
  		scalarBar = vtkScalarBarActor::New();

		/* get lut from first grid */
		scalarBar->SetLookupTable(ugrids[0]->GetLookupTable());

		const StringT& var_name = (body_data->NodeLabels())[body_data->CurrentVariableNumber()];
		scalarBar->SetTitle(var_name);
		
		scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
		scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
		scalarBar->SetOrientationToHorizontal();
		scalarBar->SetWidth(0.8); 
		scalarBar->SetHeight(0.17);
		
		fRenderer->AddActor(scalarBar);
		Render();
		return true;
  	}
  } 
  else if (command.Name() == "HideColorBar")
  {
	if (!scalarBar)
	{
		cout << "show scalar bar first" << endl;
		return false;
	}
	else
  	{
		fRenderer->RemoveActor(scalarBar);		
		scalarBar->Delete();
		scalarBar = NULL;
		Render();
		return true;
	}  
  }      
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
	  if (RemoveBody(bodies[body]->BodyData()))
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
    		CommandSpecT* show = bodies[0]->iCommand("ShowNodeNumbers");
    		if (!show)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels ON */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i]->iDoCommand(*show, tmp);

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
    		CommandSpecT* hide = bodies[0]->iCommand("HideNodeNumbers");
    		if (!hide)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels OFF */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i]->iDoCommand(*hide, tmp);

			Render();
			return true;
    	}    
	}
  	else if (command.Name() == "ShowAxes")
    {
    	if (bodies.Length() == 0)
    		return false;
    	else
    	{
    		/* command spec */
    		CommandSpecT* show = bodies[0]->iCommand("ShowAxes");
    		if (!show)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels ON */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i]->iDoCommand(*show, tmp);

			Render();
			return true;
    	}    
    }
	else if (command.Name() == "HideAxes")
    {
    	if (bodies.Length() == 0)
    		return false;
    	else
    	{
    		/* command spec */
    		CommandSpecT* hide = bodies[0]->iCommand("HideAxes");
    		if (!hide)
    		{
    			cout << "command not found" << endl;
    			return false;
    		}
    	
    		/* labels OFF */
    		StringT tmp;
    		for (int i = 0; i < bodies.Length(); i++)
    			bodies[i]->iDoCommand(*hide, tmp);

			Render();
			return true;
    	}    
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
		if (fabs(y) > 1.0e-6) fRenderer->GetActiveCamera()->Azimuth(-y);
		if (fabs(z) > 1.0e-6) fRenderer->GetActiveCamera()->Roll(z);

		fRenderer->GetActiveCamera()->OrthogonalizeViewUp();
		fRenderer->ResetCameraClippingRange(); //to avoid clipping with rotation
		Render();
		return true;
	}
  else if (command.Name() == "Zoom")
    {
    	double factor;
    	command.Argument(0).GetValue(factor);
      fRenderer->GetActiveCamera()->Zoom(factor);
      Render();
      return true;
    }
    /* based on vtkInteractorStyle::PanCamera */
	else if (command.Name() == "Pan")
	{
		double x, y;
		command.Argument("x").GetValue(x);
		command.Argument("y").GetValue(y);
		
		/* renderer viewport size */
		int* size = fRenderer->GetSize();
  
		/* calculate the focal depth */
		double ViewFocus[4];
		fRenderer->GetActiveCamera()->GetFocalPoint(ViewFocus);
		ComputeWorldToDisplay(ViewFocus[0], ViewFocus[1], ViewFocus[2], ViewFocus);
		double focalDepth = ViewFocus[2];

		double NewPickPoint[4];
		ComputeDisplayToWorld(ViewFocus[0] + size[0]*x/2.0, 
		                      ViewFocus[1] + size[1]*y/2.0, 
		                      focalDepth, NewPickPoint);

		/* get the current focal point and position */
		fRenderer->GetActiveCamera()->GetFocalPoint(ViewFocus);
		double *ViewPoint = fRenderer->GetActiveCamera()->GetPosition();

  		/* compute a translation vector */
		double MotionVector[3];
		double a = 1.0;
		MotionVector[0] = a*(ViewFocus[0] - NewPickPoint[0]);
		MotionVector[1] = a*(ViewFocus[1] - NewPickPoint[1]);
		MotionVector[2] = a*(ViewFocus[2] - NewPickPoint[2]);

		fRenderer->GetActiveCamera()->SetFocalPoint(
			MotionVector[0] + ViewFocus[0],
			MotionVector[1] + ViewFocus[1],
			MotionVector[2] + ViewFocus[2]);
		fRenderer->GetActiveCamera()->SetPosition(
			MotionVector[0] + ViewPoint[0],
			MotionVector[1] + ViewPoint[1],
			MotionVector[2] + ViewPoint[2]);

//TEMP - don't know what this does			
#if 0
		/* get render window interactor */
		vtkRenderWindow* rw = fRenderer->GetRenderWindow();
  		vtkRenderWindowInteractor *rwi = rw->GetInteractor();
		if (rwi->GetLightFollowCamera())
		{
			/* get the first light */
			fRenderer->CurrentLight->SetPosition(this->CurrentCamera->GetPosition());
			fRenderer->CurrentLight->SetFocalPoint(this->CurrentCamera->GetFocalPoint());
		}
#endif
		Render();
		return true;
	}
	else if (command.Name() =="ChangeBackgroundColor")
	{
		int bg_color;
		command.Argument().GetValue(bg_color);
            
		switch (bg_color) {
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
			case 5:
				fRenderer->SetBackground(0,0,1);
				break;
      		default:
      			return false;
		}
		Render();
		return true;
	}
  else if (command.Name() == "ChooseVariable")
	{
		int varNum;
		command.Argument(0).GetValue(varNum);
		const StringT& var = (bodies[0]->BodyData()->NodeLabels())[varNum];
		for (int i = 0; i < bodies.Length(); i++)
			bodies[i]->ChangeVars(var); // will not return true if body does not have the var
		
		/* reset color bar name */
		if (scalarBar)	
		{
  			VTKBodyDataT* body_data = bodies[0]->BodyData();
			const StringT& var_name = (body_data->NodeLabels())[body_data->CurrentVariableNumber()];
			scalarBar->SetTitle(var_name);
		}
			
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
  else if (command.Name() == "ChooseVariable")
    {
		const ArrayT<StringT>& labels = bodies[0]->BodyData()->NodeLabels();
		for (int i = 0; i < labels.Length(); i++)
			out << setw(5) << i << ": " << labels[i] << '\n';
    }
  else if (command.Name() == "ChangeBackgroundColor")
  {
	out << "choose background color:\n 1: black\n 2: white\n 3: red\n 4: green\n 5: blue\n";
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

// Description:
// transform from display to world coordinates.
// WorldPt has to be allocated as 4 vector
void VTKFrameT::ComputeDisplayToWorld(double x, double y, double z,
	double *worldPt)
{
	fRenderer->SetDisplayPoint(x, y, z);
	fRenderer->DisplayToWorld();
	fRenderer->GetWorldPoint(worldPt);
	if (worldPt[3])
    {
		worldPt[0] /= worldPt[3];
		worldPt[1] /= worldPt[3];
		worldPt[2] /= worldPt[3];
		worldPt[3] = 1.0;
	}
}


// Description:
// transform from world to display coordinates.
// displayPt has to be allocated as 3 vector
void VTKFrameT::ComputeWorldToDisplay(double x, double y, double z,
	double *displayPt)
{
	fRenderer->SetWorldPoint(x, y, z, 1.0);
	fRenderer->WorldToDisplay();
	fRenderer->GetDisplayPoint(displayPt);
}
