/* $Id: VTKConsoleT.cpp,v 1.66 2004-04-16 01:09:32 paklein Exp $ */
#include "VTKConsoleT.h"

/* ANSI headers */
//#include <iostream.h>
//#include <iomanip.h>
//#include <cstdio>

/* tahoe toolbox headers */
#include "ExceptionCodes.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"
#include "IOBaseT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"

/* VTK headers */
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyle.h"
#include "vtkRendererSource.h"
#include "vtkTIFFWriter.h"
#include "vtkJPEGWriter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPointPicker.h"
#include "vtkFloatArray.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSphereSource.h"
#include "vtkCellPicker.h"
#include "vtkTextSource.h"
#include "vtkActorCollection.h"
#include "vtkRendererCollection.h"
#include "vtkProperty.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"

/* VTK console headers */
#include "VTKFrameT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

/* generating PostScript from OpenGL */
#include "gl2ps.h"

using namespace Tahoe;

AutoArrayT<vtkActor*> VTKConsoleT::pickedPoints;

VTKConsoleT::VTKConsoleT(const ArrayT<StringT>& arguments):
  fArguments(arguments),
  fBodyCount(0),
  fRenderHold(false)
{
  /* set console name */
  iSetName("vtk");

  /* add console commands */
  iAddCommand(CommandSpecT("Interactive"));
  iAddCommand(CommandSpecT("Update"));

  CommandSpecT addbody("AddBody");
  ArgSpecT file(ArgSpecT::string_);
  file.SetPrompt("path to database file");
  ArgSpecT file_format(ArgSpecT::string_);
  file_format.SetPrompt("database file format");
  file_format.SetDefault("auto");
  addbody.AddArgument(file);
  addbody.AddArgument(file_format);
  iAddCommand(addbody);

  CommandSpecT removebody("RemoveBody");
  removebody.SetPrompter(this);
  ArgSpecT bodynum(ArgSpecT::int_);
  bodynum.SetPrompt("body number to remove");
  removebody.AddArgument(bodynum);
  iAddCommand(removebody);
  
  CommandSpecT layout("Layout", false);
  ArgSpecT nx(ArgSpecT::int_, "nx");
  nx.SetPrompt("number of horizontal frames");
  ArgSpecT ny(ArgSpecT::int_, "ny");
  ny.SetPrompt("number of vertical frames");
  layout.AddArgument(nx);
  layout.AddArgument(ny);
  iAddCommand(layout);

  CommandSpecT window_size("WindowSize", false);
  ArgSpecT width(ArgSpecT::int_, "w");
  width.SetPrompt("width");
  ArgSpecT height(ArgSpecT::int_, "h");
  height.SetPrompt("height");
  window_size.AddArgument(width);
  window_size.AddArgument(height);
  iAddCommand(window_size);
  
  iAddCommand(CommandSpecT("ShowFrameNumbers"));
  iAddCommand(CommandSpecT("HideFrameNumbers"));

	CommandSpecT save("Save");
	ArgSpecT save_file(ArgSpecT::string_);
	save_file.SetPrompt("image file name");
	ArgSpecT save_format(ArgSpecT::string_);
	save_format.SetPrompt("image file format (TIFF|JPG|PS)");
	save_format.SetDefault("JPG");
	save.AddArgument(save_file);
	save.AddArgument(save_format);
	iAddCommand(save);
	
	CommandSpecT flipbook("FlipBook", false);
	ArgSpecT delay(ArgSpecT::double_, "delay");
	delay.SetDefault(0.0);
	delay.SetPrompt("frame delay in seconds");
	ArgSpecT step(ArgSpecT::int_, "step");
	step.SetDefault(1);
	step.SetPrompt("time step increment");
	flipbook.AddArgument(delay);
	flipbook.AddArgument(step);
	iAddCommand(flipbook);

	CommandSpecT save_flip("SaveFlipBook");
	save_flip.AddArgument(save_file);
	save_flip.AddArgument(save_format);
	save_flip.AddArgument(step);
	iAddCommand(save_flip);

	CommandSpecT select_step("SelectTimeStep");
	select_step.SetPrompter(this);
	ArgSpecT stepnum(ArgSpecT::int_);
	stepnum.SetPrompt("time step");
	select_step.AddArgument(stepnum);
	iAddCommand(select_step);

	CommandSpecT next_step("NextTimeStep");
	next_step.SetPrompter(this);
	ArgSpecT steps_ahead(ArgSpecT::int_);
	steps_ahead.SetDefault(1);
	steps_ahead.SetPrompt("number steps forwards");
	next_step.AddArgument(steps_ahead);
	iAddCommand(next_step);

	CommandSpecT prev_step("PreviousTimeStep");
	prev_step.SetPrompter(this);
	ArgSpecT steps_back(ArgSpecT::int_);
	steps_back.SetDefault(1);
	steps_back.SetPrompt("number steps backwards");
	prev_step.AddArgument(steps_back);
	iAddCommand(prev_step);

	/* display objects */
	renWin = vtkRenderWindow::New();
	renWin->SetWindowName("VTK for Tahoe");
	int x_position = 10;
	int y_position = 10;
	int index, start = 0;
	if (CommandLineOption("-x", index, start) && fArguments.Length() >= index+2) {
		StringT pos = fArguments[index+1];
		x_position = atoi(pos);
	}
	if (CommandLineOption("-y", index, start) && fArguments.Length() >= index+2) {
		StringT pos = fArguments[index+1];
		y_position = atoi(pos);
	}
	renWin->SetPosition(x_position, y_position);

	//renWin->SetDesiredUpdateRate(30.0);
	renWin->SetSize(500,500);
	iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	pointPicker = vtkPointPicker::New();
	pointPicker->SetTolerance(0.01);

	cellPicker = vtkCellPicker::New();
	cellPicker->SetTolerance(0.01);

	
	//iren->SetPicker(cellPicker);
	//iren->SetEndPickMethod(PickCells,(void *)iren);  

	iren->SetPicker(pointPicker);
	iren->SetEndPickMethod(PickPoints,(void *)iren);  
	//iren->SetDesiredUpdateRate(30.0);
	
	/* set interator style to trackball instead of joystick in the
	 * same way it occurs from the window */
	//iren->GetInteractorStyle()->OnChar(0, 0, 't', 1);
	//NOTE: doesn't work with VTK 4.1.1

	/* set up single frame */
	SetFrameLayout(1,1);

	/* borrow some console commands */
	CommandSpecT* command;
	command = fFrames[0]->iCommand("Rotate");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);

	command = fFrames[0]->iCommand("Zoom");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);

	command = fFrames[0]->iCommand("Pan");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);

	command = fFrames[0]->iCommand("ResetView");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	
	command = fFrames[0]->iCommand("ShowNodeNumbers");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	command = fFrames[0]->iCommand("HideNodeNumbers");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	command = fFrames[0]->iCommand("ShowElementNumbers");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	command = fFrames[0]->iCommand("HideElementNumbers");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	command = fFrames[0]->iCommand("ShowColorBar");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	command = fFrames[0]->iCommand("HideColorBar");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);

	/* representation */
  	iAddCommand(CommandSpecT("Wire"));
  	iAddCommand(CommandSpecT("Surface"));
  	iAddCommand(CommandSpecT("Point"));
// 	iAddCommand(CommandSpecT("ShowContours"));
// 	iAddCommand(CommandSpecT("HideContours"));
// 	iAddCommand(CommandSpecT("ShowCuttingPlane"));
// 	iAddCommand(CommandSpecT("HideCuttingPlane"));

	/* look for command line options */
	index, start = 0;
	while (CommandLineOption("-f", index, start))
	{
		StringT path = fArguments[index+1];
		path.ToNativePathName();
		if (!AddBody("auto", path))
			cout << "could not add body: " << path << endl;
		start = index+1;
	}

	/* draw */
	if (fBodies.Length() > 0) renWin->Render();
}

/* destructor*/
VTKConsoleT::~VTKConsoleT(void)
{
  /* free data for remaining bodies */
  for (int i = 0; i < fBodies.Length(); i++)
	{
	  delete fBodies[i];
	  fBodies[i] = NULL;
	}

  /* free data for remaining bodies */
  for (int i = 0; i < fFrames.Length(); i++)
	{
	  delete fFrames[i];
	  fFrames[i] = NULL;
	}
}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const CommandSpecT& command, StringT& line)
{
  if (command.Name() == "Interactive")
    {
      if (pickedPoints.Length()>0){
	vtkRendererCollection* renderers = renWin->GetRenderers();
	renderers->InitTraversal();
	vtkRenderer* temp;

	for (int j = 0; j < renderers->GetNumberOfItems(); j++){
	  temp = renderers->GetNextItem();
	  for (int i =0; i< pickedPoints.Length(); i++){
	    temp->RemoveActor(pickedPoints[i]);  
	  }  
	  
	}      
      }
      


      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      cout << "pick point with mouse and type 'p' to display scalar value" << endl;
      iren->Start();
      return true;
    }  
  else if (command.Name() == "Update")
	{
		if (!fRenderHold)
		{
			/* update bodies */
			for (int i = 0; i < fBodies.Length(); i++)
				fBodies[i]->UpdateData();

			/* update frames */
			for (int i = 0; i < fFrames.Length(); i++)
				fFrames[i]->UpdateData();

			renWin->Render();
		}
		return true;
    }
  else if (command.Name() == "Rotate" || 
           command.Name() == "Zoom"   || 
           command.Name() == "Pan"    || 
           command.Name() == "ResetView" ||
		   command.Name() == "ShowNodeNumbers" ||
	       command.Name() == "HideNodeNumbers" ||
		   command.Name() == "ShowElementNumbers" ||
		   command.Name() == "HideElementNumbers" ||
		   command.Name() == "ShowColorBar" ||
		   command.Name() == "HideColorBar")
	{
		/* hold rendering */
		fRenderHold = true;
	
		/* rotate all frames */
		bool OK = true;
		for (int i = 0; OK && i < fFrames.Length(); i++)
			OK = fFrames[i]->iDoCommand(command, line);
			
		/* render */
		fRenderHold = false;
		StringT tmp;
		iDoCommand(*iCommand("Update"), tmp);
		return OK;
	}
	else if (command.Name() == "SelectTimeStep")	
	{
		/* no bodies */
		if (fBodies.Length() == 0)
		{
			cout << "no bodies" << endl;
			return false;
		}
		else
		{
 			fRenderHold = true;
 	 		int old_step = fBodies[0]->CurrentStepNumber();
 			int step;
 			command.Argument(0).GetValue(step);
 			bool OK = true;
			for (int i = 0; OK && i < fBodies.Length(); i++)
				OK = fBodies[i]->SelectTimeStep(step);
 			fRenderHold = false;
	
			if (OK)
			{
				StringT tmp;
				iDoCommand(*iCommand("Update"), tmp);
				return true;
			}
  			else
  			{
  				/* reset */
				for (int i = 0; i < fBodies.Length(); i++)
					fBodies[i]->SelectTimeStep(old_step);
  				return false;
  			}
  		}
  } 
  else
    if (command.Name() == "NextTimeStep")	
	{
		/* no bodies */
		if (fBodies.Length() == 0)
		{
			cout << "no bodies" << endl;
			return false;
		}
		else
		{
 			fRenderHold = true;
 	 		int old_step = fBodies[0]->CurrentStepNumber();
 			int step;
 			command.Argument(0).GetValue(step);
			step += old_step;
 			bool OK = true;
			for (int i = 0; OK && i < fBodies.Length(); i++)
				OK = fBodies[i]->SelectTimeStep(step);
 			fRenderHold = false;
	
			if (OK)
			{
				StringT tmp;
				iDoCommand(*iCommand("Update"), tmp);
				return true;
			}
  			else
  			{
  				/* reset */
				for (int i = 0; i < fBodies.Length(); i++)
					fBodies[i]->SelectTimeStep(old_step);
  				return false;
  			}
  		}
  }
  else if (command.Name() == "PreviousTimeStep")	
	{
		/* no bodies */
		if (fBodies.Length() == 0)
		{
			cout << "no bodies" << endl;
			return false;
		}
		else
		{
 			fRenderHold = true;
 	 		int old_step = fBodies[0]->CurrentStepNumber();
 			int step;
 			command.Argument(0).GetValue(step);
			step = old_step - step;
 			bool OK = true;
			for (int i = 0; OK && i < fBodies.Length(); i++)
				OK = fBodies[i]->SelectTimeStep(step);
 			fRenderHold = false;
	
			if (OK)
			{
				StringT tmp;
				iDoCommand(*iCommand("Update"), tmp);
				return true;
			}
  			else
  			{
  				/* reset */
				for (int i = 0; i < fBodies.Length(); i++)
					fBodies[i]->SelectTimeStep(old_step);
  				return false;
  			}
  		}
  }
  else 
    if (command.Name() == "Wire")
  {
	bool OK = true;
  	for (int i = 0; OK && i < fBodies.Length(); i++)
  	{
  		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
  		if (!comm) return false;
  		OK = fBodies[i]->iDoCommand(*comm, line);
  	}
  	return OK;
  }  
  else if (command.Name() == "Surface")
  {
	bool OK = true;
  	for (int i = 0; OK && i < fBodies.Length(); i++)
  	{
  		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
  		if (!comm) return false;
  		OK = fBodies[i]->iDoCommand(*comm, line);
  	}
  	return OK;
  }
  else if (command.Name() == "Point")
  {
	bool OK = true;
  	for (int i = 0; OK && i < fBodies.Length(); i++)
  	{
  		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
  		if (!comm) return false;
  		OK = fBodies[i]->iDoCommand(*comm, line);
  	}
  	return OK;
  }

//   else if (command.Name() == "ShowContours")
//   {
// 	bool OK = true;
//   	for (int i = 0; OK && i < fBodies.Length(); i++)
//   	{
//   		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
//   		if (!comm) return false;
//   		OK = fBodies[i]->iDoCommand(*comm, line);
//   	}
//   	return OK;
//   }

//   else if (command.Name() == "HideContours")
//   {
// 	bool OK = true;
//   	for (int i = 0; OK && i < fBodies.Length(); i++)
//   	{
//   		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
//   		if (!comm) return false;
//   		OK = fBodies[i]->iDoCommand(*comm, line);
//   	}
//   	return OK;
//   }

//   else if (command.Name() == "ShowCuttingPlane")
//   {
// 	bool OK = true;
//   	for (int i = 0; OK && i < fBodies.Length(); i++)
//   	{
//   		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
//   		if (!comm) return false;
//   		OK = fBodies[i]->iDoCommand(*comm, line);
//   	}
//   	return OK;
//   }

//   else if (command.Name() == "HideCuttingPlane")
//   {
// 	bool OK = true;
//   	for (int i = 0; OK && i < fBodies.Length(); i++)
//   	{
//   		const CommandSpecT* comm = fBodies[i]->iResolveCommand(command.Name(), line);
//   		if (!comm) return false;
//   		OK = fBodies[i]->iDoCommand(*comm, line);
//   	}
//   	return OK;
//   }


  else if (command.Name() == "Layout")
	{
		int nx, ny;
		command.Argument("nx").GetValue(nx);	
		command.Argument("ny").GetValue(ny);	
	  	SetFrameLayout(nx, ny);
	  	return true;
	}
  else if (command.Name() == "AddBody")
    {
      StringT path, format;
      command.Argument(0).GetValue(path);
      command.Argument(1).GetValue(format);
      path.ToNativePathName();	  
      return AddBody(format, path);
    }
  else if (command.Name() == "RemoveBody")
    {
      int index;
      command.Argument(0).GetValue(index);
      if (index < 0 && index >= fBodies.Length())
	return false;
      else {
	
	/* remove body from all frames */
	int count = 0;
	for (int i = 0; i < fFrames.Length(); i++)
	  if (fFrames[i]->RemoveBody(fBodies[index])) count++;
	cout << "body " << index << " removed from " << count << " frames" << endl;
	
	/* remove from base console */
	iDeleteSub(*fBodies[index]);

	/* free memory */
	delete fBodies[index];
		
	/* resize array */
	fBodies.DeleteAt(index);

	/* free */
	return true;
      }
    }
  else if (command.Name() == "FlipBook")
	{
    	if (fBodies.Length() == 0) return false;
    
		double delay;
		int step;
		command.Argument("delay").GetValue(delay);
		command.Argument("step").GetValue(step);
		
		/* the select time step command */
		StringT empty_line;
		CommandSpecT* set_step = iCommand("SelectTimeStep");
		if (!set_step) throw eGeneralFail;
      
		/* assume all the bodies have the same number of steps as body 0 */
		int j;
		for (j = 0; j<fBodies[0]->NumTimeSteps(); j += step){

			/* time delay */
			clock_t start_time, cur_time;
			start_time = clock();
			while((clock() - start_time) < delay * CLOCKS_PER_SEC) { }

			/* set step */
			set_step->Argument(0).SetValue(j);
			iDoCommand(*set_step, empty_line);
      	}
      	
      	/* make sure to hit last frame */
      	if (j != fBodies[0]->NumTimeSteps())
      	{
			/* time delay */
			clock_t start_time, cur_time;
			start_time = clock();
			while((clock() - start_time) < delay * CLOCKS_PER_SEC) { }
      	
			/* get last step number */
      		j = fBodies[0]->NumTimeSteps() - 1;

			/* set step */
			set_step->Argument(0).SetValue(j);
			iDoCommand(*set_step, empty_line);
      	}
      	
      return true;
    }
  else if (command.Name() == "SaveFlipBook")
    {
    	if (fBodies.Length() == 0) return false;
    
      StringT name, format;
      command.Argument(0).GetValue(name);
      command.Argument(1).GetValue(format);
      int step;
      command.Argument(2).GetValue(step);
      step = (step < 1) ? 1 : step;
            
      vtkRendererSource* image = vtkRendererSource::New();
      image->SetInput(fFrames[0]->Renderer());
      image->WholeWindowOn();

		/* set writer (set file name extension) */
		StringT ext;
		ext.Suffix(name, '.');
		ext.ToUpper();
		vtkImageWriter* writer;
		StringT suffix;
		if (format == "TIFF" || format == "TIF")
		{
			/* construct TIFF writer */
			writer = vtkTIFFWriter::New();

			/* add extension */
			if (ext != ".TIF" && ext != ".TIFF") suffix.Append(".tif");
		}
		else /* use jpg for all others */
		{
			/* construct JPEG writer */
			writer = vtkJPEGWriter::New();
			
			/* add extension */
			if (ext != ".JPG" && ext != ".JPEG") suffix.Append(".jpg");
		}
		writer->SetInput(image->GetOutput());
      
		/* assume all the bodies have the same number of steps as body 0 */
		for (int j = 0; j < fBodies[0]->NumTimeSteps(); j += step) {
	
			/* update bodies */
			for (int i = 0; i < fBodies.Length(); i++)
				fBodies[i]->SelectTimeStep(j);

			/* update frames */
			for (int i = 0; i < fFrames.Length(); i++)
				fFrames[i]->UpdateData();
	
			renWin->Render();  
			StringT frame_name = name;
			frame_name.Append(j,3); // pad to a width of 3 digits
			frame_name.Append(suffix);
			writer->SetFileName(frame_name);
			writer->Write();
			cout << frame_name << " has been saved" << endl;
      }
      
      /* clean up */
      writer->Delete();
      image->Delete();
      
      cout << "Flip book images have been saved." << endl;
      renWin->Render();
      return true;
    }
	else if (command.Name() == "ShowFrameNumbers")
	{
		fRenderHold = true;
		for (int j = 0; j < fFrames.Length(); j++)	  
			fFrames[j]->ShowFrameLabel(fFrames[j]->iName());

		fRenderHold = false;
		StringT tmp;
		iDoCommand(*iCommand("Update"), tmp);
		return true;
	}
	else if (command.Name() == "HideFrameNumbers")
	{
		fRenderHold = true;
		for (int j = 0; j < fFrames.Length(); j++)	  
			fFrames[j]->HideFrameLabel();

		fRenderHold = false;
		StringT tmp;
		iDoCommand(*iCommand("Update"), tmp);
		return true;
	}
	else if (command.Name() == "Save")
	{
		/* arguments */
		StringT name, format;
		command.Argument(0).GetValue(name);
		command.Argument(1).GetValue(format);

		/* PostScipt uses gl2ps library */
		format.ToUpper();
		if (format == "PS")
		{
			/* file name */
			StringT ext;
			ext.Suffix(name, '.');
			if (ext != ".ps" && ext != ".PS")
				name.Append(".ps");
				
			/* open stream */
			FILE* fp = fopen(name, "w");
			int buffsize = 0;
			int state = GL2PS_OVERFLOW;
			
			/* loop on buffer size */
			int iter = 0;
			while (state == GL2PS_OVERFLOW && iter++ < 20)
			{
				buffsize += 1024*1024;
			
				/* start redirect */
				gl2psBeginPage(name, "VTK_for_Tahoe", GL2PS_PS, GL2PS_BSP_SORT,
					GL2PS_SIMPLE_LINE_OFFSET, GL_RGBA, 0, NULL, 
		      		buffsize, fp, NULL);
		      		
		      	/* draw */
		      	renWin->Render();
		      	
		      	/* end redirect */
		      	state = gl2psEndPage();
			}
			
			/* close stream */
			fclose(fp);
		}
		else /* VTK image writers */
		{
			/* window to image filter */
			vtkRendererSource* image = vtkRendererSource::New();
			image->SetInput(fFrames[0]->Renderer());
			image->Update();
			image->WholeWindowOn();
			
			/* file name extension */
			StringT ext;
			ext.Suffix(name, '.');
			ext.ToUpper();
			vtkImageWriter* writer;
			if (format == "TIFF" || format == "TIF")
			{
				/* construct TIFF writer */
				writer = vtkTIFFWriter::New();

				/* add extension */
      			if (ext != ".TIF" && ext != ".TIFF")
      				name.Append(".tif");
			}
			else if (format == "JPG" || format == "JPEG")
			{
				/* construct JPEG writer */
				writer = vtkJPEGWriter::New();
			
				/* add extension */
      			if (ext != ".JPG" && ext != ".JPEG")
      				name.Append(".jpg");
			}
			else
			{
				cout << "unrecognized file format: " << format << endl;
				return false;
			}
      
      		/* write image */
			writer->SetInput(image->GetOutput());      
			writer->SetFileName(name);
			writer->Write();
      
            /* clean up */
			writer->Delete();
			image->Delete();
			//renWin->Render();
		}
      
      cout << name << " has been saved" << endl;   
      return true;
    }
	else if (command.Name() == "WindowSize") {
	
		int w, h;
		command.Argument("w").GetValue(w);
		command.Argument("h").GetValue(h);
		
		/* limit minimum window size */
		h = (h < 100) ? 100 : h;
		w = (w < 100) ? 100 : w; 
		
		/* set size */	
		renWin->SetSize(w, h);

		/* update */
		StringT tmp;
		return iDoCommand(*iCommand("Update"), tmp);
	}


  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}

/**********************************************************************
* Protected
**********************************************************************/

/* write prompt for the specific argument of the command */
void VTKConsoleT::ValuePrompt(const CommandSpecT& command, int index, 
	ostream& out) const
{
	if (command.Name() == "SelectTimeStep")
	{
		if (fBodies.Length() > 0)
			out << "range 0 to " << fBodies[0]->NumTimeSteps() - 1 << '\n';
	}
	else if (command.Name() == "RemoveBody")
	{
		out << "range 0 to " << fBodies.Length() - 1 << '\n';
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* construct body from the given file path */
bool VTKConsoleT::AddBody(const StringT& format, const StringT& file)
{
  /* temp */
  VTKBodyDataT* body;

  /* try to construct body */
  try {
  
  	IOBaseT::FileTypeT file_format = IOBaseT::kExodusII;
  	if (format == "auto")
  		file_format = IOBaseT::name_to_FileTypeT(file);
  
	body = new VTKBodyDataT(file_format, file);
	StringT name;
	name.Append(fBodyCount++);
	name.Append(".body");
	body->iSetName(name);
	fBodies.Append(body);

	/* add bodies to frame 0 by default */
	fFrames[0]->AddBody(body);
	
	/* add to base console */
	iAddSub(*body);
  }
  catch (int) {
	cout << "\n exception constructing body from file: " << file << endl;
//	delete body;
	return false;
  }
  /* OK */
  return true;
}

/* reset the frame layout */
void VTKConsoleT::SetFrameLayout(int num_x, int num_y)
{
  /* remove all frames from console and window*/
  for (int i = 0; i < fFrames.Length(); i++)
	{
	  iDeleteSub(*fFrames[i]);
	  renWin->RemoveRenderer(fFrames[i]->Renderer());
	}

  /* no less than one in each direction */
  num_x = (num_x < 1) ? 1 : num_x;
  num_y = (num_y < 1) ? 1 : num_y;

  /* temp space for new layout */
  Array2DT<VTKFrameT*> new_frames(num_y, num_x);
  new_frames = NULL;

  /* copy in old frames */
  for (int i = 0; i < new_frames.MajorDim() && i < fFrames.MajorDim(); i++)
	for (int j = 0; j < new_frames.MinorDim() && j < fFrames.MinorDim(); j++)
	  {
		new_frames(i,j) = fFrames(i,j);
		fFrames(i,j) = NULL;
	  }

  /* delete any extra frames */
  for (int i = 0; i < fFrames.Length(); i++)
	delete fFrames[i];

  /* swap */
  new_frames.Swap(fFrames);

  /* set up frames */
  double dx = 1.0/fFrames.MinorDim();
  double dy = 1.0/fFrames.MajorDim();
  for (int i = 0; i < fFrames.MajorDim(); i++)
	for (int j = 0; j < fFrames.MinorDim(); j++)
	  {
		/* need a new one */
		if (fFrames(i,j) == NULL) fFrames(i,j) = new VTKFrameT(*this);

		/* name */
		StringT name;
		name.Append(j);
		name.Append(".", i);
		name.Append(".frame");
		fFrames(i,j)->iSetName(name);

		/* set port location/size */
		int ii = fFrames.MajorDim() - 1 - i;
		fFrames(i,j)->Renderer()->SetViewport(j*dx, ii*dy, (j+1)*dx, (ii+1)*dy);

		/* add to window */
		renWin->AddRenderer(fFrames(i,j)->Renderer());

		/* add to console */
		iAddSub(*fFrames(i,j));
	  }
}

void VTKConsoleT::PickCells(void *arg)
{
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor *)arg;
  vtkCellPicker *cellPicker = (vtkCellPicker *)iren->GetPicker(); 
  if (cellPicker->GetCellId() != -1){
    cout << "Cell: " << cellPicker->GetCellId()+1 << endl;
    //  cout << "Value: " << (cellPicker->GetDataSet()->GetCellData()->GetScalars()->GetComponent(cellPicker->GetCellId(),0)) << endl;
  }
  else
    cout << "Invalid Cell" << endl;
  
}

/* prints out picked point and scalar value */
void VTKConsoleT::PickPoints(void *arg)
{
  vtkRenderWindowInteractor *iren = (vtkRenderWindowInteractor *)arg;
  vtkPointPicker *pointPicker = (vtkPointPicker *)iren->GetPicker();
  vtkPolyDataMapper* sphereMapper = vtkPolyDataMapper::New();
  vtkActor* sphereActor = vtkActor::New();  
  vtkSphereSource *sphere = vtkSphereSource::New();
  sphere->SetThetaResolution(8); sphere->SetPhiResolution(8);
  float* bounds = pointPicker->GetDataSet()->GetBounds();
  sphere->SetRadius(.01*(bounds[1]-bounds[0]));
  sphereMapper->SetInput(sphere->GetOutput());
  sphereActor->SetMapper(sphereMapper);
  sphereActor->GetProperty()->SetColor(1,1,1);
  sphereActor->VisibilityOn();
  sphereActor->PickableOff();

  	if (pointPicker->GetPointId() != -1) {
  		
  		vtkDataSet* data_set = pointPicker->GetDataSet();
		if (data_set != NULL ) 
		{
      		float* coords = data_set->GetPoint(pointPicker->GetPointId());
			sphereActor->SetPosition(coords);
			pointPicker->GetRenderer()->AddActor(sphereActor);
			pickedPoints.Append(sphereActor);
			iren->GetRenderWindow()->Render();
    
 			/* report id and coordinates */
			cout << "      Point: " << pointPicker->GetPointId()+1 << '\n';
			cout << "Coordinates: " << "(" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << '\n';

			vtkDataArray* scalars = data_set->GetPointData()->GetScalars();
			if (scalars != NULL)
			{
      			float* values = scalars->GetTuple(pointPicker->GetPointId()); 
				int num_values = scalars->GetNumberOfComponents(); 
				cout << "      Value: " << (scalars->GetComponent(pointPicker->GetPointId(), 0)) << '\n';
			}
		}
		else cout << "      Point: " << pointPicker->GetPointId()+1 << '\n';
	}
	else cout << "Invalid Point" << '\n';  

	cout << endl;
}

/* returns the index of the requested option */
bool VTKConsoleT::CommandLineOption(const char* str, int& index, int start) const
{
	int dex = -1;
	for (int i = start; dex == -1 && i < fArguments.Length(); i++)
		if (fArguments[i] == str)
			dex = i;

	if (dex == -1)
		return false;
	else
	{
		index = dex;
		return true;
	}
}
