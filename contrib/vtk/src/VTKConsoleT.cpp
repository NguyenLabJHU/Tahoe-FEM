/* $Id: VTKConsoleT.cpp,v 1.31 2001-12-03 21:59:22 paklein Exp $ */

#include "VTKConsoleT.h"
#include "VTKFrameT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRendererSource.h"
#include "vtkTIFFWriter.h"
#include "vtkScalarBarActor.h"
#include "vtkDataSetMapper.h"
#include "vtkRenderLargeImage.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"
#include "VTKBodyT.h"
#include "VTKBodyDataT.h"

#include "CommandSpecT.h"
#include "ArgSpecT.h"

VTKConsoleT::VTKConsoleT(const ArrayT<StringT>& arguments):
  fArguments(arguments)
{
  /* set console name */
  iSetName("vtk");

  /* add console commands */
  iAddCommand(CommandSpecT("Interactive"));
  iAddCommand(CommandSpecT("Update"));
  iAddCommand(CommandSpecT("AddBody"));
  iAddCommand(CommandSpecT("RemoveBody"));
  iAddCommand(CommandSpecT("ResetView"));
  iAddCommand(CommandSpecT("Layout"));
  iAddCommand(CommandSpecT("Show_Frame_Numbers"));
  iAddCommand(CommandSpecT("Save"));
  iAddCommand(CommandSpecT("Large_Save"));
  iAddCommand(CommandSpecT("Save_flip_book_images"));
  iAddCommand(CommandSpecT("Flip_book"));

  /* display objects */
  renWin = vtkRenderWindow::New();
  renWin->SetWindowName("VTK for Tahoe");
//  renWin->SetPosition(668, 0);
//  renWin->SetSize(600,700);
  iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  /* set up single frame */
  SetFrameLayout(1,1);

  /* draw */
  renWin->Render();
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
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }  
  else if (command.Name() == "Layout")
	{
	  int num_x = -1, num_y = -1;
	  cout << "horizontal: ";
	  cin >> num_x;
	  cout << "vertical: ";
	  cin >> num_y;
	  Clean(cin);
	  SetFrameLayout(num_x, num_y);
	  return true;
	}
  else if (command.Name() == "Update")
    {
      for (int i = 0; i < fBodies.Length(); i++)
	fBodies[i]->UpdateData();
      renWin->Render();
      return true;
    }
  else if (command.Name() == "ResetView")
    {
      for (int i = 0; i < fFrames.Length(); i++)
	fFrames[i]->ResetView();
      renWin->Render();
      return true;
	}
  else if (command.Name() == "AddBody")
    {
      StringT path;
      cout << "path to data file: ";
      cin >> path;
      Clean(cin);
      path.ToNativePathName();
	  
      return AddBody(path);
    }
  else if (command.Name() == "RemoveBody")
    {
      int index;
      cout << "select body to remove (0," << fBodies.Length()-1 << "): ";
      cin >> index;
      Clean(cin);
      if (index < 0 && index >= fBodies.Length())
	return false;
      else {
	
	/* remove body from all frames */
	int count = 0;
	for (int i = 0; i < fFrames.Length(); i++)
	  if (fFrames[i]->RemoveBody(fBodies[index])) count++;
	cout << "body " << index << " removed from " << count << " frames" << endl;
	
	/* free */
	delete fBodies[index];
	
	/* resize array */
	fBodies.DeleteAt(index);
	return true;
      }
    }
  else if (command.Name() == "Flip_book")
    {
      double timeStep;
      cout << "Enter time step in seconds: ";
      cin >> timeStep;
	  Clean(cin);
      cout << "Show images at: \n 1: current view\n 2: default view: ";
     
      int sfbTest;
      cin >> sfbTest;
      Clean(cin);
      
      /* if default camera desired */
      if (sfbTest == 2) {
	for (int i = 0; i < fFrames.Length(); i++)
	  fFrames[i]->ResetView();
	renWin->Render();
      }	
      
      /* assume all the bodies have the same number of steps as body 0 */
      for (int j = 0; j<fBodies[0]->num_time_steps; j++){
        /* time delay */
	clock_t start_time, cur_time;
	start_time = clock();
	while((clock() - start_time) < timeStep * CLOCKS_PER_SEC)
	  {
	  }
	for (int i = 0; i < fBodies.Length(); i++)
	  fBodies[i]->SelectTimeStep(j);
	renWin->Render();
      }
      return true;
    }
  else if (command.Name() == "Save_flip_book_images")
    {
      StringT fbName;
      cout << "Enter name for flipbook to be saved (without .tif extension): ";
      cin >> fbName;
      Clean(cin);
      cout << "Save images at: \n 1: current view\n 2: default view: ";
      int sfbTest;
      cin >> sfbTest;
      Clean(cin);
      
      /* if default camera desired */
      if (sfbTest == 2) {
	for (int i = 0; i < fFrames.Length(); i++)
	  fFrames[i]->ResetView();
	renWin->Render();
      }	
      
      vtkRendererSource* image = vtkRendererSource::New();
      image->SetInput(fFrames[0]->Renderer());
      image->WholeWindowOn();
      
      /* construct TIFF writer */
      vtkTIFFWriter* writer = vtkTIFFWriter::New();
      writer->SetInput(image->GetOutput());
      
      /* assume all the bodies have the same number of steps as body 0 */
      for (int j = 0; j<fBodies[0]->num_time_steps; j++){
	
	for (int i = 0; i < fBodies.Length(); i++)
	  fBodies[i]->SelectTimeStep(j);
	
	renWin->Render();  
	StringT name = fbName;
	name.Append(j,3); // pad to a width of 3 digits
	name.Append(".tif");
	writer->SetFileName(name);
	writer->Write();
	cout << name << " has been saved" << endl;
      }
      
      /* clean up */
      writer->Delete();
      image->Delete();
      
      cout << "Flip book images have been saved." << endl;
      renWin->Render();
      return true;
    }
  
  else if (command.Name() == "Show_Frame_Numbers")
  {

    for (int i = 0; i < fFrames.MajorDim(); i++)
      for (int j = 0; j < fFrames.MinorDim(); j++)
	{
	  /* name */
	  StringT name = "frame";
	  name.Append(".",i);
	  name.Append(".",j);
	  
	  fFrames(i,j)->ShowFrameNum(name);
	}
    return true;
  }
  //  else if (command.Name() == "Reset_to_Default_Values")
  //{
  //  fBodies[0]->DefaultValues();
  //  fBodies[0]->UpdateData();
  //  renWin->Render();
  //  // iren->Start();
  //  return true;
  //}
  else if (command.Name() == "Save")
    {
      StringT fbName;
      cout << "Enter name for image to be saved (without .tif extension): ";
      cin >> fbName;
      Clean(cin);
      cout << "Save image at: \n 1: current view\n 2: default view: ";
      int sfbTest;
      cin >> sfbTest;
      Clean(cin);
      
      /* if default camera desired */
      if (sfbTest == 2) {
	for (int i = 0; i < fFrames.Length(); i++)
	  fFrames[i]->ResetView();
	renWin->Render();
      }	
      
      /* window to image filter */
      vtkRendererSource* image = vtkRendererSource::New();
      image->SetInput(fFrames[0]->Renderer());
	  image->WholeWindowOn();
      
      /* construct TIFF writer */
      vtkTIFFWriter* writer = vtkTIFFWriter::New();
      writer->SetInput(image->GetOutput());
      
      StringT name = fbName;
      name.Append(".tif");
      writer->SetFileName(name);
      writer->Write();
      cout << name << " has been saved" << endl;
   
      /* clean up */
      writer->Delete();
      image->Delete();
      renWin->Render();
      return true;
    }

  else if (command.Name() == "Large_Save")
    {

      StringT Name;
      cout << "Enter name for image to be saved (without .tif extension): ";
      cin >> Name;
      Clean(cin);
      Name.Append(".tif");
      vtkRenderLargeImage* renderLarge = vtkRenderLargeImage::New();
      renderLarge->SetInput(fFrames[0]->Renderer());
      renderLarge->SetMagnification(5);
      renderLarge->Update();
      vtkTIFFWriter* writer = vtkTIFFWriter::New();
      writer->SetInput(renderLarge->GetOutput());
      writer->SetFileName(Name);
      writer->Write();
      return true;

    }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}

/**********************************************************************
* Private
**********************************************************************/

/* construct body from the given file path */
bool VTKConsoleT::AddBody(const StringT& file)
{
  /* temp */
  VTKBodyDataT* body;

  /* try to construct body */
  try {
	body = new VTKBodyDataT(file);
	fBodies.Append(body);

	/* add bodies to frame 0 by default */
	fFrames[0]->AddBody(body);
  }
  catch (int) {
	cout << "\n exception constructing body from file: " << file << endl;
	delete body;
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
		if (fFrames(i,j) == NULL) fFrames(i,j) = new VTKFrameT;

		/* name */
		StringT name = "frame";
		name.Append(".",i);
		name.Append(".",j);
		fFrames(i,j)->iSetName(name);

		/* connect */
		fFrames(i,j)->setRenWin(renWin);
		fFrames(i,j)->setIren(iren);
		fFrames(i,j)->setConsole(this);

		/* set port location/size */
		fFrames(i,j)->Renderer()->SetViewport(j*dx, i*dy, (j+1)*dx, (i+1)*dy);

		/* add to window */
		renWin->AddRenderer(fFrames(i,j)->Renderer());

		/* add to console */
		iAddSub(*fFrames(i,j));
	  }
}
