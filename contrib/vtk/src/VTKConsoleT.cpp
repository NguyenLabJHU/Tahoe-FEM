/* $Id: VTKConsoleT.cpp,v 1.26 2001-11-08 00:42:35 paklein Exp $ */

#include "VTKConsoleT.h"
#include "VTKFrameT.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRendererSource.h"
#include "vtkTIFFWriter.h"
#include "vtkScalarBarActor.h"
#include "vtkDataSetMapper.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"
#include "VTKBodyT.h"

VTKConsoleT::VTKConsoleT(const ArrayT<StringT>& arguments):
  fArguments(arguments)
{
  /* set console name */
  iSetName("vtk");

  /* add console commands */
  iAddCommand("Interactive");
  iAddCommand("Update");
  iAddCommand("AddBody");
  iAddCommand("RemoveBody");
  iAddCommand("ResetView");
  iAddCommand("Layout");
  iAddCommand("Show_Frame_Numbers");
 
  iAddCommand("Save_flip_book_images");
  iAddCommand("Flip_book");

  /* display objects */
  renWin = vtkRenderWindow::New();
  renWin->SetPosition(668, 0);
  renWin->SetSize(600,700);
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
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{
  if (command == "Interactive")
    {
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }  
  else if (command == "Layout")
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
  else if (command == "Update")
    {
	  for (int i = 0; i < fBodies.Length(); i++)
		fBodies[i]->UpdateData();
      renWin->Render();
      return true;
    }
  else if (command == "ResetView")
	{
	  for (int i = 0; i < fFrames.Length(); i++)
		fFrames[i]->ResetView();
	  renWin->Render();
	  return true;
	}
  else if (command == "AddBody")
	{
	  StringT path;
	  cout << "path to data file: ";
	  cin >> path;
	  Clean(cin);
	  path.ToNativePathName();
	  
	  return AddBody(path);
	}
  else if (command == "RemoveBody")
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
  else if (command == "Flip_book")
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
  else if (command== "Save_flip_book_images")
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
      
      /* window to image filter */
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
  
  else if (command == "Show_Frame_Numbers")
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
  //  else if (command == "Reset_to_Default_Values")
  //{
  //  fBodies[0]->DefaultValues();
  //  fBodies[0]->UpdateData();
  //  renWin->Render();
  //  // iren->Start();
  //  return true;
  //}
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
//   else if (command == "Choose_variable")
//     {
//       char line[255];
//       cout << "choose variable number from 0 to " << num_node_variables-1 <<" to be displayed\n" << varList;      
//       cin >> currentVarNum;
//       cin.getline(line, 254);
//       ugrid->GetPointData()->SetScalars(scalars[frameNum][currentVarNum]);
//       ugridMapper->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
//       if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
// 	ugrid->GetPointData()->SetVectors(vectors[frameNum][currentVarNum]);
//       sbTitle = "";
//       sbTitle.Append(node_labels[currentVarNum]); 
//       sbTitle.Append(" for frame ");
//       sbTitle.Append(frameNum,3);
//       scalarBar->SetTitle(sbTitle);
//       renWin->Render();
//       return true;
//     }
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
  VTKBodyT* body;

  /* try to construct body */
  try {
	body = new VTKBodyT(file);
	fBodies.Append(body);
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
