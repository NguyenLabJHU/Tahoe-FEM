/* $Id: VTKConsoleT.cpp,v 1.19 2001-10-25 21:40:19 recampb Exp $ */

#include "VTKConsoleT.h"
#include "VTKFrameT.h"
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

VTKConsoleT::VTKConsoleT(void)
{
  /* set console name */
  iSetName("vtk");

  /* add console commands */
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
 
  iAddCommand("Save");
  iAddCommand("Save_flip_book_images");
  iAddCommand("Show_Node_Numbers");
  iAddCommand("Hide_Node_Numbers");
  iAddCommand("Color_bar_on");
  iAddCommand("Color_bar_off");
  iAddCommand("X_axis_rotation");
  iAddCommand("Y_axis_rotation");
  iAddCommand("Z_axis_rotation");
  iAddCommand("Flip_book");
  iAddCommand("Change_background_color");
  iAddCommand("Select_frame_number");
  iAddCommand("Show_axes");
  iAddCommand("Hide_axes");
  iAddCommand("Choose_variable");
 
 


 /* prompt for input file */
  // StringT file;
  char line[255];
    cout << "Choose number of File\n 1: heat.io0.exo\n 2: test.io0.exo\n 3: test2.io0.exo\n 4: big.exo\n 5: enter different .exo: ";
    cin >> test;
    cin.getline(line, 254);
    
    if (test == 1) inFile = "../../example_files/heat/heat.io0.exo";
    else if (test ==2) inFile ="test.io0.exo";
    else if (test == 3) inFile = "test2.io0.exo";
    else if (test == 4) inFile = "big.exo";
    else if (test == 5) {
      cout << "Enter file name with .exo: ";
      cin >> inFile;
      cin.getline(line, 254);
    }
    else cout << "bad entry";
  
 
  //TEMP - adding sub-scopes to the console
  fFrames.Allocate(4);
  for (int i = 0; i < 4; i++)
    {
      StringT temp;
      temp.Append("frame",i,2);
      fFrames[i].iSetName(temp);
      iAddSub(fFrames[i]);
    }

	  
	  fFrames[0].renderer = vtkRenderer::New();

// 	  scalarBar = vtkScalarBarActor::New();
	  //writer = vtkTIFFWriter::New();
	  fFrames[0].renSrc = vtkRendererSource::New();
// 	  ids = vtkIdFilter::New();
// 	  visPts = vtkSelectVisiblePoints::New();
// 	  ldm = vtkLabeledDataMapper::New();
// 	  pointLabels = vtkActor2D::New();
// 	  cam = vtkCamera::New();
// 	  axes = vtkCubeAxesActor2D::New();
	 
	  	  
// 	  fFrames[0].renderer->AddActor(ugridActor);
// 	  fFrames[0].renderer->SetBackground(0,0,0);

// 	  scalarBar->SetLookupTable(ugridMapper->GetLookupTable());
	  
// 	  sbTitle.Append(node_labels[currentVarNum]); 
// 	  sbTitle.Append(" for frame 000");
// 	  // sbTitle.Append(time,5);
// 	  scalarBar->SetTitle(sbTitle);
// 	  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
// 	  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
// 	  scalarBar->SetOrientationToHorizontal();
// 	  scalarBar->SetWidth(0.8);
// 	  scalarBar->SetHeight(0.17);
	//   fFrames[0].renderer->AddActor(scalarBar);
	  
	  

  cout << "How many plots in window (1 or 4)?";
  cin >> numRen;
  cin.getline(line, 254);

  renWin = vtkRenderWindow::New();
  iren = vtkRenderWindowInteractor::New();

  writer = vtkTIFFWriter::New();
  //  renSrc = vtkRendererSource::New();
  /* divide window into 4 parts */
  //     if (numRen ==4){
	         fFrames[0].renderer->SetViewport(0,0,.5,.5);
	         fFrames[1].renderer->SetViewport(.5,0,1,.5);
	         fFrames[2].renderer->SetViewport(0,.5,.5,1);
	         fFrames[3].renderer->SetViewport(.5,.5,1,1);
	         fFrames[0].renderer->GetActiveCamera()->Zoom(0.85);
	         fFrames[1].renderer->GetActiveCamera()->Zoom(0.85);
	         fFrames[2].renderer->GetActiveCamera()->Zoom(0.85);
		 fFrames[3].renderer->GetActiveCamera()->Zoom(0.85);
	  //       }
	  
	  
	  //   renWin->SetPosition(668, 0);
	  //   renWin->SetSize(600,700);
	


  renWin->AddRenderer(fFrames[0].renderer);
  // if (numRen==4){   
    renWin->AddRenderer(fFrames[1].renderer);
    renWin->AddRenderer(fFrames[2].renderer);
    renWin->AddRenderer(fFrames[3].renderer);
  // }
  iren->SetRenderWindow(renWin);
  

  renWin->SetPosition(668, 0);
  renWin->SetSize(600,700);


}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{


  if (command == "Start_Rendering")
  {
    renWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;
    iren->Start();
    return true;
  }



  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
