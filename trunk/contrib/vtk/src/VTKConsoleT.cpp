/* $Id: VTKConsoleT.cpp,v 1.9 2001-10-02 18:40:30 recampb Exp $ */

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
#include "vtkIdFilter.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkLabeledDataMapper.h"
#include "vtkActor2D.h"
#include "vtkFieldData.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"

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
  //iAddVariable("output_file", output_file);

  /* add console commands */
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Save");
  iAddCommand("Show_Node_Numbers");
  //iAddCommand("Hide_Node_Numbers");
  iAddCommand("Color_bar_on");
  iAddCommand("Color_bar_off");
  iAddCommand("X_axis_rotation");
  iAddCommand("Y_axis_rotation");
  iAddCommand("Z_axis_rotation");
  iAddCommand("Flip_book");
  iAddCommand("Change_background_color");
  iAddCommand("Select_frame_number");
 

   StringT file = "../../example_files/heat/heat.io0.exo";

   ExodusT exo(cout);
  if (!exo.OpenRead(file))
	{
	  cout << " ERROR: could not open file: " << file << endl;
	  throw;
	}
  else
	cout << "read database file: " << file << endl;


  /* read coordinates */
  int num_nodes = exo.NumNodes();
  int num_dim   = exo.NumDimensions();
  dArray2DT coordinates(num_nodes, num_dim);
 
  //ArrayT<dArray2DT> coordinates(num_nodes);
  exo.ReadCoordinates(coordinates);

  if (coordinates.MinorDim() == 2) 
    { 
      /* temp space */ 
      dArray2DT tmp(coordinates.MajorDim(), 3); 
      
      /* write in */ 
      tmp.BlockColumnCopyAt(coordinates, 0);    
      tmp.SetColumn(2, 0.0); 
      
      /* swap memory */ 
      tmp.Swap(coordinates); 
    } 


  /* read element block ID's */
  int num_elem_blocks = exo.NumElementBlocks();
  iArrayT element_ID(num_elem_blocks);
  exo.ElementBlockID(element_ID);

  /* read element connectivities */
  ArrayT<iArray2DT> connectivities(num_elem_blocks);
  for (int i = 0 ; i < num_elem_blocks; i++)
	{
	  /* read dimensions */
	  int num_elements, num_element_nodes;
	  exo.ReadElementBlockDims(element_ID[i], num_elements, num_element_nodes);

	  /* read connectivities */
	  connectivities[i].Allocate(num_elements, num_element_nodes);
	  GeometryT::CodeT geometry;
	  exo.ReadConnectivities(element_ID[i], geometry, connectivities[i]);


	}

/* look for results data */
  //int num_time_steps = exo.NumTimeSteps();
  num_time_steps = exo.NumTimeSteps();
  double time;
  //if (num_time_steps > 0)
  //	{
	  /* variables defined at the nodes */
	  int num_node_variables = exo.NumNodeVariables();
	  ArrayT<StringT> node_labels(num_node_variables);
	  exo.ReadNodeLabels(node_labels);
// 	  cout << " nodal variables:\n";
// 	  for (int i = 0; i < node_labels.Length(); i++)
// 		cout << node_labels[i] << '\n';

	  /* variables defined over the elements */
	  int num_element_variables = exo.NumElementVariables();
	  ArrayT<StringT> element_labels;
	  exo.ReadElementLabels(element_labels);
	 //  cout << " element variables:\n" << endl;
// 	  for (int i = 0; i < element_labels.Length(); i++)
// 		cout << element_labels[i] << '\n';
// 	  cout.flush();

	  /* read nodal data */
	  dArray2DT nodal_data(num_nodes, num_node_variables);
	  dArrayT ndata(num_nodes);

	 //  vtkFieldData *fd[100];


	  if (num_time_steps > 0)
	    {
	      for (int i = 0; i < num_time_steps; i++)
		{
		  
		  exo.ReadTime(i+1, time);
		  scalars[i] =  vtkScalars::New(VTK_DOUBLE);
		  
		  /* loop over variables */
		  for (int j = 0; j < num_node_variables; j++)
		    {
		      exo.ReadNodalVariable(i+1, j+1, ndata);
		      nodal_data.SetColumn(j, ndata);
		    }

// 	      fd[i] = vtkFieldData::New();
// 	      fd[i]->SetArray(num_nodes, nodal_data);
		  for (int k = 0; k<num_nodes; k++) scalars[i]->InsertScalar(k+1, nodal_data[k]);
// 		  cout << " time: " << time << endl;
// 		  cout << " nodal data:\n" << nodal_data << endl;
		}
	
	    }



//   vtkPoints *points = vtkPoints::New();
//  for (int i=0; i<num_nodes; i++) points->InsertPoint(i,coordinates[i]);
		  
  vtkPoints *points = vtkPoints::New();
 for (int i=0; i<num_nodes; i++) points->InsertPoint(i+1,coordinates(i));

 //NOTE: this code is only for a single block of cells
 //      the data for visualization will be provided one
 //      group of cells at a time.
 iArray2DT& connects = connectivities[0];

 /* create array of VTK-style connectivities */
 iArray2DT vtk_connects(connects.MajorDim(), connects.MinorDim()+1); //has 1 extra entry!!!
 vtk_connects.BlockColumnCopyAt(connects, 1);
 vtk_connects.SetColumn(0, connects.MinorDim()); //first value in each row is row size 

 /* release "ownership" of memory */
 int* p_vtk_connects;
 vtk_connects.ReleasePointer(&p_vtk_connects);

 /* create VTK integer array */
 vtkIntArray* intArray = vtkIntArray::New();
 intArray->SetNumberOfComponents(vtk_connects.MinorDim()); //is this needed???
 intArray->SetArray(p_vtk_connects, vtk_connects.Length(), 0);

  
 ugrid = vtkUnstructuredGrid::New();


 /* create VTK array of cells */
 vtkCellArray* vtk_cell_array = vtkCellArray::New();
 vtk_cell_array->SetCells(vtk_connects.MajorDim(), intArray);

 //NOTE: do all at once for higher efficiency 
 //ugrid->Allocate(num_elem_blocks);
 // for (int i=0; i<num_elem_blocks; i++) 
 //   ugrid->InsertNextCell(VTK_QUAD, 4, connectivities[i]); 
 

 //NOTE: the example database has triangles. Generally, you would need
 //      to determine the cell type from the database, or require that
 //      the cell type be specified when the data is sent to the visualizer.
 ArrayT<int> cell_types(vtk_connects.MajorDim());
 cell_types = VTK_QUAD; // all the cells are the same

 /* insert cells in the grid */
 //  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
 ugrid->SetCells(cell_types.Pointer(), vtk_cell_array);
 


  ugrid->SetPoints(points);
  points->Delete();
   ugrid->GetPointData()->SetScalars(scalars[0]);
  //ugrid->SetFieldData(fd[0]);

  // source_file = "../../example_files/heat/data_file1.vtk";
  output_file = "A.tif";
  numColors = 256;
  hueRange1 = 0; hueRange2 = 0.6667;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
  scalarRange1 = 6; scalarRange2 = 17;
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
  ids = vtkIdFilter::New();
  visPts = vtkSelectVisiblePoints::New();
  ldm = vtkLabeledDataMapper::New();
  pointLabels = vtkActor2D::New();
  
  renWin->AddRenderer(renderer);
 
  iren->SetRenderWindow(renWin);

  lut->SetHueRange(hueRange1, hueRange2);
  lut->SetSaturationRange(satRange1,satRange2);
  lut->SetValueRange(valRange1,valRange2);
  lut->SetAlphaRange(alphaRange1,alphaRange2);
  lut->SetNumberOfColors(numColors);
  lut->Build();

//   // x,y,z axes
//       vtkCubeAxesActor2D *axes = vtkCubeAxesActor2D::New();
//       axes->SetInput(ugrid);
//       axes->SetCamera(renderer->GetActiveCamera());
//     //  axes->SetLabelFormat("%6.4g");
//       axes->ShadowOn();
//       axes->SetFlyModeToOuterEdges();
//       axes->SetFontFactor(0.8);
//       axes->GetProperty()->SetColor(1,0,1);
  
  ugridMapper->SetInput(ugrid);
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
  scalarBar->SetTitle("Temperature for time .5000");
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
  scalarBar->SetOrientationToHorizontal();
  scalarBar->SetWidth(0.8);
  scalarBar->SetHeight(0.17);
  renderer->AddActor(scalarBar);
 
}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{

int xDir, yDir, zDir;
double xRot, yRot, zRot;

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

    renWin->Render();
    iren->Start();
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
      renWin->Render();
      iren->Start();
      return true;
    }

  else if (command == "Save")
    {
      cout << "Enter name for file to be saved to: ";
      cin >> output_file;
      char line[255];
      cin.getline(line, 254);
	
      renSrc->SetInput(renderer);
       renSrc->WholeWindowOn();
      writer->SetInput(renSrc->GetOutput());
      writer->SetFileName(output_file);
      writer->Write();
      return true;
    }

  else if (command == "Show_Node_Numbers")
    {
       ids->SetInput(ugrid);
      ids->PointIdsOn();
      ids->CellIdsOn();
      ids->FieldDataOn();
      visPts->SetInput(ids->GetOutput());
      visPts->SetRenderer(renderer);
      //visPts->SelectionWindowOn();
      //  visPts->SetSelection();
      ldm->SetInput(visPts->GetOutput());
      // ldm->SetlabelFormat();
      ldm->SetLabelModeToLabelFieldData();
      pointLabels->SetMapper(ldm);
      renderer->AddActor2D(pointLabels);
      renWin->Render();
      iren->Start();
      return true;      

    }


  else if (command == "Color_bar_off")
    {
      renderer->RemoveActor(scalarBar);
      renWin->Render();
      iren->Start();
      return true;
    }

  else if (command == "Color_bar_on")
    {
      renderer->AddActor(scalarBar);
      renWin->Render();
      iren->Start();
      return true;
    }


  else if (command == "X_axis_rotation")
    {
      cout << "Rotate:\n 1: Counter-clockwise\n 2: Clockwise: ";
      cin >> xDir;
      char line[255];
      cin.getline(line, 254);
      cout << "Rotate by how many degrees?: ";
      cin >> xRot;
      cin.getline(line, 254);
      if (xDir == 1) renderer->GetActiveCamera()->Elevation(xRot);
      else if (xDir == 2) renderer->GetActiveCamera()->Elevation(-xRot);
      renWin->Render();
      iren->Start();
      return true;
    }


  else if (command == "Y_axis_rotation")
    {
      cout << "Rotate:\n 1: Counter-clockwise\n 2: Clockwise: ";
      cin >> yDir;
      char line[255];
      cin.getline(line, 254);
      cout << "Rotate by how many degrees?: ";
      cin >> yRot;
      cin.getline(line, 254);
      if (yDir == 1) renderer->GetActiveCamera()->Azimuth(yRot);
      else if (yDir == 2) renderer->GetActiveCamera()->Azimuth(-yRot);
      renWin->Render();
      iren->Start();
      return true;
    }

  else if (command == "Z_axis_rotation")
    {
      cout << "Rotate:\n 1: Counter-clockwise\n 2: Clockwise: ";
      cin >> zDir;
      char line[255];
      cin.getline(line, 254);
      cout << "Rotate by how many degrees?: ";
      cin >> zRot;
      cin.getline(line, 254);
      if (zDir == 1) renderer->GetActiveCamera()->Roll(zRot);
      else if (zDir == 2) renderer->GetActiveCamera()->Roll(-zRot);
      renWin->Render();
      iren->Start();
      return true;
    }

  else if (command == "Flip_book")
    {
      time_t start_time, cur_time;
      double timeStep;
      cout << "Enter time step in seconds: ";
      cin >> timeStep;
      char line[255];
      cin.getline(line, 254);
      

      for (int j = 0; j<num_time_steps; j++){
// 	time(&start_time);
// 	  do {
// 	    time(&cur_time);
// 	    ugrid->GetPointData()->SetScalars(scalars[j]);
// 	    renWin->Render();
// 	  } while((cur_time - start_time) < timeStep);
      
      
        clock_t start_time, cur_time;
         start_time = clock();
         while((clock() - start_time) < timeStep * CLOCKS_PER_SEC)
         {
         }
// 	 renSrc->SetInput(renderer);
// 	 renSrc->WholeWindowOn();
// 	 writer->SetInput(renSrc->GetOutput());
// 	 outFileName = j + ".tif";
// 	 writer->SetFileName(outFileName);
// 	 writer->Write();
	 ugrid->GetPointData()->SetScalars(scalars[j]);
	 renWin->Render();

      }

//       for (int j=2; j<12; j++){
// 	time(&start_time);    
// 	if (j % 2 == 0) {
// 	  do {
// 	    time(&cur_time);
// 	    ugrid->GetPointData()->SetScalars(scalars[0]);
// 	    renWin->Render();
// 	  } while((cur_time - start_time) < timeStep);
// 	}
// 	else {
// 	  do {
// 	    time(&cur_time);
// 	    ugrid->GetPointData()->SetScalars(scalars[1]);
// 	    renWin->Render();
// 	  } while((cur_time - start_time) < timeStep);
// 	}
//       }

      renWin->Render();
      iren->Start();
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
      renWin->Render();
      iren->Start();
      return true;
    }

else if (command == "Select_frame_number")
  {
    int frameNum;
    cout << "choose frame number from 0 to " << num_time_steps-1 <<" to be displayed: ";
    cin >> frameNum;
    char line[255];
    cin.getline(line, 254);
    ugrid->GetPointData()->SetScalars(scalars[frameNum]);
    renWin->Render();
      iren->Start();
      return true;

  }

//    else if (command == "Hide_Node_Numbers")

//     {

// //       ids->PointIdsOff();
// //       ids->Update();

//       renderer->removeActor((vtkActor)pointLabels);
//       renWin->Render();
//       iren->Start();
//       return true;   
 
//     }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
