/* $Id: VTKConsoleT.cpp,v 1.16 2001-10-23 00:23:08 recampb Exp $ */

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
#include "vtkCamera.h"
#include "vtkWarpVector.h"
#include "vtkVectors.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include <GeometryT.h> 

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
  // iAddVariable("min_Scalar_Range", scalarRange1);
  //iAddVariable("max_Scalar_Range", scalarRange2);
  iAddVariable("numColors", numColors);
  iAddVariable("source_file", source_file);
  //iAddVariable("output_file", output_file);
   iAddVariable("scale_factor", scale_factor);
   
  /* add console commands */
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Reset_view");
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
  StringT file;
  cout << "Choose number of File\n 1: heat.io0.exo\n 2: test.io0.exo\n 3: test2.io0.exo\n 4: big.exo\n 5: enter different .exo: ";
  cin >> test;
  char line[255];
  cin.getline(line, 254);
  if (test == 1) file = "../../example_files/heat/heat.io0.exo";
  else if (test ==2) file ="test.io0.exo";
  else if (test == 3) file = "test2.io0.exo";
  else if (test == 4) file = "big.exo";
  else if (test == 5) {
    cout << "Enter file name with .exo: ";
    cin >> file;
    cin.getline(line, 254);
  }
  else cout << "bad entry";

  cout << "How many plots in window (1 or 4)?";
  cin >> numRen;
  cin.getline(line, 254);

  /* read exodus file */
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

  // ArrayT<dArray2DT> coordinates(num_time_steps);
  exo.ReadCoordinates(coordinates);
  // exo.ReadCoordinates(coordinates(0));

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
  GeometryT::CodeT geometry;
  for (int i = 0 ; i < num_elem_blocks; i++)
	{
	  /* read dimensions */
	  int num_elements, num_element_nodes;
	  exo.ReadElementBlockDims(element_ID[i], num_elements, num_element_nodes);

	  /* read connectivities */
	  connectivities[i].Allocate(num_elements, num_element_nodes);
	  
	  exo.ReadConnectivities(element_ID[i], geometry, connectivities[i]);
	
	}

/* look for results data */
  //int num_time_steps = exo.NumTimeSteps();
  num_time_steps = exo.NumTimeSteps();
  // double time;
  //if (num_time_steps > 0)
  //	{
	  /* variables defined at the nodes */
	  num_node_variables = exo.NumNodeVariables();
	  // ArrayT<StringT> node_labels(num_node_variables);
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
	  // scalarRange1 = 10000;
	  // scalarRange2 = -10000;
	  if (num_time_steps > 0)
	    {
	      for (int i = 0; i < num_time_steps; i++)
		{
		  exo.ReadTime(i+1, time);
		  for (int j = 0; j < num_node_variables; j++)
		    {
		      exo.ReadNodalVariable(i+1, j+1, ndata);
		      nodal_data.SetColumn(j, ndata);
		      scalars[i][j] =  vtkScalars::New(VTK_DOUBLE);
		    
		      /* instantiate displacement vector if needed */
		      if (node_labels[0] == "D_X" || node_labels[0] == "D_Y" || node_labels[0] == "D_Z")
			vectors[i][j] = vtkVectors::New(VTK_DOUBLE);
		      
		      /* initialize min and max scalar range values */
		      scalarRange1[j] = 10000;
		      scalarRange2[j] = -10000;
		      
		      for (int k = 0; k<num_nodes; k++) {
			
			//	scalars[i][j]->InsertScalar(k+1, nodal_data[k][j]);
			/* determine min and max scalar range values */
			if (nodal_data(k,j) < scalarRange1[j]) scalarRange1[j] = nodal_data(k,j);
			if (nodal_data(k,j) > scalarRange2[j]) scalarRange2[j] = nodal_data(k,j);
			/* insert scalar value at each node for each variable and time step */
			scalars[i][j]->InsertScalar(k+1, nodal_data(k,j));
			//InsertVector(k,...)?????
			/* if displacement vector needed then insert vector at each node for each time step */
			if (node_labels[0] == "D_X" && node_labels[1] == "D_Y" && node_labels[2] == "D_Z")              
			  vectors[i][j]->InsertVector(k+1, nodal_data(k,0),nodal_data(k,1),nodal_data(k,2));
			else if (node_labels[0] == "D_X" && node_labels[1] == "D_Y") 
			   vectors[i][j]->InsertVector(k+1, nodal_data(k,0), nodal_data(k,1),0);
			else if (node_labels[0] == "D_X")
			  vectors[i][j]->InsertVector(k+1, nodal_data(k,0),0,0);
			else;
			
		      }
		    }
// 		  cout << " time: " << time << endl;
// 		  cout << " nodal data:\n" << nodal_data << endl;
		}
	
	    }

//   vtkPoints *points = vtkPoints::New();
//  for (int i=0; i<num_nodes; i++) points->InsertPoint(i,coordinates[i]);


		  
 points = vtkPoints::New();
 for (int i=0; i<num_nodes; i++) points->InsertPoint(i+1,coordinates(i));

//  warpGrid = vtkUnstructuredGrid::New();
//  dPoints = vtkPoints::New();

//  for (int i = 0; i<num_nodes; i++) dPoints->InsertPoint(i+1, nodal_data(i,0), nodal_data(i,1), nodal_data(i,2));
//  warpGrid->SetPoints(dPoints);

 
 

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

 /* create displacement "warp" vector */
 warp = vtkWarpVector::New();
 ugrid = vtkUnstructuredGrid::New();
 warp->SetInput(ugrid);
 scale_factor = 5;
  warp->SetScaleFactor(scale_factor);

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
//  if (coordinates.MinorDim() == 2) 
//    cell_types =VTK_QUAD;
 // else 
//    cell_types = VTK_HEXAHEDRON; // all the cells are the same

//  if (geometry ==GeometryT::kHexahedron) 
//    cell_types = VTK_HEXAHEDRON;
// if (geometry == GeometryT::kNone)

 /* convert exodus geometry into appropriate vtk geometry */
 if (geometry == GeometryT::kPoint) cell_types = VTK_VERTEX;  
 else if (geometry == GeometryT::kLine ) cell_types =  VTK_LINE;
 else if (geometry == GeometryT::kQuadrilateral) cell_types= VTK_QUAD; 
 else if (geometry == GeometryT::kTriangle) cell_types = VTK_TRIANGLE;
 else if (geometry == GeometryT::kHexahedron) cell_types  = VTK_HEXAHEDRON;
 else if (geometry == GeometryT::kTetrahedron) cell_types = VTK_TETRA; 
 else if (geometry == GeometryT::kPentahedron) cell_types = VTK_WEDGE;
 else cout << "Bad geometry";
 
 /* insert cells in the grid */
 //  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
 ugrid->SetCells(cell_types.Pointer(), vtk_cell_array);
 /* set default variable to be displayed */ 
 currentVarNum = num_node_variables-1;
 
  ugrid->SetPoints(points);
  points->Delete();
  ugrid->GetPointData()->SetScalars(scalars[frameNum][currentVarNum]);
  if (node_labels[0] == "D_X" || node_labels[0] == "D_Y" || node_labels[0] == "D_Z")
    ugrid->GetPointData()->SetVectors(vectors[frameNum][currentVarNum]);
  // ugrid->GetPointData()->SetScalars(scalars[0]);
  //ugrid->SetFieldData(fd[0]);

  // source_file = "../../example_files/heat/data_file1.vtk";
  output_file = "A.tif";
  /* color mapping variables */
  numColors = 256;
  hueRange1 = 0; hueRange2 = 0.6667;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
  // numRen = 1;
  // scalarRange1 = 6; scalarRange2 = 17;
  //for (int i=0; i<1; i++)
  // renderer[i] = vtkRenderer::New();
  renderer = vtkRenderer::New();
  renderer2 = vtkRenderer::New();
  renderer3 = vtkRenderer::New();
  renderer4 = vtkRenderer::New();
  
  
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
  cam = vtkCamera::New();
  axes = vtkCubeAxesActor2D::New();
  
  
//   renWin->AddRenderer(renderer);
//    renWin->AddRenderer(renderer2);
  // for (int i=0; i<4; i++)
  renWin->AddRenderer(renderer);
  if (numRen==4){   
    renWin->AddRenderer(renderer2);
    renWin->AddRenderer(renderer3);
    renWin->AddRenderer(renderer4);
  }
  iren->SetRenderWindow(renWin);
  
  /* color mapping stuff */
  lut->SetHueRange(hueRange1, hueRange2);
  lut->SetSaturationRange(satRange1,satRange2);
  lut->SetValueRange(valRange1,valRange2);
  lut->SetAlphaRange(alphaRange1,alphaRange2);
  lut->SetNumberOfColors(numColors);
  lut->Build();
    
  /* set warping vector if needed */
  if (node_labels[0] == "D_X" || node_labels[0] == "D_Y" || node_labels[0] == "D_Z")
    ugridMapper->SetInput(warp->GetOutput());
  else
    ugridMapper->SetInput(ugrid);
  /* if min and max scalar values are equal, decrease min slightly to provide for useful color mapping */
  for (int j =0; j< num_node_variables; j++)
    {
      if (scalarRange1[j] == scalarRange2[j])
	scalarRange1[j]-=.001;
    }
  ugridMapper->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
  ugridMapper->SetLookupTable(lut);
  //ugridMapper->ImmediateModeRenderingOn();

  ugridActor->SetMapper(ugridMapper);
  ugridActor->AddPosition(0,0.001,0);
  
//   wireActor->SetMapper(ugridMapper);
//   wireActor->GetProperty()->SetRepresentationToWireframe();
//   wireActor->GetProperty()->SetColor(0,0,0);

//   for (int i=0; i<4; i++){
    renderer->AddActor(ugridActor);
    renderer->SetBackground(0,0,0);
    renderer4->AddActor(ugridActor);
    renderer4->SetBackground(0,0,0);
    renderer2->AddActor(ugridActor);
    renderer2->SetBackground(0,0,0);
    renderer3->AddActor(ugridActor);
    renderer3->SetBackground(0,0,0);


//   }
    /* divide window into 4 parts */
    if (numRen ==4){
      renderer->SetViewport(0,0,.5,.5);
      renderer2->SetViewport(.5,0,1,.5);
      renderer3->SetViewport(0,.5,.5,1);
      renderer4->SetViewport(.5,.5,1,1);
      renderer->GetActiveCamera()->Zoom(0.85);
      renderer2->GetActiveCamera()->Zoom(0.85);
      renderer3->GetActiveCamera()->Zoom(0.85);
      renderer4->GetActiveCamera()->Zoom(0.85);
      }

//   renderer->AddActor(ugridActor);
//   // renderer->AddActor(wireActor);
//   renderer->SetBackground(0,0,0);
//    renderer->SetViewport(0,0,0.5,1);

//   renderer2->AddActor(ugridActor);
//   renderer2->SetBackground(1,1,1);
//   renderer2->SetViewport(0.5,0,1,1);

  renWin->SetPosition(668, 0);
  renWin->SetSize(600,700);
  
  scalarBar->SetLookupTable(ugridMapper->GetLookupTable());
  // sbTitle = "Temperature for frame 000";
  sbTitle.Append(node_labels[currentVarNum]); 
  sbTitle.Append(" for frame 000");
  // sbTitle.Append(time,5);
  scalarBar->SetTitle(sbTitle);
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
  scalarBar->SetOrientationToHorizontal();
  scalarBar->SetWidth(0.8);
  scalarBar->SetHeight(0.17);
   renderer->AddActor(scalarBar);

   /* assemble a string that is a list of all the variables */
  varList = "";
  for (int i = 0; i<num_node_variables; i++){
    varList.Append(i);
    varList.Append(":");
    varList.Append(" ");
    varList.Append(node_labels[i]);
    varList.Append("\n");
  }

  //TEMP - adding sub-scopes to the console
  fFrames.Allocate(numRen);
  for (int i = 0; i < numRen; i++)
    {
      fFrames[i].Append("frame",i,2);
      fFrames[i].iSetName(fFrames[i]);
      iAddSub(fFrames[i]);
    }
}

/* execute given command - returns false on fail */
bool VTKConsoleT::iDoCommand(const StringT& command, StringT& line)
{

  int sfbTest;
  int  varNum;
  double xRot, yRot, zRot;
  double timeStep;

  if (command == "Start_Rendering")
  {
    renWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;
    iren->Start();
    return true;
  }

  else if (command == "Update_Rendering")
  {

    ugridMapper->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
    lut->SetHueRange(hueRange1, hueRange2);
    lut->SetSaturationRange(satRange1,satRange2);
    lut->SetValueRange(valRange1,valRange2);
    lut->SetAlphaRange(alphaRange1,alphaRange2);
    lut->SetNumberOfColors(numColors);
    warp->SetScaleFactor(scale_factor);
    renWin->Render();
    cout << "type 'e' in the graphics window to exit interactive mode" << endl;   
    iren->Start();
    return true;
  }
    
  else if (command == "Reset_to_Default_Values")
    {
      // source_file = "../../example_files/heat/data_file1.vtk";
      numColors = 256;
      hueRange1 = 0; hueRange2 = 0.6667;
      valRange1 = 1; valRange2 = 1;
      satRange1 = 1; satRange2 = 1;
      alphaRange1 = 1; alphaRange2 = 1;
      // scalarRange1 = 6; scalarRange2 = 17;
      //  ugridMapper->SetScalarRange(scalarRange1,scalarRange2);
      lut->SetHueRange(hueRange1, hueRange2);
      lut->SetSaturationRange(satRange1,satRange2);
      lut->SetValueRange(valRange1,valRange2);
      lut->SetAlphaRange(alphaRange1,alphaRange2);
      lut->SetNumberOfColors(numColors);
      
      renWin->Render();
      // iren->Start();
      return true;
    }

  else if (command == "Reset_view")
    {
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(0,0,1);
	  cam->ComputeViewPlaneNormal();
	  cam->SetViewUp(0,1,0);
	  cam->OrthogonalizeViewUp();
	  renderer->SetActiveCamera(cam);
	  renderer->ResetCamera();
	  renWin->Render();
	  cout << "type 'e' in the graphics window to exit interactive mode" << endl;
	  iren->Start();
	  return true;
    }

  else if (command == "Save")
    {
      cout << "Enter name for file to be saved to: ";
      cin >> output_file;
      char line[255];
      cin.getline(line, 254);
      cout << "Save image at: \n 1: current view\n 2: default view: ";
      cin >> sfbTest;
      cin.getline(line, 254);
      /* if default camera desired */
      if (sfbTest == 2) {
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(0,0,1);
	  cam->ComputeViewPlaneNormal();
	  cam->SetViewUp(0,1,0);
	  cam->OrthogonalizeViewUp();
	  renderer->SetActiveCamera(cam);
	  renderer->ResetCamera();
	  renWin->Render();
	}
      int saveOpt;
      renSrc->SetInput(renderer);
      if (numRen == 4){ 
	cout << "Save options:\n 1: single plot\n 2: all 4 plots: ";
	cin >> saveOpt;
	cin.getline(line, 254);
	if (saveOpt == 2)
	  renSrc->WholeWindowOn();
	
      }
      writer->SetInput(renSrc->GetOutput());
      writer->SetFileName(output_file);
      writer->Write();
      renWin->Render();
      cout << "File " << output_file << " has been saved." << endl;
      //  iren->Start();
      return true;
    }

  else if (command == "Show_Node_Numbers")
    {
       ids->SetInput(ugrid);
      ids->PointIdsOn();
      //ids->CellIdsOn();
      //ids->FieldDataOn();
      //visPts->SetInput(ids->GetOutput());
      visPts->SetInput(ugrid);
      visPts->SetRenderer(renderer);

      ldm->SetInput(visPts->GetOutput());
      //ldm->SetInput(ids->GetOutput());
      ldm->SetLabelModeToLabelIds();
      ldm->ShadowOff();
      // ldm->SetLabelModeToLabelFieldData();
      pointLabels->SetMapper(ldm);
      pointLabels->VisibilityOn();
      renderer->AddActor2D(pointLabels);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;      

    }


  else if (command == "Color_bar_off")
    {
      renderer->RemoveActor(scalarBar);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command == "Color_bar_on")
    {
      renderer->AddActor(scalarBar);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }


  else if (command == "X_axis_rotation")
    {
 
      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> xRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Elevation(xRot);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }


  else if (command == "Y_axis_rotation")
    {

      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> yRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Azimuth(yRot);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command == "Z_axis_rotation")
    {

      cout << "Using the right-hand rule, rotate by how many degrees?: ";
      cin >> zRot;
      char line[255];
      cin.getline(line, 254);
      renderer->GetActiveCamera()->Roll(zRot);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command == "Flip_book")
    {
 
      cout << "Enter time step in seconds: ";
      cin >> timeStep;
      char line[255];
      cin.getline(line, 254);
      cout << "Show images at: \n 1: current view\n 2: default view: ";
      cin >> sfbTest;
      cin.getline(line, 254);
	
      /* if default camera angle desired */
	if (sfbTest == 2) {
	  cam->SetFocalPoint(0,0,0);
	  cam->SetPosition(0,0,1);
	  cam->ComputeViewPlaneNormal();
	  cam->SetViewUp(0,1,0);
	  cam->OrthogonalizeViewUp();
	  renderer->SetActiveCamera(cam);
	  renderer->ResetCamera();
	  renWin->Render();
	}

      for (int j = 0; j<num_time_steps; j++){
	/* time delay */
        clock_t start_time, cur_time;
         start_time = clock();
         while((clock() - start_time) < timeStep * CLOCKS_PER_SEC)
         {
         }
	 // sbTitle = "";
	 // sbTitle.Append(node_labels(0)); 
	 // sbTitle.Append(" for frame ");
	 // sbTitle = "Temperature for frame ";
	 sbTitle.Drop(-3);
	  sbTitle.Append(j,3);
	 scalarBar->SetTitle(sbTitle);
	 ugrid->GetPointData()->SetScalars(scalars[j][currentVarNum]);
	 if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	   ugrid->GetPointData()->SetVectors(vectors[j][currentVarNum]);
	 // ugrid->GetPointData()->SetScalars(scalars[j]);
	 renWin->Render();

      }

      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command== "Save_flip_book_images")
    {
      StringT fbName;
      cout << "Enter name for flipbook to be saved (without .tif extension): ";
      cin >> fbName;
      char line[255];
      cin.getline(line,254);
      cout << "Save images at: \n 1: current view\n 2: default view: ";
      cin >> sfbTest;
      cin.getline(line, 254);
      /* if default camera desired */
      if (sfbTest == 2) {
	cam->SetFocalPoint(0,0,0);
	cam->SetPosition(0,0,1);
	cam->ComputeViewPlaneNormal();
	cam->SetViewUp(0,1,0);
	cam->OrthogonalizeViewUp();
	renderer->SetActiveCamera(cam);
	renderer->ResetCamera();
	renWin->Render();
      }	
      
      for (int j = 0; j<num_time_steps; j++){

	sbTitle.Drop(-3);
	sbTitle.Append(j,3);
	scalarBar->SetTitle(sbTitle);	 
	ugrid->GetPointData()->SetScalars(scalars[j][currentVarNum]);
	
	if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	  ugrid->GetPointData()->SetVectors(vectors[j][currentVarNum]);
	
	renWin->Render();  
	renSrc->SetInput(renderer);
	renSrc->WholeWindowOn();
	writer->SetInput(renSrc->GetOutput());
	outFileName = fbName;
	outFileName.Append(j,3); // pad to a width of 3 digits
	outFileName.Append(".tif");
	writer->SetFileName(outFileName);
	writer->Write();
	cout << outFileName << " has been saved" << endl;
	renWin->Render();
      }
      cout << "Flip book images have been saved." << endl;
      renWin->Render();
      //    iren->Start();
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
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command == "Select_frame_number")
    {

      cout << "choose frame number from 0 to " << num_time_steps-1 <<" to be displayed: ";
      cin >> frameNum;
      char line[255];
      cin.getline(line, 254);
      // sbTitle = "Temperature for frame ";
      // sbTitle = "";
      // sbTitle.Append(node_labels(0)); 
      // sbTitle.Append(" for frame ");
      sbTitle.Drop(-3);
      sbTitle.Append(frameNum,3);
      // sbTitle.Append(j,3);
      scalarBar->SetTitle(sbTitle);
      // ugrid->GetPointData()->SetScalars(scalars[frameNum]);
      ugrid->GetPointData()->SetScalars(scalars[frameNum][currentVarNum]);
      if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	ugrid->GetPointData()->SetVectors(vectors[frameNum][currentVarNum]);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
      
    }

   else if (command == "Hide_Node_Numbers")

    {

// //       ids->PointIdsOff();
// //       ids->Update();
      pointLabels->VisibilityOff();
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;   
 
    }
  
  else if (command == "Show_axes")
  {

  // x,y,z axes
      
      axes->SetInput(ugrid);
     axes->SetCamera(renderer->GetActiveCamera());
    //  axes->SetLabelFormat("%6.4g");
     // axes->ShadowOn();
     // axes->SetFlyModeToOuterEdges();
      axes->SetFlyModeToClosestTriad();
      axes->SetFontFactor(3.8);
      axes->GetProperty()->SetColor(0,1,1);
      // axes->SetBounds(0,1,0,1,0,1);
      //axes->ZAxisVisibilityOff();
      axes->VisibilityOn();
      renderer->AddActor2D(axes);
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
  }

  else if (command == "Hide_axes")
    {
      axes->VisibilityOff();
      renWin->Render();
      cout << "type 'e' in the graphics window to exit interactive mode" << endl;
      iren->Start();
      return true;
    }

  else if (command == "Choose_variable")
    {
 //      char line[255];
//       if (numRen == 1){
// 	cout << "choose variable number from 0 to " << num_node_variables-1 <<" to be displayed\n" << varList;      
// 	cin >> currentVarNum;
// 	cin.getline(line, 254);
	//  }
 //      else if (numRen == 4){
// 	for (int i = 0; i<4; i++){
// 	  cout << "choose variable number from 0 to " << num_node_variables-1 <<" to be displayed for renderer " << i << endl << varList; 
// 	cin >> currentVarNum;
// 	cin.getline(line, 254);
	
// 	}
//       }

      char line[255];
      cout << "choose variable number from 0 to " << num_node_variables-1 <<" to be displayed\n" << varList;      
      cin >> currentVarNum;
      cin.getline(line, 254);
      ugrid->GetPointData()->SetScalars(scalars[frameNum][currentVarNum]);
      ugridMapper->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
      if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	ugrid->GetPointData()->SetVectors(vectors[frameNum][currentVarNum]);
      sbTitle = "";
      sbTitle.Append(node_labels[currentVarNum]); 
      sbTitle.Append(" for frame ");
      sbTitle.Append(frameNum,3);
      scalarBar->SetTitle(sbTitle);
      renWin->Render();
      return true;

    }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
