/* $Id: VTKConsoleT.cpp,v 1.18 2001-10-24 18:20:36 paklein Exp $ */

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
  for (int p=0; p<4; p++){
    cout << "Choose number of File\n 1: heat.io0.exo\n 2: test.io0.exo\n 3: test2.io0.exo\n 4: big.exo\n 5: enter different .exo: ";
    cin >> test;
    cin.getline(line, 254);
    
    if (test == 1) fFrames[p].inFile = "../../example_files/heat/heat.io0.exo";
    else if (test ==2) fFrames[p].inFile ="test.io0.exo";
    else if (test == 3) fFrames[p].inFile = "test2.io0.exo";
    else if (test == 4) fFrames[p].inFile = "big.exo";
    else if (test == 5) {
      cout << "Enter file name with .exo: ";
      cin >> fFrames[p].inFile;
      cin.getline(line, 254);
    }
    else cout << "bad entry";
  
  /* read exodus file */
    ExodusT exo(cout);
    if (!exo.OpenRead(fFrames[p].inFile))
      {
	  cout << " ERROR: could not open file: " << fFrames[p].inFile << endl;
	  throw;
      }
    else
      cout << "read database file: " << fFrames[p].inFile << endl;
    

  //TEMP - adding sub-scopes to the console
  fFrames.Allocate(4);
  for (int i = 0; i < 4; i++)
    {
      
      fFrames[i].Append("frame",i,2);
      fFrames[i].iSetName(fFrames[i]);
      iAddSub(fFrames[i]);
    }

  /* read coordinates */
  dArray2DT coordinates;
  fFrames[p].num_nodes = exo.NumNodes();
  fFrames[p].num_dim   = exo.NumDimensions();
  coordinates(fFrames[p].num_nodes, fFrames[p].num_dim);

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
  fFrames[p].num_time_steps = exo.NumTimeSteps();
  // double time;
  //if (num_time_steps > 0)
  //	{
	  /* variables defined at the nodes */
	  fFrames[p].num_node_variables = exo.NumNodeVariables();
	  // ArrayT<StringT> node_labels(num_node_variables);
	  exo.ReadNodeLabels(fFrames[p].node_labels);
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
	  dArray2DT nodal_data(fFrames[p].num_nodes, fFrames[p].num_node_variables);
	  dArrayT ndata(fFrames[p].num_nodes);
	  // scalarRange1 = 10000;
	  // scalarRange2 = -10000;
	  if (fFrames[p].num_time_steps > 0)
	    {
	      for (int i = 0; i < fFrames[p].num_time_steps; i++)
		{
		  exo.ReadTime(i+1, time);
		  for (int j = 0; j < fFrames[p].num_node_variables; j++)
		    {
		     exo.ReadNodalVariable(i+1, j+1, ndata);
		      nodal_data.SetColumn(j, ndata);
		      fFrames[p].scalars[i][j] =  vtkScalars::New(VTK_DOUBLE);
		    
		      /* instantiate displacement vector if needed */
		      if (fFrames[p].node_labels[0] == "D_X" || fFrames[p].node_labels[0] == "D_Y" || fFrames[p].node_labels[0] == "D_Z")
			fFrames[p].vectors[i][j] = vtkVectors::New(VTK_DOUBLE);
		      
		      /* initialize min and max scalar range values */
		      fFrames[p].scalarRange1[j] = 10000;
		      fFrames[p].scalarRange2[j] = -10000;
		      
		      for (int k = 0; k<fFrames[p].num_nodes; k++) {
			
			//	scalars[i][j]->InsertScalar(k+1, nodal_data[k][j]);
			/* determine min and max scalar range values */
			if (nodal_data(k,j) < fFrames[p].scalarRange1[j]) fFrames[p].scalarRange1[j] = nodal_data(k,j);
			if (nodal_data(k,j) > fFrames[p].scalarRange2[j]) fFrames[p].scalarRange2[j] = nodal_data(k,j);
			/* insert scalar value at each node for each variable and time step */
			fFrames[p].scalars[i][j]->InsertScalar(k+1, nodal_data(k,j));
			//InsertVector(k,...)?????
			/* if displacement vector needed then insert vector at each node for each time step */
			if (fFrames[p].node_labels[0] == "D_X" && fFrames[p].node_labels[1] == "D_Y" && fFrames[p].node_labels[2] == "D_Z")              
			  fFrames[p].vectors[i][j]->InsertVector(k+1, nodal_data(k,0),nodal_data(k,1),nodal_data(k,2));
			else if (fFrames[p].node_labels[0] == "D_X" && fFrames[p].node_labels[1] == "D_Y") 
			   fFrames[p].vectors[i][j]->InsertVector(k+1, nodal_data(k,0), nodal_data(k,1),0);
			else if (fFrames[p].node_labels[0] == "D_X")
			  fFrames[p].vectors[i][j]->InsertVector(k+1, nodal_data(k,0),0,0);
			else;
			
		      }
		    }
// 		  cout << " time: " << time << endl;
// 		  cout << " nodal data:\n" << nodal_data << endl;
		}
	
	    }
 

	  fFrames[p].points = vtkPoints::New();
	  for (int i=0; i<fFrames[p].num_nodes; i++) fFrames[p].points->InsertPoint(i+1,coordinates(i));


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
	  fFrames[p].warp = vtkWarpVector::New();
	  fFrames[p].ugrid = vtkUnstructuredGrid::New();
	  fFrames[p].warp->SetInput(fFrames[p].ugrid);
	  fFrames[p].scale_factor = 5;
	  fFrames[p].warp->SetScaleFactor(fFrames[p].scale_factor);

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
	  fFrames[p].ugrid->SetCells(cell_types.Pointer(), vtk_cell_array);
	  /* set default variable to be displayed */ 
	  fFrames[p].currentVarNum = fFrames[p].num_node_variables-1;
	  
	  fFrames[p].ugrid->SetPoints(fFrames[p].points);
	  fFrames[p].points->Delete();
	  fFrames[p].ugrid->GetPointData()->SetScalars(fFrames[p].scalars[fFrames[p].frameNum][fFrames[p].currentVarNum]);
	  if (fFrames[p].node_labels[0] == "D_X" || fFrames[p].node_labels[0] == "D_Y" || fFrames[p].node_labels[0] == "D_Z")
	    fFrames[p].ugrid->GetPointData()->SetVectors(fFrames[p].vectors[fFrames[p].frameNum][fFrames[p].currentVarNum]);

	  // source_file = "../../example_files/heat/data_file1.vtk";
	  fFrames[p].output_file = "A.tif";
	  /* color mapping variables */
	  fFrames[p].numColors = 256;
	  fFrames[p].hueRange1 = 0; fFrames[p].hueRange2 = 0.6667;
	  fFrames[p].valRange1 = 1; fFrames[p].valRange2 = 1;
	  fFrames[p].satRange1 = 1; fFrames[p].satRange2 = 1;
	  fFrames[p].alphaRange1 = 1; fFrames[p].alphaRange2 = 1;
	  
	  fFrames[p].renderer = vtkRenderer::New();
	  
	  // renWin = vtkRenderWindow::New();
	  //iren = vtkRenderWindowInteractor::New();
	  fFrames[p].lut = vtkLookupTable::New();
	  fFrames[p].ugridMapper = vtkDataSetMapper::New();
	  fFrames[p].ugridActor = vtkActor::New();
	  fFrames[p].wireActor = vtkActor::New();
	  fFrames[p].scalarBar = vtkScalarBarActor::New();
	  //writer = vtkTIFFWriter::New();
	  fFrames[p].renSrc = vtkRendererSource::New();
	  fFrames[p].ids = vtkIdFilter::New();
	  fFrames[p].visPts = vtkSelectVisiblePoints::New();
	  fFrames[p].ldm = vtkLabeledDataMapper::New();
	  fFrames[p].pointLabels = vtkActor2D::New();
	  fFrames[p].cam = vtkCamera::New();
	  fFrames[p].axes = vtkCubeAxesActor2D::New();

	  /* color mapping stuff */
	  fFrames[p].lut->SetHueRange(fFrames[p].hueRange1,fFrames[p].hueRange2);
	  fFrames[p].lut->SetSaturationRange(fFrames[p].satRange1,fFrames[p].satRange2);
	  fFrames[p].lut->SetValueRange(fFrames[p].valRange1,fFrames[p].valRange2);
	  fFrames[p].lut->SetAlphaRange(fFrames[p].alphaRange1,fFrames[p].alphaRange2);
	  fFrames[p].lut->SetNumberOfColors(fFrames[p].numColors);
	  fFrames[p].lut->Build();
	  
	  /* set warping vector if needed */
	  if (fFrames[p].node_labels[0] == "D_X" ||fFrames[p].node_labels[0] == "D_Y" || fFrames[p].node_labels[0] == "D_Z")
	    fFrames[p].ugridMapper->SetInput(fFrames[p].warp->GetOutput());
	  else
	    fFrames[p].ugridMapper->SetInput(fFrames[p].ugrid);
	  /* if min and max scalar values are equal, decrease min slightly to provide for useful color mapping */
	  for (int j =0; j< fFrames[p].num_node_variables; j++)
	    {
	      if (fFrames[p].scalarRange1[j] == fFrames[p].scalarRange2[j])
		fFrames[p].scalarRange1[j]-=.001;
	    }
	  fFrames[p].ugridMapper->SetScalarRange(fFrames[p].scalarRange1[fFrames[p].currentVarNum],fFrames[p].scalarRange2[fFrames[p].currentVarNum]);
	  fFrames[p].ugridMapper->SetLookupTable(fFrames[p].lut);
	  //ugridMapper->ImmediateModeRenderingOn();
	  
	  fFrames[p].ugridActor->SetMapper(fFrames[p].ugridMapper);
	  fFrames[p].ugridActor->AddPosition(0,0.001,0);
	  	  
	  fFrames[p].renderer->AddActor(fFrames[p].ugridActor);
	  fFrames[p].renderer->SetBackground(0,0,0);

	  fFrames[p].scalarBar->SetLookupTable(fFrames[p].ugridMapper->GetLookupTable());
	  // sbTitle = "Temperature for frame 000";
	  fFrames[p].sbTitle.Append(fFrames[p].node_labels[fFrames[p].currentVarNum]); 
	  fFrames[p].sbTitle.Append(" for frame 000");
	  // sbTitle.Append(time,5);
	  fFrames[p].scalarBar->SetTitle(fFrames[p].sbTitle);
	  fFrames[p].scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	  fFrames[p].scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
	  fFrames[p].scalarBar->SetOrientationToHorizontal();
	  fFrames[p].scalarBar->SetWidth(0.8);
	  fFrames[p].scalarBar->SetHeight(0.17);
	  fFrames[p].renderer->AddActor(fFrames[p].scalarBar);
	  
	  /* assemble a string that is a list of all the variables */
	  fFrames[p].varList = "";
	  for (int i = 0; i<fFrames[p].num_node_variables; i++){
	    fFrames[p].varList.Append(i);
	    fFrames[p].varList.Append(":");
	    fFrames[p].varList.Append(" ");
	    fFrames[p].varList.Append(fFrames[p].node_labels[i]);
	    fFrames[p].varList.Append("\n");
	  }
	  
  }

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


  //TEMP - adding sub-scopes to the console
//   fFrames.Allocate(numRen);
//   for (int i = 0; i < numRen; i++)
//     {
//       fFrames[i].Append("frame",i,2);
//       fFrames[i].iSetName(fFrames[i]);
//       iAddSub(fFrames[i]);
//     }
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
