/* $Id: VTKConsoleT.cpp,v 1.4 2001-09-28 00:18:06 recampb Exp $ */

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
  iAddVariable("output_file", output_file);

  /* add console commands */
  iAddCommand("Start_Rendering");
  iAddCommand("Update_Rendering");
  iAddCommand("Reset_to_Default_Values");
  iAddCommand("Save");
  iAddCommand("Show_Node_Numbers");
  //  iAddCommand("Hide_Node_Numbers");

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
  exo.ReadCoordinates(coordinates);

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
  int num_time_steps = exo.NumTimeSteps();
  //if (num_time_steps > 0)
  //	{
	  /* variables defined at the nodes */
	  int num_node_variables = exo.NumNodeVariables();
	  ArrayT<StringT> node_labels(num_node_variables);
	  exo.ReadNodeLabels(node_labels);
	  cout << " nodal variables:\n";
	  for (int i = 0; i < node_labels.Length(); i++)
		cout << node_labels[i] << '\n';

	  /* variables defined over the elements */
	  int num_element_variables = exo.NumElementVariables();
	  ArrayT<StringT> element_labels;
	  exo.ReadElementLabels(element_labels);
	  cout << " element variables:\n" << endl;
	  for (int i = 0; i < element_labels.Length(); i++)
		cout << element_labels[i] << '\n';
	  cout.flush();

	  /* read nodal data */
	  dArray2DT nodal_data(num_nodes, num_node_variables);
	  dArrayT ndata(num_nodes);

	  if (num_time_steps > 0)
	    {
	      for (int i = 0; i < num_time_steps; i++)
		{
		  double time;
		  exo.ReadTime(i+1, time);
		  
		  /* loop over variables */
		  for (int j = 0; j < num_node_variables; j++)
		    {
		      exo.ReadNodalVariable(i+1, j+1, ndata);
		      nodal_data.SetColumn(j, ndata);
		    }
		  
		  cout << " time: " << time << endl;
		  cout << " nodal data:\n" << nodal_data << endl;
		}
	
	    }


  vtkPoints *points = vtkPoints::New();
 for (int i=0; i<num_nodes; i++) points->InsertPoint(i,coordinates[i]);
  
  vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
    ugrid->Allocate(num_elem_blocks);
  for (int i=0; i<num_elem_blocks; i++) 
    ugrid->InsertNextCell(VTK_TETRA, 4, connectivities[i]);

  vtkScalars *scalars = vtkScalars::New(VTK_DOUBLE);
  for (int i=0; i<num_nodes; i++) scalars->InsertScalar(i,nodal_data[i]);

  ugrid->SetPoints(points);
  points->Delete();
  ugrid->GetPointData()->SetScalars(scalars);


  source_file = "../../example_files/heat/data_file1.vtk";
  output_file = "A.tif";
  numColors = 256;
  hueRange1 = 0; hueRange2 = 0.6667;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
  scalarRange1 = 6; scalarRange2 = 17;
  ugr = vtkUnstructuredGridReader::New();
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


  ugr->SetFileName(source_file);

  renWin->AddRenderer(renderer);
 
  iren->SetRenderWindow(renWin);

  lut->SetHueRange(hueRange1, hueRange2);
  lut->SetSaturationRange(satRange1,satRange2);
  lut->SetValueRange(valRange1,valRange2);
  lut->SetAlphaRange(alphaRange1,alphaRange2);
  lut->SetNumberOfColors(numColors);
  lut->Build();
  
  ugridMapper->SetInput(ugr->GetOutput());
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
  scalarBar->SetTitle("Temperature for time .50000");
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
//   if (command == "Integer_Print")
//     {
//       cout << "int = " << fInteger << endl;
//       return true;
//     }


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
    ugr->SetFileName(source_file);
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
      return true;
    }

  else if (command == "Save")
    {
      //  cout << "Enter name for file to be saved to: ";
      //  cin >> output_file;
    
      renSrc->SetInput(renderer);
      renSrc->WholeWindowOn();
      writer->SetInput(renSrc->GetOutput());
      writer->SetFileName(output_file);
      writer->Write();
      return true;
    }

  else if (command == "Show_Node_Numbers")
    {
    
      ids->SetInput(ugr->GetOutput());
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

//   else if (command == "Hide_Node_Numbers")
//     {
//      ids->SetInput(ugr->GetOutput());
      
//       ids->CellIdsOff();
//       ids->FieldDataOff();
//       ids->PointIdsOff();
//       ids->Update();

//       visPts->SetInput(ids->GetOutput());
//       visPts->SetRenderer(renderer);
//       //visPts->SelectionWindowOn();
//       //  visPts->SetSelection();
//       ldm->SetInput(visPts->GetOutput());
//       // ldm->SetlabelFormat();
//       ldm->SetLabelModeToLabelFieldData();
//       pointLabels->SetMapper(ldm);
//       renderer->AddActor2D(pointLabels);
//       renWin->Render();
//       iren->Start();
//       return true;   
 
//     }

  else
    /* drop through to inherited */
    return iConsoleObjectT::iDoCommand(command, line);
}
