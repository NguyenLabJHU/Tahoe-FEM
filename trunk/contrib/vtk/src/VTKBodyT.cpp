/* $Id: VTKBodyT.cpp,v 1.5 2001-11-01 19:16:43 recampb Exp $ */

#include "VTKBodyT.h"

#include "vtkPoints.h"

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkRendererSource.h"
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
#include "vtkScalarBarActor.h"
#include "StringT.h"

#include <iostream.h>
#include <iomanip.h>
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "GeometryT.h"

/* array behavior */
const bool ArrayT<VTKBodyT*>::fByteCopy = true;

/* constructor */
VTKBodyT::VTKBodyT(const StringT& file_name): 
  inFile(file_name)
{

 /* read exodus file */
  ExodusT exo(cout);
  if (!exo.OpenRead(inFile))
    {
      cout << " ERROR: could not open file: " << inFile << endl;
      throw;
    }
  else
    cout << "read database file: " << inFile << endl;
  
  /* read coordinates */
  //dArray2DT coordinates;
  num_nodes = exo.NumNodes();
  num_dim   = exo.NumDimensions();
  dArray2DT coordinates(num_nodes, num_dim);

  
//   // ArrayT<dArray2DT> coordinates(num_time_steps);
  exo.ReadCoordinates(coordinates); //******???
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
  //num_time_steps =1; //num_node_variables =1; // used to test data
  if (num_time_steps > 0)
    {
      for (int i = 0; i < num_time_steps; i++)
      //for (int i = 0; i<1; i++)
      	{
	  exo.ReadTime(i+1, time);
	   for (int j = 0; j < num_node_variables; j++)
	  //for(int j=0; j<1; j++)  
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
	}
      
    }
  

  /* allocate points */
  points = vtkPoints::New();
  for (int i=0; i<num_nodes; i++) points->InsertPoint(i+1,coordinates(i));
  
  /* allocate cells */
  vtk_cell_array = vtkCellArray::New();

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
  
 
  
  /* create VTK array of cells */
  //vtkCellArray* vtk_cell_array = vtkCellArray::New();
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
  

  /* set default variable to be displayed */ 
  currentVarNum = num_node_variables-1;
  
  
  //read model data
  //set up grid
  ugrid = vtkUnstructuredGrid::New();
  currentStepNum = 0;
  /* insert cells in the grid */
  ugrid->SetCells(cell_types.Pointer(), vtk_cell_array);
  ugrid->SetPoints(points);
  ugrid->GetPointData()->SetScalars(scalars[currentStepNum][currentVarNum]); 
  if (node_labels[0] == "D_X" || node_labels[0] == "D_Y" || node_labels[0] == "D_Z")
    ugrid->GetPointData()->SetVectors(vectors[currentStepNum][currentVarNum]);
  
 
  scalarBar = vtkScalarBarActor::New();
  warp = vtkWarpVector::New();
  warp->SetInput(ugrid);
  scale_factor = 1;

  //set up mapper
  ugridMapper = vtkDataSetMapper::New();
  ugridMapper->SetInput(warp->GetOutput());
 /* set warping vector if needed */  
  if (node_labels[0] == "D_X" ||node_labels[0] == "D_Y" || node_labels[0] == "D_Z")
    ugridMapper->SetInput(warp->GetOutput());
  else
    ugridMapper->SetInput(ugrid);
  
  //set up actor
  ugridActor = vtkActor::New();
  ugridActor->SetMapper(ugridMapper);
  ugridActor->AddPosition(0,0.001,0);


  /* color mapping variables */
  DefaultValues();
  
  /* color mapping stuff */  
  lut = vtkLookupTable::New();
  UpdateData();
  lut->Build();
  SetLookupTable();

  /* assemble a string that is a list of all the variables */
  varList = "";
  for (int i = 0; i<num_node_variables; i++){
    varList.Append(i);
    varList.Append(":");
    varList.Append(" ");
    varList.Append(node_labels[i]);
    varList.Append("\n");
  }

  
}

/* destructor */
VTKBodyT::~VTKBodyT(void)
{
  /* clean up */
  points->Delete();
  vtk_cell_array->Delete();
  ugrid->Delete();
  warp->Delete();
  ugridMapper->Delete();
  ugridActor->Delete();
  lut->Delete();
  scalars[0][0]->Delete();
  vectors[0][0]->Delete();
  
}

void VTKBodyT::SetLookupTable(void)
{
  scalarBar->SetLookupTable(ugridMapper->GetLookupTable());
  sbTitle.Append(node_labels[currentVarNum]); 
  sbTitle.Append(" for frame 000");
  scalarBar->SetTitle(sbTitle);
  scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  scalarBar->GetPositionCoordinate()->SetValue(0.1,0.01);
  scalarBar->SetOrientationToHorizontal();
  scalarBar->SetWidth(0.8);
  scalarBar->SetHeight(0.17);
  //renderer->AddActor(scalarBar);
}

void VTKBodyT::UpdateData(void)
{
  lut->SetHueRange(hueRange1, hueRange2);
  lut->SetSaturationRange(satRange1, satRange2);
  lut->SetValueRange(valRange1, valRange2);
  lut->SetAlphaRange(alphaRange1, alphaRange2);
  lut->SetNumberOfColors(numColors);
  warp->SetScaleFactor(scale_factor);
  ugridMapper->SetScalarRange(scalarRange1[currentVarNum], scalarRange2[currentVarNum]);
  ugridMapper->SetLookupTable(lut);
}

void VTKBodyT::DefaultValues(void)
{
  numColors = 256;
  hueRange1 = 0; hueRange2 = 0.6667;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
}

void VTKBodyT::ChangeVars(const int varNum)
{
      ugrid->GetPointData()->SetScalars(scalars[currentStepNum][varNum]);
      ugridMapper->SetScalarRange(scalarRange1[varNum],scalarRange2[varNum]);
      if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
	ugrid->GetPointData()->SetVectors(vectors[currentStepNum][varNum]);
      sbTitle = "";
      sbTitle.Append(node_labels[varNum]); 
      sbTitle.Append(" for frame 000 ");
      //sbTitle.Append(currentStepNum,3);
      scalarBar->SetTitle(sbTitle);
      currentVarNum = varNum;

}

void VTKBodyT::SelectTimeStep(const int stepNum)
{
  sbTitle.Drop(-3);
  sbTitle.Append(stepNum,3);
  scalarBar->SetTitle(sbTitle);
  // ugrid->GetPointData()->SetScalars(scalars[currentStepNum]);
  ugrid->GetPointData()->SetScalars(scalars[stepNum][currentVarNum]);
  if (node_labels[0] == "D_X" || node_labels[1] == "D_Y" || node_labels[2] == "D_Z")
    ugrid->GetPointData()->SetVectors(vectors[stepNum][currentVarNum]);
  currentStepNum= stepNum;
}
