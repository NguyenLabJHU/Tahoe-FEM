/* $Id: VTKUGridT.cpp,v 1.1 2001-12-10 12:44:08 paklein Exp $ */
#include "VTKUGridT.h"

#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkWarpVector.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"

#include "iArray2DT.h"

/* array behavior */
const bool ArrayT<VTKUGridT*>::fByteCopy = true;

/* constructor */
VTKUGridT::VTKUGridT(TypeT my_type, int id, int nsd):
	fType(my_type),
	fID(id),
	fNumSD(nsd),
	fCellArray(NULL),
	fConnects(NULL),
	fUGrid(NULL),
	fMapper(NULL),
	fLookUpTable(NULL),
	fActor(NULL),
	fWarp(NULL)
{
	/* initialize grid */
	fUGrid = vtkUnstructuredGrid::New();
	
	/* the mapper */
	fMapper = vtkDataSetMapper::New();
	fMapper->SetInput(fUGrid); /* just map grid by default */

	/* the actor */
	fActor = vtkActor::New();

//TEMP
//	fActor->GetProperty()->SetOpacity(0.1);
	
	/* line color */
	if (fType == kElementSet)
		fActor->GetProperty()->SetColor(1,0,0);
	else if (fType == kNodeSet)
	{
		fActor->GetProperty()->SetColor(0,0,1);
		fActor->GetProperty()->SetPointSize(3.0);
	}
	else
		fActor->GetProperty()->SetColor(0,1,0);
	fActor->SetMapper(fMapper);
	fActor->AddPosition(0,0.001,0);
}

/* destructor */
VTKUGridT::~VTKUGridT(void)
{
	/* clean up */
	if (fCellArray) fCellArray->Delete();
	if (fConnects) fConnects->Delete();
	if (fUGrid) fUGrid->Delete();
	if (fMapper) fMapper->Delete();
	if (fActor) fActor->Delete();
	if (fWarp) fWarp->Delete();
	if (fLookUpTable) fLookUpTable->Delete();
}

/* set the point data */
void VTKUGridT::SetPoints(vtkPoints* points)
{
	/* insert points */
	fUGrid->SetPoints(points);
}

/* set the connectivities */
void VTKUGridT::SetConnectivities(GeometryT::CodeT code, const iArray2DT& connects)
{
	/* create array of VTK-style connectivities */
	iArray2DT vtk_connects(connects.MajorDim(), connects.MinorDim()+1); //has 1 extra entry!!!
	vtk_connects.BlockColumnCopyAt(connects, 1);
	vtk_connects.SetColumn(0, connects.MinorDim()); //first value in each row is row size 

	/* release memory */
	int* p_vtk_connects;
	vtk_connects.ReleasePointer(&p_vtk_connects);

	/* construct id array */
	if (!fConnects) fConnects = vtkIdTypeArray::New();
	fConnects->SetNumberOfComponents(vtk_connects.MinorDim());
	fConnects->SetArray(p_vtk_connects, vtk_connects.Length(), 0);

	/* construct cell array */
	if (!fCellArray) fCellArray = vtkCellArray::New();
	fCellArray->SetCells(vtk_connects.MajorDim(), fConnects);

	/* convert Tahoe geometry into appropriate vtk geometry */
	fCellType.Allocate(vtk_connects.MajorDim());
	if (code == GeometryT::kPoint) fCellType = VTK_VERTEX;  
	else if (code == GeometryT::kLine ) fCellType =  VTK_LINE;
	else if (code == GeometryT::kQuadrilateral) fCellType = VTK_QUAD; 
	else if (code == GeometryT::kTriangle) fCellType = VTK_TRIANGLE;
	else if (code == GeometryT::kHexahedron) fCellType = VTK_HEXAHEDRON;
	else if (code == GeometryT::kTetrahedron) fCellType = VTK_TETRA; 
	else if (code == GeometryT::kPentahedron) fCellType = VTK_WEDGE;
	else {
		cout << "\n VTKUGridT::SetConnectivities: unknown geometry code: " << code << endl;
		throw eGeneralFail;
	}

	/* insert cells in the grid */
	fUGrid->SetCells(fCellType.Pointer(), fCellArray);
}
  
/* set the scalar data */
void VTKUGridT::SetScalars(vtkFloatArray* scalars)
{
	/* insert in grid */
	fUGrid->GetPointData()->SetScalars(scalars); 
}

/* set the scalar data range */
void VTKUGridT::SetScalarRange(double min, double max)
{
	fMapper->SetScalarRange(min, max);
}

/* set the vector data */
void VTKUGridT::SetVectors(vtkFloatArray* vectors)
{
	/* insert in grid */
	fUGrid->GetPointData()->SetVectors(vectors);
}

/* set vectors that warp */
void VTKUGridT::SetWarpVectors(vtkFloatArray* vectors)
{
	/* insert in grid */
	SetVectors(vectors);
	
	/* set up warp vector */
	if (!fWarp) fWarp = vtkWarpVector::New();
	fWarp->SetInput(fUGrid);
	fMapper->SetInput(fWarp->GetOutput());
}
