/* $Id: VTKUGridT.cpp,v 1.6 2002-04-18 00:27:49 paklein Exp $ */
#include "VTKUGridT.h"

#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkWarpVector.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"

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
	
	/* change color range from blue to red */
	fLookUpTable = vtkLookupTable::New();
	fLookUpTable->SetHueRange(0.6667, 0);
	fMapper->SetLookupTable(fLookUpTable);

	/* the actor */
	fActor = vtkActor::New();
	
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
	iArray2DT vtk_connects;

	/* quads with mid-side nodes - don't display properly */
	if (code == GeometryT::kQuadrilateral && connects.MinorDim() > 4)
	{
#if 1
//NOTE: display quad8's as quad4's

		cout << "\n VTKUGridT::SetConnectivities: quad8's reduced to quad4's\n" << endl;

		/* allocate */
		vtk_connects.Dimension(connects.MajorDim(), 4 + 1); //has 1 extra entry!!!
		vtk_connects.SetColumn(0,4); //first value in each row is row size 

		/* copy 1st 4 nodes */
		for (int i = 0; i < connects.MajorDim(); i++)
		{
			int* a = connects(i);
			int* b = vtk_connects(i) + 1; /* first value in each row is row size */
			
			memcpy(b, a, 4*sizeof(int));
		}
#endif

#if 0
//NOTE: display quad8's

		/* allocate */
		vtk_connects.Dimension(connects.MajorDim(), connects.MinorDim()+1); //has 1 extra entry!!!
		vtk_connects.SetColumn(0, connects.MinorDim()); //first value in each row is row size 

		/* reorder around the element edge */
		int n_mid = connects.MinorDim() - 4;
		for (int i = 0; i < connects.MajorDim(); i++)
		{
			int* a = connects(i);
			int* b = vtk_connects(i) + 1; /* first value in each row is row size */
			int* a_mid = a + 4;
			for (int j = 0; j < 4; j++)
			{
				*b++ = *a++;
				
				/* interleave mid-side nodes */
				if (j < n_mid)
					*b++ = *a_mid++;
			}
			
			/* rotate to make first node a mid-side. For some reason, fields over
			 * 8-noded quads were not displayed correctly with the first node in
			 * a corner. This makes uniform Y gradients look correct, but uniform
			 * x gradients still weren't right */
			a = vtk_connects(i) + 1;
			int tmp = a[connects.MinorDim() - 1];
			memmove(a+1, a, (connects.MinorDim() - 1)*sizeof(int));
			a[0] = tmp;
		}
#endif
	}
	/* else just copy in */
	else
	{
		vtk_connects.Dimension(connects.MajorDim(), connects.MinorDim()+1); //has 1 extra entry!!!
		vtk_connects.SetColumn(0, connects.MinorDim()); //first value in each row is row size 
		vtk_connects.BlockColumnCopyAt(connects, 1);
	}

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
	else if (code == GeometryT::kQuadrilateral) {
		if (vtk_connects.MinorDim() == 5)
			fCellType = VTK_QUAD;
		else
			fCellType = VTK_POLYGON;
	}
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

/* set the number of color levels */
void VTKUGridT::SetNumberOfColors(int num)
{
	fLookUpTable->SetNumberOfColors(num);
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

/* set the wrap displacement scale factor */
void VTKUGridT::SetScaleFactor(float factor)
{
	if (fWarp) fWarp->SetScaleFactor(factor);
}

/* set grid representation.
 * \param code grid representation, either VTK_SURFACE, VTK_WIRE, or VTK_POINTS */
bool VTKUGridT::SetRepresentation(RepresentationT rep)
{
	vtkProperty* property = fActor->GetProperty();
	switch (rep)
	{
		case kWire:
			property->SetRepresentation(VTK_WIREFRAME);	
			break;
		case kSurface:
			property->SetRepresentation(VTK_SURFACE);	
			break;
		case kPoint:
			property->SetRepresentation(VTK_POINTS);	
			break;		
		default:
			cout << "VTKUGridT::SetRepresentation: not a valid representation: " << rep << endl;
			return false;
	}
	return true;
}

/* set the grid opacity */
void VTKUGridT::SetOpacity(double opacity)
{
	/* bounds */
	if (opacity > 1) opacity = 1;
	else if (opacity < 0 ) opacity = 0;

	/* set */
	fActor->GetProperty()->SetOpacity(opacity);
}

/* return the look up table for the specified ugrid */
vtkScalarsToColors* VTKUGridT::GetLookupTable(void)
{
	return fMapper->GetLookupTable();
}
