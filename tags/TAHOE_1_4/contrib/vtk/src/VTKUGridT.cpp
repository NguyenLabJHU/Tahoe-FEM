/* $Id: VTKUGridT.cpp,v 1.30 2004-01-02 04:26:34 paklein Exp $ */
#include "VTKUGridT.h"

/* Tahoe toolbox headers */
#include "ExceptionCodes.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

/* VTK headers */
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkWarpVector.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkContourFilter.h"
#include "vtkOutlineFilter.h"
#include "vtkExtractEdges.h"
#include "vtkLODActor.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkDataSetToPolyDataFilter.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkGlyph3D.h"
#include "vtkConeSource.h"
#include "vtkArrowSource.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkRenderer.h"
#include "vtkTahoeGlyph3D.h"
#include "vtkActorCollection.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"

using namespace Tahoe;

/* array behavior */
const bool ArrayT<VTKUGridT*>::fByteCopy = true;
const bool ArrayT<vtkActor*>::fByteCopy = true;
const bool ArrayT<vtkPlane*>::fByteCopy = true;
const bool ArrayT<vtkCutter*>::fByteCopy = true;
const bool ArrayT<vtkPolyDataMapper*>::fByteCopy = true;
const bool ArrayT<vtkContourFilter*>::fByteCopy = true;


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

	fContour = vtkContourFilter::New();
	fContourMapper = vtkPolyDataMapper::New();
	fContourActor = vtkActor::New();
	
// 	dsToPd = vtkDataSetToPolyDataFilter::New();
// 	dsToPd->SetInput(fUGrid);
// 	smoother = vtkSmoothPolyDataFilter::New();
// 	smoother->SetInput(dsToPd->GetOutput());
// 	smoother->SetNumberOfIterations(1000);
// 	fMapper->SetInput(smoother->GetOutput());

	outline = vtkOutlineFilter::New();
	outline->SetInput(fUGrid);
	outlineMapper = vtkDataSetMapper::New();
	outlineMapper->SetInput(outline->GetOutput());

	edges = vtkExtractEdges::New();
	edgesMapper = vtkDataSetMapper::New();
	edges->SetInput(fUGrid);
	edgesMapper->SetInput(edges->GetOutput());
	
	edgesActor = vtkActor::New();
	edgesActor->GetProperty()->SetColor(1,1,1);
	edgesActor->SetMapper(edgesMapper);

	boundBoxMapper = vtkDataSetMapper::New();
	boundBoxMapper->SetInput(fUGrid);
	boundBoxMapper->ScalarVisibilityOff();

	boundBoxActor = vtkActor::New();
	boundBoxActor->SetMapper(boundBoxMapper);
	boundBoxActor->GetProperty()->SetOpacity(.27);
	boundBoxActor->SetVisibility(false);
	boundBoxActor->GetProperty()->SetColor(1,1,1);
	
	// fWarp = vtkWarpVector::New();

	glyph = vtkTahoeGlyph3D::New();
	cone = vtkArrowSource::New();
	
	visPoints = vtkSelectVisiblePoints::New();
	visPoints->SetInput(fUGrid);
	glyph->SetInput(fUGrid);
	glyph->SetSource(cone->GetOutput());
	glyph->SetVectorModeToUseVector();
	glyph->SetScaleModeToScaleByVector();
	spikeMapper = vtkPolyDataMapper::New();
	spikeMapper->SetInput(glyph->GetOutput());
	spikeActor = vtkActor::New();
	spikeActor->SetMapper(spikeMapper);
	spikeActor->GetProperty()->SetColor(0,1,0);
	spikeActor->SetVisibility(false);
	spikeActor->PickableOff();
	glyphFilter = true;
	warpBool = false;
	warpArrows = false;
	contours = false;
	cutting = false;
	
	
	/* change color range from blue to red */
	fLookUpTable = vtkLookupTable::New();
	fLookUpTable->SetHueRange(0.6667, 0);
	fMapper->SetLookupTable(fLookUpTable);
	fContourMapper->SetLookupTable(fLookUpTable);


	/* the actor */
	fActor = vtkActor::New();
	fActor->GetProperty()->SetPointSize(3.0);
	//fActor = vtkLODActor::New();
	//fActor->GetProperty()->SetInterpolationToGouraud();

	/* line color */
	if (fType == kElementSet)
		fActor->GetProperty()->SetColor(1,0,0);
	else if (fType == kNodeSet)
		fActor->GetProperty()->SetColor(0,0,1);
	else
		fActor->GetProperty()->SetColor(0,1,0);
	fActor->SetMapper(fMapper);

	fOutlineActor = vtkActor::New();
	fOutlineActor->SetMapper(outlineMapper);
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
	if (fContour) fContour->Delete();
	if (fContourMapper) fContourMapper->Delete();
	if (outlineMapper) outlineMapper->Delete();
	if (outline) outline->Delete();
	if (fOutlineActor) fOutlineActor->Delete();
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
			const int* a = connects(i);
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

/* show contour surfaces for 3D or contour lines for 2D */
void VTKUGridT::ShowContours(vtkFloatArray* scalars, int numContours, double min, double max, vtkRenderer* renderer)
{
  vtkActorCollection* temp = renderer->GetActors();
  
    if (!warpBool)
      fContour->SetInput(fUGrid);
    else
      fContour->SetInput(fWarp->GetOutput());
    fContourMapper->SetInput(fContour->GetOutput());  
    fContour->GenerateValues(numContours+2, min, max);
    fContourMapper->SetScalarRange(min,max);
   
    //fActor->SetMapper(fContourMapper);
    fContourActor->SetMapper(fContourMapper);
    if (!cutting){    
      boundBoxActor->SetVisibility(true);
      boundBoxActor->PickableOff();
      cout << "Contour Values:" << endl;
      for (int i=0; i<numContours+2; i++) 
	cout << i <<"  " << fContour->GetValue(i) << endl;
      contours = true;
      renderer->RemoveActor(fActor);
      if (temp->IsItemPresent(fContourActor) == 0)
	renderer->AddActor(fContourActor);
      if (temp->IsItemPresent(boundBoxActor) == 0)
	renderer->AddActor(boundBoxActor);
    }
  
  else
    {
  
    for (int i=0; i<cutter.Length(); i++)
      {
	vtkContourFilter* tContour = vtkContourFilter::New();
	vtkPolyDataMapper* tContourMapper = vtkPolyDataMapper::New();

	tContourMapper->SetLookupTable(fLookUpTable);
		
	contourA.Append(tContour);
	contourMapperA.Append(tContourMapper);
	tContour->SetInput(cutter[i]->GetOutput());
	tContour->GenerateValues(numContours+2, min, max);
	tContourMapper->SetInput(tContour->GetOutput());
	tContourMapper->SetScalarRange(min, max);
	vtkActor* tActor = vtkActor::New();
	boundPlane.Append(tActor);
	boundPlane[i]->SetMapper(cutterMapper[i]);

	cutterMapper[i]->ScalarVisibilityOff();

	boundPlane[i]->GetProperty()->SetOpacity(.20);
	boundPlane[i]->GetProperty()->SetColor(1,1,1);
	cut[i]->SetMapper(tContourMapper);
	
	
	if (temp->IsItemPresent(cut[i]) == 0)
	  renderer->AddActor(cut[i]);
	if (temp->IsItemPresent(boundPlane[i]) == 0)
	  renderer->AddActor(boundPlane[i]);
	
	cout << "Contour Values:" << endl;
	for (int i=0; i<numContours+2; i++) 
	  cout << i <<"  " << tContour->GetValue(i) << endl;
      }	
    
    contours = true;
    renderer->RemoveActor(fActor);
    
    
    }
  
}

/* hide contour surfaces */
void VTKUGridT::HideContours(vtkRenderer* renderer)
{
  contours = false;
  
  if (!cutting)
    {
      renderer->AddActor(fActor);
      fActor->SetVisibility(true);
      renderer->RemoveActor(fContourActor);
      renderer->RemoveActor(boundBoxActor);

      //boundBoxActor->SetVisibility(false);
    }
  
  else 
    {

      for (int i = 0; i < cut.Length(); i++)
	{
	  if (warpBool)
	      cutter[i]->SetInput(fWarp->GetOutput());
	  else
	    cutter[i]->SetInput(fUGrid);
	  cutterMapper[i]->ScalarVisibilityOn();
	  cut[i]->SetMapper(cutterMapper[i]);
	  renderer->RemoveActor(boundPlane[i]);

	}
      

    }
  
}


void VTKUGridT::CuttingPlane(vtkRenderer* renderer, double oX, double oY, double oZ,double nX, double nY, double nZ, bool warp, double min, double max)
{
  cutting = true;
  if (!contours)
    {
      vtkActorCollection* temp = renderer->GetActors();
      
      warpBool = warp;
      vtkPlane* tplane = vtkPlane::New();
      vtkCutter* tcutter = vtkCutter::New();
      vtkPolyDataMapper* tcutterMapper = vtkPolyDataMapper::New();
      vtkActor* tcut = vtkActor::New();
      tcut->SetMapper(tcutterMapper);
      tcutter->SetInput(fUGrid);
      tcutter->SetCutFunction(tplane);
      tcutterMapper->SetInput(tcutter->GetOutput());
      tcutterMapper->SetLookupTable(fLookUpTable);
      tcutterMapper->SetScalarRange(min, max);
      
      if (temp->IsItemPresent(boundBoxActor) == 0)
	renderer->AddActor(boundBoxActor);  
      
      boundBoxActor->SetVisibility(true);
      boundBoxActor->PickableOff();
      tplane->SetOrigin(oX, oY, oZ);
      tplane->SetNormal(nX, nY, nZ);
      plane.Append(tplane);
      cutter.Append(tcutter);
      cutterMapper.Append(tcutterMapper);
      cut.Append(tcut);
      
      //fActor->SetMapper(cutterMapper);
      
      if (warp)
	for (int i=0; i<cutter.Length(); i++)
	  {
	    cutter[i]->SetInput(fWarp->GetOutput());
	    cutterMapper[i]->SetInput(cutter[i]->GetOutput());
	  }  
      
      renderer->RemoveActor(fActor);
      //fActor->SetVisibility(false);
      renderer->AddActor(tcut);
      
    }

else
  {
    vtkActorCollection* temp = renderer->GetActors();
    cutting = true;
    warpBool = warp;
    vtkPlane* tplane = vtkPlane::New();
    vtkCutter* tcutter = vtkCutter::New();
    vtkCutter* tcutter2 = vtkCutter::New();
    vtkPolyDataMapper* tcutterMapper = vtkPolyDataMapper::New();
    vtkPolyDataMapper* tboundPlaneMapper = vtkPolyDataMapper::New();
    vtkActor* tcut = vtkActor::New();
    tcutter->SetInput(fContour->GetOutput());
    tcutter->SetCutFunction(tplane);
    tcutter2->SetCutFunction(tplane);
    if (!warpBool)
      tcutter2->SetInput(fUGrid);
    else
      tcutter2->SetInput(fWarp->GetOutput());
    tcutterMapper->SetInput(tcutter->GetOutput());
    tcutterMapper->SetLookupTable(fLookUpTable);
    tcutterMapper->SetScalarRange(min, max);
    tboundPlaneMapper->SetInput(tcutter2->GetOutput());
    tboundPlaneMapper->SetLookupTable(fLookUpTable);
    tboundPlaneMapper->SetScalarRange(min, max);
    tboundPlaneMapper->ScalarVisibilityOff();
    tcutterMapper->ScalarVisibilityOn();
    tcut->SetMapper(tcutterMapper);
    
    if (temp->IsItemPresent(boundBoxActor) == 0)
      renderer->AddActor(boundBoxActor);  

    vtkActor* tActor = vtkActor::New();
    tActor->SetMapper(tboundPlaneMapper);
    tActor->GetProperty()->SetOpacity(.20);
    tActor->GetProperty()->SetColor(1,1,1);    
    boundPlane.Append(tActor);
    renderer->AddActor(tActor);
    
    boundBoxActor->SetVisibility(true);
    boundBoxActor->PickableOff();
    tplane->SetOrigin(oX, oY, oZ);
    tplane->SetNormal(nX, nY, nZ);
    plane.Append(tplane);
    cutter.Append(tcutter);
    cutter2.Append(tcutter2);
    boundPlaneMapper.Append(tboundPlaneMapper);
    cutterMapper.Append(tcutterMapper);
    cut.Append(tcut);
    
    
    
      renderer->RemoveActor(fActor);
      //fActor->SetVisibility(false);
      renderer->AddActor(tcut);   
      renderer->RemoveActor(fContourActor);
 
  }

  
}

void VTKUGridT::HideCuttingPlane(vtkRenderer* renderer)
{
  cutting = false;
  if (!contours)
    {
      fActor->SetMapper(fMapper);
      renderer->AddActor(fActor);
      fActor->SetVisibility(true);
      boundBoxActor->SetVisibility(false);
      renderer->RemoveActor(boundBoxActor);
      for (int i = 0; i < cut.Length(); i++)
	{
	  renderer->RemoveActor(cut[i]);

	}
    }
  else
    {
      for (int i = 0; i < cut.Length(); i++)
	{
	  renderer->RemoveActor(cut[i]);
	  renderer->RemoveActor(boundPlane[i]);

	}
     renderer->AddActor(fContourActor);
     //renderer->RemoveActor(boundBoxActor);
     fContourActor->SetVisibility(true);
     
     cut.Free();
     cutter.Free();
     cutter2.Free();
     boundPlaneMapper.Free();
     cutterMapper.Free();
     plane.Free();
      
    }




}

void VTKUGridT::Glyphing(VTKBodyDataT* body, const Tahoe::StringT& field, vtkRenderer* renderer, bool filter, bool warpA, bool scale, bool color) 
{
  glyphFilter = filter;
  warpArrows = warpA;
  visPoints->SetRenderer(renderer);
  if (filter)
    {
      if (warpBool && warpArrows)
	visPoints->SetInput(fWarp->GetOutput());
	//visPoints->SetInput(fUGrid);
      else  
	visPoints->SetInput(fUGrid);
	  
      glyph->SetInput(visPoints->GetOutput());     
    }
  
  else
    {
      if (warpBool && warpArrows)
	glyph->SetInput(fWarp->GetOutput());

      else
	glyph->SetInput(fUGrid);
      
    }
  

  if (scale)
    glyph->SetScaleModeToScaleByVector();
  else
    glyph->SetScaleModeToDataScalingOff();

  if (color){
    glyph->SetColorModeToColorByVector();
    spikeMapper->ScalarVisibilityOn();
  }  
  else
    spikeMapper->ScalarVisibilityOff();
  
  glyph->SetVectors(body, field);
  
  spikeActor->SetVisibility(true);
  spikeActor->PickableOff();
}

void VTKUGridT::SetGlyphScale(double scale)
{
	glyph->SetScaleFactor(scale);
}

void VTKUGridT::HideGlyphing(void)
{
  spikeActor->SetVisibility(false);
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

  if (!fWarp)
    fWarp = vtkWarpVector::New();
  /* insert in grid */
  SetVectors(vectors);
  
  
    /* set up warp vector */
    warpBool = true;
    fWarp->SetInput(fUGrid);

    fMapper->SetInput(fWarp->GetOutput());
    outline->SetInput(fWarp->GetOutput());
    outlineMapper->SetInput(outline->GetOutput());
    visPoints->SetInput(fWarp->GetOutput());


  if (glyphFilter)
    {
      if (warpArrows)
	visPoints->SetInput(fWarp->GetOutput());
      else  
	visPoints->SetInput(fUGrid);
      
      glyph->SetInput(visPoints->GetOutput());
  
    }

  else
    {
      if (warpArrows)
	glyph->SetInput(fWarp->GetOutput());
      else
	glyph->SetInput(fUGrid);
      
    }
    
    for (int i=0; i<cutter.Length(); i++)
      {
	cutter[i]->SetInput(fWarp->GetOutput());
	cutterMapper[i]->SetInput(cutter[i]->GetOutput());
      }
      
      fContour->SetInput(fWarp->GetOutput());
      fContourMapper->SetInput(fContour->GetOutput());
      edges->SetInput(fWarp->GetOutput());
      edgesMapper->SetInput(edges->GetOutput());
      boundBoxMapper->SetInput(fWarp->GetOutput());
 
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
	{
	  property->SetRepresentation(VTK_POINTS);	
	  property->SetPointSize(3.0);	
	  break;
	} 
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

void VTKUGridT::SetBoundingOpacity(double opacity)
{
  	/* bounds */
	if (opacity > 1) opacity = 1;
	else if (opacity < 0 ) opacity = 0;

	/* set */
	boundBoxActor->GetProperty()->SetOpacity(opacity);
}


/* return the look up table for the specified ugrid */
vtkScalarsToColors* VTKUGridT::GetLookupTable(void)
{
	return fMapper->GetLookupTable();
}

