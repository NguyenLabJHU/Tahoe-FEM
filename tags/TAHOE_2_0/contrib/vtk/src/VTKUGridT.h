/* $Id: VTKUGridT.h,v 1.21 2003-01-08 22:20:15 rjones Exp $ */
#ifndef _VTK_U_GRID_T_H_
#define _VTK_U_GRID_T_H_

/* direct members */
#include "ArrayT.h"
#include "GeometryT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "AutoArrayT.h"

using namespace Tahoe; 

/* VTK forward declarations */
class vtkPoints;
class vtkCellArray;
class vtkUnstructuredGrid;
class vtkDataSetMapper;
class vtkDataSet;
class vtkActor;
class vtkWarpVector;
class vtkFloatArray;
class vtkIdTypeArray;
class vtkLookupTable;
class vtkScalarsToColors;
class vtkContourFilter;
class vtkPolyDataMapper;
class vtkOutlineFilter;
class vtkExtractEdges;
class vtkLODActor;
class vtkPlane;
class vtkCutter;
class vtkSmoothPolyDataFilter;
class vtkDataSetToPolyDataFilter;
class vtkTahoeGlyph3D;
class vtkArrowSource;
class vtkSelectVisiblePoints;
class vtkRenderer;

class VTKBodyDataT;
namespace Tahoe {
	class StringT;
}

/* toolbox forward declarations */
namespace Tahoe {
class iArray2DT;
}

/** interface for display of unstructured grid data. */
class VTKUGridT
{
 public:
 
	/** ugrid types */
	enum TypeT {kElementSet = 0,
                   kNodeSet = 1,
                   kSideSet = 2};

	/** representation */
	enum RepresentationT {kWire, kSurface, kPoint};

	/** constructor */
	VTKUGridT(TypeT my_type, int id, int nsd); 

	/** destructor */
	~VTKUGridT(void);
 
	/** return the grid type */
	TypeT Type(void) const { return fType; }; 
 	
	/** return the grid ID */
	int ID(void) const { return fID; };

	/** return the dimensionality of the grid */
	int NumSD(void) const { return fNumSD; };
  
	/** set the point data */
	void SetPoints(vtkPoints* points);

	/** set the connectivities */
	void SetConnectivities(GeometryT::CodeT code, const iArray2DT& connects);
  
	/** set the scalar data */
	void SetScalars(vtkFloatArray* scalars);

	/** display data as contour surfaces in 3D or contour lines in 2D */
	void ShowContours(vtkFloatArray* scalars, int numContours, double min, double max, vtkRenderer* renderer);
	
	/** remove contours */
	void HideContours(vtkRenderer* renderer);

	void CuttingPlane(vtkRenderer* renderer, double oX, double oY, double oZ, double nX, double nY, double nZ, bool warp, double min, double max);

	void HideCuttingPlane(vtkRenderer* renderer);

	void Glyphing(VTKBodyDataT* body, const Tahoe::StringT& field, vtkRenderer* renderer, bool filter, bool warpA, bool scale, bool color);
	void SetGlyphScale(double scale);
	void HideGlyphing(void);

	/** set the scalar data range */
	void SetScalarRange(double min, double max);
	
	/** set the number of color levels */
	void SetNumberOfColors(int num);

	/** set the vector data */
	void SetVectors(vtkFloatArray* vectors);

	/** set vectors that warp */
	void SetWarpVectors(vtkFloatArray* vectors);
  
	/** return the grid actor */
  	vtkActor* Actor(void) { return fActor; };

	/** return the box outline actor */
	vtkActor* OutlineActor(void) { return fOutlineActor;};
	
	/** return the bounding wire-frame actor */
	vtkActor* EdgesActor(void) {return edgesActor;};
	
	/** return the semi-transparent bounding volume actor */
	vtkActor* BoundBoxActor(void) {return boundBoxActor;};

	vtkActor* SpikeActor(void) { return spikeActor;};
  	
  	/** return the grid wrap vector */
  	vtkWarpVector* Warp(void) { return fWarp; };
  	
  	/** set the wrap displacement scale factor */
  	void SetScaleFactor(float factor);
  	
  	/** return the unstructured grid */
  	vtkUnstructuredGrid* UGrid(void) { return fUGrid; };
  
   	/** return the look up table for the specified ugrid */
 	vtkScalarsToColors* GetLookupTable(void);
 
 	/** set grid representation.
 	 * \param code grid representation */
 	bool SetRepresentation(RepresentationT rep);
 
 	/** set the grid opacity.
 	 * \param opacity ranges from 0 to 1 for transparent to opaque */
	void SetOpacity(double opacity);

	/** set bounding volume opacity
	 * \param boundingOpacity ranges from 0 to 1 for transparent to opaque */
	void SetBoundingOpacity(double boundingOpacity);
 
 	/** return a reference to the cell numbering map */
	const iArrayT& CellNumberMap(void) const { return fCellNumberMap; };
	
	/** set the cell number map */
	void SetCellNumberMap(const iArrayT& map) { fCellNumberMap = map; };

	bool GetContoursBool(void) { return contours;};

 private:

 	/** type of the unstructured grid set */
 	TypeT fType;
 
 	/** set ID */
 	int fID;

	/** dimensionality of the cells */
	int fNumSD; 

	/** cell numbering map */
	iArrayT fCellNumberMap;
	
	/** connectivities in VTK format */
  	vtkCellArray* fCellArray;
  	vtkIdTypeArray* fConnects; /**< connectivities */
  	ArrayT<int>     fCellType; /**< VTK type of each cell */

	/** grid object */
	vtkUnstructuredGrid* fUGrid;

	/** grid mapper */
	vtkDataSetMapper* fMapper; //or does this belong in VTKBodyT?
	                           //if you want different colors in different frames

	/** color look-up table */
	vtkLookupTable* fLookUpTable; //or does this belong in VTKBodyT?
	                           //if you want different colors in different frames

	/** actor */
	vtkActor* fActor;

	/** displaces grid */
	vtkWarpVector* fWarp;

	/** contour lines/surfaces */
	vtkContourFilter* fContour;
	vtkPolyDataMapper* fContourMapper;
	vtkActor* fContourActor;
	bool contours;
	vtkSmoothPolyDataFilter* smoother;
	vtkDataSetToPolyDataFilter* dsToPd;
	
	/** bounding box outline */
	vtkOutlineFilter* outline;
	vtkDataSetMapper* outlineMapper;
	vtkActor* fOutlineActor;

	/** bounding volume wireframe */
	/* not used at the moment */
	vtkExtractEdges* edges;
	vtkDataSetMapper* edgesMapper;
	vtkActor* edgesActor;

	/** bounding semi-transparent volume */
	vtkDataSetMapper* boundBoxMapper;
	vtkActor* boundBoxActor;

	/** cutting plane variables */
	AutoArrayT<vtkPlane*> plane;
	AutoArrayT<vtkCutter*> cutter;
	AutoArrayT<vtkCutter*> cutter2;
	AutoArrayT<vtkActor*> cut;
	AutoArrayT<vtkPolyDataMapper*> cutterMapper;
	AutoArrayT<vtkPolyDataMapper*> boundPlaneMapper;
	AutoArrayT<vtkContourFilter*> contourA;
	AutoArrayT<vtkPolyDataMapper*> contourMapperA;
	AutoArrayT<vtkActor*> boundPlane;

	/** glyphing variables */
	vtkTahoeGlyph3D* glyph;
	vtkPolyDataMapper* spikeMapper;
	vtkActor* spikeActor;
	vtkArrowSource* cone;
	vtkSelectVisiblePoints* visPoints;
	bool glyphFilter;
	bool warpBool;
	bool warpArrows;
	bool cutting;
	
      
};

#endif /* _VTK_U_GRID_T_H_ */
