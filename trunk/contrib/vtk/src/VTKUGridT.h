/* $Id: VTKUGridT.h,v 1.4 2002-02-01 18:11:41 paklein Exp $ */
#ifndef _VTK_U_GRID_T_H_
#define _VTK_U_GRID_T_H_

/* direct members */
#include "ArrayT.h"
#include "GeometryT.h"
#include "iArrayT.h"

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

/* toolbox forward declarations */
class iArray2DT;

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
 
 	/** return a reference to the cell numbering map */
	const iArrayT& CellNumberMap(void) const { return fCellNumberMap; };
	
	/** set the cell number map */
	void SetCellNumberMap(const iArrayT& map) { fCellNumberMap = map; };

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
};

#endif /* _VTK_U_GRID_T_H_ */
