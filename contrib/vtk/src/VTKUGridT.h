/* $Id: VTKUGridT.h,v 1.1 2001-12-10 12:44:09 paklein Exp $ */
#ifndef _VTK_U_GRID_T_H_
#define _VTK_U_GRID_T_H_

/* direct members */
#include "ArrayT.h"
#include "GeometryT.h"

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

	/** set the vector data */
	void SetVectors(vtkFloatArray* vectors);

	/** set vectors that warp */
	void SetWarpVectors(vtkFloatArray* vectors);
  
	/** return the grid actor */
  	vtkActor* Actor(void) { return fActor; };
  	
  	/** return the grid wrap vector */
  	vtkWarpVector* Warp(void) { return fWarp; };
  	
  	/** return the unstructured grid */
  	vtkUnstructuredGrid* UGrid(void) { return fUGrid; };
  
 private:

 	/** type of the unstructured grid set */
 	TypeT fType;
 
 	/** set ID */
 	int fID;

	/** dimensionality of the cells */
	int fNumSD; 
	
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
