/* $Id: VTKBodyT.h,v 1.2 2001-10-25 21:40:19 recampb Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"

/* forward declarations */
class vtkPoints;
class vtkCellArray;
class vtkUnstructuredGrid;
class vtkDataSetMapper;
class vtkActor;
class vtkLookupTable;
class vtkIdFilter;
class vtkSelectVisiblePoints;
class vtkLabeledDataMapper;
class vtkActor2D;
class vtkScalars;
class vtkCamera;
class vtkWarpVector;
class vtkVectors;

class VTKBodyT: public iConsoleObjectT
{
 public:

  /** constructor */
   VTKBodyT(const StringT& file_name); 

  /** destructor */
  ~VTKBodyT(void);

  /** return pointer to actor for the body */
  vtkActor* Actor(void) { return ugridActor; };

  //private:

  /* source file */
  const StringT inFile;

  /* model data */
  int num_nodes;
  int num_dim;
  vtkPoints *points;
 /*  vtkCellArray *cells; */
  vtkCellArray *vtk_cell_array;
  int num_node_variables;
  double hueRange1, hueRange2;
  double satRange1, satRange2;
  double valRange1, valRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1[100], scalarRange2[100];
  double scale_factor;
  double time;
  int numColors;
  int num_time_steps;
  ArrayT<StringT> node_labels;
  int currentVarNum;
  int frameNum;
  StringT output_file;
  StringT outFileName;
  StringT sbTitle;
  StringT varList;
  vtkLookupTable *lut;
  vtkDataSetMapper *ugridMapper;
  vtkActor *ugridActor;
  vtkIdFilter *ids;
  vtkSelectVisiblePoints *visPts;
  vtkLabeledDataMapper *ldm;
  vtkActor2D *pointLabels;
  vtkUnstructuredGrid *ugrid;
  vtkScalars *scalars [100][20];
  vtkVectors *vectors [100][20];
  vtkWarpVector *warp;
};

#endif
