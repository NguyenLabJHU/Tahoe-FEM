/* $Id: VTKBodyT.h,v 1.1 2001-10-24 18:20:36 paklein Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* direct members */
#include "StringT.h"

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

class VTKBodyT
{
 public:

  /** constructor */
  VTKBodyT(const StringT& file_name);

  /** destructor */
  ~VTKBodyT(void);

  /** return pointer to actor for the body */
  vtkActor* Actor(void) { return ugridActor; };

  private:

  /* source file */
  const StringT inFile;

  /* model data */
  int num_nodes;
  int num_dim;
  vtkPoints *points;
  vtkCellArray *cells;
  int num_node_variables;

  double scalarRange1[100], scalarRange2[100];
  double scale_factor;
  int numColors;
  int num_time_steps;
  ArrayT<StringT> node_labels;
  int currentVarNum;
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
  vtkScalars *scalars [1000][100];
  vtkVectors *vectors [1000][100];
  vtkWarpVector *warp;
};

#endif
