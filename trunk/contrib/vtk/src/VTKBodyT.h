/* $Id: VTKBodyT.h,v 1.5 2001-11-06 02:39:51 recampb Exp $ */

#ifndef _VTK_BODY_T_H_
#define _VTK_BODY_T_H_

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"
#include "ExodusT.h"


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
class vtkScalarBarActor;
class ExodusT;
class VTKBodyT
{
 public:

  /** constructor */
   VTKBodyT(const StringT& file_name); 

  /** destructor */
  ~VTKBodyT(void);
  
  void SetLookupTable(void);
  void UpdateData(void);
  void DefaultValues(void);
  void ChangeVars(const int);
  void SelectTimeStep(const int);

  /** return pointer to actor for the body */
  vtkActor* Actor(void) { return ugridActor; };

  vtkScalarBarActor* SBActor(void) {return scalarBar;}; 
  
/*   vtkDataSetMapper *Mapper(void) {return ugridMapper;}; */
/*   StringT SBTitle(void) {return sbTitle;}; */
/*   int CurrentVar (void) {return currentVarNum;}; */
/*   StringT NodeLabels(const int i) {return node_labels[i];}; */
  int num_node_variables;
 const StringT inFile;
  StringT varList;
int num_time_steps;
  int currentVarNum;
  int currentStepNum;

 private:
  
  /* source file */


  /* model data */
  int num_nodes;
  int num_dim;
  vtkPoints *points;
 /*  vtkCellArray *cells; */
  vtkCellArray *vtk_cell_array;
 
  double hueRange1, hueRange2;
  double satRange1, satRange2;
  double valRange1, valRange2;
  double alphaRange1, alphaRange2;
  double scalarRange1[100], scalarRange2[100];
  double scale_factor;
  double time;
  int numColors;
  
  ArrayT<StringT> node_labels;

  StringT output_file;
  
  StringT sbTitle;

  vtkScalarBarActor *scalarBar;
  vtkLookupTable *lut;
  vtkDataSetMapper *ugridMapper;
  vtkActor *ugridActor;
  vtkIdFilter *ids;//     //TEMP - construct one body
//     try {
//       VTKBodyT* body2 = new VTKBodyT(inFile);
//       fBodies.Append(body2);
//     }
//     catch (int) {
//       cout << "\n exception constructing body" << endl;
//     }

  vtkSelectVisiblePoints *visPts;
  vtkLabeledDataMapper *ldm;
  vtkActor2D *pointLabels;
  vtkUnstructuredGrid *ugrid;
  vtkScalars *scalars [100][20];
  vtkVectors *vectors [100][20];
  vtkWarpVector *warp;
/*   ExodusT exo; */
};

#endif
