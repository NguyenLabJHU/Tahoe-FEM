#ifndef _VTK_BODY_DATA_T_H_
#define _VTK_BODY_DATA_T_H_

/* direct members */
#include "StringT.h"
#include "iConsoleObjectT.h"
#include "ExodusT.h"
#include "VTKBodyT.h"

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
class vtkWarpVector;
class vtkVectors;
class vtkScalarBarActor;
class ExodusT;
class VTKBodyT;
class vtkDataArray;

class VTKBodyDataT: public iConsoleObjectT
{
 public:

  /** constructor */
   VTKBodyDataT(const StringT& file_name); 

  /** destructor */
  ~VTKBodyDataT(void);
  
  void SetLookupTable(void);
  void UpdateData(void);
  void DefaultValues(void);
  bool ChangeVars(const StringT& var);
  void SelectTimeStep(int);
  void ChangeDataColor(int);
  void AddVars(void);
 
  const ArrayT<StringT>& NodeLabels(void) const { return node_labels; };

  /** return pointer to actor for the body */
  vtkActor* Actor(void) { return ugridActor; };

  vtkActor* WireActor(void) { return wireActor;};

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
  vtkWarpVector *warp;
  vtkUnstructuredGrid *ugrid;
  
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
  vtkActor *wireActor;

#ifdef __VTK_NEW__
  vtkDataArray *scalars [100][20];
  vtkDataArray *vectors [100][20];
#else 
  vtkScalars *scalars [100][20];
  vtkVectors *vectors [100][20];
#endif 
/*   ExodusT exo; */
};

#endif
