/* $Id: VTKBodyDataT.h,v 1.5 2001-12-10 12:44:08 paklein Exp $ */
#ifndef _VTK_BODY_DATA_T_H_
#define _VTK_BODY_DATA_T_H_

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "StringT.h"
#include "Array2DT.h"
#include "dArrayT.h"

/* forward declarations */
class vtkPoints;
class vtkRenderer;
class VTKUGridT;
class vtkFloatArray;

/** interface for model data. A single VTKBodyDataT is constructed
 * for each body (unless loaded more than once). This data may appear
 * in multiple frames, each instance with its own VTKBodyT. */
class VTKBodyDataT: public iConsoleObjectT
{
public:

	/** constructor */
	VTKBodyDataT(const StringT& file_name); 

	/** destructor */
	~VTKBodyDataT(void);

 	/** return the number of spatial dimensions */
 	int NumSD(void);
 
	/** return the source file for the body data */
	const StringT& SourceFile(void) const { return fInFile; };
  
#if 0
	void SetLookupTable(void);
	void ChangeDataColor(int);
#endif 

	/** update the data state */
	void UpdateData(void);

	/** change the plot variable */
	bool ChangeVars(const StringT& var);

	/** set the current time step */
	bool SelectTimeStep(int);

 	/** add actors in self to the given renderer */
 	void AddToRenderer(vtkRenderer* renderer) const;

 	/** add actors in self to the given renderer */
 	void RemoveFromRenderer(vtkRenderer* renderer) const;
 
	/** return the number of time steps of results data */
	int NumTimeSteps(void) const { return fScalars.MajorDim(); };

	/** return current step number */
	int CurrentStepNumber(void) const { return currentStepNum; };
	
	/** return tbe number of nodal variables */
	int NumNodeVariables(void) const { return fScalars.MinorDim(); };
	
	/** return a reference to the nodal labels */
	const ArrayT<StringT>& NodeLabels(void) const { return fNodeLabels; };
 
 	/** show node numbers */
// 	void ShowNodeNumbers(vtkRenderer* renderer);

 	/** show node numbers */
//	void HideNodeNumbers(vtkRenderer* renderer);

	/** return array of unstructured grid displays */
	const ArrayT<VTKUGridT*>& UGrids(void) { return fUGrids; };
 
 private:
 
	/** array type conversion */
	void double_to_float(const ArrayT<double>& d, float* f) const;
	
	/** set run time variables to defaults */
	void DefaultValues(void);
	  
 private:

	/** source file */
	const StringT fInFile;
  
	/** point coordinates */
	vtkPoints* fPoints;

	/* scalar data per node */
	ArrayT<StringT> fNodeLabels; /**< labels for the nodal output variables */
	Array2DT<vtkFloatArray*> fScalars; /**< dimension: [time_steps] x [num_vars] : [num_nodes] x [1] */
	
	/** vector data per node. 
	 * dimension: [time_steps] : [num_nodes] x [ndof] */
	ArrayT<vtkFloatArray*> fVectors; 

	/** array of unstructured grid displays */
	ArrayT<VTKUGridT*> fUGrids;

	/* runtime data */
  	int currentStepNum; /**< current time step on display */
  	int currentVarNum;  /**< current display variable */
	dArrayT scalarRange1, scalarRange2;
	double hueRange1, hueRange2;
	double satRange1, satRange2;
	double valRange1, valRange2;
	double alphaRange1, alphaRange2;
	double scale_factor;
	int numColors;
};

/* type conversion */
inline void VTKBodyDataT::double_to_float(const ArrayT<double>& d, float* f) const
{
	int len = d.Length();
	double* pd = d.Pointer(); 
 	for (int i = 0; i < len; i++)
 		*f++ = float(*pd++);
}

#endif
