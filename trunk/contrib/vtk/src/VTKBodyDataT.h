/* $Id: VTKBodyDataT.h,v 1.23 2003-11-25 19:55:48 paklein Exp $ */
#ifndef _VTK_BODY_DATA_T_H_
#define _VTK_BODY_DATA_T_H_

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "StringT.h"
#include "Array2DT.h"
#include "dArrayT.h"
#include "IOBaseT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

/* forward declarations */
class vtkPoints;
class vtkRenderer;
class VTKUGridT;
class vtkFloatArray;
class vtkRenderer;

using namespace Tahoe;

/** interface for model data. A single VTKBodyDataT is constructed
 * for each body (unless loaded more than once). This data may appear
 * in multiple frames, each instance with its own VTKBodyT. */
class VTKBodyDataT: public iConsoleObjectT
{
public:

	/** constructor */
	VTKBodyDataT(IOBaseT::FileTypeT format, const StringT& file_name); 

	/** destructor */
	~VTKBodyDataT(void);

 	/** return the number of spatial dimensions */
 	int NumSD(void);
 
	/** return the source file for the body data */
	const StringT& SourceFile(void) const { return fInFile; };

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

	/** return current time */
	double CurrentTime(void) const;

	/** return current step number */
	int CurrentStepNumber(void) const { return currentStepNum; };

	/** return current variable number */
	int CurrentVariableNumber(void) const { return currentVarNum; };
	
	/** return tbe number of nodal variables */
	int NumNodeVariables(void) const { return fScalars.MinorDim(); };

 	double CurrentScalarRange1(void) { return scalarRange1[currentVarNum];}; 
	double CurrentScalarRange2(void) { return scalarRange2[currentVarNum];}; 

	/** return a reference to the nodal labels */
	const ArrayT<StringT>& NodeLabels(void) const { return fNodeLabels; };
 
 	/** \name unstructured grid information */
 	/*@{*/
	/** return array of unstructured grid displays */
	const ArrayT<VTKUGridT*>& UGrids(void) { return fUGrids; };

	/** return array of unstructured grid names */
	const ArrayT<StringT>& UGridNames(void) { return fUGridNames; };
 
	/** visibility flags of grids */
	ArrayT<bool>& UGridVisible(void) { return fUGridVisible; };
	
	/** update visibility */
	void UpdateVisibility(void);
 	/*@}*/
 
 	/** execute console command. \return true is executed normally */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

	/** return a reference to the point numbering map */
	const iArrayT& PointNumberMap(void) const { return fPointNumberMap; };

	const dArray2DT& Coordinates(void) const { return fCoords; };

	/** specified vector field for the current step number. Returns NULL if
	 * the vector field does not exist */
	vtkFloatArray* VectorField(const char* name);

	Array2DT<vtkFloatArray*> getScalars(void) const {return fScalars;};

	vtkPoints* GetPoints(void) {return fPoints;};

	void ShowContours(vtkRenderer* renderer);
	void HideContours(vtkRenderer* renderer);

	/** return the index of the given node number in the node number map
	 * or -1 if the node is not present in the map */
	int NodeMapIndex(int node) const;

protected:

	/** return the index of specified vector field, -1 if not present */
	int VectorFieldNumber(const char* name) const;

 private:
 
	/** array type conversion */
	void double_to_float(const ArrayT<double>& d, float* f) const;
	
	/** set run time variables to defaults */
	void DefaultValues(void);

	/** load data for the current time step into the VTKBodyDataT::fScalars
	 * and the VTKBodyDataT::fVectors arrays. If the data is already stored,
	 * nothing is done. If the data is not already stored, it is read from
	 * the database file. */
	void LoadData(int step);
	Array2DT<vtkFloatArray*> fScalars; /**< dimension: [time_steps] x [num_vars] : [num_nodes] x [1] */

	/** determine the number and dimension of the vector fields. Scans labels
	 * looking for vector variables. Vector variables are any variables with
	 * labels ending following:
	 *     [root]_[suffix] 
	 * The root and suffix can be anything, but it is assumed vector variables
	 * with the same root will be grouped together in the list of labels. */
	void SetVectorFields(const ArrayT<StringT>& labels, ArrayT<StringT>& field, 
		iArrayT& dimension, iArrayT& index) const;
	  
 private:
 
 	/** source file format */
 	IOBaseT::FileTypeT fFormat;

	/** source file */
	const StringT fInFile;

	/** body coordinates */
	dArray2DT fCoords;
	
	/** point coordinates */
	vtkPoints* fPoints;
	
	/** point numbering map */
	iArrayT fPointNumberMap;

	/** array of output times */
	dArrayT fTimeList;

	/** scalar data per node */
	ArrayT<StringT> fNodeLabels;

	/** \name vector fields information */
	/*@{*/
	/** vector field data for each (loaded) time step.
	 * dimensions: [time_steps] x [num_fields] : [num_nodes x field_dim] */
	Array2DT<vtkFloatArray*> fVectorFields; 
	
	/** dimension of each of the vector fields */
	iArrayT fVectorFieldDim;

	/** index of the first component of each field in VTKBodyDataT::fNodeLabels */
	iArrayT fVectorFieldDex;

	/** root of the label for each vector field */
	ArrayT<StringT> fVectorFieldLabel;
	/*@}*/

	/** \name unstructured grids */
	/*@{*/
	/** array of unstructured grid displays */
	ArrayT<VTKUGridT*> fUGrids;

	/** array of names in VTKBodyDataT::fUGrids */
	ArrayT<StringT> fUGridNames;

	/** visibility of each member of VTKBodyDataT::fUGrids */
	ArrayT<bool> fUGridVisible;
	/*@}*/
	
	vtkRenderer* bodyDataRenderer;

	/** \name runtime data */
	/*@{*/
  	int currentStepNum; /**< current time step on display */
  	int currentVarNum;  /**< current display variable */
	dArrayT scalarRange1, scalarRange2;
	double scale_factor;
	int numColors;
	int numContours;
	int numColorBarLabels;
	double opacity;
	double boundingOpacity;
	/*@}*/
};

/* type conversion */
inline void VTKBodyDataT::double_to_float(const ArrayT<double>& d, float* f) const
{
	int len = d.Length();
	const double* pd = d.Pointer(); 
 	for (int i = 0; i < len; i++)
 		*f++ = float(*pd++);
}

/* specified vector field for the current step number */
inline vtkFloatArray* VTKBodyDataT::VectorField(const char* name)
{
	int dex = VectorFieldNumber(name);
	if (dex > -1)
		return fVectorFields(currentStepNum, dex);
	else
		return NULL;
}

/* return current time */
inline double VTKBodyDataT::CurrentTime(void) const
{
	if (fTimeList.Length() == 0) 
		return 0.0;
	else
		return fTimeList[currentStepNum];
}

#endif
