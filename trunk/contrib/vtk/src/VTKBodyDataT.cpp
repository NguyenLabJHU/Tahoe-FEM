/* $Id: VTKBodyDataT.cpp,v 1.6 2001-12-13 02:57:59 paklein Exp $ */
#include "VTKBodyDataT.h"

#include "VTKUGridT.h"

#include <iostream.h>
#include <iomanip.h>
#include <float.h>

#include "vtkIdTypeArray.h"
#include "vtkPoints.h"
#include "vtkRenderer.h"
#include "vtkFloatArray.h"

#include "iArray2DT.h"
#include "ExodusT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "GeometryT.h"
#include "StringT.h"
#include "CommandSpecT.h"

/* array behavior */
const bool ArrayT<VTKBodyDataT*>::fByteCopy = true;

/* constructor */
VTKBodyDataT::VTKBodyDataT(const StringT& file_name): 
	fInFile(file_name),
	fPoints(NULL),
	currentStepNum(0)
{
	/* read exodus file */
	ExodusT exo(cout);
	if (!exo.OpenRead(fInFile))
    {
		cout << " ERROR: could not open file: " << fInFile << endl;
		throw eDatabaseFail;;
	}
	else
		cout << "read database file: " << fInFile << endl;
  
	/* read coordinates */
	int num_nodes = exo.NumNodes();
	int num_dim   = exo.NumDimensions();
  	dArray2DT coords(num_nodes, num_dim);
	exo.ReadCoordinates(coords); 
	if (coords.MinorDim() == 2) /* fill to 3D */
    { 
		/* temp space */ 
		dArray2DT tmp(coords.MajorDim(), 3); 
      
		/* write in */ 
		tmp.BlockColumnCopyAt(coords, 0);    
		tmp.SetColumn(2, 0.0); 
      
		/* swap memory */ 
		tmp.Swap(coords); 
    }

#if 1
	/* set up points */
  	fPoints = vtkPoints::New();
  	for (int i=0; i < num_nodes; i++) 
		fPoints->InsertPoint(i+1, coords(i));
//  	fPoints->InsertPoint(i, coords(i)); //SHIFT
#endif

#if 0
//NOTE: not the most efficient way to do things, but the code below
//      doesn't work properly.  		

	/* convert to float */
	nArray2DT<float> coords_float(num_nodes+1, 3);
	double_to_float(coords, coords_float(1));

	coordinates = vtkFloatArray::New();
	coordinates->SetNumberOfComponents(3);
	float* pcoords;
	coords_float.ReleasePointer(&pcoords);
	coordinates->SetArray(pcoords, coords.Length(), 0);
	fPoints = vtkPoints::New();
	fPoints->SetData(coordinates);
#endif

	/* dimensions */
	int num_elem_blocks = exo.NumElementBlocks();
	int num_node_sets = exo.NumNodeSets();
	fUGrids.Allocate(num_elem_blocks + num_node_sets);	
  
	/* load element connectivities */
  	iArrayT element_ID(num_elem_blocks);
  	exo.ElementBlockID(element_ID);
	for (int i = 0 ; i < element_ID.Length(); i++)
    {
		/* read dimensions */
		int num_elements, num_element_nodes;
		exo.ReadElementBlockDims(element_ID[i], num_elements, num_element_nodes);

#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading element block: " 
		     << num_elements << " x " << num_element_nodes << endl;
#endif

		/* read connectivities */
		iArray2DT connectivities(num_elements, num_element_nodes);
		GeometryT::CodeT geom_code;
		exo.ReadConnectivities(element_ID[i], geom_code, connectivities);      
		//connectivities[i]--; //SHIFT
		
		/* construct VTK grid */
		fUGrids[i] = new VTKUGridT(VTKUGridT::kElementSet, element_ID[i], num_dim);
		fUGrids[i]->SetPoints(fPoints);
		fUGrids[i]->SetConnectivities(geom_code, connectivities);
	}    
    cout << "read element blocks" << endl;

	/* load node sets */
	iArrayT node_ID(num_node_sets);
  	exo.NodeSetID(node_ID);
	for (int i = 0; i < node_ID.Length(); i++)
    {
		/* read dimensions */
		int num_nodes = exo.NumNodesInSet(node_ID[i]);

#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading node set: " 
		     << num_nodes << endl;
#endif

		/* read nodes */
		iArray2DT connectivities(num_nodes, 1);
		GeometryT::CodeT geom_code = GeometryT::kPoint;
		exo.ReadNodeSet(node_ID[i], connectivities);
		
		/* construct VTK grid */
		int ii = i + num_elem_blocks;
		fUGrids[ii] = new VTKUGridT(VTKUGridT::kNodeSet, node_ID[i], num_dim);
		fUGrids[ii]->SetPoints(fPoints);
		fUGrids[ii]->SetConnectivities(geom_code, connectivities);
	}    
    cout << "read node sets" << endl;
  
	/* number of results sets */
	int num_time_steps = exo.NumTimeSteps();

    cout << "read time steps" << endl;

	/* variables defined at the nodes */
	int num_node_variables = exo.NumNodeVariables();
	exo.ReadNodeLabels(fNodeLabels);
	vec_dim = num_dim;
	if (fNodeLabels.Length() >= vec_dim)
	{
		const char *d[] = {"D_X", "D_Y", "D_Z"};
		for (int i = 0; vec_dim > 0 && i < vec_dim; i++)
			if (fNodeLabels[i] != d[i])
				vec_dim = 0;

		//TEMP - other displacement variable names
		if (vec_dim == 0)
		{
			vec_dim = num_dim;
			const char *d[] = {"DISX", "DISY", "DISZ"};
			for (int i = 0; vec_dim > 0 && i < vec_dim; i++)
				if (fNodeLabels[i] != d[i])
					vec_dim = 0;	
		}
	}
	else vec_dim = 0;
	
	/* close file */
	exo.Close();
	
	/* results history */
	fScalars.Allocate(num_time_steps, num_node_variables);
	fScalars = NULL;
	if (vec_dim > 0)
	{
		fVectors.Allocate(num_time_steps);
		fVectors = NULL;
	}

//NOT USED YET
#if 0
	/* variables defined over the elements */
	int num_element_variables = exo.NumElementVariables();
	ArrayT<StringT> element_labels;
	exo.ReadElementLabels(element_labels);
#endif
  
	/* read nodal data */
	cout << "num node vars: "<< num_node_variables << endl;
	cout << "num nodes: " << num_nodes << endl;
	cout << "num time steps: " << num_time_steps << endl;
	
	/* initialize variable ranges */
	scalarRange1.Allocate(num_node_variables);
	scalarRange2.Allocate(num_node_variables);
	scalarRange1 = DBL_MAX;
	scalarRange2 = DBL_MIN;
	
	/* set default variable to be displayed */ 
	if (num_node_variables > 0)  
		currentVarNum = num_node_variables-1;
	else
		currentVarNum = -1;
  
  	/* load step zero */
	if (num_time_steps > 0) SelectTimeStep(0);  

	/* color mapping variables */
	DefaultValues();
  	if (num_node_variables > 0) UpdateData();

	/* add variables to the console */
	iAddVariable("min_Hue_Range", hueRange1);
	iAddVariable("max_Hue_Range", hueRange2);
	iAddVariable("min_Value_Range", valRange1);
	iAddVariable("max_Value_Range", valRange2);
	iAddVariable("min_Saturation_Range", satRange1);
	iAddVariable("max_Saturation_Range", satRange2);
	iAddVariable("min_Alpha_Range", alphaRange1);
	iAddVariable("max_Alpha_Range", alphaRange2);
	if (currentVarNum > 0)
	{
		iAddVariable("min_Scalar_Range", scalarRange1[currentVarNum]);
		iAddVariable("max_Scalar_Range", scalarRange2[currentVarNum]);
	}
	iAddVariable("numColors", numColors);
	iAddVariable("scale_factor", scale_factor);
	iAddVariable("opacity", opacity);
  	
  	/* commands */
  	iAddCommand(CommandSpecT("Wire"));
  	iAddCommand(CommandSpecT("Surface"));
  	iAddCommand(CommandSpecT("Point"));
}

/* destructor */
VTKBodyDataT::~VTKBodyDataT(void)
{
	/* coordinate point data */
  	if (fPoints) fPoints->Delete();

	/* unstructured grids */
	for (int i = 0; i < fUGrids.Length(); i++)
		if (fUGrids[i])
			delete fUGrids[i];

	/* free memory for all stored results */
	for (int i = 0; i <fVectors.Length(); i++)
	{
		if (fVectors[i])fVectors[i]->Delete();
    	for (int j = 0; j < fScalars.MinorDim(); j++)
    		if (fScalars(i,j)) fScalars(i,j)->Delete();
	}
}

/* return the number of spatial dimensions */
int VTKBodyDataT::NumSD(void)
{
	if (fUGrids.Length() == 0)
		return 0;
	else
		return fUGrids[0]->NumSD();
}

/* add actors in self to the given renderer */
void VTKBodyDataT::AddToRenderer(vtkRenderer* renderer) const
{
	/* add all actors */
	for (int i = 0; i < fUGrids.Length(); i++)
		renderer->AddActor(fUGrids[i]->Actor());
}

/** remove actors in self to the given renderer */
void VTKBodyDataT::RemoveFromRenderer(vtkRenderer* renderer) const
{
	/* remove all actors */
	for (int i = 0; i < fUGrids.Length(); i++)
		renderer->RemoveActor(fUGrids[i]->Actor());
}


void VTKBodyDataT::UpdateData(void)
{
	if (fScalars.MinorDim() > 0)
	{
  		/* update range */
  		for (int i = 0; i < fUGrids.Length(); i++)
  			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
  			{
				fUGrids[i]->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
				fUGrids[i]->SetOpacity(opacity);
			}
	}	
}

bool VTKBodyDataT::ChangeVars(const StringT& var)
{
	/* find variable number */
	int varNum = -1;
	for (int i = 0; varNum == -1 && i < fNodeLabels.Length(); i++)
		if (fNodeLabels[i] == var)
			varNum = i;

	/* change if found */
	if (varNum == -1)
		return false;
  	else 
  	{
		currentVarNum = varNum;
  		for (int i = 0; i < fUGrids.Length(); i++)
  		{
  			/* change scalar */
  			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
  			{
  				fUGrids[i]->SetScalars(fScalars(currentStepNum, currentVarNum));
				fUGrids[i]->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
				
				/* reset references in console variables */
				iDeleteVariable("min_Scalar_Range");
				iDeleteVariable("max_Scalar_Range");
				iAddVariable("min_Scalar_Range", scalarRange1[currentVarNum]);
				iAddVariable("max_Scalar_Range", scalarRange2[currentVarNum]);
  			}
  		}
		return true;
  	}
}

bool VTKBodyDataT::SelectTimeStep(int stepNum)
{
	if (fScalars.MinorDim() > 0) {

		if (stepNum >= 0 && stepNum < fScalars.MajorDim())
		{
			/* load data into fScalars and fVectors */
			LoadData(stepNum);
  			for (int i = 0; i < fUGrids.Length(); i++)
  			{
  				/* color */
  				if (!fScalars(stepNum, currentVarNum)) throw eGeneralFail;
  				if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
  				{
  					fUGrids[i]->SetScalars(fScalars(stepNum, currentVarNum));
					fUGrids[i]->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
				}
	
	  			/* displaced shape */
				if (fVectors.Length() > 0){
					if (!fVectors[stepNum]) throw eGeneralFail;
		  			fUGrids[i]->SetWarpVectors(fVectors[stepNum]);
		  		}
	  		}
			currentStepNum = stepNum;
			return true;
		}
		else
		{
			cout << "step number out of range: " << stepNum << endl;
			return false;
		}
  	}
  	else /* no history */
  		return true;
}

#if 0
void VTKBodyDataT::ChangeDataColor(int color)
{
  ugridMapper->ScalarVisibilityOff();
  if (color ==1)
    ugridActor->GetProperty()->SetColor(1,0,0);
  else if (color==2)
    ugridActor->GetProperty()->SetColor(0,1,0);
  else if (color==3)
    ugridActor->GetProperty()->SetColor(0,0,1);
  else
    cout << "invalid color";
}
#endif

/* execute console command. \return true is executed normally */
bool VTKBodyDataT::iDoCommand(const CommandSpecT& command, StringT& line)
{
	if (command.Name() == "Wire")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kWire);
		return true;
	}
	else if (command.Name() == "Surface")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kSurface);
		return true;
	}
	else if (command.Name() == "Point")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kPoint);
		return true;
	}
	else /* inherited */
		return iConsoleObjectT::iDoCommand(command, line);
}

/*************************************************************************
* private
*************************************************************************/

void VTKBodyDataT::DefaultValues(void)
{
  numColors = 256;
  hueRange1 = 0.6667; hueRange2 = 0;
  valRange1 = 1; valRange2 = 1;
  satRange1 = 1; satRange2 = 1;
  alphaRange1 = 1; alphaRange2 = 1;
  opacity = 1;
}

/* load data for the current time step */
void VTKBodyDataT::LoadData(int step)
{
	if (step < 0 || step >= fScalars.MajorDim())
	{
		cout << "VTKBodyDataT::LoadData: step is out of range: " << step << endl;
		throw eOutOfRange;
	}
	
	/* database file */
	ExodusT exo(cout);
	if (!exo.OpenRead(fInFile)) {
		cout << "VTKBodyDataT::LoadData: unable to open file: " << fInFile << endl;
		throw eGeneralFail;
	}

	//not used
	//double time;
	//exo.ReadTime(i+1, time);
	
	/* dimensions */
	int num_nodes = exo.NumNodes();
	int num_node_variables = fScalars.MinorDim();

	/* load variable data in scalar */
	bool did_read = false;
	for (int j = 0; j < num_node_variables; j++)
	{
		/* not read yet */
		if (!fScalars(step,j))
		{
			did_read = true;
		
			/* allocate scalars */
			fScalars(step,j) = vtkFloatArray::New();
			fScalars(step,j)->SetNumberOfComponents(1);

			/* read variable */
			dArrayT ndata(num_nodes);
			exo.ReadNodalVariable(step+1, j+1, ndata);

			/* range over steps that have been loaded */
			double min, max;
			ndata.MinMax(min, max);
			scalarRange1[j] = (min < scalarRange1[j]) ? min : scalarRange1[j];
			scalarRange2[j] = (max > scalarRange2[j]) ? max : scalarRange2[j];

#if 0
			/* set one tuple at a time */
			for (int k = 0; k < num_nodes; k++)
				fScalars(step,j)->InsertTuple1(k+1, ndata[k]);
#endif
					
#if 1			
			/* translate to float */
			ArrayT<float> fdata(num_nodes + 1); /* zeroth tuple is ignored */
			double_to_float(ndata, fdata.Pointer(1));
//			double_to_float(ndata, fdata.Pointer(0)); //SHIFT
				
			/* load in */
			float* p;
			fdata.ReleasePointer(&p);
			fScalars(step, j)->SetArray(p, ndata.Length(), 0);
#endif
		}
	}
			
	/* instantiate displacement vector if needed */
	if (vec_dim > 0 && !fVectors[step])
	{
		did_read = true;
	
		/* allocate vectors */
		fVectors[step] = vtkFloatArray::New();
		fVectors[step]->SetNumberOfComponents(3);

#if 0
		/* set one tuple at a time */
		for (int j = 0; j < num_nodes; j++)
		{	
			double* p = disp_tmp(j);
			fVectors[step]->InsertTuple3(j+1, p[0], p[1], p[2]); //?????????do tuple numbers start at 1?????????
		}
#endif

#if 1
		/* temp space */
		nArray2DT<float> disp(num_nodes+1, 3);
		disp = 0.0;
		for (int j = 0; j < vec_dim; j++)
			disp.SetColumn(j, fScalars(step, j)->GetPointer(0));
					
		/* load into vectors */
		float* p;
		disp.ReleasePointer(&p);
		fVectors[step]->SetArray(p, disp.Length()-3, 0);
#endif
	}		

	/* message */
	if (did_read) cout << "read data for step: " << step << endl;
}
