/* $Id: VTKBodyDataT.cpp,v 1.13 2002-02-28 16:27:58 sawimme Exp $ */
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
#include "ModelManagerT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "GeometryT.h"
#include "StringT.h"
#include "CommandSpecT.h"

/* array behavior */
const bool ArrayT<VTKBodyDataT*>::fByteCopy = true;

/* constructor */
VTKBodyDataT::VTKBodyDataT(IOBaseT::FileTypeT format, const StringT& file_name): 
	fFormat(format),
	fInFile(file_name),
	fPoints(NULL),
	currentStepNum(0)
{
	/* data reader */
	ModelManagerT model(cout);

	/* read exodus file */
	try { model.Initialize(fFormat, fInFile, true); }
	catch (int error) {
		cout << " EXCEPTION: caught exception " << error << " reading file: " << fInFile << endl;
		throw eDatabaseFail;
	}
	cout << "initialized database file: " << fInFile << endl;
  
	/* read coordinates */
	dArray2DT coords;
	coords.Alias(model.Coordinates());
	if (coords.MinorDim() == 2) /* fill to 3D */
    {
		/* temp space */ 
		dArray2DT tmp(coords.MajorDim(), 3); 
      
		/* write in */ 
		tmp.BlockColumnCopyAt(coords, 0);    
		tmp.SetColumn(2, 0.0); 
      
		/* swap memory */
		coords.Free();
		tmp.Swap(coords); 
    }
    
    /* read the node numbering map */
    fPointNumberMap.Allocate(coords.MajorDim());
	model.AllNodeMap(fPointNumberMap);
	fPointNumberMap++; // modelmanagerT offsets map to zero

	//TEMP
	//cout << "node number map:\n" << fPointNumberMap.wrap(10) << endl;

//TEMP
	if (fPointNumberMap.Length() != coords.MajorDim()) {
		cout << "VTKBodyDataT: no node number map?" << endl;
		throw eGeneralFail;
	}
    

#if 1
	/* set up points */
	int num_nodes = model.NumNodes();
	fPoints = vtkPoints::New();
  	//fPoints->SetNumberOfPoints(num_nodes + 1);
  	for (int i=0; i < num_nodes; i++) 
//		fPoints->InsertPoint(i+1, coords(i));
  	fPoints->InsertPoint(i, coords(i)); //SHIFT
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
	int num_elem_blocks = model.NumElementGroups();
	int num_node_sets = model.NumNodeSets();
	fUGrids.Allocate(num_elem_blocks + num_node_sets);	
  
	/* load element connectivities */
	const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
	for (int i = 0 ; i < num_elem_blocks; i++)
    {
		/* read connectivities */
		GeometryT::CodeT geom_code = model.ElementGroupGeometry(elem_ID[i]);
		const iArray2DT& connectivities = model.ElementGroup(elem_ID[i]);
		//connectivities[i]--; //SHIFT

#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading element block: " 
		     << connectivities.MajorDim() << " x " << connectivities.MinorDim() << endl;
#endif

		/* element numbering map */
		iArrayT map(connectivities.MajorDim());
		model.ElementMap(elem_ID[i], map);
		
		//TEMP
		//cout << "element number map:\n" << map.wrap(10) << endl;
	
		/* construct VTK grid */
		fUGrids[i] = new VTKUGridT(VTKUGridT::kElementSet, i, model.NumDimensions());
		fUGrids[i]->SetPoints(fPoints);
		fUGrids[i]->SetConnectivities(geom_code, connectivities);
		fUGrids[i]->SetCellNumberMap(map);
	}    
    cout << "read element blocks" << endl;

	/* load node sets */
	const ArrayT<StringT>& node_ID = model.NodeSetIDs();
	for (int i = 0; i < num_node_sets; i++)
    {
		/* read nodes */
		GeometryT::CodeT geom_code = GeometryT::kPoint;
		const iArrayT& nodes = model.NodeSet(node_ID[i]);
		
		iArray2DT connectivities(nodes.Length(), 1, nodes.Pointer());
	
#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading node set: " 
		     << nodes.Length() << endl;
#endif
	
		/* construct VTK grid */
		int ii = i + num_elem_blocks;
		fUGrids[ii] = new VTKUGridT(VTKUGridT::kNodeSet, i, model.NumDimensions());
		fUGrids[ii]->SetPoints(fPoints);
		fUGrids[ii]->SetConnectivities(geom_code, connectivities);
	}    
    cout << "read node sets" << endl;
  
	/* number of results sets */
	int num_time_steps = model.NumTimeSteps();

    cout << "read time steps" << endl;

	/* variables defined at the nodes */
	int num_node_variables = model.NumNodeVariables();
	model.NodeLabels(fNodeLabels);
	vec_dim = model.NumDimensions();
	if (fNodeLabels.Length() >= vec_dim)
	{
		const char *d[] = {"D_X", "D_Y", "D_Z"};
		for (int i = 0; vec_dim > 0 && i < vec_dim; i++)
			if (fNodeLabels[i] != d[i])
				vec_dim = 0;

		//TEMP - other displacement variable names
		if (vec_dim == 0)
		{
			vec_dim = model.NumDimensions();
			const char *d[] = {"DISX", "DISY", "DISZ"};
			for (int i = 0; vec_dim > 0 && i < vec_dim; i++)
				if (fNodeLabels[i] != d[i])
					vec_dim = 0;	
		}
	}
	else vec_dim = 0;
	
	/* close file */
	model.CloseModel();
	
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
  		{
  			/* set displacement scale factor */
  			fUGrids[i]->SetScaleFactor(scale_factor);
  		
  			/* rendering properties */
  			if (fUGrids[i]->Type() == VTKUGridT::kElementSet)
  			{
				fUGrids[i]->SetScalarRange(scalarRange1[currentVarNum],scalarRange2[currentVarNum]);
				fUGrids[i]->SetOpacity(opacity);
				fUGrids[i]->SetNumberOfColors(numColors);
			}
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
//  hueRange1 = 0.6667; hueRange2 = 0;
//  valRange1 = 1; valRange2 = 1;
//  satRange1 = 1; satRange2 = 1;
//  alphaRange1 = 1; alphaRange2 = 1;
	numColors = 256;
	opacity = 1;
	scale_factor = 1.0;
}

/* load data for the current time step */
void VTKBodyDataT::LoadData(int step)
{
	if (step < 0 || step >= fScalars.MajorDim())
	{
		cout << "VTKBodyDataT::LoadData: step is out of range: " << step << endl;
		throw eOutOfRange;
	}
	
	//not used
	//double time;
	//exo.ReadTime(i+1, time);
	
	/* dimensions */
	int num_node_variables = fScalars.MinorDim();

	/* load variable data in scalar */
	dArray2DT nodal_data;
	dArrayT ndata;
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
			if (nodal_data.MinorDim() == 0) /* data not read yet */
			{
				/* data reader */
				ModelManagerT model(cout);

				/* read exodus file */
				try { model.Initialize(fFormat, fInFile, true); }
				catch (int error) {
					cout << " EXCEPTION: caught exception " << error << " reading file: " << fInFile << endl;
					throw eDatabaseFail;;
				}

				int num_nodes = model.NumNodes();
				nodal_data.Allocate(num_nodes, num_node_variables);
				model.AllNodeVariables(step, nodal_data);
				
				/* column */
				ndata.Allocate(num_nodes);
			}
			
			/* get data */
			nodal_data.ColumnCopy(j, ndata);

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
			ArrayT<float> fdata(ndata.Length());
//			double_to_float(ndata, fdata.Pointer(1)); //allocation must match
			double_to_float(ndata, fdata.Pointer(0)); //SHIFT
				
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
		int num_nodes = fScalars(step,0)->GetNumberOfTuples();
		nArray2DT<float> disp(num_nodes+1, 3);
		disp = 0.0;
		for (int j = 0; j < vec_dim; j++)
			disp.SetColumn(j, fScalars(step,j)->GetPointer(0));
					
		/* load into vectors */
		float* p;
		disp.ReleasePointer(&p);
		fVectors[step]->SetArray(p, disp.Length()-3, 0);
#endif
	}		

	/* message */
	if (did_read) cout << "read data for step: " << step << endl;
}
