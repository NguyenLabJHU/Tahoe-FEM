/* $Id: VTKBodyDataT.cpp,v 1.34 2003-11-25 19:54:45 paklein Exp $ */
#include "VTKBodyDataT.h"

#include "VTKUGridT.h"

#include <iostream.h>
#include <iomanip.h>
#include <float.h>

#include "vtkIdTypeArray.h"
#include "vtkPoints.h"
#include "vtkRenderer.h"
#include "vtkFloatArray.h"

#include "ExceptionCodes.h"
#include "iArray2DT.h"
#include "ModelManagerT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "GeometryT.h"
#include "StringT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"

using namespace Tahoe;

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
  

	/* read and store coordinates */
	fCoords = model.Coordinates();
	if (fCoords.MinorDim() < 3) /* fill to 3D */
    {
		/* temp space */ 
		dArray2DT tmp(fCoords.MajorDim(), 3); 
      
		/* write in */ 
		tmp.BlockColumnCopyAt(fCoords, 0);
		for (int i = fCoords.MinorDim(); i < 3; i++)
			tmp.SetColumn(i, 0.0); 
      
		/* swap memory */
		fCoords.Free();
		tmp.Swap(fCoords); 
    }
    
    /* read the node numbering map */
    fPointNumberMap.Allocate(fCoords.MajorDim());
	model.AllNodeIDs(fPointNumberMap);

//TEMP
	if (fPointNumberMap.Length() != fCoords.MajorDim()) {
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
		fPoints->InsertPoint(i, fCoords(i)); //SHIFT
#endif

#if 0
//NOTE: not the most efficient way to do things, but the code below
//      doesn't work properly.  		

	/* convert to float */
	nArray2DT<float> coords_float(num_nodes+1, 3);
	double_to_float(fCoords, coords_float(1));

	coordinates = vtkFloatArray::New();
	coordinates->SetNumberOfComponents(3);
	float* pcoords;
	coords_float.ReleasePointer(&pcoords);
	coordinates->SetArray(pcoords, fCoords.Length(), 0);
	fPoints = vtkPoints::New();
	fPoints->SetData(coordinates);
#endif

	/* dimensions */
	int num_elem_blocks = model.NumElementGroups();
	int num_node_sets = model.NumNodeSets();
	int num_side_sets = model.NumSideSets();
	fUGrids.Dimension(num_elem_blocks + num_node_sets + num_side_sets);
	fUGridNames.Dimension(fUGrids.Length());
	fUGridVisible.Dimension(fUGrids.Length());
	fUGridVisible = true;
  
	/* load element connectivities */
	const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
	for (int i = 0 ; i < num_elem_blocks; i++)
    {
		/* read connectivities */
		GeometryT::CodeT geom_code = model.ElementGroupGeometry(elem_ID[i]);
		const iArray2DT& connectivities = model.ElementGroup(elem_ID[i]);
		//connectivities[i]--; // shifted by the model manager

#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading element block: " 
		     << connectivities.MajorDim() << " x " << connectivities.MinorDim() << endl;
#endif

		/* element numbering map */
		iArrayT map(connectivities.MajorDim());
		model.ElementIDs(elem_ID[i], map);
		
		//TEMP
		//cout << "element number map:\n" << map.wrap(10) << endl;
	
		/* construct VTK grid */
		fUGrids[i] = new VTKUGridT(VTKUGridT::kElementSet, i, model.NumDimensions());
		fUGrids[i]->SetPoints(fPoints);
		fUGrids[i]->SetConnectivities(geom_code, connectivities);
		fUGrids[i]->SetCellNumberMap(map);
		
		/* name */
		fUGridNames[i].Append("elem_set_", elem_ID[i]);
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

		/* name */
		fUGridNames[ii].Append("node_set_", node_ID[i]);
	}    
    cout << "read node sets" << endl;

	/* load side sets */
	const ArrayT<StringT>& side_ID = model.SideSetIDs();
	for (int i = 0; i < num_side_sets; i++)
    {
		/* read side set */
		ArrayT<GeometryT::CodeT> facet_geom;
		iArrayT facet_nodes;
		iArray2DT faces;
		model.SideSet(side_ID[i], facet_geom, facet_nodes, faces);
	
#if __option(extended_errorcheck)
		cout << "VTKBodyDataT::VTKBodyDataT: reading side set: " 
		     << faces.MajorDim() << endl;
#endif
		
		GeometryT::CodeT geom = GeometryT::kPoint;
		if (faces.MajorDim() == 0) {
			cout << "VTKBodyDataT::VTKBodyDataT: side set " << side_ID[i] << " is empty" << endl; 
		} else geom = facet_geom[0];
	
		/* construct VTK grid */
		int ii = i + num_elem_blocks + num_node_sets;
		fUGrids[ii] = new VTKUGridT(VTKUGridT::kSideSet, i, model.NumDimensions());
		fUGrids[ii]->SetPoints(fPoints);
		fUGrids[ii]->SetConnectivities(geom, faces);

		/* name */
		fUGridNames[ii].Append("side_set_", side_ID[i]);
	}    
    cout << "read side sets" << endl;
  
	/* number of results sets */
	int num_time_steps = model.NumTimeSteps();
	fTimeList.Dimension(num_time_steps);
	model.TimeSteps(fTimeList);
	
    cout << "read time steps" << endl;

	/* variables defined at the nodes */
	int num_node_variables = model.NumNodeVariables();
	model.NodeLabels(fNodeLabels);
	
	/* set the vector fields and dimensions */
	SetVectorFields(fNodeLabels, fVectorFieldLabel, fVectorFieldDim, fVectorFieldDex);
	cout << "vector fields: " << fVectorFieldLabel.Length() << '\n';
	for (int i = 0; i < fVectorFieldLabel.Length(); i++)
	{
		/* labels */
		if (i == 0)
			cout << setw(10) << "field" << setw(kIntWidth) << "dim" << '\n';
		cout << setw(10) << fVectorFieldLabel[i] 
		     << setw(kIntWidth) << fVectorFieldDim[i] << '\n';
	}
	cout << endl;
	
	/* close file */
	model.CloseModel();
	
	/* results history */
	fScalars.Allocate(num_time_steps, num_node_variables);
	fScalars = NULL;
	if (fVectorFieldDim.Length() > 0)
	{
		fVectorFields.Allocate(num_time_steps, fVectorFieldDim.Length());
		fVectorFields = NULL;
	}

//NOT USED YET
#if 0
	/* variables defined over the elements */
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
	if (currentVarNum > -1)
	{
		iAddVariable("min_Scalar_Range", scalarRange1[currentVarNum]);
		iAddVariable("max_Scalar_Range", scalarRange2[currentVarNum]);
	}
	iAddVariable("numColors", numColors);
	iAddVariable("scale_factor", scale_factor);
	iAddVariable("opacity", opacity);
	iAddVariable("numContours", numContours);
	iAddVariable("boundingOpacity", boundingOpacity);
  
	

  	/* commands */
  	iAddCommand(CommandSpecT("Wire"));
  	iAddCommand(CommandSpecT("Surface"));
  	iAddCommand(CommandSpecT("Point"));
// 	iAddCommand(CommandSpecT("ShowContours"));
// 	iAddCommand(CommandSpecT("HideContours"));
// 	iAddCommand(CommandSpecT("ShowGlyphs"));
	iAddCommand(CommandSpecT("HideGlyphs"));


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

	/* free vector fields */
	for (int i = 0; i < fVectorFields.Length(); i++)
		if (fVectorFields[i]) fVectorFields[i]->Delete();

	/* free scalars */
	for (int i = 0; i < fScalars.MajorDim(); i++)
    	for (int j = 0; j < fScalars.MinorDim(); j++)
    		if (fScalars(i,j)) fScalars(i,j)->Delete();
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
  for (int i = 0; i < fUGrids.Length(); i++){
		renderer->AddActor(fUGrids[i]->Actor());
		//renderer->AddActor(fUGrids[i]->OutlineActor());
		//renderer->AddActor(fUGrids[i]->BoundBoxActor());
		renderer->AddActor(fUGrids[i]->SpikeActor());
  }
}

/** remove actors in self to the given renderer */
void VTKBodyDataT::RemoveFromRenderer(vtkRenderer* renderer) const
{
	/* remove all actors */
  for (int i = 0; i < fUGrids.Length(); i++){
		renderer->RemoveActor(fUGrids[i]->Actor());
		//renderer->RemoveActor(fUGrids[i]->OutlineActor());
		//renderer->RemoveActor(fUGrids[i]->BoundBoxActor());
		renderer->RemoveActor(fUGrids[i]->SpikeActor());
  }
		
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
			  
				if (fUGrids[i]->GetContoursBool()){
				  fUGrids[i]->ShowContours(fScalars(currentStepNum, currentVarNum), numContours,scalarRange1[currentVarNum],scalarRange2[currentVarNum], bodyDataRenderer);
				  fUGrids[i]->SetBoundingOpacity(boundingOpacity);
				}
				//fUGrids[i]->SetNumberOfColorBarLabels(numColorBarLabels);
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

	  			int disp_dex = VectorFieldNumber("D");
				if (disp_dex > -1) {
					if (!fVectorFields(stepNum,disp_dex)) throw eGeneralFail;
		  			fUGrids[i]->SetWarpVectors(fVectorFields(stepNum,disp_dex));

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


void VTKBodyDataT::ShowContours(vtkRenderer* renderer)
{
  if (fScalars.MinorDim() > 0){
    bodyDataRenderer = renderer;
    for (int i = 0; i < fUGrids.Length(); i++)
      {
	fUGrids[i]->ShowContours(fScalars(currentStepNum, currentVarNum), numContours, scalarRange1[currentVarNum], scalarRange2[currentVarNum], renderer);
	
      }
  }
  
}

void VTKBodyDataT::HideContours(vtkRenderer* renderer)
{
  for (int i = 0; i < fUGrids.Length(); i++)
    {
      fUGrids[i]->HideContours(renderer);
      
    }
  
}

/* return the index of the given node number in the node number map */
int VTKBodyDataT::NodeMapIndex(int node) const
{
  for (int i = 0; i < fPointNumberMap.Length(); i++)
    if (fPointNumberMap[i] == node)
      return i;
  return -1;
}

/* update visibility */
void VTKBodyDataT::UpdateVisibility(void)
{
	for (int i = 0; i < fUGrids.Length(); i++)
	{
		if (fUGridVisible[i])
			fUGrids[i]->Actor()->SetVisibility(1);
		else
			fUGrids[i]->Actor()->SetVisibility(0);
	}
}

/* execute console command. \return true is executed normally */
bool VTKBodyDataT::iDoCommand(const CommandSpecT& command, StringT& line)
{
	if (command.Name() == "Wire")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet ||
			    fUGrids[i]->Type() == VTKUGridT::kSideSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kWire);
		return true;
	}
	else if (command.Name() == "Surface")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet ||
			    fUGrids[i]->Type() == VTKUGridT::kSideSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kSurface);
		return true;
	}
	else if (command.Name() == "Point")
	{
		for (int i = 0; i < fUGrids.Length(); i++)
			if (fUGrids[i]->Type() == VTKUGridT::kElementSet ||
			    fUGrids[i]->Type() == VTKUGridT::kSideSet)
				fUGrids[i]->SetRepresentation(VTKUGridT::kPoint);
		return true;
	}
	else if (command.Name() == "HideGlyphs")
	  {
	     if (fScalars.MinorDim() > 0){
	       for (int i = 0; i < fUGrids.Length(); i++)
		 {
		   fUGrids[i]->HideGlyphing();

		 }
	     }
	     return true;
	  }
	

	else /* inherited */
		return iConsoleObjectT::iDoCommand(command, line);
}

/* return the index of specified vector field, -1 if not present */
int VTKBodyDataT::VectorFieldNumber(const char* name) const
{
	for (int i = 0; i < fVectorFieldLabel.Length(); i++)
		if (fVectorFieldLabel[i] == name)
			return i;
	return -1;
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
	numContours = 10;
	boundingOpacity = .27;
}

/* load data for the current time step */
void VTKBodyDataT::LoadData(int step)
{
	if (step < 0 || step >= fScalars.MajorDim())
	{
		cout << "VTKBodyDataT::LoadData: step is out of range: " << step << endl;
		throw eOutOfRange;
	}
	
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
	
	/* copy scalars data into vector fields */
	for (int i = 0; i < fVectorFieldLabel.Length(); i++) {
	
		/* not read yet */
		if (!fVectorFields(step,i)) {
			did_read = true;
			
			/* make vector field at least dimension 3 */
			int field_dim = fVectorFieldDim[i]; 
			int vec_dim = (field_dim < 3) ? 3 : field_dim;
	
			/* copy into temp space */
			int num_nodes = fScalars(step,0)->GetNumberOfTuples();
			nArray2DT<float> tmp(num_nodes+1, vec_dim);
			tmp = 0.0;
			int dex = fVectorFieldDex[i];
			for (int j = 0; j < field_dim; j++)
				tmp.SetColumn(j, fScalars(step, dex+j)->GetPointer(0));
					
			/* load into vectors */
			vtkFloatArray* new_vector = vtkFloatArray::New();
			new_vector->SetNumberOfComponents(vec_dim);
			float* p;
			tmp.ReleasePointer(&p);
			new_vector->SetArray(p, tmp.Length()-vec_dim, 0);

			/* store */
			fVectorFields(step,i) = new_vector;
		}
	}
			
#if 0
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
#endif		

	/* message */
	if (did_read) cout << "read data for step: " << step << endl;
}

/** determine the number and dimension of the vector fields */
void VTKBodyDataT::SetVectorFields(const ArrayT<StringT>& labels, ArrayT<StringT>& field, 
	iArrayT& dimension, iArrayT& index) const
{
	/* initialize */
	field.Dimension(0);
	dimension.Dimension(0);
	index.Dimension(0);

	/* scan */
	for (int i = 0; i < labels.Length(); i++)
	{
		StringT suffix;
		suffix.Suffix(labels[i], '_');
		bool is_123 = false;
		bool is_xyz = false;
		if (suffix.StringLength() > 0) {
			is_123 = (suffix[1] == '1' || suffix[1] == '2' || suffix[1] == '3');
			if (!is_123)		
				is_xyz = (suffix[1] == 'x' || suffix[1] == 'y' || suffix[1] == 'y' ||
				          suffix[1] == 'X' || suffix[1] == 'Y' || suffix[1] == 'Z');
		}
		
		if (is_123 || is_xyz)
		{
			StringT root, next_root;
			root.Root(labels[i], '_');
			next_root = root;
			int count = 0;
			while (next_root.StringLength() > 0 && 
				root == next_root &&
				suffix.StringLength() > 0)
			{
				count++;
				if (count == 1) {
					field.Resize(field.Length() + 1);
					field.Last() = root;
					dimension.Resize(field.Length());
					index.Resize(field.Length());
					index.Last() = i;
				}
				
				i++;
				if (i < labels.Length()) {
					suffix.Suffix(labels[i], '_');
					next_root.Root(labels[i], '_');
				} else {
					suffix.Clear();
					next_root.Clear();
				}
			}
			i--; /* rewind one */
			
			/* dimension field */
			dimension.Last() = count;
		}
	}
}
