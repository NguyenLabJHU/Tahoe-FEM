/* $Id: VTKBodyT.cpp,v 1.15 2001-12-10 22:57:58 paklein Exp $ */

#include "VTKBodyT.h"
#include "VTKBodyDataT.h"
#include "VTKFrameT.h"
#include "VTKUGridT.h"
#include "CommandSpecT.h"

#include "vtkCubeAxesActor2D.h"
#include "vtkRenderer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkWarpVector.h"
#include "vtkIdFilter.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkLabeledDataMapper.h"
#include "vtkActor2D.h"

/* array behavior */
const bool ArrayT<VTKBodyT>::fByteCopy = true;
const bool ArrayT<vtkCubeAxesActor2D*>::fByteCopy = true;

const bool ArrayT<vtkIdFilter*>::fByteCopy = true;
const bool ArrayT<vtkSelectVisiblePoints*>::fByteCopy = true;
const bool ArrayT<vtkLabeledDataMapper*>::fByteCopy = true;
const bool ArrayT<vtkActor2D*>::fByteCopy = true;

/* constructor */
VTKBodyT::VTKBodyT(VTKFrameT* frame, VTKBodyDataT* body_data):
	fFrame(frame),
	fBodyData(body_data)
{
	/* set name */
	iSetName(fBodyData->iName());

	/* add all variables */
	AddVariables(*fBodyData);

	/* add some frame commands */
	CommandSpecT* command;
	command = fFrame->iCommand("Update");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);

	command = fFrame->iCommand("Interactive");
	if (!command) throw eGeneralFail;
	iAddCommand(*command);
	
	/* other commands */
	iAddCommand(CommandSpecT("ShowNodeNumbers"));
	iAddCommand(CommandSpecT("HideNodeNumbers"));
	iAddCommand(CommandSpecT("ShowAxes"));
	iAddCommand(CommandSpecT("HideAxes"));
}

/* destructor */
VTKBodyT::~VTKBodyT(void)
{
	vtkRenderer* renderer = fFrame->Renderer();

	/* free axis actors */
	for (int i = 0; i < fAxes.Length(); i++)
		if (fAxes[i])
		{
			renderer->RemoveActor(fAxes[i]);
			fAxes[i]->Delete();
		}
		
	/* node number actors */
	for (int i = 0; i < fIDFilter.Length(); i++)
		if (fIDFilter[i])
		{	
			renderer->RemoveActor(fNodeLabelActor[i]);
			if (fNodeLabelActor[i]) fNodeLabelActor[i]->Delete();
			if (fNodeLabelMapper[i]) fNodeLabelMapper[i]->Delete();
			if (fVisPoints[i]) fVisPoints[i]->Delete();
			if (fIDFilter[i]) fIDFilter[i]->Delete();
		}
}

/* execute console command. \return true is executed normally */
bool VTKBodyT::iDoCommand(const CommandSpecT& command, StringT& line)
{
	/* check */
	if (!fFrame) {
		cout << "\n VTKBodyT::iDoCommand: frame pointer not set" << endl;
		throw eGeneralFail;
	}

	/* resolve command */
	if (command.Name() == "Update")
		return fFrame->iDoCommand(command, line);
	else if (command.Name() == "Interactive")
		return fFrame->iDoCommand(command, line);
	else if (command.Name() == "ShowNodeNumbers")
	{
		if (fIDFilter.Length() > 0)
		{
			cout << "hide numbers first" << endl;
			return false;
		}
		else
		{
			/* unstructured grids */
			const ArrayT<VTKUGridT*>& ugrids = fBodyData->UGrids();
		
			/* coordinate axis list */
			fIDFilter.Allocate(ugrids.Length());
			fVisPoints.Allocate(ugrids.Length());
			fNodeLabelMapper.Allocate(ugrids.Length());
			fNodeLabelActor.Allocate(ugrids.Length());

			fIDFilter = NULL;
			fVisPoints = NULL;
			fNodeLabelMapper = NULL;
			fNodeLabelActor = NULL;

			for (int i = 0; i < ugrids.Length(); i++)
				if (ugrids[i]->Type() == VTKUGridT::kElementSet)
				{
					VTKUGridT* ugrid = ugrids[i];
				
					/* generate id's */
					vtkIdFilter* idFilter = vtkIdFilter::New();
					idFilter->PointIdsOn();
					idFilter->FieldDataOff();
					if (ugrid->Warp())
						idFilter->SetInput(ugrid->Warp()->GetOutput());
					else
						idFilter->SetInput(ugrid->UGrid());

					/* label mapper */
					vtkLabeledDataMapper* nodeLabelMapper = vtkLabeledDataMapper::New();
					//nodeLabelMapper->SetInput(idFilter->GetOutput());
					//nodeLabelMapper->SetLabelModeToLabelIds();
					//nodeLabelMapper->SetLabelModeToLabelFieldData();
					nodeLabelMapper->SetLabelModeToLabelScalars(); /* idFilter output's id's as scalars */
					nodeLabelMapper->ShadowOff();

					/* visibility */
					if (ugrid->NumSD() == 3)
					{
						/* visibility filter */
						vtkSelectVisiblePoints* visPoints = vtkSelectVisiblePoints::New();
						visPoints->SetInput(idFilter->GetOutput());
						visPoints->SetRenderer(fFrame->Renderer());
						//visPoints->SelectionWindowOn(); // this slows things down considerably
						fVisPoints[i] = visPoints;
			
						/* label mapper */
						nodeLabelMapper->SetInput(visPoints->GetOutput());
					}
					/* assume ALL visible in 2D */
					else
						/* label mapper */
						nodeLabelMapper->SetInput(idFilter->GetOutput());

					/* labels */
					vtkActor2D* nodeLabelActor = vtkActor2D::New();
					nodeLabelActor->SetMapper(nodeLabelMapper);
					nodeLabelActor->VisibilityOn();		
					fFrame->Renderer()->AddActor(nodeLabelActor);
					
					/* add to lists */
					fIDFilter[i] = idFilter;
					fNodeLabelMapper[i] = nodeLabelMapper;
					fNodeLabelActor[i] = nodeLabelActor;
				}
		
			return true;
		}
	}
	else if (command.Name() == "HideNodeNumbers")
	{
		if (fIDFilter.Length() == 0)
		{
			cout << "numbers not showing" << endl;
			return false;
		}
		else
		{
			vtkRenderer* renderer = fFrame->Renderer();
			for (int i = 0; i < fIDFilter.Length(); i++)
			{
				/* clean-up */
				if (fNodeLabelActor[i])
				{
					renderer->RemoveActor(fNodeLabelActor[i]);
					fNodeLabelActor[i]->Delete();
				}
				if (fIDFilter[i]) fIDFilter[i]->Delete();
				if (fVisPoints[i]) fVisPoints[i] ->Delete();
				if (fNodeLabelMapper[i]) fNodeLabelMapper[i]->Delete();
			}
			
			fIDFilter.Allocate(0);
			fVisPoints.Allocate(0);
			fNodeLabelMapper.Allocate(0);
			fNodeLabelActor.Allocate(0);
			return true;
		}	
	}
	else if (command.Name() == "ShowAxes")
	{
		if (fAxes.Length() != 0) 
		{
			cout << "hide axes first" << endl;
			return false;
		} 
		else 
		{
			/* unstructured grids */
			const ArrayT<VTKUGridT*>& ugrids = fBodyData->UGrids();
		
			/* coordinate axis list */
			fAxes.Allocate(ugrids.Length());
			fAxes = NULL;
		
			vtkRenderer* renderer = fFrame->Renderer();
			for (int i = 0; i < ugrids.Length(); i++)
				if (ugrids[i]->Type() == VTKUGridT::kElementSet)
				{
					vtkCubeAxesActor2D* axes = vtkCubeAxesActor2D::New();
					axes->SetInput(ugrids[i]->UGrid());
					axes->SetCamera(renderer->GetActiveCamera());
					axes->SetLabelFormat("%6.4g");
					//axes->SetCornerOffset(.2);
					//axes->ShadowOn();
					//axes->SetFlyModeToOuterEdges();
					axes->SetFlyModeToClosestTriad();
					//axes->SetFontFactor(1.8);
					axes->GetProperty()->SetColor(0,1,1);
					// axes->SetBounds(0,1,0,1,0,1);
					//axes->ZAxisVisibilityOff();
					axes->VisibilityOn();
					renderer->AddActor(axes);
	
					fAxes[i] = axes;
				}
				
			return true;
		}	
	}
	else if (command.Name() == "HideAxes")
	{
		if (fAxes.Length() == 0)
		{
			cout << "no axes not showing" << endl;
			return false;
		}
		else
		{
			/* free all axes actors */
			vtkRenderer* renderer = fFrame->Renderer();
			for (int i = 0; i < fAxes.Length(); i++)
				if (fAxes[i])
				{
					fAxes[i]->VisibilityOff();
					renderer->RemoveActor(fAxes[i]);
					fAxes[i]->Delete();
				}
			fAxes.Allocate(0);
			return true;
		}
	}
	else
		/* inherited */
		return iConsoleObjectT::iDoCommand(command, line);
}

/* add actors in self to the given renderer */
void VTKBodyT::AddToFrame(void)
{
	/* the body */
	fBodyData->AddToRenderer(fFrame->Renderer());

	/* axes */
	if (fAxes.Length())
	{
		StringT tmp;	
		iDoCommand(*iCommand("ShowAxes"), tmp);
	}
}

/* add actors in self to the given renderer */
void VTKBodyT::RemoveFromFrame(void)
{
	/* the body */
	fBodyData->RemoveFromRenderer(fFrame->Renderer());

	/* axes */
	if (fAxes.Length())
	{
		StringT tmp;	
		iDoCommand(*iCommand("HideAxes"), tmp);
	}
}
