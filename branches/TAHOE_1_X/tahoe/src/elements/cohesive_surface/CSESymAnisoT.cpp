/* $Id: CSESymAnisoT.cpp,v 1.6 2003-12-01 23:53:15 cjkimme Exp $ */
/* created: paklein (11/19/1997) */
#include "CSESymAnisoT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "toolboxConstants.h"
#include "SurfaceShapeT.h"
#include "SurfacePotentialT.h"
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "eIntegratorT.h"
#include "NodeManagerT.h"
#include "OutputSetT.h"
#endif
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "dSymMatrixT.h"
#include "LocalArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

/* potential functions */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "TiedPotentialBaseT.h"
#endif

using namespace Tahoe;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* constructor */
CSESymAnisoT::CSESymAnisoT(const ElementSupportT& support, const FieldT& field, bool rotate):
	CSEAnisoT(support, field, rotate)
{
	SetName("modeI_CSE");

	/* reset format for the element stiffness matrix */
	if (fRotate) fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

/* constructor */
CSESymAnisoT::CSESymAnisoT(const ElementSupportT& support):
	CSEAnisoT(support)
{
	SetName("modeI_CSE");
}
#else
CSESymAnisoT::CSESymAnisoT(ElementSupportT& support, bool rotate):
	CSEAnisoT(support)
{
	SetName("modeI_CSE");
}
#endif

CSESymAnisoT::~CSESymAnisoT(void)
{

}

/* writing output */
void CSESymAnisoT::RegisterOutput(void)
{
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"

	/* "deformed" geometry */
	GeometryT::CodeT geo_code;
	switch (fGeometryCode)
	{
		case GeometryT::kLine:		
			geo_code = GeometryT::kLine;
			break;

		case GeometryT::kQuadrilateral:
			geo_code = GeometryT::kQuadrilateral;
			break;

		case GeometryT::kTriangle:
			geo_code = GeometryT::kTriangle;
			break;

		default:	
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		cout << "\n CSEBaseT::RegisterOutput: could not translate\n";
		cout << "     geometry code " << fGeometryCode
			 << " to a pseudo-geometry code for the volume." << endl;
#endif
		throw ExceptionT::kGeneralFail;	
	}	

	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	
	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();
	
	/* set output specifier */
	OutputSetT output_set(geo_code, block_ID, sideSet_ID, fOutput_Connectivities, n_labels, e_labels, false);

	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
#else
	ElementSupport().RegisterOutput(n_labels,e_labels);
#endif
}


/***********************************************************************
 * Protected
 ***********************************************************************/

void CSESymAnisoT::LHSDriver(GlobalT::SystemTypeT)
{
	const char caller[] = "CSESymAnisoT::LHSDriver";

		/* matrix format */
	dMatrixT::SymmetryFlagT format = (fRotate) ?
		dMatrixT::kWhole : dMatrixT::kUpperOnly;

	/* time-integration parameters */
	double constK = 1.0;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;
#endif

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

	/* work space for collecting element variables */
	LocalArrayT nodal_values(LocalArrayT::kUnspecified);
	dArray2DT elementVals;
	dArrayT localFrameIP;
	iArrayT ndIndices;

	/* Next three lines just for transposed displacements */
	dArray2DT disps_global;
	LocalArrayT disps(LocalArrayT::kDisp);
	iArrayT disps_inds;
	
	int  nsd = NumSD();
	
	AutoArrayT<double> state2;
	dArrayT state;
	Top();
	while (NextElement())
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
	
		/* surface potential */
		SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
		int num_state = fNumStateVariables[element.MaterialNumber()];
		state2.Dimension(num_state);
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		TiedPotentialBaseT* tiedpot = fTiedPots[element.MaterialNumber()];
#endif

		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(element.NodesX());
		if (!disps_global.IsAllocated())
		{
			disps_global.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
			disps.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
			disps.SetGlobal(disps_global);
			disps_inds.Dimension(element.NodesX().Length());
			disps_inds.SetValueToPosition();
		}
		disps_global.RowCollect(element.NodesX().Pointer(), ElementSupport().CurrentCoordinates());
	  			

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
		
		disps.SetLocal(disps_inds);
		
		disps -= fLocInitCoords1;
		
		/* initialize */
		fLHS = 0.0;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
		/* Get local bulk values for CSE*/
		if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
			SurfaceValuesFromBulk(element, ndIndices, elementVals, nodal_values);
#endif

		/* loop over integration points */
		double* pstate = fStateVariables(CurrElementNumber());
		fShapes->TopIP();
		while (fShapes->NextIP())
		{  
			/* set state variables */
			state.Set(num_state, pstate);
			pstate += num_state;
		
			/* integration weights */
			double w = fShapes->IPWeight();		

			/* coordinate transformations */
			double j0, j;
			if (fRotate)
			{
				j0 = fShapes->Jacobian();
				j  = fCurrShapes->Jacobian(fQ, fdQ);
			}
			else
				j0 = j = fShapes->Jacobian(fQ);
			fIPArea = w*j0;

			/* check */
			if (j0 <= 0.0 || j <= 0.0) ExceptionT::BadJacobianDet(caller);
		
			/* gap vector and gradient (facet1 to facet2) */
			const dArrayT&    delta = fShapes->InterpolateJumpU(disps);
			const dMatrixT& d_delta = fShapes->Grad_d();

			/* gap vector in local frame */
			fQ.MultTx(delta, fdelta);	
			
			/* enforce symmetry */
			if (nsd == 2)
			{
				fdelta[0] = 0.;
				fdelta[1] *= 2.;
			}
			else
			{
				fdelta[0] = 0.;
				fdelta[1] = 0.;
				fdelta[2] *= 2.;
			}
					
			
#ifndef _FRACTURE_INTERFACE_LIBRARY_
			/* Interpolate nodal info to IPs using stress tensor in local frame*/
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
				FromNodesToIPs(tiedpot->RotateNodalQuantity(), localFrameIP, nodal_values);
#endif

			/* stiffness in local frame */
			fK = surfpot->Stiffness(fdelta, state, localFrameIP);
			
			/* enforce symmetry */
			if (nsd == 2)
			{
				fK(0,0) = fK(1,0) = fK(0,1) = 0.;
				fK(1,1) *= 2.;
			}
			else
			{
				fK(0,0) = fK(1,0) = fK(0,1) = fK(1,1) = fK(2,0) = fK(0,2) = fK(1,2) = fK(2,1) = 0.;
				fK(2,2) *= 2.;
			} 
				
			/* rotation */
			if (fRotate)
			{
				/* traction in local frame */
				const dArrayT& T = surfpot->Traction(fdelta, state, localFrameIP, false);

				/* 1st term */
				fT.SetToScaled(j0*w*constK, T);
				Q_ijk__u_j(fdQ, fT, fnsd_nee_1);
				fNEEmat.MultATB(d_delta, fnsd_nee_1);
				fLHS += fNEEmat;

				/* 2nd term */
				fddU.SetToScaled(j0*w*constK, fK);
				fnsd_nee_1.MultATB(fQ, d_delta);
				fnsd_nee_2.MultATB(fddU, fnsd_nee_1);
				u_i__Q_ijk(delta, fdQ, fnsd_nee_1);
				fNEEmat.MultATB(fnsd_nee_2, fnsd_nee_1);
				fLHS += fNEEmat;
			}
			
			/* 3rd term */
			fddU.MultQBQT(fQ, fK);
			fddU *= j0*w*constK;
			fLHS.MultQTBQ(d_delta, fddU, format, dMatrixT::kAccumulate);
		}

		/* assemble */
		AssembleLHS();
	}
}

void CSESymAnisoT::RHSDriver(void)
{
	const char caller[] = "CSESymAnisoT::RHSDriver";

	/* time-integration parameters */
	double constKd = 1.0;

#ifndef _FRACTURE_INTERFACE_LIBRARY_

	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* heat source if needed */
	const FieldT* temperature = ElementSupport().Field("temperature");

	/* initialize sources */
	if (temperature) 
		InitializeTemperature(temperature);

#else // _FRACTURE_INTERFACE_LIBRARY_ defined

    /*Read in SIERRA's new state variables. We need their memory. */	
	fStateVariables.Set(fStateVariables.MajorDim(),fStateVariables.MinorDim(),
		ElementSupport().StateVariableArray());

#endif

	/* set state to start of current step */
	fStateVariables = fStateVariables_last;
	if (freeNodeQ.IsAllocated())
		freeNodeQ = freeNodeQ_last;

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);
	
	/* If potential needs info from nodes, start to gather it now */
	if (fCalcNodalInfo)
		StoreBulkOutput();

	/* fracture surface area */
	fFractureArea = 0.0;

	/* work space for collecting element variables */
	LocalArrayT nodal_values(LocalArrayT::kUnspecified);
	dArray2DT elementVals;
	dArrayT localFrameIP;
	iArrayT ndIndices;

	int block_count = 0, block_dex = 0;
	
	/* Next three lines just for tranposed displacements */
	dArray2DT disps_global;
	LocalArrayT disps(LocalArrayT::kDisp);
	iArrayT disps_inds;
	
	int  nsd = NumSD();
	
	dArrayT state;
	Top();
	while (NextElement())
	{
		/* advance to block (skip empty blocks) */
		while (block_count == fBlockData[block_dex].Dimension()) {
			block_count = 0;
			block_dex++;		
		}

		/* current element */
		ElementCardT& element = CurrentElement();
	
		/* get ref geometry (1st facet only) */
		fLocInitCoords1.SetLocal(element.NodesX());
		if (!disps_global.IsAllocated())
		{
			disps_global.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
			disps.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
			disps.SetGlobal(disps_global);
			disps_inds.Dimension(element.NodesX().Length());
			disps_inds.SetValueToPosition();
		}
		disps_global.RowCollect(element.NodesX().Pointer(), ElementSupport().CurrentCoordinates());
	  			
		if (element.Flag() != kOFF)
		{
			/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
			int num_state = fNumStateVariables[element.MaterialNumber()];
	
			/* get current geometry */
			SetLocalX(fLocCurrCoords);
			disps.SetLocal(disps_inds);
			
			disps -= fLocInitCoords1;
			
	  		/* initialize */
	  		fRHS = 0.0;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
			TiedPotentialBaseT* tiedpot = fTiedPots[element.MaterialNumber()];

			/* Get local bulk values for CSE*/
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 	
				SurfaceValuesFromBulk(element, ndIndices, elementVals, nodal_values);
#endif
			
			/* loop over integration points */
			double* pstate = fStateVariables(CurrElementNumber());
			int all_failed = 1;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* set state variables */
				state.Set(num_state, pstate);
				pstate += num_state;
			
				/* integration weights */
				double w = fShapes->IPWeight();		

				/* coordinate transformations */
				double j0, j;
				if (fRotate)
				{
					j0 = fShapes->Jacobian();
					j  = fCurrShapes->Jacobian(fQ);
				}
				else
					j0 = j = fShapes->Jacobian(fQ);
				fIPArea = w*j0;

				/* check */
				if (j0 <= 0.0 || j <= 0.0) ExceptionT::BadJacobianDet(caller);
	
				/* gap vector from facet1 to facet2 */
				const dArrayT& delta = fShapes->InterpolateJumpU(disps);
	
				/* gap vector in local frame */
				fQ.MultTx(delta, fdelta);
				
				/* enforce symmetry */
				if (nsd == 2)
				{
					fdelta[0] = 0.;
					fdelta[1] *= 2.;
				}
				else
				{
					fdelta[0] = 0.;
					fdelta[1] = 0.;
					fdelta[2] *= 2.;
				}

#ifndef _FRACTURE_INTERFACE_LIBRARY_					
				/* Interpolate nodal info to IPs using stress tensor in local frame*/
				if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
				{
					FromNodesToIPs(tiedpot->RotateNodalQuantity(), localFrameIP, nodal_values);
					UntieOrRetieNodes(CurrElementNumber(), ndIndices.Length(), 
									tiedpot, state, localFrameIP);
				}
#endif
				/* traction vector in/out of local frame */
				fQ.Multx(surfpot->Traction(fdelta, state, localFrameIP, true), fT);
				
				/* expand */
				fShapes->Grad_d().MultTx(fT, fNEEvec);
	
				/* accumulate */
				fRHS.AddScaled(-j0*w*constKd, fNEEvec);
				
				/* check status */
				SurfacePotentialT::StatusT status = surfpot->Status(fdelta, state);
				if (status != SurfacePotentialT::Failed) all_failed = 0;
				
				/* fracture area */
				if (fOutputArea && status != SurfacePotentialT::Precritical)
					fFractureArea += j0*w;
					
#ifndef _FRACTURE_INTERFACE_LIBRARY_
				/* incremental heat */
				if (temperature) 
					fIncrementalHeat[block_dex](block_count, fShapes->CurrIP()) = 
						surfpot->IncrementalHeat(fdelta, state);
#endif
			}

			/* assemble */
			AssembleRHS();
			
			/* mark elements */
			if (all_failed)
			{
				int& flag = element.Flag();
				if (flag == kON) flag = kMarked;
			}
		}
		else if (fOutputArea)
		{
			/* integrate fracture area */
			fShapes->TopIP();
			while (fShapes->NextIP())
				fFractureArea += (fShapes->Jacobian())*(fShapes->IPWeight());
		}

		/* next in block */
		block_count++;
	}
}

/* extrapolate the integration point stresses and strains and extrapolate */
void CSESymAnisoT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{      

	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	/* dimensions */
	int  nsd = NumSD();
	int ndof = NumDOF();
	int  nen = NumElementNodes();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* allocate element results space */
	e_values.Dimension(NumElements(), e_out);
	e_values = 0.0;

	/* work arrays */
	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT coords, disp;
	dArray2DT jump, T;
	dArray2DT matdat;	

	/* ip values */
	LocalArrayT loc_init_coords(LocalArrayT::kInitCoords, nen, nsd);
	LocalArrayT loc_disp(LocalArrayT::kDisp, nen, ndof);
	ElementSupport().RegisterCoordinates(loc_init_coords);
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	Field().RegisterLocal(loc_disp);
#else
	loc_disp.SetGlobal(ElementSupport().CurrentCoordinates());
#endif
	dArrayT ipmat(n_codes[MaterialData]);
	
	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[NodalCoord], pall) ; pall += coords.Length();
	disp.Set(nen, n_codes[NodalDisp], pall)    ; pall += disp.Length();
	jump.Set(nen, n_codes[NodalDispJump], pall); pall += jump.Length();
	T.Set(nen, n_codes[NodalTraction], pall)   ; pall += T.Length();
	matdat.Set(nen, n_codes[MaterialData], pall);

	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT centroid;
	if (e_codes[Centroid])
	{
		centroid.Set(nsd, pall); 
		pall += nsd;
	}
	double phi_tmp, area;
	double& phi = (e_codes[CohesiveEnergy]) ? *pall++ : phi_tmp;
	dArrayT traction;
	if (e_codes[Traction])
	{
		traction.Set(ndof, pall); 
		pall += ndof;
	}

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

	dArray2DT elementVals;
	LocalArrayT nodal_values(LocalArrayT::kUnspecified);
	dArrayT localFrameIP;
	iArrayT ndIndices;

	/* Next three lines just for tranposed displacements */
	dArray2DT disps_global;
	LocalArrayT disps(LocalArrayT::kDisp);
	iArrayT disps_inds;
	
	AutoArrayT<double> state;
	Top();
	while (NextElement())
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
		
		/* initialize */
		nodal_space = 0.0;
		element_values = 0.0;

		/* coordinates for whole element */
		if (n_codes[NodalCoord])
		{
			SetLocalX(loc_init_coords);
			loc_init_coords.ReturnTranspose(coords);
		}
		
		/* displacements for whole element */
		if (n_codes[NodalDisp])
		{
			SetLocalU(loc_disp);
			loc_disp.ReturnTranspose(disp);
		}
		
		//NOTE: will not get any element output if the element not kON
		//      although quantities like the reference element centroid could
		//      still be safely calculated.
		/* compute output */
		if (element.Flag() == kON)
		{
	  		/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
			int num_state = fNumStateVariables[element.MaterialNumber()];
			state.Dimension(num_state);

			/* get ref geometry (1st facet only) */
			fNodes1.Collect(facet1, element.NodesX());
			fLocInitCoords1.SetLocal(element.NodesX());
			
			if (!disps_global.IsAllocated())
			{
				disps_global.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
				disps.Dimension(fLocInitCoords1.NumberOfNodes(),fLocInitCoords1.MinorDim());
				disps.SetGlobal(disps_global);
				disps_inds.Dimension(element.NodesX().Length());
				disps_inds.SetValueToPosition();
			}
			disps_global.RowCollect(element.NodesX().Pointer(), ElementSupport().CurrentCoordinates());
	  			
			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
			
			disps.SetLocal(disps_inds);
			
			disps -= fLocInitCoords1;

			/* initialize element values */
			phi = area = 0.0;
			if (e_codes[Centroid]) centroid = 0.0;
			if (e_codes[Traction]) traction = 0.0;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
			TiedPotentialBaseT* tiedpot = fTiedPots[element.MaterialNumber()];
			
			/* Get local bulk values for CSE*/
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 	
				SurfaceValuesFromBulk(element, ndIndices, elementVals, nodal_values);
#endif

			/* integrate */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				double* pstate = fStateVariables(CurrElementNumber()) + 
					fShapes->CurrIP()*num_state;
			
				/* element integration weight */
				double ip_w = fShapes->Jacobian()*fShapes->IPWeight();
				area += ip_w;

				/* gap */
				const dArrayT& gap = fShapes->InterpolateJumpU(disps);

				/* coordinate transformation */
				double j = fCurrShapes->Jacobian(fQ);
				fQ.MultTx(gap, fdelta);
				
				/* enforce symmetry */
				if (nsd == 2)
				{
					fdelta[0] = 0.;
					fdelta[1] *= 2.;
				}
				else
				{
					fdelta[0] = 0.;
					fdelta[1] = 0.;
					fdelta[2] *= 2.;
				}
				
				/* gap */				
				if (n_codes[NodalDispJump])
					fShapes->Extrapolate(fdelta, jump);
	     
				/* traction */
				if (n_codes[NodalTraction] || e_codes[Traction])
				{
					/* copy state variables (not integrated) */
					state.Set(num_state,pstate);

#ifndef _FRACTURE_INTERFACE_LIBRARY_
					/* Interpolate nodal info to IPs using stress tensor in local frame*/
					if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
						FromNodesToIPs(tiedpot->RotateNodalQuantity(), localFrameIP, nodal_values);			
#endif

					/* compute traction in local frame */
					const dArrayT& tract = surfpot->Traction(fdelta, state, localFrameIP, false);
				       
					/* project to nodes */
					if (n_codes[NodalTraction])
						fShapes->Extrapolate(tract, T);
					
					/* element average */
					if (e_codes[Traction])
						traction.AddScaled(ip_w, tract);
				}
					
				/* material output data */
				if (n_codes[MaterialData])
				{
					/* evaluate */
					surfpot->ComputeOutput(fdelta, state, ipmat);
					fShapes->Extrapolate(ipmat, matdat);
				}

				/* moment */
				if (e_codes[Centroid])
					centroid.AddScaled(ip_w, fShapes->IPCoords());
				
				/* cohesive energy */
				if (e_codes[CohesiveEnergy])
				{
					/* copy state variables (not integrated) */
					state.Copy(pstate);

					/* surface potential */
					double potential = surfpot->Potential(fdelta, state);

					/* integrate */
					phi += potential*ip_w;
				}
			}
			
			/* element values */
			if (e_codes[Centroid]) centroid /= area;
			if (e_codes[Traction]) traction /= area;
		}
		/* element has failed */
		else
		{
			/* can still be calculated */
			if (e_codes[Centroid] || e_codes[CohesiveEnergy])
			{
				/* get ref geometry (1st facet only) */
				fNodes1.Collect(facet1, element.NodesX());
				fLocInitCoords1.SetLocal(fNodes1);

				/* initialize element values */
				phi = area = 0.0;
				if (e_codes[Centroid]) centroid = 0.0;
		
		  		/* surface potential */
				SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
				int num_state = fNumStateVariables[element.MaterialNumber()];
				state.Dimension(num_state);
		
				/* integrate */
				fShapes->TopIP();
				while (fShapes->NextIP())
				{

					double* pstate = fStateVariables_last(CurrElementNumber()) + fShapes->CurrIP()*num_state;
					/* element integration weight */
					double ip_w = fShapes->Jacobian()*fShapes->IPWeight();
					area += ip_w;
		
					/* moment */
					if (e_codes[Centroid])
						centroid.AddScaled(ip_w, fShapes->IPCoords());

					/* cohesive energy */
					if (e_codes[CohesiveEnergy])
					  {  
						/* surface potential */
						state.Copy(pstate);
						double potential = surfpot->FractureEnergy(state);
	
						/* integrate */
						phi += potential*ip_w;
					}
				}
				
				/* element values */
				if (e_codes[Centroid]) centroid /= area;
			}		
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(jump  , colcount); colcount += jump.MinorDim();
		nodal_all.BlockColumnCopyAt(T     , colcount); colcount += T.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat, colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		ElementSupport().AssembleAverage(element.NodesX(), nodal_all);
		
		/* store results */
		e_values.SetRow(CurrElementNumber(), element_values);	      
	}

	/* get nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
}

/* read element connectivity data */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
void CSESymAnisoT::ReadConnectivity(ifstreamT& in, ostream& out)
{
#else
void CSESymAnisoT::ReadConnectivity(void)
{
#endif
	
	/* read from parameter file */
	iArrayT matnums;
	ModelManagerT& model = ElementSupport().Model();
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	bool multiDatabaseSets = false;
	model.SideSetList(in, sideSet_ID, multiDatabaseSets);
#else
	/* For Sierra, can't use input stream */
	sideSet_ID.Dimension(1);
	matnums.Dimension(1);
	sideSet_ID[0] = fSupport.BlockID();
	/*Might have to generalize this later*/
	matnums = 1; 
#endif

	/* allocate block map */
	int num_blocks = sideSet_ID.Length();
	
	if (num_blocks != 1)
	{
		ExceptionT::GeneralFail("CSESymAnisoT::ReadConnectivitiy","Multiple element blocks not implemented\n");
	}
	
	fBlockData.Dimension(num_blocks);
	fConnectivities.Allocate(num_blocks);

	/* Workspace */
	ArrayT<GeometryT::CodeT> ssArray;
	iArrayT facetNodes;
	iArray2DT faces;

	/* read from database */
	int elem_count = 0;
	int nen = 0;
	for (int b=0; b < num_blocks; b++)
	{
	    /* check number of nodes */
	    int num_elems, num_nodes;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	    model.SideSet(sideSet_ID[b], ssArray, facetNodes, faces);
	    num_elems = faces.MajorDim();
	    if (num_elems == 0) // empty side sets are possible and should be allowed
	    {
	    	ssArray.Dimension(1);
	    	ssArray[0] = NumSD() == 2 ? GeometryT::kLine : GeometryT::kQuadrilateral;
	    	faces.Dimension(0, DefaultNumElemNodes());
	    }
	    num_nodes = faces.MinorDim();
#else
		num_elems = fSupport.NumElements();
		num_nodes = model.NumNodes(); 
#endif

	    /* set if unset */
	    if (nen == 0) nen = num_nodes;
	    
	    /* consistency check */
	    if (num_nodes != 0 && nen != num_nodes)
		{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
			cout << "\n CSymAnisoT::ReadConnectivity: minor dimension "
                 << num_nodes << " of block " << b+1 << '\n';
			cout <<   "     does not match dimension of previous blocks "
                 << nen << endl;
#endif                 
			throw ExceptionT::kBadInputValue;
		}
		
		/* Make a new ID that's the last element group in the database */
		StringT new_id;
		new_id.Append(model.NumElementGroups()+1);
		new_id = model.FreeElementID(new_id);
	    
	    /* store block data  ASSUMING bth block is material number b*/
	    fBlockData[b].Set(new_id, elem_count, num_elems, b); 

	    /* increment element count */
	    elem_count += num_elems;
		
		if (!model.RegisterElementGroup(new_id, faces, ssArray[b], false))
		{
			ExceptionT::GeneralFail("CSESymAnisoT::ReadConnectivities","Cannot register element group\n");
		}
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		out << " Creating block ID "<< new_id << " for this element group\n";
#endif
	    /* set pointer to connectivity list */
	    fConnectivities[b] = model.ElementGroupPointer(new_id);
	}

	/* connectivities came back empty */
	if (nen == 0)
		ExceptionT::GeneralFail("CSESymAnisoT::ReadConnectivity","No connectivity from side sets\n");
	for (int i = 0; i < fConnectivities.Length(); i++)
		if (fConnectivities[i]->MinorDim() == 0)
		{
			/* not really violating const-ness */
			iArray2DT* connects = const_cast<iArray2DT*>(fConnectivities[i]);
			connects->Dimension(0, nen);
		}
	  
	/* set dimensions */
	fElementCards.Dimension(elem_count);
	
	/* set for output */
	fOutput_Connectivities = fConnectivities;
}

/* write all current element information to the stream */
void CSESymAnisoT::CurrElementInfo(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* inherited */
	CSEBaseT::CurrElementInfo(out);
	
	/* current element configuration */
	out << " current integration point:" << fCurrShapes->CurrIP() << '\n';
	try
	{
		dMatrixT Qtemp(fQ);
		double j = fCurrShapes->Jacobian(Qtemp);
		out << " surface jacobian = " << j << '\n';
		out << " coordinate transformation:\n";
		out << fQ << '\n';
		out << " gap vector (global frame) (not correct):\n";
		const dArrayT& delta = fShapes->InterpolateJumpU(fLocCurrCoords);
		out << delta << '\n';
		out << " gap vector (local frame):\n";
		dArrayT delta_temp(fdelta);
		fQ.MultTx(delta, delta_temp);
		
		/* enforce symmetry */
		if (delta.Length() == 2)
		{
			delta_temp[0] = 0.;
			delta_temp[1] *= 2.;
		}
		else
		{
			delta_temp[0] = 0.;
			delta_temp[1] = 0.;
			delta_temp[2] *= 2.;
		}
		
		out << delta_temp << '\n';	
	}
	
	catch (ExceptionT::CodeT error)
	{
		out << " CSESymAnisoT::CurrElementInfo: error on surface jacobian\n";
	}
#else
#pragma unused(out)
	throw ExceptionT::kGeneralFail;
#endif
}

/***********************************************************************
* Private
***********************************************************************/

/* operations with pseudo rank 3 (list in j) matrices */
void CSESymAnisoT::u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
	dMatrixT& Qu)
{
	for (int i = 0; i < u.Length(); i++)
	{	
		Q[i].MultTx(u, fNEEvec);
		Qu.SetRow(i, fNEEvec);
	}
}

void CSESymAnisoT::Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
	dMatrixT& Qu)
{
	if (Q.Length() == 2)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1]);
	else if (Q.Length() == 3)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1], u[2], Q[2]);
	else
		throw ExceptionT::kGeneralFail;
}

