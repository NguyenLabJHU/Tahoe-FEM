/* $Id: BridgingScaleManagerT.cpp,v 1.2 2004-07-26 09:40:12 paklein Exp $ */
#include "BridgingScaleManagerT.h"

#ifdef BRIDGING_ELEMENT

#include "ifstreamT.h"
#include "SolverT.h"
#include "FEManagerT_THK.h"
#include "NodeManagerT.h"
#include "TimeManagerT.h"
#include "ParticlePairT.h"
#include "FieldT.h"
#include "dSPMatrixT.h"

using namespace Tahoe;

/* constructor */
BridgingScaleManagerT::BridgingScaleManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv):
	MultiManagerT(input_file, output, comm, argv),
	fFine_THK(NULL)
{
	SetName("tahoe_bridging_scale");
}

/* (re-)set system to initial conditions */
ExceptionT::CodeT BridgingScaleManagerT::InitialCondition(void)
{
	/* shouldn't enter here. BridgingScaleManagerT::Solve overrides all */
	ExceptionT::GeneralFail("BridgingScaleManagerT::InitialCondition", "should not be called");
	return ExceptionT::kNoError;
}

/* (re-)set the equation number for the given group */
void BridgingScaleManagerT::Solve(void)
{
	const char caller[] = "FEExecutionManagerT::RunDynamicBridging";

	/* configure ghost nodes */
	NodeManagerT& fine_node_manager = *(fFine_THK->NodeManager());
	int nsd = fine_node_manager.NumSD();
	int group = 0;
	int order1 = 0;	// For InterpolateField, 3 calls to obtain displacement/velocity/acceleration
	int order2 = 1;
	int order3 = 2;
	double dissipation = 0.0;
	dArray2DT field_at_ghosts, totalu, fubig, fu, projectedu, boundghostdisp, boundghostvel, boundghostacc;
	dArray2DT thkforce, totaldisp, elecdens, embforce, wavedisp;
	dSPMatrixT ntf;
	iArrayT activefenodes, boundaryghostatoms;
	const StringT& bridging_field = fFineField->FieldName();

	//already called in MultiManagerT::TakeParameterList
	//fFine_THK->InitGhostNodes(fCoarse->ProjectImagePoints());

	/* figure out boundary atoms for use with THK boundary conditions, 
	   ghost atoms for usage with MD force calculations */
	if (nsd == 2)
		boundaryghostatoms = fFine_THK->InterpolationNodes2D();
	else if (nsd == 3)
		boundaryghostatoms = fFine_THK->InterpolationNodes3D();
	else
		ExceptionT::GeneralFail(caller, "%d dimensions not implemented", nsd);
		
	int numgatoms = (fFine_THK->GhostNodes()).Length();	// total number of ghost atoms
	int numbatoms = boundaryghostatoms.Length() - numgatoms;	// total number of boundary atoms
	dArray2DT gadisp(numgatoms,nsd), gavel(numgatoms,nsd), gaacc(numgatoms,nsd);
	dArray2DT badisp(numbatoms,nsd), bavel(numbatoms,nsd), baacc(numbatoms,nsd);
	iArrayT allatoms(boundaryghostatoms.Length()), gatoms(numgatoms), batoms(numbatoms), boundatoms(numbatoms);
	allatoms.SetValueToPosition();
	batoms.CopyPart(0, allatoms, numgatoms, numbatoms);
	gatoms.CopyPart(0, allatoms, 0, numgatoms);      
	boundatoms.CopyPart(0, boundaryghostatoms, numgatoms, numbatoms);
	//elecdens.Dimension(gatoms.Length(), 1);
	//embforce.Dimension(gatoms.Length(), 1);

	/* reset interpolation data set in MultiManagerT::TakeParameterList */
	fCoarse->InitInterpolation(bridging_field, boundaryghostatoms, fine_node_manager);

	//dArrayT mdmass;
	//fFine_THK->LumpedMass(fFine_THK->NonGhostNodes(), mdmass);	// acquire array of MD masses to pass into InitProjection, etc...

	/* reset interpolation data set in MultiManagerT::TakeParameterList with makeinactive = false */
	bool makeinactive = false;	
	fCoarse->InitProjection(bridging_field, *(fFine_THK->CommManager()), fFine_THK->NonGhostNodes(), 
		fine_node_manager, makeinactive);		

	/* nodes to include/exclude in calculation of the atomistic displacements. Dimension
	 * of these matricies is the number atom types */
	nMatrixT<int> ghostonmap(2), ghostoffmap(2);  // 3D wave propagation
	//nMatrixT<int> ghostonmap(5), ghostoffmap(5);  // for 2D fracture problem
	//nMatrixT<int> ghostonmap(4), ghostoffmap(4);    // for planar wave propagation problem
	//nMatrixT<int> ghostonmap(7), ghostoffmap(7);	// 3D fracture problem
	ghostonmap = 0;
	ghostoffmap = 0;
	ghostoffmap(0,1) = ghostoffmap(1,0) = 1;	// 3D wave propagation
	//ghostonmap(4,5) = ghostonmap(5,4) = 1;	// 3D fracture problem
	//ghostoffmap(0,6) = ghostoffmap(6,0) = ghostoffmap(1,6) = ghostoffmap(6,1) = 1;
	//ghostoffmap(2,6) = ghostoffmap(6,2) = ghostoffmap(3,6) = ghostoffmap(6,3) = 1;
	//ghostoffmap(4,6) = ghostoffmap(6,4) = ghostoffmap(5,6) = ghostoffmap(6,5) = 1;
	//ghostoffmap(1,0) = ghostoffmap(0,1) = 1;  // for wave propagation problem
	//ghostoffmap(4,0) = ghostoffmap(0,4) = ghostoffmap(4,1) = ghostoffmap(1,4) = 1;  // center MD crack
	//ghostoffmap(4,2) = ghostoffmap(2,4) = ghostoffmap(2,3) = ghostoffmap(3,2) = 1;
	//ghostoffmap(4,3) = ghostoffmap(3,4) = 1;
	//ghostonmap(2,3) = ghostonmap(3,2) = 1;
	//ghostoffmap(1,0) = ghostoffmap(0,1) = ghostoffmap(3,0) = ghostoffmap(0,3) = 1; // left edge MD crack
	//ghostoffmap(1,3) = ghostoffmap(3,1) = ghostoffmap(2,3) = ghostoffmap(3,2) = 1;
	//ghostonmap(1,0) = ghostonmap(0,1) = 1;

	if (nsd == 3)
	{
		/* set pointers to embedding force/electron density in FEManagerT_bridging atoms */
		//fFine_THK->SetExternalElecDensity(elecdens, fFine_THK->GhostNodes());
		//fFine_THK->SetExternalEmbedForce(embforce, fFine_THK->GhostNodes());
		//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
	}
	
	/* time managers */
	TimeManagerT* atom_time = fFine_THK->TimeManager();
	TimeManagerT* continuum_time = fCoarse->TimeManager();
	
	/* compute initial conditions */
	fFine->SetComputeInitialCondition(true);
	fCoarse->SetComputeInitialCondition(true);

	/* set to initial condition */
	fFine_THK->InitialCondition();
		
	/* calculate fine scale part of MD displacement and total displacement u */
	fCoarse->InitialProject(bridging_field, fine_node_manager, projectedu, order1);
	
	/* solve for initial FEM force f(u) as function of fine scale + FEM */
	/* use projected totalu instead of totalu for initial FEM displacements */
	nMatrixT<int>& promap = fFine_THK->PropertiesMap(0);   // element group for particles = 0
	ghostonmap = promap; // copy original properties map
	promap = ghostoffmap;  // turn ghost atoms off for f(u) calculation
	fubig = InternalForce(bridging_field, projectedu, *fFine_THK);
	promap = ghostonmap;   // turn ghost atoms back on for MD force calculations

	/* calculate global interpolation matrix ntf */
	fCoarse->Ntf(ntf, fFine_THK->NonGhostNodes(), activefenodes);
		
	/* compute FEM RHS force as product of ntf and fu */
//	dArrayT fx(ntf.Rows()), fy(ntf.Rows()), fz(ntf.Rows()), tempx(ntf.Cols()), tempy(ntf.Cols()), tempz(ntf.Cols());
	dArrayT fx(ntf.Rows()), tempx(ntf.Cols());
	dArray2DT ntfproduct(ntf.Rows(), nsd);
	fu.Dimension((fFine_THK->NonGhostNodes()).Length(), nsd);
	fu.RowCollect(fFine_THK->NonGhostNodes(), fubig); 
	
	for (int i = 0; i < nsd; i++) {
		fu.ColumnCopy(i, tempx);
		ntf.Multx(tempx, fx);
		ntfproduct.SetColumn(i, fx);
	}
	
//old implementation
#if 0
	fu.ColumnCopy(0,tempx);
	ntf.Multx(tempx,fx);
	ntfproduct.SetColumn(0,fx);
	fu.ColumnCopy(1,tempy);
	ntf.Multx(tempy,fy);
	ntfproduct.SetColumn(1,fy);
	if (nsd == 3)
	{
		fu.ColumnCopy(2,tempz);
		ntf.Multx(tempz,fz);
		ntfproduct.SetColumn(2,fz);
	}
#endif

	/* Add FEM RHS force to RHS using SetExternalForce */
	fCoarse->SetExternalForce(bridging_field, ntfproduct, activefenodes);
	
	/* now d0, v0 and a0 are known after InitialCondition */
	fCoarse->InitialCondition();
		
	if (nsd == 3)
	{
		/* Calculate EAM electron density/embedding terms for ghost atoms using continuum information */
		//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
	}
		
	/* interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
	fCoarse->InterpolateField(bridging_field, order1, boundghostdisp);
	fCoarse->InterpolateField(bridging_field, order2, boundghostvel);
	fCoarse->InterpolateField(bridging_field, order3, boundghostacc);
		
	/* sort boundary + ghost atom info into separate arrays */
	gadisp.RowCollect(gatoms, boundghostdisp);
	gavel.RowCollect(gatoms, boundghostvel);
	gaacc.RowCollect(gatoms, boundghostacc);
	badisp.RowCollect(batoms, boundghostdisp);
	bavel.RowCollect(batoms, boundghostvel);
	baacc.RowCollect(batoms, boundghostacc);

	/* Removed fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp); here */	
			
	if (nsd == 2)
	{
		/* store initial MD boundary displacement histories */
		thkforce = fFine_THK->THKForce(bridging_field, badisp);
		fFine_THK->SetExternalForce(bridging_field, thkforce, boundatoms);  // sets pointer to thkforce 
	}
	else
	{
		/* thkdisp = fine scale part of ghost atom displacement */
		thkforce = fFine_THK->THKDisp(bridging_field, badisp);
		fFine_THK->SetExternalForce(bridging_field, thkforce, boundatoms);
		//totaldisp+=gadisp;	// add FEM coarse scale part of displacement
		//fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, totaldisp);
	}
		
	/* figure out timestep ratio between fem and md simulations */
	int nfesteps = continuum_time->NumberOfSteps();
	double mddt = atom_time->TimeStep();
	double fedt = continuum_time->TimeStep();
	double d_ratio = fedt/mddt;		
	int ratio = int((2.0*d_ratio + 1.0)/2.0);
		
	/* running status flag */
	ExceptionT::CodeT error = ExceptionT::kNoError;	

	for (int i = 0; i < nfesteps; i++)	
	{
		for (int j = 0; j < ratio; j++)	// MD update first
		{
			atom_time->Step();	
								
			/* initialize step */
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->InitStep();
				
			/* update FEM solution interpolated at boundary atoms and ghost atoms assuming 
			constant acceleration - because of constant acceleration assumption, predictor and 
			corrector are combined into one function */
			fFine_THK->BAPredictAndCorrect(mddt, badisp, bavel, baacc);
			fFine_THK->BAPredictAndCorrect(mddt, gadisp, gavel, gaacc);

			if (nsd == 2)
			{
				/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);	
				
				/* calculate THK force on boundary atoms, update displacement histories */
				thkforce = fFine_THK->THKForce(bridging_field, badisp);  // SetExternalForce set via pointer
			}
			else
			{
				/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
				fFine_THK->SetFieldValues(bridging_field, fFine_THK->GhostNodes(), order1, gadisp);
					
				/* calculate thkforces */
				thkforce = fFine_THK->THKDisp(bridging_field, badisp);
			}
				
			/* solve MD equations of motion */
			if (1 || error == ExceptionT::kNoError) {
					fFine_THK->ResetCumulativeUpdate(group);
					error = fFine_THK->SolveStep();
			}

			/* close  md step */
			if (1 || error == ExceptionT::kNoError) error = fFine_THK->CloseStep();    
		}
		continuum_time->Step();

		/* initialize step */
		if (1 || error == ExceptionT::kNoError) error = fCoarse->InitStep();
            
		/* calculate total displacement u = FE + fine scale here using updated FEM displacement */
		fCoarse->BridgingFields(bridging_field, fine_node_manager, *(fCoarse->NodeManager()), totalu);
		
		/* calculate FE internal force as function of total displacement u here */
		promap = ghostoffmap;  // turn off ghost atoms for f(u) calculations
		fubig = InternalForce(bridging_field, totalu, *fFine_THK);
		promap = ghostonmap;   // turn on ghost atoms for MD force calculations
		fu.RowCollect(fFine_THK->NonGhostNodes(), fubig); 

		for (int i = 0; i < nsd; i++) {
			fu.ColumnCopy(i, tempx);
			ntf.Multx(tempx, fx);
			ntfproduct.SetColumn(i, fx);
		}

//old implementation
#if 0
		fu.ColumnCopy(0,tempx);
		ntf.Multx(tempx,fx);
		ntfproduct.SetColumn(0,fx);
		fu.ColumnCopy(1,tempy);
		ntf.Multx(tempy,fy);
		ntfproduct.SetColumn(1,fy);	// SetExternalForce updated via pointer
		if (nsd == 3)
		{
			fu.ColumnCopy(2,tempz);
			ntf.Multx(tempz,fz);
			ntfproduct.SetColumn(2,fz);
		}
#endif
			
		/* no need to call SetExternalForce again due to pointers to ntfproduct */
			
		/* solve FE equation of motion using internal force just calculated */
		if (1 || error == ExceptionT::kNoError) {
				fCoarse->ResetCumulativeUpdate(group);
				error = fCoarse->SolveStep();
		}

		/* Interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
		fCoarse->InterpolateField(bridging_field, order1, boundghostdisp);
		fCoarse->InterpolateField(bridging_field, order2, boundghostvel);
		fCoarse->InterpolateField(bridging_field, order3, boundghostacc);
			
		/* sort boundary + ghost atom info into separate arrays */
		gadisp.RowCollect(gatoms, boundghostdisp);
		gavel.RowCollect(gatoms, boundghostvel);
		gaacc.RowCollect(gatoms, boundghostacc);
		badisp.RowCollect(batoms, boundghostdisp);
		bavel.RowCollect(batoms, boundghostvel);
		baacc.RowCollect(batoms, boundghostacc);
			
		if (nsd == 3)
		{
			/* Calculate EAM electron density/embedding terms for ghost atoms using updated continuum information */
			//fCoarse->ElecDensity(gatoms.Length(), elecdens, embforce);
		}
			
		/* close fe step */
		if (1 || error == ExceptionT::kNoError) error = fCoarse->CloseStep();                   
	}

	/* check for error */
	if (0) ExceptionT::GeneralFail(caller, "hit error %d", error);
}

/* accept parameter list */
void BridgingScaleManagerT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "BridgingScaleManagerT::TakeParameterList";

	/* inherited */
	MultiManagerT::TakeParameterList(list);

	/* resolve atomistic solver */
	fFine_THK = TB_DYNAMIC_CAST(FEManagerT_THK*, fFine);
	if (!fFine_THK || fFine_THK->Name() != "tahoe_THK")
		ExceptionT::GeneralFail(caller, "expecting \"tahoe_THK\" fine scale not \"%s\"",
			fFine_THK->Name().Pointer());
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* calculate internal force for the given displacement u */
const dArray2DT& BridgingScaleManagerT::InternalForce(const StringT& field_name, const dArray2DT& field_values, 
	FEManagerT_bridging& bridging) const
{
	/* first obtain the MD displacement field */
	FieldT* field = bridging.NodeManager()->Field(field_name);
	if (!field) 
		ExceptionT::GeneralFail("BridgingScaleManagerT::InternalForce", "field \"%s\" not found",
			field_name.Pointer());
	
	dArray2DT disp_0 = (*field)[0];	// temporarily store current displacements
	int group = 0; // assume particle group number = 0
	
	/* obtain atom node list - can calculate once and store... */
	int nnd = field_values.MajorDim();
	iArrayT nodes(nnd);
	nodes.SetValueToPosition();
		
	/* now write total bridging scale displacement u into field */
	int order = 0;	// write displacement only
	bridging.SetFieldValues(field_name, nodes, order, field_values);
		
	/* compute RHS - ParticlePairT fForce calculated by this call */
	bridging.FormRHS(group);
	
	/* write actual MD displacements back into field */
	bridging.SetFieldValues(field_name, nodes, order, disp_0);

	/* get the internal force contribution associated with the last call to FormRHS */
	return bridging.InternalForce(group);
}

#endif /* BRIDGING_ELEMENT */
