/* $Id: CSEAnisoT.cpp,v 1.47 2003-04-22 19:02:05 cjkimme Exp $ */
/* created: paklein (11/19/1997) */
#include "CSEAnisoT.h"

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
#ifndef _SIERRA_TEST_
#include "eIntegratorT.h"
#include "NodeManagerT.h"
#endif
#include "ElementSupportT.h"
#include "dSymMatrixT.h"

/* potential functions */
#ifndef _SIERRA_TEST_
#include "XuNeedleman2DT.h"
#include "TvergHutch2DT.h"
#include "ViscTvergHutch2DT.h"
#include "Tijssens2DT.h"
#include "RateDep2DT.h"
#include "TiedPotentialT.h"
#include "TiedPotentialBaseT.h"
#include "YoonAllen2DT.h"
#endif

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
#include "InelasticDuctile2DT.h"
#include "MR2DT.h"
#include "MR_RP2DT.h"
#endif

#include "TvergHutch3DT.h"
#include "YoonAllen3DT.h"
#include "XuNeedleman3DT.h"

using namespace Tahoe;

#ifndef _SIERRA_TEST_
/* constructor */
CSEAnisoT::CSEAnisoT(const ElementSupportT& support, const FieldT& field, bool rotate):
	CSEBaseT(support, field),
	fRotate(rotate),
	fCurrShapes(NULL),
	fQ(NumSD()),
	fdelta(NumSD()),
	fT(NumSD()),
	fddU(NumSD()),
	fRunState(support.RunState())
{
	/* reset format for the element stiffness matrix */
	if (fRotate) fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}
#else
CSEAnisoT::CSEAnisoT(ElementSupportT& support, bool rotate):
	CSEBaseT(support),
	fRotate(rotate),
	fCurrShapes(NULL),
	fQ(NumSD()),
	fdelta(NumSD()),
	fT(NumSD()),
	fddU(NumSD())
{
	/* reset format for the element stiffness matrix */
	if (fRotate) fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}
#endif

/* destructor */
CSEAnisoT::~CSEAnisoT(void)
{
	if (fRotate)
	{
		delete fCurrShapes;
		fCurrShapes = NULL;
	}
}

/* form of tangent matrix */
GlobalT::SystemTypeT CSEAnisoT::TangentType(void) const
{
	if (fRotate)
		/* tangent matrix is not symmetric */
		return GlobalT::kNonSymmetric;
	else
		return GlobalT::kSymmetric;
}

void CSEAnisoT::Initialize(void)
{
	const char caller[] = "CSEAnisoT::Initialize";

	/* inherited */
	CSEBaseT::Initialize();
	
	/* rotating local frame */
	if (fRotate)
	{
		/* shape functions wrt. current coordinates (linked parent domains) */
		fCurrShapes = new SurfaceShapeT(*fShapes, fLocCurrCoords);
		if (!fCurrShapes) throw ExceptionT::kOutOfMemory;
		fCurrShapes->Initialize();
 		
		/* allocate work space */
		int nee = NumElementNodes()*NumDOF();
		fnsd_nee_1.Dimension(NumSD(), nee);
		fnsd_nee_2.Dimension(NumSD(), nee);
		fdQ.Dimension(NumSD());
		for (int k = 0; k < NumSD(); k++)
			fdQ[k].Dimension(NumSD(), nee);
	}
	else
		fCurrShapes = fShapes;

	/* streams */
	ifstreamT& in = ElementSupport().Input();
	ostream&   out = ElementSupport().Output();
#ifndef _SIERRA_TEST_		
	fCalcNodalInfo = false;

	/* construct props */
	int numprops;
	in >> numprops;
	
	fTiedPots.Dimension(numprops);
#else
	int numprops;
	numprops = 1;
	fCalcNodalInfo = false;
#endif

	fSurfPots.Dimension(numprops);
	fNumStateVariables.Dimension(numprops);

	for (int i = 0; i < fSurfPots.Length(); i++)
	{
		int num, code;
#ifndef _SIERRA_TEST_		
		in >> num >> code;
#else
		num = 1; 
		code = ElementSupport().ReturnInputInt(ElementSupportT::kMaterialCode);
#endif
		num--;

		/* check for repeated number */
		if (fSurfPots[num] != NULL) throw ExceptionT::kBadInputValue;

		switch (code)
		{
			case SurfacePotentialT::kXuNeedleman:
			{			
				if (NumDOF() == 2)
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new XuNeedleman2DT(in);			
#else
					throw ExceptionT::kBadInputValue;
#endif
				}
				else
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new XuNeedleman3DT(in);
#else
					dArrayT *params = ElementSupport().FloatInput();
					fSurfPots[num] = new XuNeedleman3DT(*params);
#endif
				}
				break;
			}
			case SurfacePotentialT::kTvergaardHutchinson:
			{
				if (NumDOF() == 2)
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new TvergHutch2DT(in);
#else
					throw ExceptionT::kBadInputValue;
#endif
				}
				else
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new TvergHutch3DT(in);
#else
					dArrayT *params = ElementSupport().FloatInput();
					fSurfPots[num] = new TvergHutch3DT(*params);				
#endif
				}
				break;
			}
#ifndef _SIERRA_TEST_
			case SurfacePotentialT::kViscTvergaardHutchinson:
			{
				if (NumDOF() == 2)
					fSurfPots[num] = new ViscTvergHutch2DT(in, ElementSupport().TimeStep());
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				break;
			}
			case SurfacePotentialT::kTijssens:
			{	
				if (NumDOF() == 2)
					fSurfPots[num] = new Tijssens2DT(in, ElementSupport().TimeStep());
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				break;
			}
			case SurfacePotentialT::kRateDep:
			{	
				if (NumDOF() == 2)
					fSurfPots[num] = new RateDep2DT(in, ElementSupport().TimeStep());
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				break;
			}
#endif
			case SurfacePotentialT::kYoonAllen:
			{	
				if (NumDOF() == 2)
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new YoonAllen2DT(in, ElementSupport().TimeStep());
#else
					throw ExceptionT::kBadInputValue;
#endif
				}
				else
				{
#ifndef _SIERRA_TEST_
					fSurfPots[num] = new YoonAllen3DT(in, ElementSupport().TimeStep());
#else
					dArrayT *fparams = ElementSupport().FloatInput();
					iArrayT *iparams = ElementSupport().IntInput();
					fSurfPots[num] = new YoonAllen3DT(*fparams,*iparams,ElementSupport().TimeStep());
#endif
				}	
				break;
			}
			case SurfacePotentialT::kInelasticDuctile:
			{
#ifdef COHESIVE_SURFACE_ELEMENT_DEV
				if (NumDOF() == 2)
					fSurfPots[num] = new InelasticDuctile2DT(in, ElementSupport().TimeStep());
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				break;
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT_DEV not enabled: %d", code);
#endif
			}
			case SurfacePotentialT::kMR:
			{
#ifdef COHESIVE_SURFACE_ELEMENT_DEV
				if (NumDOF() == 2)
					fSurfPots[num] = new MR2DT(in);
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				break;
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT_DEV not enabled: %d", code);
#endif
			}
#ifndef _SIERRA_TEST_
			// Handle all tied potentials here till the code is finalized.
			case SurfacePotentialT::kTiedPotential:
			{	
				if (NumDOF() == 2)
				{
					if (freeNodeQ.Length() != 0)
						ExceptionT::GeneralFail(caller,"only 1 TiedNodes potential can be extant");
				
					fSurfPots[num] = new TiedPotentialT(in);
				}
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
				//Fall through
			} 
			case SurfacePotentialT::kMR_RP:
			{
				if (code != SurfacePotentialT::kTiedPotential)
				{
#ifdef COHESIVE_SURFACE_ELEMENT_DEV
				if (NumDOF() == 2)
				{
					if (freeNodeQ.Length() != 0)
						ExceptionT::GeneralFail(caller,"only 1 TiedNodes potential can be extant");
				
					fSurfPots[num] = new MR_RP2DT(in);
				}
				else
					ExceptionT::BadInputValue(caller, "potential not implemented for 3D: %d", code);
#else
				ExceptionT::BadInputValue(caller, "COHESIVE_SURFACE_ELEMENT_DEV not enabled: %d", code);
#endif
				}
				
				// Common allocation for tied potentials. 
				SurfacePotentialT* surfpot = fSurfPots[num];
				fTiedPots[num] = new TiedPotentialBaseT*;
				*fTiedPots[num] = dynamic_cast<TiedPotentialBaseT*>(surfpot);
				if (!fTiedPots[num]) // bad cast
					ExceptionT::GeneralFail(caller,"Unable to case tied potential");
				iTiedFlagIndex = surfpot->NumStateVariables()-1;

				freeNodeQ.Dimension(NumElements(),NumElementNodes());
				freeNodeQ = 0.;
				freeNodeQ_last = freeNodeQ;
					 
				/* Initialize things if a potential needs more info than the gap vector */
				if ((*fTiedPots[num])->NeedsNodalInfo()) 
				{
					fCalcNodalInfo = true;
					fNodalInfoCode = (*fTiedPots[num])->NodalQuantityNeeded();
					iBulkGroups = (*fTiedPots[num])->BulkGroups();
				}
				if ((*fTiedPots[num])->NodesMayRetie())
					qRetieNodes = true;
				else 
					qRetieNodes = false;
				
				/* utility array for stress smoothing over tied node pairs */
				otherInds.Dimension(NumElementNodes());
			  	if (NumSD() == 2)
			  	{
			  		otherInds[0] = 3; otherInds[1] = 2; otherInds[2] = 1; otherInds[3] = 0;
			  	} 
			  	else // 3D
			  	{
			  		otherInds[0] = 4; otherInds[1] = 5; otherInds[2] = 6; otherInds[3] = 7;
			  		otherInds[4] = 0; otherInds[5] = 1; otherInds[6] = 2; otherInds[7] = 3;
			  	}
				
				break;
			}
#endif // ndef _SIERRA_TEST_
			default:
#ifndef _SIERRA_TEST_
				cout << "\n CSEAnisoT::Initialize: unknown potential code: " << code << endl;
#endif
				throw ExceptionT::kBadInputValue;
		}
		if (!fSurfPots[num]) throw ExceptionT::kOutOfMemory;
		
		/* get number of state variables */
		fNumStateVariables[num] = fSurfPots[num]->NumStateVariables();
		  
	}

#ifndef _SIERRA_TEST_
	/* check compatibility of constitutive outputs */
	if (fSurfPots.Length() > 1 && fNodalOutputCodes[MaterialData])
		for (int k = 0; k < fSurfPots.Length(); k++)
		{
			const SurfacePotentialT* pot_k = fSurfPots[k];
			for (int i = k+1; i < fSurfPots.Length(); i++)
			{
				const SurfacePotentialT* pot_i = fSurfPots[i];
				if (!SurfacePotentialT::CompatibleOutput(*pot_k, *pot_i))
				{
					cout << "\n CSEAnisoT::Initialize: incompatible output between potentials\n"
					     <<   "     " << k+1 << " and " << i+1 << endl;
					throw ExceptionT::kBadInputValue;
				}
			}
		}
			
	/* write */
	out << " Rotating local coordinate frame . . . . . . . . = " <<
	    ((fRotate) ? "ACTIVE" : "INACTIVE") << '\n';
	out << "\n Cohesive surface potentials:\n";
	out << " Number of potentials. . . . . . . . . . . . . . = ";
	out << fSurfPots.Length() << '\n';
	for (int j = 0; j < fSurfPots.Length(); j++)
	{
		out << "\n Potential number. . . . . . . . . . . . . . . . = " << j + 1 << '\n';
		out << " Potential name:\n";
		fSurfPots[j]->PrintName(out);
		fSurfPots[j]->Print(out);
	}
#endif
	
	/* initialize state variable space */
	if (fNumStateVariables.Min() > 0)
	{
		/* number of integration points */
		int num_ip = fCurrShapes->NumIP();
	
		/* get state variables per element */
		int num_elements = fElementCards.Length();
		iArrayT num_elem_state(num_elements);
		for (int i = 0; i < num_elements; i++)
			num_elem_state[i] = num_ip*fNumStateVariables[fElementCards[i].MaterialNumber()];

#ifndef _SIERRA_TEST_
		/* allocate space */
		fStateVariables.Configure(num_elem_state);
#else
		fStateVariables.Set(1,num_elem_state[0],ElementSupport().StateVariableArray());
#endif

		/* initialize state variable space */
		dArrayT state;
		for (int i = 0; i < num_elements; i++)
		{
			/* material number */
			int mat_num = fElementCards[i].MaterialNumber();
			int num_var = fNumStateVariables[mat_num];
			
			/* loop over integration points */
			double* pstate = fStateVariables(i);
			for (int j = 0; j < num_ip; j++)
			{
				state.Set(num_var, pstate);
				fSurfPots[mat_num]->InitStateVariables(state);
				pstate += num_var;
			}
		}		
	}
	else /* set dimensions to zero */
		fStateVariables.Dimension(fElementCards.Length(), 0);

#ifndef _SIERRA_TEST_
	/* set history */
	fStateVariables_last = fStateVariables;
	/* For SIERRA, don't do anything. Wait until InitStep. */
#endif

}

#ifdef _SIERRA_TEST_	
/* Get state variables from ElementSupportT here */
void CSEAnisoT::InitStep(void) 
{
	
	fStateVariables_last.Set(fStateVariables.MajorDim(),fStateVariables.MinorDim(),
		ElementSupport().StateVariableArray());
};
#endif

/* close current time increment */
void CSEAnisoT::CloseStep(void)
{
	/* inherited */
	CSEBaseT::CloseStep();

#ifndef _SIERRA_TEST_
	/* reset state variables from history */
	fStateVariables_last = fStateVariables;

	if (freeNodeQ.IsAllocated())
		freeNodeQ_last = freeNodeQ;
#endif
}

#ifndef _SIERRA_TEST_
/* write restart data to the output stream. */
void CSEAnisoT::WriteRestart(ostream& out) const
{
	/* inherited */
	CSEBaseT::WriteRestart(out);
	
	/* write state variable data */
	fStateVariables.WriteData(out);
	out << '\n';
	
	out << freeNodeQ.Length() << '\n';
	out << freeNodeQ.wrap_tight(10) << endl;
}

/* read restart data to the output stream */
void CSEAnisoT::ReadRestart(istream& in)
{
	/* inherited */
	CSEBaseT::ReadRestart(in);

	/* read state variable data */
	fStateVariables.ReadData(in);

	/* set history */
	fStateVariables_last = fStateVariables;
	
	int freeNode_length;
	in >> freeNode_length;
	if (freeNodeQ.Length() != freeNode_length)
		ExceptionT::GeneralFail("CSEAnisoT::ReadRestart","Length mismatch for freeNodeQ");
	in >> freeNodeQ;
	freeNodeQ_last = freeNodeQ;
}

#else

void CSEAnisoT::WriteRestart(double* outgoingData) const
{
	/* inherited */
	CSEBaseT::WriteRestart(outgoingData);

	// Nothing to do here right now since Sierra controls state variables	
	/* write state variable data */
//	fStateVariables.WriteData(out);

}

/* read restart data to the output stream */
void CSEAnisoT::ReadRestart(double* incomingData)
{
	/* inherited */
	CSEBaseT::ReadRestart(incomingData);

	// Nothing to do here right now since Sierra controls state variables	
	
	/* read state variable data */
//	fStateVariables.ReadData(in);

	/* set history */
//	fStateVariables_last = fStateVariables;
//	if (freeNodeQ.IsAllocated()) //This is useless
//		freeNodeQ_last = freeNodeQ;
}
#endif

/***********************************************************************
* Protected
***********************************************************************/

void CSEAnisoT::LHSDriver(GlobalT::SystemTypeT)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format = (fRotate) ?
		dMatrixT::kWhole : dMatrixT::kUpperOnly;

	/* time-integration parameters */
	double constK = 1.0;
#ifndef _SIERRA_TEST_
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;
#endif

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

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
#ifndef _SIERRA_TEST_
		TiedPotentialBaseT* tiedpot;
		if (fTiedPots[element.MaterialNumber()] != NULL)
			tiedpot = *fTiedPots[element.MaterialNumber()];
		else
			tiedpot = NULL;
#endif

		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
		
		/* initialize */
		fLHS = 0.0;

		bool nodalReleaseQ = false;
		LocalArrayT fNodalValues(LocalArrayT::kUnspecified);
#ifndef _SIERRA_TEST_
		if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
		{
		  	iArrayT ndIndices = element.NodesX();
		  	int numElemNodes = ndIndices.Length();
		  	dArray2DT elementVals(numElemNodes,fNodalQuantities.MinorDim());
		  	fNodalValues.Dimension(numElemNodes,fNodalQuantities.MinorDim());
		  	for (int iIndex = 0; iIndex < numElemNodes; iIndex++) 
			{
				elementVals.SetRow(iIndex,0.);
			 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[iIndex]));
			 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[otherInds[iIndex]]));
			}
		  	fNodalValues.SetGlobal(elementVals);
		  	ndIndices.SetValueToPosition();
		  	fNodalValues.SetLocal(ndIndices);
		}
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

			/* check */
			if (j0 <= 0.0 || j <= 0.0) throw ExceptionT::kBadJacobianDet;
		
			/* gap vector and gradient (facet1 to facet2) */
			const dArrayT&    delta = fShapes->InterpolateJumpU(fLocCurrCoords);
			const dMatrixT& d_delta = fShapes->Grad_d();

			/* gap vector in local frame */
			fQ.MultTx(delta, fdelta);			
			
			/* Interpolate nodal info to IPs */
			/* stress tensor in local frame */
			dSymMatrixT localFrameIP(NumSD());
#ifndef _SIERRA_TEST_
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
			{
				dArrayT tensorIP(fNodalValues.MinorDim());
			    fShapes->Interpolate(fNodalValues,tensorIP);			    
				localFrameIP.MultQBQT(fQ,dSymMatrixT(NumSD(),tensorIP.Pointer()));
			}
#endif

			/* stiffness in local frame */
			const dMatrixT& K = surfpot->Stiffness(fdelta, state, localFrameIP);
			
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
				fddU.SetToScaled(j0*w*constK, K);
				fnsd_nee_1.MultATB(fQ, d_delta);
				fnsd_nee_2.MultATB(fddU, fnsd_nee_1);
				u_i__Q_ijk(delta, fdQ, fnsd_nee_1);
				fNEEmat.MultATB(fnsd_nee_2, fnsd_nee_1);
				fLHS += fNEEmat;
			}
			
			/* 3rd term */
			fddU.MultQBQT(fQ, K);
			fddU *= j0*w*constK;
			fLHS.MultQTBQ(d_delta, fddU, format, dMatrixT::kAccumulate);
		}

		/* assemble */
		AssembleLHS();
	}
}

void CSEAnisoT::RHSDriver(void)
{
	/* time-integration parameters */
	double constKd = 1.0;
#ifndef _SIERRA_TEST_
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* heat source if needed */
	const FieldT* temperature = ElementSupport().Field("temperature");

	/* initialize sources */
	if (temperature) {
	
		if (fIncrementalHeat.Length() == 0) {

			/* initialize heat source arrays */
			fIncrementalHeat.Dimension(fBlockData.Length());
			for (int i = 0; i < fIncrementalHeat.Length(); i++)
			{
				/* dimension */
				fIncrementalHeat[i].Dimension(fBlockData[i].Dimension(), fShapes->NumIP());
	
				/* register */
				temperature->RegisterSource(fBlockData[i].ID(), fIncrementalHeat[i]);
			}
		}
		
		/* clear sources */
		for (int i = 0; i < fIncrementalHeat.Length(); i++)
			fIncrementalHeat[i] = 0.0;
	}
#else // _SIERRA_TEST_ defined
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
	
#ifndef _SIERRA_TEST_
	/* If the potential needs info from the nodes, start to gather it now */
	if (fCalcNodalInfo) 
	{
		ElementBaseT& surroundingGroup = ElementSupport().ElementGroup(iBulkGroups[0]);
		surroundingGroup.SendOutput(fNodalInfoCode);
		if (fNodalQuantities.Length() > 0) 
		{
			fNodalQuantities.Free();
		}
		fNodalQuantities = ElementSupport().OutputAverage();
	}
#endif

	/* fracture surface area */
	fFractureArea = 0.0;

	int block_count = 0, block_dex = 0;
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
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);
	  			
		if (element.Flag() != kOFF)
		{
			/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
			int num_state = fNumStateVariables[element.MaterialNumber()];
	
			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
	
	  		/* initialize */
	  		fRHS = 0.0;

#ifndef _SIERRA_TEST_
			bool nodalReleaseQ = false;
			LocalArrayT fNodalValues(LocalArrayT::kUnspecified);
			TiedPotentialBaseT* tiedpot;
			if (fTiedPots[element.MaterialNumber()] != NULL)
				tiedpot = *fTiedPots[element.MaterialNumber()];
			else
				tiedpot = NULL;

			int currElNum = CurrElementNumber();
			iArrayT ndIndices;
			int numElemNodes;
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
			{
				ndIndices = element.NodesX();
			  	numElemNodes = ndIndices.Length();
			  	dArray2DT elementVals(numElemNodes,fNodalQuantities.MinorDim()); 	
			  	fNodalValues.Dimension(numElemNodes,fNodalQuantities.MinorDim());
			  	for (int iIndex = 0; iIndex < numElemNodes; iIndex++) 
				{
					elementVals.SetRow(iIndex,0.);
				 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[iIndex]));
				 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[otherInds[iIndex]]));
				}
			  	fNodalValues.SetGlobal(elementVals);
			  	ndIndices.SetValueToPosition();
			  	fNodalValues.SetLocal(ndIndices);
			}
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
				
				/* check */
				if (j0 <= 0.0 || j <= 0.0)
				{
					cout << "\n CSEAnisoT::RHSDriver: jacobian error" << endl;
					throw ExceptionT::kBadJacobianDet;
				}
	
				/* gap vector from facet1 to facet2 */
				const dArrayT& delta = fShapes->InterpolateJumpU(fLocCurrCoords);
	
				/* gap vector in local frame */
				fQ.MultTx(delta, fdelta);
					
				/* Interpolate nodal info to IPs */
				/* stress tensor in local frame */
				dSymMatrixT localFrameIP(NumSD());
#ifndef _SIERRA_TEST_
				if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
				{
					dArrayT tensorIP(fNodalValues.MinorDim());
				    fShapes->Interpolate(fNodalValues,tensorIP);			    
					localFrameIP.MultQBQT(fQ,dSymMatrixT(NumSD(),tensorIP.Pointer()));
					
					if (state[iTiedFlagIndex] == kTiedNode && 
						tiedpot->InitiationQ(localFrameIP.Pointer()))
					{
						for (int i = 0; i < numElemNodes; i++) 
							freeNodeQ(currElNum,i) = 1.;
						nodalReleaseQ = true;
						state[iTiedFlagIndex] = kReleaseNextStep;
					}
					else
						/* see if nodes need to be retied */
						if (qRetieNodes && state[iTiedFlagIndex] == kFreeNode && 
							tiedpot->RetieQ(localFrameIP.Pointer(), state, fdelta))
						{
							state[iTiedFlagIndex] = kTieNextStep;
							for (int i = 0; i < numElemNodes; i++)
								freeNodeQ(currElNum,i) = 0.;
						}
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
					
#ifndef _SIERRA_TEST_
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

/* nodal value calculations */
void CSEAnisoT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* inherited */
	CSEBaseT::SetNodalOutputCodes(mode, flags, counts);

	/* resize for vectors not magnitudes */
	if (flags[NodalDispJump] == mode)
		counts[NodalDispJump] = NumDOF();	
	if (flags[NodalTraction] == mode)
		counts[NodalTraction] = NumDOF();
	if (flags[MaterialData] == mode)
		counts[MaterialData] = fSurfPots[0]->NumOutputVariables();
}

void CSEAnisoT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* inherited */
	CSEBaseT::SetElementOutputCodes(mode, flags, counts);
	
	/* resize for vectors not magnitudes */
	if (flags[Traction] == mode) counts[Traction] = NumDOF();
}

void CSEAnisoT::SendOutput(int kincode)
{
	if (kincode != InternalData)
		CSEBaseT::SendOutput(kincode);
	else // TiedNodesT wants its freeNode info
		ComputeFreeNodesForOutput();
}

/* extrapolate the integration point stresses and strains and extrapolate */
void CSEAnisoT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
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
#ifndef _SIERRA_TEST_
	Field().RegisterLocal(loc_disp);
#else
#pragma message("This routine needs displacements from SIERRA.")
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

	AutoArrayT<double> state;
	Top();
	while (NextElement())
	{
		/* current element */
		ElementCardT& element = CurrentElement();
		
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
			fLocInitCoords1.SetLocal(fNodes1);

			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

			/* initialize element values */
			phi = area = 0.0;
			if (e_codes[Centroid]) centroid = 0.0;
			if (e_codes[Traction]) traction = 0.0;

			LocalArrayT fNodalValues(LocalArrayT::kUnspecified);
			bool nodalReleaseQ = false;
#ifndef _SIERRA_TEST_
			TiedPotentialBaseT* tiedpot;
			if (fTiedPots[element.MaterialNumber()] != NULL)
				tiedpot = *fTiedPots[element.MaterialNumber()];
			else
				tiedpot = NULL;
			if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
			{
			  	int numElemNodes = element.NodesX().Length();
			  	dArray2DT elementVals(numElemNodes,fNodalQuantities.MinorDim());
			  	iArrayT ndIndices = element.NodesX();
			  	fNodalValues.Dimension(numElemNodes,fNodalQuantities.MinorDim());
			  	for (int iIndex = 0; iIndex < numElemNodes; iIndex++) 
				{
					elementVals.SetRow(iIndex,0.);
				 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[iIndex]));
				 	elementVals.AddToRowScaled(iIndex,.5,fNodalQuantities(ndIndices[otherInds[iIndex]]));
				}
				fNodalValues.SetGlobal(elementVals);
				ndIndices.SetValueToPosition();
			  	fNodalValues.SetLocal(ndIndices);
			}
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
				const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);

				/* coordinate transformation */
				double j = fCurrShapes->Jacobian(fQ);
				fQ.MultTx(gap, fdelta);
				
				/* gap */				
				if (n_codes[NodalDispJump])
					fShapes->Extrapolate(fdelta, jump);
	     
				/* traction */
				if (n_codes[NodalTraction] || e_codes[Traction])
				{
					/* copy state variables (not integrated) */
					state.Set(num_state,pstate);

					/* Interpolate nodal info to IPs */
					/* stress tensor in local frame */
					dSymMatrixT localFrameIP(NumSD());
#ifndef _SIERRA_TEST_
					if (tiedpot != NULL && tiedpot->NeedsNodalInfo()) 
					{
						dArrayT tensorIP(fNodalValues.MinorDim());
					    fShapes->Interpolate(fNodalValues,tensorIP);			    
						localFrameIP.MultQBQT(fQ,dSymMatrixT(NumSD(),tensorIP.Pointer()));
					}
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

void CSEAnisoT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* inherited */
	CSEBaseT::GenerateOutputLabels(n_codes, n_labels, e_codes, e_labels);

	/* overwrite nodal labels */
	n_labels.Dimension(n_codes.Sum());
	int count = 0;
	if (n_codes[NodalDisp])
	{
#ifndef _SIERRA_TEST_
		/* labels from the field */
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
#else
		const char* labels[] = {"D_1", "D_2", "D_3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = labels[i];
#endif
	}

	if (n_codes[NodalCoord])
	{
		const char* xlabels[3] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[NodalDispJump])
	{
		const char* d_2D[2] = {"d_t", "d_n"};
		const char* d_3D[3] = {"d_t1", "d_t2", "d_n"};
		const char** dlabels;
		if (NumDOF() == 2)
			dlabels = d_2D;
		else if (NumDOF() == 3)
			dlabels = d_3D;
		else
			throw ExceptionT::kGeneralFail;

		for (int i = 0; i < NumDOF(); i++)
			n_labels[count++] = dlabels[i];
	}

	if (n_codes[NodalTraction])
	{
		const char* t_2D[2] = {"T_t", "T_n"};
		const char* t_3D[3] = {"T_t1", "T_t2", "T_n"};
		const char** tlabels;
		if (NumDOF() == 2)
			tlabels = t_2D;
		else if (NumDOF() == 3)
			tlabels = t_3D;
		else
			throw ExceptionT::kGeneralFail;

		for (int i = 0; i < NumDOF(); i++)
			n_labels[count++] = tlabels[i];
	}

	/* material output labels */
	if (n_codes[MaterialData])
	{
		ArrayT<StringT> matlabels;
		fSurfPots[0]->OutputLabels(matlabels);
		for (int i = 0; i < n_codes[MaterialData]; i++)
			n_labels[count++] = matlabels[i];
	}
	
	/* allocate nodal output labels */
	e_labels.Dimension(e_codes.Sum());
	count = 0;
	if (e_codes[Centroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < NumSD(); i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[CohesiveEnergy]) e_labels[count++] = "phi";
	if (e_codes[Traction])
	{
		const char* t_2D[2] = {"T_t", "T_n"};
		const char* t_3D[3] = {"T_t1", "T_t2", "T_n"};
		const char** tlabels;
		if (NumDOF() == 2)
			tlabels = t_2D;
		else if (NumDOF() == 3)
			tlabels = t_3D;
		else
			throw ExceptionT::kGeneralFail;

		for (int i = 0; i < NumDOF(); i++)
			e_labels[count++] = tlabels[i];
	}
}

/* write all current element information to the stream */
void CSEAnisoT::CurrElementInfo(ostream& out) const
{
#pragma unused(out)
#ifndef _SIERRA_TEST_
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
		out << " gap vector (global frame):\n";
		const dArrayT& delta = fShapes->InterpolateJumpU(fLocCurrCoords);
		out << delta << '\n';
		out << " gap vector (local frame):\n";
		dArrayT delta_temp(fdelta);
		fQ.MultTx(delta, delta_temp);
		out << delta_temp << '\n';	
	}
	
	catch (ExceptionT::CodeT error)
	{
		out << " CSEAnisoT::CurrElementInfo: error on surface jacobian\n";
	}
#else
	throw ExceptionT::kGeneralFail;
#endif
}

/***********************************************************************
* Private
***********************************************************************/

/* operations with pseudo rank 3 (list in j) matrices */
void CSEAnisoT::u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
	dMatrixT& Qu)
{
	for (int i = 0; i < u.Length(); i++)
	{	
		Q[i].MultTx(u, fNEEvec);
		Qu.SetRow(i, fNEEvec);
	}
}

void CSEAnisoT::Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
	dMatrixT& Qu)
{
	if (Q.Length() == 2)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1]);
	else if (Q.Length() == 3)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1], u[2], Q[2]);
	else
		throw ExceptionT::kGeneralFail;
}

void CSEAnisoT::ComputeFreeNodesForOutput()
{
	if (!freeNodeQ.IsAllocated())
		ExceptionT::GeneralFail("CSEAnisoT::ComputeFreeNodesForOutput","No TiedNodes Data!");

	ElementSupport().ResetAverage(1);
	
	int nen = NumElementNodes();
	dArray2DT oneElement(nen,1);
	Top();
	while (NextElement())
	{
		oneElement.Set(nen,1,freeNodeQ(CurrElementNumber()));
		ElementSupport().AssembleAverage(CurrentElement().NodesX(),oneElement);
	}
}
