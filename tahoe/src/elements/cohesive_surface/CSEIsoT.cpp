/* $Id: CSEIsoT.cpp,v 1.1.1.1 2001-01-29 08:20:37 paklein Exp $ */
/* created: paklein (11/19/1997)                                          */
/* Cohesive surface elements with scalar traction potentials,             */
/* i.e., the traction potential is a function of the gap magnitude,       */
/* or effective gap magnitude only.                                       */

#include "CSEIsoT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "SurfaceShapeT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "C1FunctionT.h"
#include "eControllerT.h"

/* potential functions */
#include "LennardJones612.h"
#include "SmithFerrante.h"

/* constructor */
CSEIsoT::CSEIsoT(FEManagerT& fe_manager):
	CSEBaseT(fe_manager)
{

}

/* form of tangent matrix */
GlobalT::SystemTypeT CSEIsoT::TangentType(void) const
{
	/* tangent matrix is not symmetric */
	return GlobalT::kSymmetric;
}

void CSEIsoT::Initialize(void)
{
	/* inherited */
	CSEBaseT::Initialize();

	/* check output codes */
	if (fOutputCodes[MaterialData])
	{
		cout << "\n CSEIsoT::Initialize: material outputs not supported, overriding" << endl;
		fOutputCodes[MaterialData] = IOBaseT::kAtNever;
	}

	/* streams */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();
	
	/* construct props */
	int numpots;
	in >> numpots;
	fSurfPots.Allocate(numpots);
	for (int i = 0; i < fSurfPots.Length(); i++)
	{
		int num, code;
		in >> num >> code;
		
		/* check for repeated number */
		if (fSurfPots[--num] != NULL) throw eBadInputValue;
	
		/* construct surface potential function */
		switch (code)
		{
			case C1FunctionT::kLennardJones:
			{	
				double A;
				in >> A;
				
				fSurfPots[num] = new LennardJones612(A);
				break;
			}	
			case C1FunctionT::kSmithFerrante:
			{
				double A, B;
				in >> A >> B;
			
				fSurfPots[num] = new SmithFerrante(A,B,0.0);
				break;
			}
			default:
			
				throw eBadInputValue;	
		}
	
		if (!fSurfPots[num]) throw eOutOfMemory;
	}

	/* echo */
	out << "\n Cohesive surface potentials:\n";
	out << " Number of potentials. . . . . . . . . . . . . . = ";
	out << fSurfPots.Length() << '\n';
	for (int j = 0; j < fSurfPots.Length(); j++)
	{
		out << "\n Potential number. . . . . . . . . . . . . . . . = "
		    << j + 1 << '\n';
		out << " Potential name:\n";
		fSurfPots[j]->PrintName(out);
		fSurfPots[j]->Print(out);
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* called by FormRHS and FormLHS */
void CSEIsoT::LHSDriver(void)
{
	/* algorithmic constants */
	double constK = 0.0;
	int     formK = fController->FormK(constK);
	if (!formK) return;

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

	/* loop over elements */
	Top();
	while ( NextElement() )
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
	
		/* surface potential */
		C1FunctionT* surfpot = fSurfPots[element.MaterialNumber()];

		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
	  	
	  	/* initialize */
		fLHS = 0.0;
	
		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* gap vector from facet1 to facet2 */
			const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);
			double           d = gap.Magnitude();
			double           j = fShapes->Jacobian();
			double           w = fShapes->IPWeight();
					
			if (fabs(d) > kSmall)
			{
				double k  = j*w*constK;			
				double k2 = k*(surfpot->DFunction(d))/d;
				double k1 = (k*surfpot->DDFunction(d) - k2)/(d*d);
				
				/* accumulate */
				fShapes->Grad_d().MultTx(gap,fNEEvec);
				fNEEmat.Outer(fNEEvec,fNEEvec);
				fLHS.AddScaled(k1,fNEEmat);
				
				fLHS.AddScaled(k2,fShapes->Grad_dTGrad_d());
			}
			else /* initial stiffness */
			{				
				double ddphi = j*w*constK*(surfpot->DDFunction(0.0));
				fLHS.AddScaled(ddphi,fShapes->Grad_dTGrad_d());
			}
		}
									
		/* assemble */
		AssembleLHS();
	}
}

void CSEIsoT::RHSDriver(void)
{
	/* time-stepping parameters */
	double constKd = 0.0;
	int     formKd = fController->FormKd(constKd);
	if (!formKd) return;
	
	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);
	
	Top();
	while ( NextElement() )
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
	
		/* surface potential */
		C1FunctionT* surfpot = fSurfPots[element.MaterialNumber()];

		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

		/* initialize */
		fRHS = 0.0;
	
		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* gap vector from facet1 to facet2 */
			const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);
			double           d = gap.Magnitude();
			
			if (fabs(d) > kSmall)
			{
				double    j = fShapes->Jacobian();
				double    w = fShapes->IPWeight();
				double dphi =-j*w*constKd*(surfpot->DFunction(d));
				
				/* accumulate */
				fShapes->Grad_d().MultTx(gap,fNEEvec);
				fRHS.AddScaled(dphi/d,fNEEvec);
			}
		}
									
		/* assemble */
		AssembleRHS();
	}
}

/* extrapolate the integration point stresses and strains and extrapolate */
void CSEIsoT::ComputeNodalValues(const iArrayT& codes)
{
/* number of nodally smoothed values */
int num_out = codes.Sum();

	/* nothing to output */
	if (num_out == 0) return;

	/* work arrays */
	dArray2DT nodal_space(fNumElemNodes, num_out);
	dArray2DT nodal_all(fNumElemNodes, num_out);
	dArray2DT coords, disp;
	dArray2DT jump, Tmag;

	/* ip values */
	dArrayT ipjump(1), ipTmag(1);
	LocalArrayT loc_init_coords(LocalArrayT::kInitCoords, fNumElemNodes, fNumSD);
	LocalArrayT loc_disp(LocalArrayT::kDisp, fNumElemNodes, fNumDOF);
	fFEManager.RegisterLocal(loc_init_coords);
	fFEManager.RegisterLocal(loc_disp);

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(fNumElemNodes, codes[NodalCoord], pall);
	pall += coords.Length();
	disp.Set(fNumElemNodes, codes[NodalDisp], pall);
	pall += disp.Length();
	jump.Set(fNumElemNodes, codes[NodalDispJump], pall);
	pall += jump.Length();
	Tmag.Set(fNumElemNodes, codes[NodalTraction], pall);

	Top();
	while (NextElement())
	{
		/* initialize */
	    nodal_space = 0.0;

		/* coordinates for whole element */
		if (codes[NodalCoord])
		{
			SetLocalX(loc_init_coords);
			loc_init_coords.ReturnTranspose(coords);
		}
		
		/* displacements for whole element */
		if (codes[NodalDisp])
		{
			SetLocalU(loc_disp);
			loc_disp.ReturnTranspose(disp);
		}

		/* gap and/or traction magnitude */
		if (codes[NodalDispJump] ||
		    codes[NodalTraction])
		{
	  		/* surface potential */
			C1FunctionT* surfpot = fSurfPots[CurrentElement().MaterialNumber()];

			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

			/* integrate */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* gap */
				const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);
				double d = gap.Magnitude();

				if (codes[NodalDispJump])
				{
					/* extrapolate ip values to nodes */
					ipjump[0] = d;
					fShapes->Extrapolate(ipjump, jump);				
				}

				/* traction */
				if (codes[NodalTraction])
				{
					/* traction magnitude */
					ipTmag[0] = surfpot->DFunction(d);

					/* extrapolate to nodes */
					fShapes->Extrapolate(ipTmag, Tmag);
				}			
			}
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(jump  , colcount); colcount += jump.MinorDim();
		nodal_all.BlockColumnCopyAt(Tmag  , colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		fNodes->AssembleAverage(CurrentElement().NodesX(), nodal_all);
	}
}
