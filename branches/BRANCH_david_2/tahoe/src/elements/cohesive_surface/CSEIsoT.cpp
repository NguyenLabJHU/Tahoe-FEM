/* $Id: CSEIsoT.cpp,v 1.21 2004-07-15 08:25:57 paklein Exp $ */
/* created: paklein (11/19/1997) */
#include "CSEIsoT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ElementSupportT.h"

#include "toolboxConstants.h"
#include "SurfaceShapeT.h"
#include "C1FunctionT.h"
#ifndef _FRACTURE_INTERFACE_LIBRARY_
#include "eIntegratorT.h"
#endif
#include "ParameterContainerT.h"

/* potential functions */
#include "LennardJones612.h"
#include "SmithFerrante.h"

using namespace Tahoe;

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* constructor */
CSEIsoT::CSEIsoT(const ElementSupportT& support):
	CSEBaseT(support)
{
	SetName("isotropic_CSE");
}
#else
CSEIsoT::CSEIsoT(ElementSupportT& support):
	CSEBaseT(support)
{
	SetName("isotropic_CSE");
}
#endif

/* form of tangent matrix */
GlobalT::SystemTypeT CSEIsoT::TangentType(void) const
{
	/* tangent matrix is symmetric */
	return GlobalT::kSymmetric;
}

/* information about subordinate parameter lists */
void CSEIsoT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	CSEBaseT::DefineSubs(sub_list);

	/* element block/material specification */
	sub_list.AddSub("isotropic_CSE_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
void CSEIsoT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "surface_potential")
	{
		/* choice */
		order = ParameterListT::Choice;
		
		/* function types */
		sub_lists.AddSub("Lennard-Jones_6-12");
		sub_lists.AddSub("Smith-Ferrante");
	}
	else /* inherited */
		CSEBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* CSEIsoT::NewSub(const StringT& name) const
{
	/* try surface potential */
	C1FunctionT* surf_pot = C1FunctionT::New(name);
	if (surf_pot)
		return surf_pot;

	/* other lists */
	if (name == "isotropic_CSE_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);

		/* choice of materials lists (inline) */
		block->AddSub("surface_potential", ParameterListT::Once, true);

		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}		
	else /* inherited */
		return CSEBaseT::NewSub(name);
}

/* accept parameter list */
void CSEIsoT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CSEIsoT::TakeParameterList";

	/* inherited */
	CSEBaseT::TakeParameterList(list);

	/* check output codes */
	if (fNodalOutputCodes[MaterialData])
		fNodalOutputCodes[MaterialData] = IOBaseT::kAtNever; /* not supported */

	/* construct list of potentials - one per block */
	int num_block = list.NumLists("isotropic_CSE_element_block");
	fSurfPots.Dimension(num_block);
	for (int i = 0; i < fSurfPots.Length(); i++) {

		/* block information */
		const ParameterListT& block = list.GetList("isotropic_CSE_element_block", i);

		/* resolve material choice */
		const ParameterListT& surf_pot_params = block.GetListChoice(*this, "surface_potential");

		/* construct material */
		C1FunctionT* surf_pot = C1FunctionT::New(surf_pot_params.Name());
		if (!surf_pot)
			ExceptionT::BadInputValue(caller, "could not construct \"%s\"", surf_pot_params.Name().Pointer());
		surf_pot->TakeParameterList(surf_pot_params);

		/* keep */
		fSurfPots[i] = surf_pot;
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* called by FormRHS and FormLHS */
void CSEIsoT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time-stepping parameters */
	double constK = 0.0;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int     formK = fIntegrator->FormK(constK);
	if (!formK) return;
#endif

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
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;
#endif
	
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
void CSEIsoT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* number of nodally smoothed values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();
#else
	int n_out = 0;
	int e_out = 0;
#endif

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	int nsd = NumSD();
	int ndof = NumDOF();
	int nen = NumElementNodes();

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_out);

	/* allocate element results space */
	e_values.Dimension(NumElements(), e_out);
	e_values = 0.0;

	/* work arrays */
	dArray2DT nodal_space(nen, n_out);
	dArray2DT nodal_all(nen, n_out);
	dArray2DT coords, disp;
	dArray2DT jump, Tmag;

	/* ip values */
	dArrayT ipjump(1), ipTmag(1);
	LocalArrayT loc_init_coords(LocalArrayT::kInitCoords, nen, nsd);
	LocalArrayT loc_disp(LocalArrayT::kDisp, nen, ndof);
	ElementSupport().RegisterCoordinates(loc_init_coords);
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	Field().RegisterLocal(loc_disp);
#else
#pragma message("What to do with loc_disp?")
#endif

	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(nen, n_codes[NodalCoord], pall);
	pall += coords.Length();
	disp.Set(nen, n_codes[NodalDisp], pall);
	pall += disp.Length();
	jump.Set(nen, n_codes[NodalDispJump], pall);
	pall += jump.Length();
	Tmag.Set(nen, n_codes[NodalTraction], pall);

	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT centroid;
	if (e_codes[Centroid])
	{
		centroid.Set(nsd, pall); 
		pall += nsd;
	}
	double e_tmp, area;
	double& phi = (e_codes[CohesiveEnergy]) ? *pall++ : e_tmp;
	double& T_e = (e_codes[Traction]) ? *pall++ : e_tmp;

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

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

		/* compute output */
		if (element.Flag() == kON)
		{
	  		/* surface potential */
			C1FunctionT* surfpot = fSurfPots[CurrentElement().MaterialNumber()];

			/* get ref geometry (1st facet only) */
			fNodes1.Collect(facet1, element.NodesX());
			fLocInitCoords1.SetLocal(fNodes1);

			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

			/* initialize element values */
			phi = area = 0.0;
			if (e_codes[Centroid]) centroid = 0.0;

			/* integrate */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* element integration weight */
				double ip_w = fShapes->Jacobian()*fShapes->IPWeight();
				area += ip_w;

				/* gap */
				const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);
				double d = gap.Magnitude();

				if (n_codes[NodalDispJump])
				{
					/* extrapolate ip values to nodes */
					ipjump[0] = d;
					fShapes->Extrapolate(ipjump, jump);				
				}

				/* traction */
				if (n_codes[NodalTraction] || e_codes[Traction])
				{
					/* compute traction */
					double t = surfpot->DFunction(d);
					
					/* extrapolate to nodes */
					if (n_codes[NodalTraction])
					{
						ipTmag[0] = t;
						fShapes->Extrapolate(ipTmag, Tmag);
					}
					
					/* element average */
					if (e_codes[Traction]) T_e += ip_w*t;
				}	
				
				/* moment */
				if (e_codes[Centroid])
					centroid.AddScaled(ip_w, fShapes->IPCoords());
				
				/* cohesive energy */
				if (e_codes[CohesiveEnergy])
				{
					/* surface potential */
					double potential = surfpot->Function(d);

					/* integrate */
					phi += potential*ip_w;
				}
			}
			
			/* element values */
			if (e_codes[Centroid]) centroid /= area;
			if (e_codes[Traction]) T_e /= area;
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
				C1FunctionT* surfpot = fSurfPots[CurrentElement().MaterialNumber()];
		
				/* integrate */
				fShapes->TopIP();
				while (fShapes->NextIP())
				{
					/* element integration weight */
					double ip_w = fShapes->Jacobian()*fShapes->IPWeight();
					area += ip_w;
		
					/* moment */
					if (e_codes[Centroid])
						centroid.AddScaled(ip_w, fShapes->IPCoords());
				}
				
				/* element values */
				if (e_codes[Centroid]) centroid /= area;
				
				/* cohesive energy */
				if (e_codes[CohesiveEnergy])
					phi = area*surfpot->Function(sqrt(area)*1000);
					//NOTE: assume "large" opening is related to the element size
			}		
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(jump  , colcount); colcount += jump.MinorDim();
		nodal_all.BlockColumnCopyAt(Tmag  , colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);
		
		/* store results */
		e_values.SetRow(CurrElementNumber(), element_values);		
	}

	/* get nodally averaged values */
	ElementSupport().OutputUsedAverage(n_values);
}
