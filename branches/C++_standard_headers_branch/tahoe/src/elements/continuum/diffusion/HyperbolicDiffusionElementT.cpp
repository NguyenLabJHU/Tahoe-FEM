/* $Id: HyperbolicDiffusionElementT.cpp,v 1.2 2005-01-05 01:25:07 paklein Exp $ */
#include "HyperbolicDiffusionElementT.h"

#include "ElementCardT.h"
#include "eIntegratorT.h"
#include "DiffusionMaterialT.h"
#include "Traction_CardT.h"

using namespace Tahoe;

/* constructor */
HyperbolicDiffusionElementT::HyperbolicDiffusionElementT(const ElementSupportT& support):
	DiffusionElementT(support),
	fTau(0.0),
	fLocAcc(LocalArrayT::kAcc)
{
	SetName("hyperbolic_diffusion");
}

/* describe the parameters needed by the interface */
void HyperbolicDiffusionElementT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	DiffusionElementT::DefineParameters(list);
	
	ParameterT tau(fTau, "tau");
	tau.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(tau);
}

/* accept parameter list */
void HyperbolicDiffusionElementT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	DiffusionElementT::TakeParameterList(list);

	/* check */
	int order = fIntegrator->Order();
	if (order != 0 && order != 2)
		ExceptionT::GeneralFail("HyperbolicDiffusionElementT::TakeParameterList",
			"expecting integrator order 0 or 2, not %d", order);

	/* relaxation time */
	fTau = list.GetParameter("tau");	

	/* set relaxation time in traction cards */
	for (int i = 0; i < fTractionList.Length(); i++)
		fTractionList[i].SetRelaxationTime(fTau);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initialize local arrays */
void HyperbolicDiffusionElementT::SetLocalArrays(void)
{
	/* inherited */
	DiffusionElementT::SetLocalArrays();

	/* allocate */
	fLocAcc.Dimension(NumElementNodes(), NumDOF());

	/* nodal accelerations */
	if (fIntegrator->Order() > 1)
		Field().RegisterLocal(fLocAcc);
}

/* construct the effective mass matrix */
void HyperbolicDiffusionElementT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	/* inherited */
	ContinuumElementT::LHSDriver(sys_type);

	/* set components and weights */
	double constM = 0.0;
	double constC = 0.0;
	double constK = 0.0;
	
	int formM = fIntegrator->FormM(constM);
	int formC = fIntegrator->FormC(constC);
	int formK = fIntegrator->FormK(constK);

	/* loop over elements */
	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
	{
		/* initialize */
		fLHS = 0.0;
		
		/* set shape function derivatives */
		SetGlobalShape();

		/* element mass */
		if (formM || formC) {
			double k = fCurrMaterial->Capacity()*(fTau*constM + constC);
			FormMass(kConsistentMass, k, axisymmetric, NULL);
		}

		/* element stiffness */
		if (formK) FormStiffness(constK);
	
		/* add to global equations */
		AssembleLHS();		
	}
}

void HyperbolicDiffusionElementT::RHSDriver(void)
{
	/* inherited */
	ContinuumElementT::RHSDriver();

	/* set components and weights */
	double constMa = 0.0;
	double constCv = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormCv(constMa);
	int formCv = fIntegrator->FormCv(constCv);
	int formKd = fIntegrator->FormKd(constKd);

	/* block info - needed for source terms */
	int block_dex = 0;
	int block_count = 0;
	dArray2DT ip_source;
	const ElementBlockDataT* block_data = fBlockData.Pointer(block_dex);
	const dArray2DT* block_source = Field().Source(block_data->ID());
	if (block_source) ip_source.Dimension(NumIP(), 1);

	/* body forces */
	int formBody = 0;
	if ((fBodySchedule && fBody.Magnitude() > kSmall) || block_source) {	
		formBody = 1;
		if (!formCv) constCv = 1.0; // correct value ??
	}

	bool axisymmetric = Axisymmetric();
	double dt = ElementSupport().TimeStep();
	double by_dt = (fabs(dt) > kSmall) ? 1.0/dt: 0.0; /* for dt -> 0 */
	Top();
	while (NextElement())
	{
		/* capacity */
		double pc = fCurrMaterial->Capacity();

		/* reset block info (skip empty) */
		while (block_count == block_data->Dimension()) {
			block_data = fBlockData.Pointer(++block_dex);
			block_source = Field().Source(block_data->ID());
			block_count = 0;
		}
		
		/* convert heat increment/volume to unit of fLocVel (T/s) */
		if (block_source) {
			block_source->RowCopy(block_count, ip_source);
			ip_source *= by_dt/pc;
		}
		block_count++;
		
		/* initialize */
		fRHS = 0.0;

		/* global shape function values */
		SetGlobalShape();

		/* conduction term */
		if (formKd) 
		{
			SetLocalU(fLocDisp);
			FormKd(-constKd);
		}

		/* capacity term */
		if (formMa || formCv || formBody)
		{
			if (formMa)
				SetLocalU(fLocAcc);
			else
				fLocAcc = 0.0;
			
			if (formCv) {
				SetLocalU(fLocVel);
				fLocAcc.AddScaled(pc, fLocVel);
			}

			if (formBody) AddBodyForce(fLocVel);

			/* add internal contribution */
			FormMa(kConsistentMass, -constCv, axisymmetric,
				&fLocAcc,
				(block_source) ? &ip_source : NULL,
				NULL);			  		
		}
				
		/* assemble */
		AssembleRHS();
	}
}
