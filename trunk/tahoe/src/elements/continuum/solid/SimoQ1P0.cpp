/* $Id: SimoQ1P0.cpp,v 1.2 2002-09-23 06:58:25 paklein Exp $ */
#include "SimoQ1P0.h"

#include "ShapeFunctionT.h"
#include "StructuralMaterialT.h"

using namespace Tahoe;

/* constructor */
SimoQ1P0::SimoQ1P0(const ElementSupportT& support, const FieldT& field):
	UpdatedLagrangianT(support, field),
	fF_tmp(NumSD()),
	fLastVolumeInit(false)
{

}

/* destructor */
SimoQ1P0::~SimoQ1P0(void)
{

}

/* data initialization */
void SimoQ1P0::Initialize(void)
{
	/* inherited */
	UpdatedLagrangianT::Initialize();

	/* check geometry code and number of element nodes -> Q1 */
	if (GeometryCode() == GeometryT::kQuadrilateral) {
		if (NumElementNodes() != 4) {
			cout << "\n SimoQ1P0::Initialize: expecting 4 node quad: " 
			     << NumElementNodes() << endl;
			throw eBadInputValue;
		}	
	}
	else if (GeometryCode() == GeometryT::kHexahedron) {
		if (NumElementNodes() != 8) {
			cout << "\n SimoQ1P0::Initialize: expecting 8 node hex: " 
			     << NumElementNodes() << endl;
			throw eBadInputValue;
		}	
	}
	else {
		cout << "\n SimoQ1P0::Initialize: expecting hex or quad geometry: "
		     << GeometryCode() << endl;
		throw eBadInputValue;
	}
	
	/* need to store last deformed element volume */
	fElementVolume.Dimension(NumElements());	
	fElementVolume = 0.0;
	fElementVolume_last.Dimension(NumElements());
	fElementVolume_last = 0.0;
	
	/* dimension work space */
	fMeanGradient.Dimension(NumSD(), NumElementNodes());
}

/* finalize current step - step is solved */
void SimoQ1P0::CloseStep(void)
{
	/* inherited */
	UpdatedLagrangianT::CloseStep();
	
	/* store converged solution */
	fElementVolume_last = fElementVolume;
	fLastVolumeInit = true;
}
	
/* restore last converged state */
void SimoQ1P0::ResetStep(void)
{
	/* inherited */
	UpdatedLagrangianT::ResetStep();
	
	/* store converged solution */
	fElementVolume = fElementVolume_last;
}

/* read restart information from stream */
void SimoQ1P0::ReadRestart(istream& in)
{
	/* inherited */
	UpdatedLagrangianT::ReadRestart(in);
	
	/* read restart data */
	in >> fElementVolume;
	
	/* reset last state */
	fElementVolume_last = fElementVolume;
	fLastVolumeInit = true;
}

/* write restart information from stream */
void SimoQ1P0::WriteRestart(ostream& out) const
{
	/* inherited */
	UpdatedLagrangianT::WriteRestart(out);
	
	/* read restart data */
	out << fElementVolume;
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form shape functions and derivatives */
void SimoQ1P0::SetGlobalShape(void)
{
	/* inherited - computes gradients and standard 
	 * deformation gradients */
	UpdatedLagrangianT::SetGlobalShape();

	/* compute mean of shape function gradients */
	double H; /* reference volume */
	double& v = fElementVolume[CurrElementNumber()];
	SetMeanGradient(fMeanGradient, H, v);
	
	/* last deformed volume */
	double& v_last = fElementVolume_last[CurrElementNumber()];
	if (!fLastVolumeInit) v_last = v;

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);

	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_List[i];
			double J = F.Det();
			F *= pow(v/(H*J), 1.0/3.0);
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_last_List[i];
			double J = F.Det();
			F *= pow(v_last/(H*J), 1.0/3.0);
		}
	}
}

/* form the element stiffness matrix */
void SimoQ1P0::FormStiffness(double constK)
{		
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* double scale factor */
		double scale = constK*(*Det++)*(*Weight++);
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fCauchyStress);
	
		/* integration constants */		
		fCauchyStress *= scale;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
	
		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fCauchyStress,
			format, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF());
}

/* calculate the internal force contribution ("-k*d") */
void SimoQ1P0::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* compute mean shape function gradient, Hughes (4.5.23) */
void SimoQ1P0::SetMeanGradient(dArray2DT& mean_gradient, double& H, double& v) const
{
	/* assume same integration rule defined for current and references
	 * shape functions */
	int nip = NumIP();
	const double*   det = fCurrShapes->IPDets();
	const double* det_0 = fShapes->IPDets();
	const double*     w = fShapes->IPWeights();

	/* H and current volume */
	H = 0.0;
	v = 0.0;
	for (int i = 0; i < nip; i++)
	{
		H += w[i]*det_0[i];
		v += w[i]*det[i];
	}

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_gradient.AddScaled(w[i]*det[i]/H, fCurrShapes->Derivatives_U(i));
}
