/* $Id: SimoQ1P0.cpp,v 1.3 2002-10-05 20:10:45 paklein Exp $ */
#include "SimoQ1P0.h"

#include "ShapeFunctionT.h"
#include "StructuralMaterialT.h"
#include "StructuralMatListT.h"

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
	
	/* element pressure */
	fPressure.Dimension(NumElements());
	fPressure = 0.0;
	
	/* dimension work space */
	fMeanGradient.Dimension(NumSD(), NumElementNodes());
	fNEEmat.Dimension(fLHS);
	fdiff_b.Dimension(fGradNa);
	fb_bar.Dimension(fGradNa);
	fb_sig.Dimension(fGradNa);
	fStressStiff2.Dimension(fStressStiff);
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
			
			double tmp = v/(H*J);
			
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

	/* current element info */
	int el = CurrElementNumber();
	double v = fElementVolume[el];
	double p_bar = fPressure[el];

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	
	//TEMP
	dMatrixT One(NumSD());
	
	//TEMP
#if 0	
	dMatrixT tmp_1(fNEEmat);
	dMatrixT tmp_2(fNEEmat);
	dMatrixT tmp_3(fNEEmat);
	dMatrixT tmp_4(fStressStiff);
	tmp_1 = 0.0;
	tmp_2 = 0.0;
	tmp_3 = 0.0;
	tmp_4 = 0.0;
#endif

	fCurrShapes->GradNa(fMeanGradient, fb_bar);	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		/* double scale factor */
		double scale = constK*(*Det++)*(*Weight++);
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		cauchy.ToMatrix(fCauchyStress);
		double p = fCurrMaterial->Pressure();

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
		
		fb_sig.MultAB(fCauchyStress, fGradNa);

		/* integration constants */		
		fCauchyStress *= scale;
	
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
		
		/* $div div$ term */	
		fNEEmat.Outer(fGradNa, fGradNa);
		fLHS.AddScaled(p_bar*scale, fNEEmat);
		
//		tmp_1.AddScaled(p_bar*scale, fNEEmat);
		
		fdiff_b.DiffOf(fGradNa, fb_bar);
		fNEEmat.Outer(fdiff_b, fdiff_b);
		fLHS.AddScaled(scale*2.0*p/3.0, fNEEmat);
		
//		tmp_2.AddScaled(scale*2.0*p/3.0, fNEEmat);
		
		fNEEmat.Outer(fb_sig, fdiff_b);
		fNEEmat.Symmetrize();
		fLHS.AddScaled(-scale*4.0/3.0, fNEEmat);

//		tmp_3.AddScaled(-scale*4.0/3.0, fNEEmat);
		
//		One.Identity(scale*(p - p_bar));
//		fStressStiff.MultQTBQ(fGradNa, One,
//			format, dMatrixT::kAccumulate);
		bSp_bRq_to_KSqRp(fGradNa, fNEEmat);
		fLHS.AddScaled(scale*(p - p_bar), fNEEmat);

//		tmp_4.MultQTBQ(fGradNa, One,
//			format, dMatrixT::kAccumulate);
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF());
	
	/* $\bar{div}\bar{div}$ term */
	fNEEmat.Outer(fb_bar, fb_bar);
	fLHS.AddScaled(-p_bar*v, fNEEmat);

//	tmp_1.AddScaled(-p_bar*v, fNEEmat);

#if 0	
	cout << "||tmp_1|| = " << tmp_1.ScalarProduct() << '\n';
	cout << "||tmp_2|| = " << tmp_2.ScalarProduct() << '\n';
	cout << "||tmp_3|| = " << tmp_3.ScalarProduct() << '\n';
	cout << "||tmp_4|| = " << tmp_4.ScalarProduct() << '\n';
	cout << endl;
#endif
}

/* calculate the internal force contribution ("-k*d") */
void SimoQ1P0::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	/* constant pressure */
	double& p_bar = fPressure[CurrElementNumber()];
	p_bar = 0.0;

	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		/* strain displacement matrix */
		Set_B_bar(fCurrShapes->Derivatives_U(), fMeanGradient, fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);
		p_bar += (*Weight)*(*Det)*fCurrMaterial->Pressure();

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}
	
	/* volume averaged */
	p_bar /= fElementVolume[CurrElementNumber()];
}

/* read materials data */
void SimoQ1P0::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	UpdatedLagrangianT::ReadMaterialData(in);

	/* make sure 2D materials are plane strain */
	if (StructuralMaterialList().HasPlaneStress()) {
		cout << "\n SimoQ1P0::ReadMaterialData: 2D materials must be plane strain" << endl;
		throw eBadInputValue;					
	}	
}

/***********************************************************************
 * Private
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
		mean_gradient.AddScaled(w[i]*det[i]/v, fCurrShapes->Derivatives_U(i));
}

void SimoQ1P0::bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const
{
	//dim check
	int dim = K.Rows();
	int sub_dim = b.Rows();
	int S = 0;
	int p = 0;
	for (int i = 0; i < dim; i++)
	{
		int R = 0;
		int q = 0;
		for (int j = 0; j < dim; j++)
		{
			K(i,j) = b(q,S)*b(p,R);
		
			q++;
			if (q == sub_dim) {
				R++;
				q = 0;
			}
		}
		p++;
		if (p == sub_dim) {
			S++;
			p = 0;
		}
	}	
}
