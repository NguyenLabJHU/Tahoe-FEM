/* $Id: SimoFiniteStrainT.cpp,v 1.3 2001-07-20 00:58:01 paklein Exp $ */
#include "SimoFiniteStrainT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "FEManagerT.h"
#include "StructuralMaterialT.h"
#include "MaterialListT.h" //NOTE - only needed for check in Initialize?
#include "SimoShapeFunctionT.h"

/* constructor */
SimoFiniteStrainT::SimoFiniteStrainT(FEManagerT& fe_manager):
	FiniteStrainT(fe_manager),
	fCurrElementModes(LocalArrayT::kUnspecified),
	fCurrElementModes_last(LocalArrayT::kUnspecified),
	fEnhancedShapes(NULL),
	fStressMat(fNumSD),
	fTempMat1(fNumSD),
	fTempMat2(fNumSD)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != 0)
	{
		cout << "\n SimoFiniteStrainT::SimoFiniteStrainT: no strain-displacement options\n" << endl;
		fStrainDispOpt = 0;
	}

	/* read parameters */
	ifstreamT& in = fe_manager.Input();
	int inc_mode = -1;
	in >> inc_mode;
	in >> fLocalIterationMax;
	in >> fAbsTol;
	in >> fRelTol;

	/* checks */
	if (inc_mode != 0 && inc_mode != 1) throw eBadInputValue;
	if (fLocalIterationMax < 1) throw eGeneralFail;
	fIncompressibleMode = (inc_mode == 1);
	if (NumSD() == 2) fIncompressibleMode = false;
	if (fAbsTol < 0 || fAbsTol > 1.0) throw eBadInputValue;
	if (fRelTol < 0 || fRelTol > 1.0) throw eBadInputValue;
	
	/* set number of mode shapes */
	if (NumSD() == 2)
		fNumModeShapes = 2;
	else if (NumSD() == 3)
		fNumModeShapes = (fIncompressibleMode) ? 4 : 3;
	else throw eGeneralFail;
}

/* destructor */
SimoFiniteStrainT::~SimoFiniteStrainT(void)
{
	/* free shape functions */
	delete fEnhancedShapes;
	fShapes = NULL; /* already gone */
}

/* data initialization */
void SimoFiniteStrainT::Initialize(void)
{
	/* inherited */
	FiniteStrainT::Initialize();

//TEMP
	if (fMaterialList->HasThermalStrains())
	{
		cout << "\n SimoFiniteStrainT::Initialize: element group ";
		cout << fFEManager.ElementGroupNumber(this) + 1 << ": total Lagrangian\n";
		cout << "     formulation has a bug associated with the multiplicative treatment\n";
		cout << "     of thermal strains. The solutions appear correct, but the convergence\n";
		cout << "     is not quadratic. Use updated Lagrangian formulation." << endl;
		throw eGeneralFail;
	}	

	/* space for enhanced part of the deformation gradient */
	if (Needs_F())
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_enh_all.Allocate(nip*nsd*nsd);
		fF_Galerkin_all.Allocate(nip*nsd*nsd);
		fF_enh_List.Allocate(NumIP());
		fF_Galerkin_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
		{
			int dex = i*nsd*nsd;
			fF_enh_List[i].Set(nsd, nsd, fF_enh_all.Pointer(dex));
			fF_Galerkin_List[i].Set(nsd, nsd, fF_Galerkin_all.Pointer(dex));
		}
	}
	
	/* space for enhanced part of the "last" deformation gradient */
	if (Needs_F_last())
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_enh_last_all.Allocate(nip*nsd*nsd);
		fF_Galerkin_last_all.Allocate(nip*nsd*nsd);
		fF_enh_last_List.Allocate(NumIP());
		fF_Galerkin_last_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
		{
			int dex = i*nsd*nsd;
			fF_enh_last_List[i].Set(nsd, nsd, fF_enh_last_all.Pointer(dex));
			fF_Galerkin_last_List[i].Set(nsd, nsd, fF_Galerkin_last_all.Pointer(dex));
		}
	}

	/* dimension */
	fGradNa.Allocate(fNumSD, fNumElemNodes);
	fStressStiff.Allocate(fNumElemNodes);
	fTemp2.Allocate(fNumElemNodes*fNumDOF);
}

/* finalize current step - step is solved */
void SimoFiniteStrainT::CloseStep(void)
{
	/* inherited */
	FiniteStrainT::CloseStep();
	
	/* store converged solution */
	fCurrElementModes_last = fCurrElementModes;
}
	
/* restore last converged state */
void SimoFiniteStrainT::ResetStep(void)
{
	/* inherited */
	FiniteStrainT::ResetStep();
	
	/* store converged solution */
	fCurrElementModes = fCurrElementModes_last;
}

/* read restart information from stream */
void SimoFiniteStrainT::ReadRestart(istream& in)
{
	/* inherited */
	FiniteStrainT::ReadRestart(in);
	
	/* read restart data */
	in >> fCurrElementModes;
	
	/* reset last state */
	fCurrElementModes_last = fCurrElementModes;
}

/* write restart information from stream */
void SimoFiniteStrainT::WriteRestart(ostream& out) const
{
	/* inherited */
	FiniteStrainT::WriteRestart(out);
	
	/* read restart data */
	out << fCurrElementModes;
}

/***********************************************************************
* Protected
***********************************************************************/

/* write element parameter to out */
void SimoFiniteStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	FiniteStrainT::PrintControlData(out);
	
	/* parameters */
	out << " Include incompressible mode . . . . . . . . . . = " << fIncompressibleMode << '\n';
	out << "    Number of enhanced mode shapes = " << fNumModeShapes << '\n';
	out << " Maximum number of local sub-iterations. . . . . = " << fLocalIterationMax << '\n';
	out << " Absolute tol. on residual of enhanced modes . . = " << fAbsTol << '\n';
	out << " Maximum number of local sub-iterations. . . . . = " << fRelTol << '\n';
}

/* construct shape function */
void SimoFiniteStrainT::SetShape(void)
{
	/* allocate memory for incompressible modes */
	fElementModes.Allocate(NumElements(), fNumModeShapes*NumSD());
	fElementModes = 0;
	
	/* dimension before sending to the shape functions */
	fCurrElementModes.Set(fNumModeShapes, NumSD(), fElementModes.Pointer());

	/* construct shape functions */
	fEnhancedShapes = new SimoShapeFunctionT(fGeometryCode, fNumIP,
		fLocInitCoords, fCurrElementModes);
	if (!fEnhancedShapes) throw eOutOfMemory;

	/* modes from the last time step */
	fElementModes_last = fElementModes;
	fCurrElementModes_last.Set(fNumModeShapes, NumSD(), fElementModes_last.Pointer());	

	/* initialize */
	fEnhancedShapes->Initialize();
	
	/* set base class pointer */
	fShapes = fEnhancedShapes;
}

/* form shape functions and derivatives */
void SimoFiniteStrainT::SetGlobalShape(void)
{
	/* set modes for the current element */
	int current_element = CurrElementNumber();
	fCurrElementModes.Set(fNumModeShapes, NumSD(), fElementModes(current_element));
	fCurrElementModes_last.Set(fNumModeShapes, NumSD(), fElementModes_last(current_element));

	/* 3D with incompressible mode */
	if (fIncompressibleMode)
		ModifiedEnhancedDeformation();
	else /* modification is strictly additive */
	{
		/* what needs to get computed */
		bool needs_F = Needs_F();
		bool needs_F_last = Needs_F_last();

		/* inherited - set Galerkin part of deformation gradient */
		FiniteStrainT::SetGlobalShape();
		
		/* store Galerkin parts of F and F_last */
		if (needs_F) fF_Galerkin_all = fF_all;
		if (needs_F_last) fF_Galerkin_last_all = fF_last_all;

		/* compute enhanced part of F and total F */
		ComputeEnhancedDeformation(needs_F, needs_F_last);
		
		/* calculate the residual from the internal force */
		FormKd_enhanced(fPK1_list, fRHS_enh);
		
		/* internal iterations */
		double res, res_0, res_rel;
		res = fRHS_enh.Magnitude();
		res_0 = res;
		res_rel = 1.0;
		int iter_enh = 0;
//TEMP
cout << "\n SimoFiniteStrainT::SetGlobalShape: solve internal modes\n";
cout << setw(kIntWidth) << "iter" << setw(kDoubleWidth) << "error" << '\n';
		while (iter_enh++ < fLocalIterationMax && res > fAbsTol && res_rel > fRelTol)
		{
//TEMP
cout << setw(kIntWidth) << iter_enh << setw(kDoubleWidth) << res << '\n';

			/* form the stiffness associated with the enhanced modes */
			FormStiffness_enhanced(fK22);

			/* update enhanced modes */
			fK22.LinearSolve(fRHS_enh);			
			
			/* recompute shape functions */
			fCurrElementModes -= fRHS_enh;
		
			/* compute enhanced part of F and total F */
			ComputeEnhancedDeformation(needs_F, needs_F_last);
		
			/* calculate the residual from the internal force */
			FormKd_enhanced(fPK1_list, fRHS_enh);
		
			/* errors */
			res = fRHS_enh.Magnitude();
			res_rel = res/res_0;
		}
//TEMP
cout << setw(kIntWidth) << iter_enh << setw(kDoubleWidth) << res << '\n';
	}
}

/* form the element stiffness matrix */
void SimoFiniteStrainT::FormStiffness(double constK)
{		
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
	/* S T R E S S   S T I F F N E S S */
				
		/* PK2 stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
	
		/* change rule shape function derivatives */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();
		fTempMat1.Inverse();
		fShapes->SetChainRule(fTempMat1, fDNa_x);

		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x, fGradNa);

		/* scale factor */
		double scale = constK*(*Det++)*(*Weight++)*J;

		/* integration constants */		
		fStressMat *= scale;
	
		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fStressMat, format,
			dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
	
		/* strain displacement matrix */
		fShapes->B(fDNa_x, fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, fNumDOF);
}

/* calculate the internal force contribution ("-k*d") */
void SimoFiniteStrainT::FormKd(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT WP(fNumSD, fStressStiff.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
//TEMP: compute 1st PK stress based on Cauchy stress since most materials
//      have not been optimized to compute PK2 directly.	
	
		/* get Cauchy stress */
		//(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F^(-1) */
		//fTempMat2 = DeformationGradient();
		double J = DeformationGradient().Det();
		if (J <= 0.0)
		{
			cout << "\n SimoFiniteStrainT::FormKd: negative jacobian determinant" << endl;
			throw eBadJacobianDet;
		}
		//else
		//	fTempMat2.Inverse();

		/* compute PK1 */
		//fStressMat.MultABT(fTempMat1, fTempMat2);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		WP.MultAB(PK1_list[CurrIP()];, fGradNa);

		/* accumulate */
		fRHS.AddScaled(J*constK*(*Weight++)*(*Det++), fNEEvec);
		
		//need to add contribution from (4.20) using the update vector
		//calculated by the last pass through the iterations -> make
		//a separate data array
	}	
}

/***********************************************************************
* Private
***********************************************************************/

/* compute modified, enhanced deformation gradient */
void SimoFiniteStrainT::ModifiedEnhancedDeformation(void)
{
	/* check */
	if (NumSD() != 3)
	{
		cout << "\n SimoFiniteStrainT::ModifiedEnhancedDeformation: for 3D only" << endl;
		throw eGeneralFail;
	}
	
	/* skip base class implementation because the deformation gradient
	 * modification are not simply additive */
	ElasticT::SetGlobalShape();
	
	cout << "\n SimoFiniteStrainT::ModifiedEnhancedDeformation: not done" << endl;
	throw;
}

/* compute enhanced part of F and total F */
void SimoFiniteStrainT::ComputeEnhancedDeformation(bool need_F, bool need_F_last)
{
	/* store Galerkin part/compute enhancement and total F */
	if (need_F)
	{
		/* compute enhanced part of the deformation gradient */
		for (int i = 0; i < NumIP(); i++)
			fEnhancedShapes->GradU_enhanced(fCurrElementModes, fF_enh_List[i], i);
			
		/* compute total, enhanced deformation gradient */
		fF_all.SumOf(fF_Galerkin_all, fF_enh_all);
	}		

	/* store Galerkin part/compute enhancement and total F from last step */
	if (need_F_last)
	{
		/* compute enhanced part of the deformation gradient */
		for (int i = 0; i < NumIP(); i++)
			fEnhancedShapes->GradU_enhanced(fCurrElementModes_last, fF_enh_last_List[i], i);
			
		/* compute total, enhanced deformation gradient */
		fF_last_all.SumOf(fF_Galerkin_last_all, fF_enh_last_all);
	}		
}

/* calculate the residual from the internal force */
void SimoFiniteStrainT::FormKd_enhanced(ArrayT<dMatrixT>& PK1_list, dArrayT& RHS_enh)
{
	/* integration rule */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* integrate over element */
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F^(-1) */
		fTempMat2 = DeformationGradient();
		double J = fTempMat2.Det();
		if (J <= 0.0)
		{
			cout << "\n SimoFiniteStrainT::FormKd_enhanced: negative jacobian determinant" << endl;
			throw eBadJacobianDet;
		}
		else
			fTempMat2.Inverse();

		/* compute PK1 */
		dMatrixT& PK1 = PK1_list[CurrIP()];
		PK1.MultABT(fTempMat1, fTempMat2);

		/* enhanced shape function gradients */
		fEnhancedShapes->GradNa_enhanced(fGradNa_enh);

		/* Wi,J PiJ */
		fWP_enh.MultAB(PK1, fGradNa_enh);

		/* accumulate */
		RHS_enh.AddScaled(J*(*Weight++)*(*Det++), fWP_enh);
	}
}

/* form the stiffness associated with the enhanced modes */
void SimoFiniteStrainT::FormStiffness_enhanced(dMatrixT& K_22)
{
	K_22 = 0.0;
}
