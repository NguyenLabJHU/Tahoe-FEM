/* $Id: SimoFiniteStrainT.cpp,v 1.9 2001-09-04 06:54:08 paklein Exp $ */
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

/* debugging flag */
//#define _SIMO_FINITE_STRAIN_DEBUG_
#undef _SIMO_FINITE_STRAIN_DEBUG_

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

	int solver_method = -1;
	in >> solver_method;
	if (solver_method == kStaticCondensation)
		fModeSolveMethod = kStaticCondensation;
	else if (solver_method == kLocalIteration)
		fModeSolveMethod = kLocalIteration;
	else
		throw eBadInputValue;
		
	/* parameters for local iteration */
	if (fModeSolveMethod == kLocalIteration)
	{
		in >> fLocalIterationMax;
		in >> fAbsTol;
		in >> fRelTol;
	}
	else
	{
		fLocalIterationMax = 1;
		fAbsTol = fRelTol = 0.5; //dummy values
	}

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

	/* dimensions */
	int nip = NumIP();
	int nsd = NumSD();
	int nst = dSymMatrixT::NumValues(nsd);

	/* allocate storage for stress and modulus at the integration
	 * points of all the elements */
	fPK1_storage.Allocate(NumElements(), nip*nsd*nsd);
	fc_ijkl_storage.Allocate(NumElements(), nip*nst*nst);
	fPK1_list.Allocate(nip);
	fc_ijkl_list.Allocate(nip);

	/* what's needed */
	bool need_F = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		need_F = need_F || Needs_F(i);		
		need_F_last = need_F_last || Needs_F_last(i);
	}	

	/* space for enhanced part of the deformation gradient */
	if (need_F)
	{
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
	if (need_F_last)
	{
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

	/* dimension work space */
	fGradNa.Allocate(fNumSD, fNumElemNodes);
	fGradNa_enh.Allocate(fNumSD, fNumModeShapes);
	fRHS_enh.Allocate(fNumSD*fNumModeShapes);
	fB_enh.Allocate(dSymMatrixT::NumValues(fNumSD), fNumSD*fNumModeShapes);
	fWP_enh.Allocate(fNumSD, fNumModeShapes);
	
	/* stiffness work space */
	fStressStiff_11.Allocate(fNumElemNodes);
	fStressStiff_12.Allocate(fNumElemNodes, fNumModeShapes);
	fStressStiff_22.Allocate(fNumModeShapes, fNumModeShapes);
	
	fK22.Allocate(fNumModeShapes*fNumDOF);	
	fK12.Allocate(fNumElemNodes*fNumDOF, fNumModeShapes*fNumDOF);
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

/* increment current element */
bool SimoFiniteStrainT::NextElement(void)
{
	/* inherited */
	bool next = FiniteStrainT::NextElement();
	
	/* set references to element arrays */
	if (next)
	{
		/* get pointer to current element data */
		int element = CurrElementNumber();
		double* s = fPK1_storage(element);
		double* c = fc_ijkl_storage(element);
	
		/* dimensions */
		int nip = NumIP();
		int nsd = NumSD();
		int nst = dSymMatrixT::NumValues(nsd);
	
		/* set references */
		int nsd_2 = nsd*nsd;
		int nst_2 = nst*nst;
		for (int i = 0; i < nip; i++)
		{
			fPK1_list[i].Set(nsd, nsd, s);
			fc_ijkl_list[i].Set(nst, nst, c);
			s += nsd_2;
			c += nst_2;
		}
	}
	return next;
}

/* write element parameter to out */
void SimoFiniteStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	FiniteStrainT::PrintControlData(out);
	
	/* parameters */
	out << " Include incompressible mode . . . . . . . . . . = " << fIncompressibleMode << '\n';
	out << "    Number of enhanced mode shapes = " << fNumModeShapes << '\n';	
	out << " Solution method for enhanced element modes. . . = " << fModeSolveMethod << '\n';
	out << "    eq." << kStaticCondensation << ", static condensation\n";
	out << "    eq." << kLocalIteration	    << ", staggered, local iteration\n";
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
	if (fShapes != NULL)
	{
		cout << "\n SimoFiniteStrainT::SetShape: deleting non-NULL fShapes" << endl;
		delete fShapes;
	}
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
		int material_number = CurrentElement().MaterialNumber();
		bool needs_F = Needs_F(material_number);
		bool needs_F_last = Needs_F_last(material_number);

		/* inherited - set Galerkin part of deformation gradient */
		FiniteStrainT::SetGlobalShape();

#ifdef _SIMO_FINITE_STRAIN_DEBUG_
		cout << "fLocInitCoords:\n" << fLocInitCoords << '\n';
		cout << "fLocDisp:\n" << fLocDisp << endl;
#endif
		
		/* store Galerkin parts of F and F_last */
		if (needs_F) fF_Galerkin_all = fF_all;
		if (needs_F_last) fF_Galerkin_last_all = fF_last_all;

		/* compute enhanced part of F and total F */
		ComputeEnhancedDeformation(needs_F, needs_F_last);
		
		/* calculate the residual from the internal force */
		if (RunState() == GlobalT::kFormRHS)
		{
			fRHS_enh = 0.0;
			FormKd_enhanced(fPK1_list, fRHS_enh);
		
			/* internal iterations */
			double res, res_0, res_rel;
			res = fRHS_enh.Magnitude();
			res_0 = res;
			res_rel = 1.0;
			int iter_enh = 0;

//#ifdef _SIMO_FINITE_STRAIN_DEBUG_
#if 1
			cout << "\n SimoFiniteStrainT::SetGlobalShape: solve internal modes\n";
			cout << setw(kIntWidth) << "iter" << setw(kDoubleWidth) << "error" << '\n';
#endif

			while (iter_enh++ < fLocalIterationMax && res > fAbsTol && res_rel > fRelTol)
			{
//#ifdef _SIMO_FINITE_STRAIN_DEBUG_
#if 1
				cout << setw(kIntWidth) << iter_enh << setw(kDoubleWidth) << res << '\n';
#endif

				/* form the stiffness associated with the enhanced modes */
				fK22 = 0.0;
				FormStiffness_enhanced(fK22, NULL);

#ifdef _SIMO_FINITE_STRAIN_DEBUG_
				cout << "\n F_int =\n" << fRHS_enh << endl;
				cout << "\n K22 =\n"   << fK22 << endl;
#endif

				/* update enhanced modes */
				fK22.LinearSolve(fRHS_enh);			

#ifdef _SIMO_FINITE_STRAIN_DEBUG_
				cout << "\n d_mode =\n" << fRHS_enh << endl;
#endif	
			
				/* recompute shape functions */
				fCurrElementModes.AddScaledTranspose(-1.0, fRHS_enh);
				//fCurrElementModes -= fRHS_enh;
		
				/* compute enhanced part of F and total F */
				ComputeEnhancedDeformation(needs_F, needs_F_last);
		
				/* calculate the residual from the internal force */
				fRHS_enh = 0.0;
				FormKd_enhanced(fPK1_list, fRHS_enh);
		
				/* errors */
				res = fRHS_enh.Magnitude();
				res_rel = res/res_0;
			}
//#ifdef _SIMO_FINITE_STRAIN_DEBUG_
#if 1
			cout << setw(kIntWidth) << iter_enh << setw(kDoubleWidth) << res << '\n';
#endif
		
			/* form the stiffness associated with the enhanced modes */
			fK22 = 0.0;
			fK12 = 0.0;
			FormStiffness_enhanced(fK22, &fK12); /* also stores c_ijkl */
		}
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
	fStressStiff_11 = 0.0;
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* deformation gradient */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();

		/* scale factor */
		double scale = constK*(*Det++)*(*Weight++)*J;

	/* S T R E S S   S T I F F N E S S */
		
		/* Cauchy stress */
		const dMatrixT& PK1 = fPK1_list[CurrIP()]; /* retrieve */
		fStressMat.MultABT(PK1, fTempMat1);
		fStressMat *= scale/J;
		
		/* PK2 stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
	
		/* chain rule shape function derivatives */
		fTempMat1.Inverse();
		fShapes->TransformDerivatives(fTempMat1, fDNa_x);

		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x, fGradNa);

		/* integration constants */		
		fStressMat *= scale;
	
		/* using the stress symmetry */
		fStressStiff_11.MultQTBQ(fGradNa, fStressMat, format, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
	
		/* strain displacement matrix */
		fShapes->B(fDNa_x, fB);

		/* get D matrix */
		fD.SetToScaled(scale, fc_ijkl_list[CurrIP()]);
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff_11, fNumDOF);
}

/* calculate the internal force contribution ("-k*d") */
void SimoFiniteStrainT::FormKd(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT WP(fNumSD, fStressStiff_11.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		WP.MultAB(fPK1_list[CurrIP()], fGradNa);

		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
	}
	
	/* apply correction for unequilibrated element modes */
	bool correct_force = false; //TEMP
	cout << "\n SimoFiniteStrainT::FormKd: SKIPPING correction for unequilibrated element modes" << endl;
	if (correct_force)
	{
		/* "project" to nodal force using (4.20) */
		fK22.LinearSolve(fRHS_enh);
		fK12.Multx(fRHS_enh, fNEEvec);
		fRHS.AddScaled(-constK, fNEEvec);
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
	throw eGeneralFail;
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

#ifdef _SIMO_FINITE_STRAIN_DEBUG_
		cout << "\n SimoFiniteStrainT::ComputeEnhancedDeformation:\n"
		     << " element modes:\n" << fCurrElementModes << endl;
		for (int i = 0; i < NumIP(); i++)
		{
			cout << " ip: " << i << '\n';
			cout << " F_Galerkin:\n" << fF_Galerkin_List[i] << '\n';
			cout << " F_enh:\n" << fF_enh_List[i] << '\n';
			cout << " F:\n" << fF_List[i] << '\n';
		}
		cout.flush();
#endif
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
#ifdef _SIMO_FINITE_STRAIN_DEBUG_
		cout << "\n SimoFiniteStrainT::FormKd_enhanced: ip: " << fShapes->CurrIP() << endl;
#endif

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
		PK1 *= J;

		/* enhanced shape function gradients */
		fEnhancedShapes->GradNa_enhanced(fGradNa_enh);

#ifdef _SIMO_FINITE_STRAIN_DEBUG_
		cout << "Cauchy stress:\n" << fTempMat1 << '\n';
		cout << "PK1:\n" << PK1 << '\n';
		cout << "GRAD[Na_enh]:\n" << fGradNa_enh << endl;
#endif

		/* Wi,J PiJ */
		fWP_enh.MultAB(PK1, fGradNa_enh);

		/* accumulate */
		RHS_enh.AddScaled((*Weight++)*(*Det++), fWP_enh);
		
#ifdef _SIMO_FINITE_STRAIN_DEBUG_
		/* gradients in the current config */
		fEnhancedShapes->TransformDerivatives_enhanced(fTempMat2, fDNa_x_enh);
		fShapes->GradNa(fDNa_x_enh, fGradNa_enh);
		cout << "grad[Na_enh]:\n" << fGradNa_enh << endl;
#endif
	}
}

/* form the stiffness associated with the enhanced modes */
void SimoFiniteStrainT::FormStiffness_enhanced(dMatrixT& K_22, dMatrixT* K_12)
{
	/* integration */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* initialize */
	fStressStiff_22 = 0.0;
	if (K_12) fStressStiff_12 = 0.0;
	
	fEnhancedShapes->TopIP();
	while (fEnhancedShapes->NextIP())
	{
		/* deformation gradient */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();

		/* integration weight (current configuration) */
		double scale = (*Det++)*(*Weight++)*J;

		/* Cauchy stress */
		const dMatrixT& PK1 = fPK1_list[CurrIP()]; /* retrieve */
		fStressMat.MultABT(PK1, fTempMat1);
		fStressMat *= scale/J;

		/* material tangent modulus- calculate and store */
		fc_ijkl_list[CurrIP()] = fCurrMaterial->c_ijkl();

		/* transform shape functions to derivatives wrt. current coords */
		fTempMat1.Inverse(); /* F^-1 */
		fEnhancedShapes->TransformDerivatives_enhanced(fTempMat1, fDNa_x_enh);

	/* S T R E S S   S T I F F N E S S */
	
		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x_enh, fGradNa_enh);
	
		/* using the stress symmetry */
		fStressStiff_22.MultQTBQ(fGradNa_enh, fStressMat, dMatrixT::kWhole, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
	
		/* strain displacement matrix */
		fShapes->B(fDNa_x_enh, fB_enh);

		/* get D matrix */
		fD.SetToScaled(scale, fc_ijkl_list[CurrIP()]);
						
		/* accumulate */
		K_22.MultQTBQ(fB_enh, fD, dMatrixT::kWhole, dMatrixT::kAccumulate);
		
		/* mixed shape function part */
		if (K_12)
		{
			/* transform Galerkin shape function derivatives */
			fShapes->TransformDerivatives(fTempMat1, fDNa_x);
		
	/* S T R E S S   S T I F F N E S S */
	
			/* get shape function gradients matrix */
			fShapes->GradNa(fDNa_x, fGradNa);
	
			/* using the stress symmetry */
			fStressStiff_12.MultATBC(fGradNa, fStressMat, fGradNa_enh, 
				dMatrixT::kWhole, dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
	
			/* strain displacement matrix */
			fShapes->B(fDNa_x, fB);

			/* accumulate */
			K_12->MultATBC(fB, fD, fB_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);
		}
	}
						
	/* expand and add in stress stiffness parts */
	K_22.Expand(fStressStiff_22, fNumDOF);
	if (K_12) K_12->Expand(fStressStiff_12, fNumDOF);
}
