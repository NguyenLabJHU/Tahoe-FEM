/* $Id: SimoFiniteStrainT.cpp,v 1.30 2004-07-15 08:26:27 paklein Exp $ */
#include "SimoFiniteStrainT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "SolidMaterialT.h"
#include "MaterialListT.h" //NOTE - only needed for check in Initialize?
#include "SimoShapeFunctionT.h"
#include "FieldT.h"
#include "XDOF_ManagerT.h"

using namespace Tahoe;

/* constructor */
SimoFiniteStrainT::SimoFiniteStrainT(const ElementSupportT& support, const FieldT& field):
	FiniteStrainT(support),
	fCurrElementModes(LocalArrayT::kUnspecified),
	fCurrElementModes_last(LocalArrayT::kUnspecified),
	fEnhancedShapes(NULL),
	fStressMat(NumSD()),
	fTempMat1(NumSD()),
	fTempMat2(NumSD())
{
ExceptionT::GeneralFail("SimoFiniteStrainT::SimoFiniteStrainT", "out of date");
#if 0
	/* read parameters */
	ifstreamT& in = ElementSupport().Input();

	int inc_mode = -1;
	in >> inc_mode;

	int solver_method = -1;
	in >> solver_method;
	if (solver_method == kStaticCondensation)
	{
		fModeSolveMethod = kStaticCondensation;
		
		//TEMP
		cout << "\n SimoFiniteStrainT::SimoFiniteStrainT: solution method not implemented: " 
		     << fModeSolveMethod << endl;
		throw ExceptionT::kBadInputValue;
	}
	else if (solver_method == kLocalIteration)
		fModeSolveMethod = kLocalIteration;
	else if (solver_method == kMonolithic)
		fModeSolveMethod = kMonolithic;
	else
	{
		cout << "\n SimoFiniteStrainT::SimoFiniteStrainT: unrecognized method to solve\n"
		     <<   "     for the enhanced element modes: " << solver_method << endl;
		throw ExceptionT::kBadInputValue;
	}
		
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
	if (inc_mode != 0 && inc_mode != 1) throw ExceptionT::kBadInputValue;
	if (fLocalIterationMax < 1) throw ExceptionT::kGeneralFail;
	fIncompressibleMode = (inc_mode == 1);
	if (NumSD() == 2) fIncompressibleMode = false;
	if (fAbsTol < 0 || fAbsTol > 1.0) throw ExceptionT::kBadInputValue;
	if (fRelTol < 0 || fRelTol > 1.0) throw ExceptionT::kBadInputValue;
	
	/* set number of mode shapes */
	if (NumSD() == 2)
		fNumModeShapes = 2;
	else if (NumSD() == 3)
		fNumModeShapes = (fIncompressibleMode) ? 4 : 3;
	else throw ExceptionT::kGeneralFail;
#endif
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
//	FiniteStrainT::Initialize();

	/* dimensions */
	int nip = NumIP();
	int nsd = NumSD();
	int nst = dSymMatrixT::NumValues(nsd);

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
		fF_enh_all.Dimension(nip*nsd*nsd);
		fF_Galerkin_all.Dimension(nip*nsd*nsd);
		fF_enh_List.Dimension(NumIP());
		fF_Galerkin_List.Dimension(NumIP());
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
		fF_enh_last_all.Dimension(nip*nsd*nsd);
		fF_Galerkin_last_all.Dimension(nip*nsd*nsd);
		fF_enh_last_List.Dimension(NumIP());
		fF_Galerkin_last_List.Dimension(NumIP());
		for (int i = 0; i < NumIP(); i++)
		{
			int dex = i*nsd*nsd;
			fF_enh_last_List[i].Set(nsd, nsd, fF_enh_last_all.Pointer(dex));
			fF_Galerkin_last_List[i].Set(nsd, nsd, fF_Galerkin_last_all.Pointer(dex));
		}
	}

	/* dimension work space */
	int nen = NumElementNodes();
	fGradNa.Dimension(NumSD(), nen);
	fGradNa_enh.Dimension(NumSD(), fNumModeShapes);
	fRHS_enh.Dimension(NumSD()*fNumModeShapes);
	fB_enh.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*fNumModeShapes);
	fWP_enh.Dimension(NumSD(), fNumModeShapes);
	
	/* stiffness work space */
	fStressStiff_11.Dimension(nen);
	fStressStiff_12.Dimension(nen, fNumModeShapes);
	fStressStiff_21.Dimension(fNumModeShapes, nen);
	fStressStiff_22.Dimension(fNumModeShapes, fNumModeShapes);
	
	fK22.Dimension(fNumModeShapes*NumDOF());	
	fK12.Dimension(nen*NumDOF(), fNumModeShapes*NumDOF());
	fK21.Dimension(fNumModeShapes*NumDOF(), nen*NumDOF());
	
	/* solve all dof's together */
	if (fModeSolveMethod == kMonolithic)
	{
		/* work space for UL quad */
		fK11.Dimension(fLHS.Rows());
	
		/* resize work arrays */
		fRHS.Dimension(nen*NumDOF() + fCurrElementModes.Length());
		fLHS.Dimension(fRHS.Length());
		
		/* XDOF support */
		XDOF_ManagerT& xdof_man = ElementSupport().XDOF_Manager();		

		/* register as XDOF group */
		iArrayT set_dimensions(1);
		set_dimensions[0] = fNumModeShapes*NumSD();
		xdof_man.XDOF_Register(this, set_dimensions);
		
		/* element nodes use space from XDOF manager */
		const dArray2DT& xdof = xdof_man.XDOF(this, 0);
		fElementModes.Alias(xdof);
		
		/* checks */
		if (fElementModes.MajorDim() != NumElements() || 
		    fElementModes.MinorDim() != fNumModeShapes*NumSD())
		{
			cout << "\n SimoFiniteStrainT::Initialize: element modes array is not\n" 
			     <<   "     the correct dimension" << endl;
			throw ExceptionT::kGeneralFail;
		}		
	}
	else /* all other solution methods split dof's */
	{
		/* allocate storage for stress and modulus at the integration
		 * points of all the elements */
		fPK1_storage.Dimension(NumElements(), nip*nsd*nsd);
		fc_ijkl_storage.Dimension(NumElements(), nip*nst*nst);
		fPK1_list.Dimension(nip);
		fc_ijkl_list.Dimension(nip);

		/* allocate memory for incompressible modes */
		fElementModes.Dimension(NumElements(), fNumModeShapes*NumSD());
		fElementModes = 0;
		
		/* element modes increment */
		//fElementModesInc.Dimension(fCurrElementModes.Length());
	}

	/* allocate memory for last incompressible modes */
	fElementModes_last.Dimension(NumElements(), fNumModeShapes*NumSD());
	fElementModes_last = 0;
}

/* finalize current step - step is solved */
void SimoFiniteStrainT::CloseStep(void)
{
	/* inherited */
	FiniteStrainT::CloseStep();
	
	/* store converged solution */
	fElementModes_last = fElementModes;
}
	
/* restore last converged state */
GlobalT::RelaxCodeT SimoFiniteStrainT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = FiniteStrainT::ResetStep();
	
	/* store converged solution */
	fElementModes = fElementModes_last;

	return relax;
}

/* read restart information from stream */
void SimoFiniteStrainT::ReadRestart(istream& in)
{
	/* inherited */
	FiniteStrainT::ReadRestart(in);
	
	/* read restart data */
	in >> fElementModes;
	
	/* reset last state */
	fElementModes_last = fElementModes;
}

/* write restart information from stream */
void SimoFiniteStrainT::WriteRestart(ostream& out) const
{
	/* inherited */
	FiniteStrainT::WriteRestart(out);
	
	/* read restart data */
	out << fElementModes;
}

/* return field connectivities. */
void SimoFiniteStrainT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	FiniteStrainT::ConnectsU(connects_1, connects_2);

	/* send connectivities with enhanced mode tags */
	if (fModeSolveMethod == kMonolithic)
		connects_1.AppendUnique(&fEnhancedConnectivities);
}

/* append element equations numbers to the list. */
void SimoFiniteStrainT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	if (fModeSolveMethod != kMonolithic)
		FiniteStrainT::Equations(eq_1, eq_2);
	else
	{
		/* equations associated with the element modes */
		const iArray2DT& xdof_eqnos_all = ElementSupport().XDOF_Manager().XDOF_Eqnos(this, 0);

		/* loop over connectivity blocks */
		int element_count = 0;
		for (int i = 0; i < fConnectivities.Length(); i++)
		{
			/* block connectivities and equations */
			const iArray2DT& connects = *fConnectivities[i];
			iArray2DT& eqnos = fEqnos[i];
		
			/* resize equations array - displacement DOF's + element modes */
			int neq = NumElementNodes()*NumDOF();
			eqnos.Dimension(connects.MajorDim(), neq + fCurrElementModes.Length());
		
			/* set displacement equations */
			Field().SetLocalEqnos(connects, eqnos);

			/* xdof equations for the block */
			const iArray2DT xdof_eqnos(connects.MajorDim(), xdof_eqnos_all.MinorDim(), xdof_eqnos_all(element_count));
			element_count += connects.MajorDim();

			/* fill columns with xdof equations */
			iArrayT tmp(xdof_eqnos.MajorDim());
			for (int i = neq; i < eqnos.MinorDim(); i++)
			{
				xdof_eqnos.ColumnCopy(i - neq, tmp);
				eqnos.SetColumn(i, tmp);
			}

			/* add to list */
			eq_1.AppendUnique(&eqnos);
		}
		
		/* reset element cards */
		SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
	}
}

/* determine number of tags needed. See DOFElementT. */
void SimoFiniteStrainT::SetDOFTags(void)
{
	/* one tag per element */
	fEnhancedModeTags.Dimension(NumElements());
}
	
/* return the array tag numbers in the specified set currently 
 * used/need by the group. See DOFElementT. */
iArrayT& SimoFiniteStrainT::DOFTags(int tag_set)
{
	/* check */
	if (tag_set != 0)
	{
		cout << "\n SimoFiniteStrainT::DOFTags: expecting tag set 0: " 
		     << tag_set << endl;
		throw ExceptionT::kOutOfRange;
	}

	/* return the tags array */
	return fEnhancedModeTags;
}

/* generate nodal connectivities. See DOFElementT. */
void SimoFiniteStrainT::GenerateElementData(void)
{
	/* just link tags with (one of) the element nodes */
	fEnhancedConnectivities.Dimension(NumElements(), 2);
	
	int *pt = fEnhancedConnectivities.Pointer();
	int *mt = fEnhancedModeTags.Pointer();
	for (int i=0; i < NumElements(); i++)
	  {
	    const iArrayT& elemnodes = fElementCards[i].NodesX();
	    *pt++ = elemnodes[0];
	    *pt++ = *mt++;
	  }
}

/* return the connectivities associated with the element generated
 * degrees of freedom. See DOFElementT. */
const iArray2DT& SimoFiniteStrainT::DOFConnects(int tag_set) const
{
	/* check */
	if (tag_set != 0)
	{
		cout << "\n SimoFiniteStrainT::DOFConnects: expecting tag set 0: " 
		     << tag_set << endl;
		throw ExceptionT::kOutOfRange;
	}

	/* partial info */
	return fEnhancedConnectivities;
}

/* restore the element degrees of freedom. See DOFElementT. */
void SimoFiniteStrainT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
	/* check */
	if (tag_set != 0)
	{
		cout << "\n SimoFiniteStrainT::ResetDOF: expecting tag set 0: " 
		     << tag_set << endl;
		throw ExceptionT::kOutOfRange;
	}

	/* restore last solution */
	DOF = fElementModes_last;
}

/* the equation group */
int SimoFiniteStrainT::Group(void) const { return FiniteStrainT::Group(); }

/***********************************************************************
* Protected
***********************************************************************/

/* increment current element */
bool SimoFiniteStrainT::NextElement(void)
{
	/* inherited */
	bool next = FiniteStrainT::NextElement();
	
	/* set references to element arrays */
	if (fModeSolveMethod != kMonolithic && next)
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

#if 0
/* write element parameter to out */
void SimoFiniteStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	FiniteStrainT::PrintControlData(out);
	
	/* parameters */
	out << " Include incompressible mode . . . . . . . . . . = " << fIncompressibleMode << '\n';
	out << "    Number of enhanced mode shapes = " << fNumModeShapes << '\n';	
	out << " Solution method for enhanced element modes. . . = " << fModeSolveMethod << '\n';
	out << "    eq." << kMonolithic         << ", monolithic\n";
	out << "    eq." << kStaticCondensation << ", static condensation\n";
	out << "    eq." << kLocalIteration	    << ", staggered, local iteration\n";

	/* solver-dependent output */
	if (fModeSolveMethod == kLocalIteration)
	{
		out << " Maximum number of local sub-iterations. . . . . = " << fLocalIterationMax << '\n';
		out << " Absolute tol. on residual of enhanced modes . . = " << fAbsTol << '\n';
		out << " Maximum number of local sub-iterations. . . . . = " << fRelTol << '\n';
	}
}
#endif

/* construct shape function */
void SimoFiniteStrainT::SetShape(void)
{
	/* dimension before sending to the shape functions */
	fCurrElementModes.Dimension(fNumModeShapes, NumSD());
	fCurrElementModes = 0;
	fCurrElementModes_last = fCurrElementModes;

	/* construct shape functions */
	fEnhancedShapes = new SimoShapeFunctionT(GeometryCode(), NumIP(),
		fLocInitCoords, fCurrElementModes);
	if (!fEnhancedShapes) throw ExceptionT::kOutOfMemory;

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
	/* 3D with incompressible mode */
	if (fIncompressibleMode)
		ModifiedEnhancedDeformation();
	else /* modification is strictly additive */
	{
		/* what needs to get computed */
		int current_element_number = CurrElementNumber();
		int material_number = CurrentElement().MaterialNumber();
		bool needs_F = Needs_F(material_number);
		bool needs_F_last = Needs_F_last(material_number);

		/* inherited - set Galerkin part of deformation gradient */
		FiniteStrainT::SetGlobalShape();

		/* store Galerkin parts of F and F_last */
		if (needs_F) fF_Galerkin_all = fF_all;
		if (needs_F_last) fF_Galerkin_last_all = fF_last_all;

		/* compute enhanced part of F and total F */		
		ComputeEnhancedDeformation(needs_F, needs_F_last);
		
		/* calculate the residual from the internal force */
		if (fModeSolveMethod == kLocalIteration && RunState() == GlobalT::kFormRHS)
		{
			/* residual of element modes */
			fRHS_enh = 0.0;
			FormKd_enhanced(fPK1_list, fRHS_enh);
		
			/* internal iterations */
			double res, res_0, res_rel;
			res = fRHS_enh.Magnitude();
			res_0 = res;
			res_rel = 1.0;
			int iter_enh = 0;
			
			//TEMP - try freezing modes for first iteration
			int global_iter = IterationNumber();
			while (global_iter > -1 && iter_enh++ < fLocalIterationMax && res > fAbsTol && res_rel > fRelTol)
			
//			while (iter_enh++ < fLocalIterationMax && res > fAbsTol && res_rel > fRelTol)
			{
				/* form the stiffness associated with the enhanced modes */
				fK22 = 0.0;
				FormStiffness_enhanced(fK22, NULL);

				/* update enhanced modes */
				fK22.LinearSolve(fRHS_enh);			

				/* update element modes */
				fElementModes.AddToRowScaled(CurrElementNumber(), -1.0, fRHS_enh);
		
				/* compute enhanced part of F and total F */
				ComputeEnhancedDeformation(needs_F, needs_F_last);
		
				/* calculate the residual from the internal force */
				fRHS_enh = 0.0;
				FormKd_enhanced(fPK1_list, fRHS_enh);
		
				/* errors */
				res = fRHS_enh.Magnitude();
				res_rel = res/res_0;
			}
			
			/* set element modes convergence flag */
			fModesConverged = res < fAbsTol || res_rel < fRelTol;
		
			/* form the stiffness associated with the enhanced modes */
			fK22 = 0.0;
			fK12 = 0.0;
			FormStiffness_enhanced(fK22, &fK12); /* also sets/stores c_ijkl */
		}
	}
}

/* form the element stiffness matrix */
void SimoFiniteStrainT::FormStiffness(double constK)
{		
	/* resolve based on method for solving element modes */
	switch (fModeSolveMethod)
	{
		case kMonolithic:
			FormStiffness_monolithic(constK);
			break;

		case kLocalIteration:
			FormStiffness_staggered(constK);
			break;
			
		default:
		
			cout << "\n SimoFiniteStrainT::FormStiffness: no implementation for solution method: "
			     << fModeSolveMethod << endl;
			throw ExceptionT::kGeneralFail;
	}
}

/* calculate the internal force contribution ("-k*d") */
void SimoFiniteStrainT::FormKd(double constK)
{
	/* resolve based on method for solving element modes */
	switch (fModeSolveMethod)
	{
		case kMonolithic:
			FormKd_monolithic(constK);
			break;

		case kLocalIteration:
			FormKd_staggered(constK);
			break;
			
		default:
		
			cout << "\n SimoFiniteStrainT::FormKd: no implementation for solution method: "
			     << fModeSolveMethod << endl;
			throw ExceptionT::kGeneralFail;
	}
}

/* write all current element information to the stream */
void SimoFiniteStrainT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	FiniteStrainT::CurrElementInfo(out);

	/* element modes */
	out << "\n element modes:\n" << fCurrElementModes << '\n';
	out << "\n element modes last:\n" << fCurrElementModes_last << '\n';

	/* write deformation gradients */
	out << "\n i.p. enhanced deformation gradients:\n";
	for (int i = 0; i < fF_enh_List.Length(); i++)
		out << " ip: " << i+1 << '\n'
		    << fF_enh_List[i] << '\n';
	out << '\n';
}

/***********************************************************************
* Private
***********************************************************************/

/* form the element stiffness matrix */
void SimoFiniteStrainT::FormStiffness_staggered(double constK)
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
	fStressStiff_22 = 0.0;
	fStressStiff_12 = 0.0;
	fK22 = 0.0;
	fK12 = 0.0;
	if (format != dMatrixT::kUpperOnly) fK21 = 0.0;
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* deformation gradient */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();

		/* scale factor */
		double scale = constK*(*Det++)*(*Weight++)*J;

		/* (retrieve) Cauchy stress */
		const dMatrixT& PK1 = fPK1_list[CurrIP()];
		fStressMat.MultABT(PK1, fTempMat1);
		fStressMat *= scale/J;

		/* (retrieve) D matrix */
		fD.SetToScaled(scale, fc_ijkl_list[CurrIP()]);

		/* transform derivatives to current coordinates */
		fTempMat1.Inverse(); /* F^-1 */
		fShapes->TransformDerivatives(fTempMat1, fDNa_x);
		fEnhancedShapes->TransformDerivatives_enhanced(fTempMat1, fDNa_x_enh);		

		/* get shape function gradients matricies */
		fShapes->GradNa(fDNa_x, fGradNa);
		fShapes->GradNa(fDNa_x_enh, fGradNa_enh);
		
		/* strain displacement matricies */
		Set_B(fDNa_x, fB);
		Set_B(fDNa_x_enh, fB_enh);

		/* stress stiffness (4.18) */
		fStressStiff_11.MultQTBQ(fGradNa, fStressMat, format, dMatrixT::kAccumulate);
		fStressStiff_22.MultQTBQ(fGradNa_enh, fStressMat, dMatrixT::kWhole, dMatrixT::kAccumulate);
		fStressStiff_12.MultATBC(fGradNa, fStressMat, fGradNa_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);

		/* material stiffness (4.14) */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
		fK22.MultQTBQ(fB_enh, fD, dMatrixT::kWhole, dMatrixT::kAccumulate);
		fK12.MultATBC(fB, fD, fB_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);
		
		/* non symmetric */
		if (format != dMatrixT::kUpperOnly)
		{
			/* material stiffness */
			fK21.MultATBC(fB_enh, fD, fB, dMatrixT::kWhole, dMatrixT::kAccumulate);
		}
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff_11, NumDOF(), dMatrixT::kAccumulate);
	fK22.Expand(fStressStiff_22, NumDOF(), dMatrixT::kAccumulate);
	fK12.Expand(fStressStiff_12, NumDOF(), dMatrixT::kAccumulate);

	/* condensation of element modes */
	fK22.Inverse();
	fK22 *= -1.0;
	if (format == dMatrixT::kUpperOnly)
		fLHS.MultQBQT(fK12, fK22, format, dMatrixT::kAccumulate);
	else
	{
		/* stress stiffness */
		fStressStiff_21.Transpose(fStressStiff_12);

		/* expand stress stiffness part */
		fK21.Expand(fStressStiff_21, NumDOF(), dMatrixT::kAccumulate);

		/* assemble */
		fLHS.MultABC(fK12, fK22, fK21, dMatrixT::kWhole, dMatrixT::kAccumulate);
	}
}

/* compute and assemble the element stiffness for the monolithic
 * solution scheme */
void SimoFiniteStrainT::FormStiffness_monolithic(double constK)
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
	fStressStiff_22 = 0.0;
	fStressStiff_12 = 0.0;
	fK11 = 0.0;
	fK22 = 0.0;
	fK12 = 0.0;
	if (format != dMatrixT::kUpperOnly) {
		fK21 = 0.0;
		fStressStiff_21 = 0.0;
	}
	
	/* integrate over element */
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* chain rule shape function derivatives */
		fTempMat1 = DeformationGradient();
		double J = fTempMat1.Det();
		fTempMat1.Inverse();
		fShapes->TransformDerivatives(fTempMat1, fDNa_x);
		fEnhancedShapes->TransformDerivatives_enhanced(fTempMat1, fDNa_x_enh);

		/* scale factor */
		double scale = constK*(*Det++)*(*Weight++)*J;

		/* Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
		fStressMat *= scale;

		/* get D matrix - material tangent modulus */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
		
		/* get shape function gradients matricies */
		fShapes->GradNa(fDNa_x, fGradNa);
		fShapes->GradNa(fDNa_x_enh, fGradNa_enh);
		
		/* strain displacement matricies */
		Set_B(fDNa_x, fB);
		Set_B(fDNa_x_enh, fB_enh);

		/* stress stiffness (4.18) */
		fStressStiff_11.MultQTBQ(fGradNa, fStressMat, format, dMatrixT::kAccumulate);
		fStressStiff_22.MultQTBQ(fGradNa_enh, fStressMat, format, dMatrixT::kAccumulate);
		fStressStiff_12.MultATBC(fGradNa, fStressMat, fGradNa_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);

		/* material stiffness (4.14) */
		fK11.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
		fK22.MultQTBQ(fB_enh, fD, format, dMatrixT::kAccumulate);
		fK12.MultATBC(fB, fD, fB_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);

		/* non symmetric */
		if (format != dMatrixT::kUpperOnly)
		{		
			/* material stiffness */
			fK21.MultATBC(fB_enh, fD, fB, dMatrixT::kWhole, dMatrixT::kAccumulate);
		}
	}
						
	/* expand/assemble stress stiffness */
	fK11.Expand(fStressStiff_11, NumDOF(), dMatrixT::kAccumulate);
	fK22.Expand(fStressStiff_22, NumDOF(), dMatrixT::kAccumulate);
	fK12.Expand(fStressStiff_12, NumDOF(), dMatrixT::kAccumulate);
	
	/* assemble into element stiffness matrix */
	fLHS.AddBlock(0          , 0          , fK11);
	fLHS.AddBlock(fK11.Rows(), fK11.Cols(), fK22);
	fLHS.AddBlock(0          , fK11.Cols(), fK12);
	
	/* non symmetric */
	if (format != dMatrixT::kUpperOnly) {

		/* stress stiffness */
		fStressStiff_21.Transpose(fStressStiff_12);
	
		/* expand stress stiffness term */
		fK21.Expand(fStressStiff_21, NumDOF(), dMatrixT::kAccumulate);

		/* assemble */
		fLHS.AddBlock(fK11.Rows(), 0, fK21);
	}
}

/* form the contribution to the the residual force associated with the 
 * Galerkin part of the deformation gradient */
void SimoFiniteStrainT::FormKd_staggered(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT WP(NumSD(), fStressStiff_11.Rows(), fNEEvec.Pointer());

	/* integration rules */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* integrate over the element */
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
	if (!fModesConverged)
	{
		/* "project" to nodal force using (4.20) */
		fK22.LinearSolve(fRHS_enh);
		fK12.Multx(fRHS_enh, fNEEvec);
		fRHS.AddScaled(-constK, fNEEvec);
	}
}

/* compute and assemble the residual force for the monolithic
 * solution scheme */
void SimoFiniteStrainT::FormKd_monolithic(double constK)
{
	const char caller[] = "SimoFiniteStrainT::FormKd_enhanced";

	/* matrix alias to fTemp */
	dMatrixT WP(NumSD(), fStressStiff_11.Rows(), fNEEvec.Pointer());
	
	/* partition residual force vector */
	int neq = NumElementNodes()*NumDOF();
	dArrayT RHS(neq, fRHS.Pointer());
	dArrayT RHS_enh(fCurrElementModes.Length(), fRHS.Pointer(neq));

	/* integration rule */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F^(-1) */
		fTempMat2 = DeformationGradient();
		double J = fTempMat2.Det();
		if (J <= 0.0) 
			ExceptionT::BadJacobianDet(caller);
		else
			fTempMat2.Inverse();

		/* compute PK1/J */
		fStressMat.MultABT(fTempMat1, fTempMat2);

		/* integration weight */
		double w = J*constK*(*Weight++)*(*Det++);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* W_i,J P_iJ */
		WP.MultAB(fStressMat, fGradNa);

		/* accumulate */
		RHS.AddScaled(w, WP);

		/* enhanced shape function gradients */
		fEnhancedShapes->GradNa_enhanced(fGradNa_enh);

		/* Wenh_i,J P_iJ */
		fWP_enh.MultAB(fStressMat, fGradNa_enh);

		/* accumulate */
		RHS_enh.AddScaled(w, fWP_enh);
	}	
}

/* compute modified, enhanced deformation gradient */
void SimoFiniteStrainT::ModifiedEnhancedDeformation(void)
{
	const char caller[] = "SimoFiniteStrainT::ModifiedEnhancedDeformation";

	/* check */
	if (NumSD() != 3) ExceptionT::GeneralFail(caller, "3D only");
	
	/* skip base class implementation because the deformation gradient
	 * modification are not simply additive */
	SolidElementT::SetGlobalShape();
	
	//TEMP
	ExceptionT::GeneralFail(caller, "not implemented");
}

/* compute enhanced part of F and total F */
void SimoFiniteStrainT::ComputeEnhancedDeformation(bool need_F, bool need_F_last)
{
	/* set modes for the current element */
	int current_element = CurrElementNumber();

	/* store Galerkin part/compute enhancement and total F */
	if (need_F)
	{
		/* collect element modes in local ordering */
		fCurrElementModes.FromTranspose(fElementModes(current_element));
	
		/* compute enhanced part of the deformation gradient */
		for (int i = 0; i < NumIP(); i++)
			fEnhancedShapes->GradU_enhanced(fCurrElementModes, fF_enh_List[i], i);
			
		/* compute total, enhanced deformation gradient */
		fF_all.SumOf(fF_Galerkin_all, fF_enh_all);
	}		

	/* store Galerkin part/compute enhancement and total F from last step */
	if (need_F_last)
	{
		/* collect element modes in local ordering */
		fCurrElementModes_last.FromTranspose(fElementModes_last(current_element));

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
	const char caller[] = "SimoFiniteStrainT::FormKd_enhanced";
	
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
			ExceptionT::BadJacobianDet(caller);
		else
			fTempMat2.Inverse();

		/* compute PK1 */
		dMatrixT& PK1 = PK1_list[CurrIP()];
		PK1.MultABT(fTempMat1, fTempMat2);
		PK1 *= J;

		/* enhanced shape function gradients */
		fEnhancedShapes->GradNa_enhanced(fGradNa_enh);

		/* Wi,J PiJ */
		fWP_enh.MultAB(PK1, fGradNa_enh);

		/* accumulate */
		RHS_enh.AddScaled((*Weight++)*(*Det++), fWP_enh);
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
		Set_B(fDNa_x_enh, fB_enh);

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
			Set_B(fDNa_x, fB);

			/* accumulate */
			K_12->MultATBC(fB, fD, fB_enh, dMatrixT::kWhole, dMatrixT::kAccumulate);
		}
	}
						
	/* expand and add in stress stiffness parts */
	K_22.Expand(fStressStiff_22, NumDOF(), dMatrixT::kAccumulate);
	if (K_12) K_12->Expand(fStressStiff_12, NumDOF(), dMatrixT::kAccumulate);
}
