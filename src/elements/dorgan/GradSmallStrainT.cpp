/* $Id: GradSmallStrainT.cpp,v 1.4 2004-04-01 22:46:54 rdorgan Exp $ */ 
#include "GradSmallStrainT.h"

#include "ShapeFunctionT.h"
#include "C1ShapeFunctionT.h"

#include <math.h>
#include "ifstreamT.h"
#include "ofstreamT.h"

/* materials */
#include "GradSSSolidMatT.h"
#include "GradSSMatSupportT.h"
#include "SolidMatListT.h"

#include "FieldT.h"
#include "FieldSupportT.h"

using namespace Tahoe;

/* constructor */
GradSmallStrainT::GradSmallStrainT(const ElementSupportT& support, 
								   const FieldT& disp, const FieldT& iso_hard):
	SmallStrainT(support, disp), //pass the displacement field to the base class

	fDisplacement(disp),
	fIsoHardening(iso_hard),

	fLocR(LocalArrayT::kDisp),
	fLocLastR(LocalArrayT::kLastDisp),

	fShapes_R (NULL),

	fK_bb(ElementMatrixT::kNonSymmetric),
	fK_bh(ElementMatrixT::kNonSymmetric),
	fK_hb(ElementMatrixT::kNonSymmetric),
	fK_hh(ElementMatrixT::kNonSymmetric),
	fK_hp(ElementMatrixT::kNonSymmetric),
	fK_hq(ElementMatrixT::kNonSymmetric),
	fK_ct(ElementMatrixT::kNonSymmetric),


	fDM_bb(dSymMatrixT::NumValues(NumSD())),
	fOM_bh(dSymMatrixT::NumValues(NumSD()), NumSD()),
	fOM_hb(NumSD(), dSymMatrixT::NumValues(NumSD())),
	fGM_hh(NumSD()),
	fGM_hp(NumSD()),
	fI(1)

{
	ifstreamT& in = ElementSupport().Input();

	/* control parameters for iso_hard mesh */
	in >> fGeometryCode_R;
	in >> fNumIP_R;

	/********DEBUG*******/
	in >> print_GlobalShape;
	in >> print_Kd;
	in >> print_Stiffness;
	in >> fKConstraintA;
//	in >> fKConstraintB;
	in >> fRConstraintA;
	in >> fRConstraintB;
	/*******************/
}

/* destructor */
GradSmallStrainT::~GradSmallStrainT(void) 
{  
	delete fShapes_R;
}

void GradSmallStrainT::Initialize(void)
{
	const char caller[] = "GradSmallStrainT::Initialize";

	/* dimensions */
	fNumIP	= NumIP();
	fNumIP_R      = NumIP_Field();
	fNumSD	= NumSD();
	fNumDOF       = fDisplacement.NumDOF();
	fNumDOF_R     = fIsoHardening.NumDOF();
	fNumDOF_Total = fNumDOF + fNumDOF_R;

	/* checks */
	if(fNumIP != fNumIP_R)
		ExceptionT::BadInputValue(caller, "NumIP != NumIP_R");
	if(fNumDOF != fNumSD)
		ExceptionT::BadInputValue(caller, "NumDOF != NumSD");
	if(fNumDOF_R != pow(2,fNumSD))
		ExceptionT::BadInputValue(caller, "NumDOF_R != 2^NumSD");
	if(fNumSD != 1)
		ExceptionT::BadInputValue(caller, "NumSD != 1");

	/* inherited */
	SmallStrainT::Initialize();
	int fNumElementNodes = NumElementNodes();

	/* allocate lists */
	fR_List.Dimension(fNumIP_R);	 // isotropic hardening
	fR_last_List.Dimension(fNumIP_R);    // "last" isotropic hardening
	fLapR_List.Dimension(fNumIP_R);      // laplacian isotropic hardening
	fLapR_last_List.Dimension(fNumIP_R); // "last" laplacian isotropic hardening
	
	/* dimension work space */
	fK_bb.Dimension(fNumElementNodes*fNumDOF,   fNumElementNodes*fNumDOF  );
	fK_bh.Dimension(fNumElementNodes*fNumDOF,   fNumElementNodes*fNumDOF_R);
	fK_hb.Dimension(fNumElementNodes*fNumDOF_R, fNumElementNodes*fNumDOF  );
	fK_hh.Dimension(fNumElementNodes*fNumDOF_R, fNumElementNodes*fNumDOF_R);
	fK_hp.Dimension(fNumElementNodes*fNumDOF_R, fNumElementNodes*fNumDOF_R);
	fK_hq.Dimension(fNumElementNodes*fNumDOF_R, fNumElementNodes*fNumDOF_R);
	fK_ct.Dimension(fNumElementNodes*fNumDOF_R, fNumElementNodes*fNumDOF_R);


	/* resize work arrays */
	fRHS.Dimension(fNumElementNodes*fNumDOF_Total);
	fLHS.Dimension(fRHS.Length());

	/* allocate shape functions */
	fh.Dimension (1, fNumElementNodes*fNumDOF_R);
	fhT.Dimension(fNumElementNodes*fNumDOF_R, 1);
	fp.Dimension (1, fNumElementNodes*fNumDOF_R);
	fpT.Dimension(fNumElementNodes*fNumDOF_R, 1);
}

void GradSmallStrainT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
				     AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* connectivities */
		const iArray2DT& connects = *(fConnectivities[i]);
		int fNumElements = connects.MajorDim();
		int fNumElementNodes = NumElementNodes();

		/* dimension */
		fEqnos[i].Dimension  (fNumElements, fNumElementNodes*fNumDOF_Total);
		iArray2DT fEqnos_Disp(fNumElements, fNumElementNodes*fNumDOF      );
		iArray2DT fEqnos_R   (fNumElements, fNumElementNodes*fNumDOF_R    );
			
		/* get equation numbers */
		fDisplacement.SetLocalEqnos(connects, fEqnos_Disp);
		fIsoHardening.SetLocalEqnos(connects, fEqnos_R   );

		/* write into one array */
		fEqnos[i].BlockColumnCopyAt(fEqnos_Disp, 0		     );
		fEqnos[i].BlockColumnCopyAt(fEqnos_R   , fEqnos_Disp.MinorDim());

		/* add to list of equation numbers */
		eq_1.Append(&fEqnos[i]);
	}
	
	/* reset pointers to element cards */
	SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* GradSmallStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate new */
	if (!p) p = new GradSSMatSupportT(NumSD(), NumDOF(), NumDOF_Field(), NumIP(), NumIP_Field());

	/* inherited initializations */
	SmallStrainT::NewMaterialSupport(p);

	/* set GradSSMatSupportT fields */
	GradSSMatSupportT* pss = dynamic_cast<GradSSMatSupportT*>(p);
	if (pss) {
		pss->SetLinearField(&fR_List);
		pss->SetLinearField_last(&fR_last_List);

		pss->SetLinearLaplacianField(&fLapR_List);
		pss->SetLinearLaplacianField_last(&fLapR_last_List);
	}

	return p;
}

/* initialize local arrays */
void GradSmallStrainT::SetLocalArrays(void)
{
	/* inherited */
	SmallStrainT::SetLocalArrays();

	int fNumElementNodes = NumElementNodes();

	/* dimension local arrays */
	fLocR.Dimension    (fNumElementNodes, fNumDOF_R);
	fLocLastR.Dimension(fNumElementNodes, fNumDOF_R);
	
	/* register local arrays */
	fIsoHardening.RegisterLocal(fLocR    );	
	fIsoHardening.RegisterLocal(fLocLastR);	
}

/* initialization functions */
void GradSmallStrainT::SetShape(void)
{
	/* inherited */
	SmallStrainT::SetShape();

	fShapes_R = new C1ShapeFunctionT(GeometryCode_Field(), NumIP_Field(), fLocInitCoords, NumDOF_Field());
	if (!fShapes_R) throw ExceptionT::kOutOfMemory;
	fShapes_R->Initialize();
}

/* form shape functions and derivatives */
void GradSmallStrainT::SetGlobalShape(void)
{
	/* inherited */
	SmallStrainT::SetGlobalShape();

	/* collect element values of R */
	SetLocalU(fLocR);
	SetLocalU(fLocLastR);

	/********DEBUG*******/
	if (fRConstraintA)
		/* enforce loading/unloading conditions */
		for (int nd = 0; nd < NumElementNodes(); nd++)
		{
			if (fLocR[nd] < 0.) fLocR[nd] = 0.; 
			if (fLocLastR[nd] < 0.) fLocLastR[nd] = 0.; 
			if (fLocR[nd] < fLocLastR[nd]) fLocR[nd] = fLocLastR[nd];
		}
	/*******************/

	/********DEBUG*******/
	if (print_GlobalShape)
	{
		if (CurrElementNumber() == 0)
			cout << "GradSmallStrainT::SetGlobalShape" << endl;
		cout << "    element: " << CurrElementNumber() << endl;
	}
	/*******************/

	fShapes_R ->TopIP();
	/* loop over integration points */
	while(fShapes_R->NextIP())
	{
		int ip = fShapes_R->CurrIP();

		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);

		/* setup workspace */
		dArrayT temp(1);

		/* compute iso_hard using h shape functions */
		fh.Multx(fLocR, temp);
		fR_List[ip] = temp[0];

		fh.Multx(fLocLastR, temp);
		fR_last_List[ip] = temp[0];

		/* compute laplacian of iso_hard using p shape functions */
		fp.Multx(fLocR, temp);
		fLapR_List[ip] = temp[0];

		fp.Multx(fLocLastR, temp);
		fLapR_last_List[ip] = temp[0];

		/********DEBUG*******/
		if (fRConstraintB)
		{
			if (fR_List[ip] < 0.) fR_List[ip] = 0.; 
			if (fR_last_List[ip] < 0.) fR_last_List[ip] = 0.; 
			if (fR_List[ip] < fR_last_List[ip]) fR_List[ip] = fR_last_List[ip];
		}
		/*******************/
	}

	/********DEBUG*******/
	if (print_GlobalShape)
	{
		for (int nd_dof = 0; nd_dof < NumElementNodes()*NumDOF_Field(); nd_dof++)
			cout << "	fLocR[" << nd_dof << "]: " << fLocR[nd_dof] << endl;
		for (int ip = 0; ip < NumIP(); ip++)
			cout << "      fR_List[" << ip << "]: " << fR_List[ip] << endl;
	}
	/*******************/
}

/* form the element stiffness matrix */
void GradSmallStrainT::FormStiffness(double constK)
{
	const char caller[] = "GradSmallStrainT::FormStiffness";
	int curr_group = ElementSupport().CurrentGroup();

	int fNumElementNodes = NumElementNodes();
	
	/* dmatrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	if(format != dMatrixT::kWhole)
		ExceptionT::BadInputValue(caller, "LHS must be NonSymmetric");

	/* integration rules (same for both fields)*/
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	if (print_Stiffness)
	{
		if (CurrElementNumber() == 0)
			cout << caller << endl;
		cout << "    element: " << CurrElementNumber() << endl;
	}
	/*******************/

	/* initialize */
	fK_bb = 0.0;
	fK_bh = 0.0;
	fK_hb = 0.0;
	fK_hh = 0.0;
	fK_ct = 0.0;
	fK_hp = 0.0;

	fShapes   ->TopIP();
	fShapes_R ->TopIP();

	/* integrate stiffness over the elements */
	while(fShapes->NextIP() && fShapes_R->NextIP())
	{
		/* integration weight */
		double scale = constK*(*Det++)*(*Weight++);

		/* accumulate strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B    (fShapes->Derivatives_U(), fB);

		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);

		/* compute elastic stiffness matrix */
		fDM_bb.SetToScaled(scale, fCurrMaterial_Grad->dm_bb_ijkl());   // D matrix - algorithmic stiffness operator
		fK_bb.MultQTBQ(fB, fDM_bb, format, dMatrixT::kAccumulate);

		/* compute off-diagonal matrices */
		fOM_bh.SetToScaled(scale, fCurrMaterial_Grad->om_bh_ij());
		fK_bh.MultATBC(fB, fOM_bh, fh, format, dMatrixT::kAccumulate);

		fOM_hb.SetToScaled(scale, fCurrMaterial_Grad->om_hb_ij());
		fK_hb.MultATBC(fh, fOM_hb, fB, format, dMatrixT::kAccumulate);

		/* compute non-symmetric, gradient dependent matrices */
		fGM_hp.SetToScaled(scale, fCurrMaterial_Grad->gm_hp());
		fK_hp.MultATBC(fh, fGM_hp, fp, format, dMatrixT::kAccumulate);

		fGM_hh.SetToScaled(scale, fCurrMaterial_Grad->gm_hh());
		fK_hh.MultQTBQ(fh, fGM_hh, format, dMatrixT::kAccumulate);

		 /* material info for constraint */
		double yield = fCurrMaterial_Grad->yc();
		double del_r = fCurrMaterial_Grad->del_Field();

		fI = 0.;

		/* apply constraint - elastic or decreasing R */
		if (yield <= 0. || del_r <= 0.0)
			fI = 1.0*scale*fKConstraintA;
		fK_ct.MultQTBQ(fh, fI, format, dMatrixT::kAccumulate);

		/********DEBUG*******/
		if (print_Stiffness)
		{
			cout<<"     ip: "             << CurrIP()                                           <<endl;
			cout<<"                  R: " << fR_List[CurrIP()]                                  <<endl;
			cout<<"           R-R_Last: " << fR_List[CurrIP()] - fR_last_List[CurrIP()]         <<endl;
			cout<<"             strain: " << (fStrain_List[CurrIP()])[0]                        <<endl;
			cout<<"             stress: " << (fCurrMaterial_Grad->s_ij())[0]                    <<endl;
			cout<<"                 yc: " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"           om_bh_ij: " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"           om_hb_ij: " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"              gm_hh: " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"              gm_hp: " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			for (int i = 0; i < NumElementNodes()*NumDOF_Field(); i++)
				cout<<"               h" << i << "]: " << fh[i]                                 <<endl;
			for (int i = 0; i < NumElementNodes(); i++)
				cout<<"               B" << i << "]: " << fB[i]                                 <<endl;
		}
		/*******************/
	}		

	/* assemble into element stiffness matrix */
	fLHS.AddBlock(0           , 0           , fK_bb);
	fLHS.AddBlock(0           , fK_bb.Cols(), fK_bh);
	fLHS.AddBlock(fK_bb.Rows(), 0           , fK_hb);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hh);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hp);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hq);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_ct);

#if 0
	// add constraint to Krr if elastic
	fK_hh = 0.;
	for (int nd = 0; nd < NumElementNodes(); nd ++)
	  if (ElasticNode())
	    fK_hh[nd][nd] = fKConstraintB * fLocR[nd];
	fLHS.AddBlock(fK_hh.Rows(), fK_hh.Cols(), fK_hh);
#endif
}

/* calculate the internal force contribution ("-k*d") */
void GradSmallStrainT::FormKd(double constK)
{
	const char caller[] = "GradSmallStrainT::FormKd";
	int curr_group = ElementSupport().CurrentGroup();

	int fNumElementNodes = NumElementNodes();

	/* partition residual force vector */
	int neq_Disp = fNumElementNodes*fNumDOF;
	int neq_R    = fNumElementNodes*fNumDOF_R;

	dArrayT RHS_Disp(neq_Disp, fRHS.Pointer()        );
	dArrayT RHS_R   (neq_R   , fRHS.Pointer(neq_Disp));

	/* integration rules (same for both fields)*/
	const double* Det     = fShapes->IPDets();
	const double* Weight  = fShapes->IPWeights();
	
	if (print_Kd)
	{
		if (CurrElementNumber() == 0)
			cout << caller << endl;
		cout<<"    element: "<<CurrElementNumber()<<endl;
	}	

	/* initialize */
	fShapes->TopIP();
	fShapes_R ->TopIP();

	/* integrate residuals over the elements */
	while(fShapes->NextIP() && fShapes_R->NextIP())
	{
		/* integration weight */
		double scale = constK*(*Det++)*(*Weight++);

		/* accumulate strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* compute internal stresses residual */
		fB.MultTx(fCurrMaterial_Grad->s_ij(), fNEEvec);

		/* accumulate in rhs dislacement equations */
		RHS_Disp.AddScaled(scale, fNEEvec);
		
		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);

		/* material info for constraint */
		double yield = fCurrMaterial_Grad->yc();
		double del_r = fCurrMaterial_Grad->del_Field();
		double r = fCurrMaterial_Grad->Field();

		/* apply constraint - elastic or decreasing R */
		double const_RHS_R = -scale*yield;

#if 1
		if (yield < kSmall || del_r < 0.0) {
		  //			const_RHS_R -= scale*fKConstraintA*del_r;
			const_RHS_R += 1.0*scale*fKConstraintA*r;
		}
#endif

		/* accumulate in rhs iso_hard equations */
		RHS_R.AddScaled(const_RHS_R, fhT.Transpose(fh));

		/********DEBUG*******/
		if (print_Stiffness)
		{
			cout<<"     ip: "             << CurrIP()                                           <<endl;
			cout<<"                  R: " << fR_List[CurrIP()]                                  <<endl;
			cout<<"           R-R_Last: " << fR_List[CurrIP()] - fR_last_List[CurrIP()]         <<endl;
			cout<<"             strain: " << (fStrain_List[CurrIP()])[0]                        <<endl;
			cout<<"             stress: " << (fCurrMaterial_Grad->s_ij())[0]                    <<endl;
			cout<<"                 yc: " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"           om_bh_ij: " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"           om_hb_ij: " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"              gm_hh: " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"              gm_hp: " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			for (int i = 0; i < NumElementNodes()*NumDOF_Field(); i++)
				cout<<"               h" << i << "]: " << fh[i]                                 <<endl;
			for (int i = 0; i < NumElementNodes(); i++)
				cout<<"               B" << i << "]: " << fB[i]                                 <<endl;
		}
		/*******************/
	}
}

/* current element operations */
bool GradSmallStrainT::NextElement(void)
{
	/* inherited */
	bool result = SolidElementT::NextElement();
	
	/* get material pointer */
	if (result)
	{
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[CurrentElement().MaterialNumber()];
	
		/* cast is safe since class contructs materials list */
		fCurrMaterial_Grad = (GradSSSolidMatT*) pcont_mat;
	}
	
	return result;
}

/* form of tangent matrix */
GlobalT::SystemTypeT GradSmallStrainT::TangentType(void) const
{
	// should be obtained from material through SolidElementT::TangentType
	return GlobalT::kNonSymmetric;
}


/* print element group data */
void GradSmallStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	SmallStrainT::PrintControlData(out);

	out << " Associated field. . . . . . . . . . . . . . . . = \"" << fIsoHardening.Name() << "\"\n";
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode_R << '\n';
	//	out << "    eq." << C1GeometryT::kPoint	 << ", point\n";
	out << "    eq." << C1GeometryT::kC1Line	<< ", C1line\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP_R    << '\n';
}

/***********************************************************************
* Private
***********************************************************************/

/* accumulate C1 shape function matrix */
void GradSmallStrainT::Set_h(dMatrixT& h) const
{
	const char caller[] = "GradSmallStrainT::Set_h";

	int  fNumElementNodes = NumElementNodes();
	const double* pShapes_R = fShapes_R->IPShapeU();
	double* ph = h.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		for (int i = 0; i < fNumElementNodes*fNumDOF_R; i++)
			*ph++ = *pShapes_R++;
	}
	else
		ExceptionT::GeneralFail(caller, "implemented for 1d only");

      double jacobian = (fLocInitCoords[1]-fLocInitCoords[0])/2;
      for (int i = fNumElementNodes; i < fNumElementNodes*NumDOF_Field(); i++)
	      h[i]=h[i]*jacobian;
}

/* accumulate Laplacian C1 shape function matrix */
void GradSmallStrainT::Set_p(dMatrixT& p) const
{
	const char caller[] = "GradSmallStrainT::Set_p";

	int  fNumElementNodes = NumElementNodes();
	const double* pShapes_LapR = fShapes_R->IPShapeDDU();
	double* pp = p.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		for (int i = 0; i < fNumElementNodes*fNumDOF_R; i++)
			*pp++ = *pShapes_LapR++;
	}
	else
		ExceptionT::GeneralFail(caller, "implemented for 1d only");

      double jacobian = (fLocInitCoords[1]-fLocInitCoords[0])/2;
      for (int i = 0; i < NumElementNodes(); i++)
	      p[i]=p[i]/jacobian/jacobian;
      for (int i = NumElementNodes(); i < NumElementNodes()*NumDOF_Field(); i++)
	      p[i]=p[i]/jacobian;
}

// CHECK COEFFICIENT 'scale' AND 'constK'

// need global loading/unloading conditions

// solves for new step of R in iterative form: inc_R_current-iteration = inc_R_previous-iteration + dR
//  where dR comes from solution of Kd=r
// this allows the dR to go positive, and then not back to zero.
// need some sort of penalty to force the following:
//   1) dR >= 0
//   2) if, during the computation of the total increment of R for the time step,
//	  f changes from greater than zero to zero, the total increment of R should go back to zero.
