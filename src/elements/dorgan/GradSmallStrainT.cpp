/* $Id: GradSmallStrainT.cpp,v 1.7 2004-06-09 00:25:53 rdorgan Exp $ */ 
#include "GradSmallStrainT.h"

/* shape functions */
#include "ShapeFunctionT.h"
#include "C1ShapeTools.h"

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
								   const FieldT& disp, const FieldT& field1, const FieldT& field2):
	SmallStrainT(support, disp), //pass the displacement field to the base class
	fDisplacement(disp),
	fField1(field1),
	fField2(field2),
	
	fLocField1(LocalArrayT::kDisp),
	fLocLastField1(LocalArrayT::kLastDisp),
	fLocField2(LocalArrayT::kDisp),
	fLocLastField2(LocalArrayT::kLastDisp),

	fShapes_Field (NULL),
	
	fK_bb(ElementMatrixT::kNonSymmetric),
	fK_bh(ElementMatrixT::kNonSymmetric),
	fK_hb(ElementMatrixT::kNonSymmetric),
	fK_hh(ElementMatrixT::kNonSymmetric),
	fK_hp(ElementMatrixT::kNonSymmetric),
	fK_hq(ElementMatrixT::kNonSymmetric),
	fK_ct(ElementMatrixT::kNonSymmetric),

	fDM_bb(dSymMatrixT::NumValues(NumSD())),
	fOM_bh(dSymMatrixT::NumValues(NumSD()), 1),
	fOM_hb(1, dSymMatrixT::NumValues(NumSD())),
	fGM_hh(1),
	fGM_hp(1, NumSD()),
	fGM_hq(1),

	fI(1)
{
	ifstreamT& in = ElementSupport().Input();

	/* control parameters for field mesh */
	in >> fGeometryCode_Field;
	in >> fNumIP_Field;
	in >> fNumElementNodes_Field;
	in >> fDegreeOfContinuity_Field;
	in >> fNodalConstraint;

	/********DEBUG*******/
	in >> print_GlobalShape;
	in >> print_Kd;
	in >> print_KdMatrix;
	in >> print_Stiffness;
	in >> print_StiffnessMatrix;
	/*******************/
}

/* destructor */
GradSmallStrainT::~GradSmallStrainT(void) 
{  
	delete fShapes_Field;
}

void GradSmallStrainT::Initialize(void)
{
	const char caller[] = "GradSmallStrainT::Initialize";

	/* dimensions */
	fNumSD        = NumSD();
	fNumIP_Disp   = NumIP();
	fNumIP_Field  = NumIP_Field();
	fNumDOF_Disp  = fDisplacement.NumDOF();
	fNumDOF_Field1 = fField1.NumDOF();
	fNumDOF_Field2 = fField2.NumDOF();
	fNumDOF_Field = fField1.NumDOF() + fField2.NumDOF();
	
	/* check spatial dimensions */
	if (fNumSD != 1)
		ExceptionT::BadInputValue(caller, "fNumSD != 1");

	/* additional dimension checks */
	if (fNumSD == 1) {
		if (fNumIP_Disp != fNumIP_Field)
			ExceptionT::BadInputValue(caller, "fNumIP_Disp != fNumIP_Field");
		if (fDegreeOfContinuity_Field != 1)
			ExceptionT::BadInputValue(caller, "fDegreeOfContinuity_Field != 1");
		if (fDegreeOfContinuity_Field == 1)
			if (fNumDOF_Field != 2)
				ExceptionT::BadInputValue(caller, "fNumDOF_Field != 2^fNumSD");
	}

	/* inherited */
	SmallStrainT::Initialize();

	int	fNumElementNodes_Disp  = NumElementNodes();
	int	fNumElementNodes_Field = NumElementNodes_Field();

	// should be changed to allow different nodes for field and for displacement
	//	if (fNumElementNodes_Disp != fNumElementNodes_Field)
	//		ExceptionT::BadInputValue(caller, "fNumElementNodes_Disp != fNumElementNodes_Field");

	/* allocate lists */
	fField_List.Dimension(fNumIP_Field);       // field
	fField_last_List.Dimension(fNumIP_Field);  // "last" field

	fGradField_List.Dimension(fNumIP_Field);       // field
	fGradField_last_List.Dimension(fNumIP_Field);  // "last" field

	fLapField_List.Dimension(fNumIP_Field);       // field
	fLapField_last_List.Dimension(fNumIP_Field);  // "last" field

	fYield_List.Dimension(fNumIP_Disp);
	
	/* dimension work space */
	fK_bb.Dimension(fNumElementNodes_Disp *fNumDOF_Disp,  fNumElementNodes_Disp *fNumDOF_Disp );
	fK_bh.Dimension(fNumElementNodes_Disp *fNumDOF_Disp,  fNumElementNodes_Field*fNumDOF_Field);
	fK_hb.Dimension(fNumElementNodes_Field*fNumDOF_Field, fNumElementNodes_Disp *fNumDOF_Disp );
	fK_hh.Dimension(fNumElementNodes_Field*fNumDOF_Field, fNumElementNodes_Field*fNumDOF_Field);
	fK_hp.Dimension(fNumElementNodes_Field*fNumDOF_Field, fNumElementNodes_Field*fNumDOF_Field);
	fK_hq.Dimension(fNumElementNodes_Field*fNumDOF_Field, fNumElementNodes_Field*fNumDOF_Field);
	fK_ct.Dimension(fNumElementNodes_Field*fNumDOF_Field, fNumElementNodes_Field*fNumDOF_Field);

	/* resize work arrays */
	fNumEQ_Total = fNumElementNodes_Disp*fNumDOF_Disp+fNumElementNodes_Field*fNumDOF_Field;
	fRHS.Dimension(fNumEQ_Total);
	fLHS.Dimension(fRHS.Length());

	/* allocate shape functions */
	fh.Dimension (1, fNumElementNodes_Field*fNumDOF_Field);
	fhT.Dimension(fNumElementNodes_Field*fNumDOF_Field, 1);
	fh1.Dimension (1, fNumElementNodes_Field*fNumDOF_Field1);
	fh1T.Dimension(fNumElementNodes_Field*fNumDOF_Field1, 1);
	fh2.Dimension (1, fNumElementNodes_Field*fNumDOF_Field2);
	fh2T.Dimension(fNumElementNodes_Field*fNumDOF_Field2, 1);

	/* allocate gradient shape functions */
	fp.Dimension (fNumSD, fNumElementNodes_Field*fNumDOF_Field);
	fp1.Dimension (1, fNumElementNodes_Field*fNumDOF_Field1);
	fp2.Dimension (1, fNumElementNodes_Field*fNumDOF_Field2);

	/* allocate Laplacian shape functions */
	fq.Dimension (1, fNumElementNodes_Field*fNumDOF_Field);
	fq1.Dimension (1, fNumElementNodes_Field*fNumDOF_Field1);
	fq2.Dimension (1, fNumElementNodes_Field*fNumDOF_Field2);
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
		int fNumElementNodes_Disp  = NumElementNodes();
		int fNumElementNodes_Field = NumElementNodes_Field();

		/* dimension */
		fEqnos[i].Dimension   (fNumElements, fNumEQ_Total);
		iArray2DT fEqnos_Disp (fNumElements, fNumElementNodes_Disp *fNumDOF_Disp );
		iArray2DT fEqnos_Field1(fNumElements, fNumElementNodes_Field*fNumDOF_Field1);
		iArray2DT fEqnos_Field2(fNumElements, fNumElementNodes_Field*fNumDOF_Field2);

		/* get equation numbers */
		fDisplacement.SetLocalEqnos(connects, fEqnos_Disp);
		fField1.SetLocalEqnos(connects, fEqnos_Field1);
		fField2.SetLocalEqnos(connects, fEqnos_Field2);

		/* write into one array */
		fEqnos[i].BlockColumnCopyAt(fEqnos_Disp, 0);
		fEqnos[i].BlockColumnCopyAt(fEqnos_Field1, fEqnos_Disp.MinorDim());
		fEqnos[i].BlockColumnCopyAt(fEqnos_Field2, fEqnos_Disp.MinorDim() + fEqnos_Field1.MinorDim());		

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
		pss->SetLinearField(&fField_List);
		pss->SetLinearField_last(&fField_last_List);

		pss->SetLinearGradField(&fGradField_List);
		pss->SetLinearGradField_last(&fGradField_last_List);

		pss->SetLinearLapField(&fLapField_List);
		pss->SetLinearLapField_last(&fLapField_last_List);
	}

	return p;
}

/* initialize local arrays */
void GradSmallStrainT::SetLocalArrays(void)
{
	/* inherited */
	SmallStrainT::SetLocalArrays();

	int fNumElementNodes_Disp = NumElementNodes();
	int fNumElementNodes_Field = NumElementNodes_Field();
	
	/* dimension local arrays */
	fLocField1.Dimension(fNumElementNodes_Field, fNumDOF_Field1);
	fLocLastField1.Dimension(fNumElementNodes_Field, fNumDOF_Field1);
	fLocField2.Dimension(fNumElementNodes_Field, fNumDOF_Field2);
	fLocLastField2.Dimension(fNumElementNodes_Field, fNumDOF_Field2);	
	
	/* register local arrays */
	fField1.RegisterLocal(fLocField1);
	fField1.RegisterLocal(fLocLastField1);
	fField2.RegisterLocal(fLocField2);
	fField2.RegisterLocal(fLocLastField2);	
}

/* initialization functions */
void GradSmallStrainT::SetShape(void)
{
	/* inherited */
	SmallStrainT::SetShape();

	fShapes_Field = new C1ShapeTools(GeometryCode_Field(), NumIP_Field(), fLocInitCoords, NumDOF_Field());
	if (!fShapes_Field) throw ExceptionT::kOutOfMemory;
	fShapes_Field->Initialize();
}

/* form shape functions and derivatives */
void GradSmallStrainT::SetGlobalShape(void)
{
	const char caller[] = "GradSmallStrainT::SetGlobalShape";

	/* inherited */
	SmallStrainT::SetGlobalShape();

	/* collect element values of Field */
	SetLocalU(fLocField1);
	SetLocalU(fLocLastField1);
	SetLocalU(fLocField2);
	SetLocalU(fLocLastField2);	

	/* compute shape function derivatives */
	fShapes_Field->SetDerivatives();

	/********DEBUG*******/
	if (print_GlobalShape)
	{
		cout << caller << endl;
		cout << "  element: " << CurrElementNumber() << endl;
	}
	/*******************/

	fShapes->TopIP();
	/* loop over integration points */
	while(fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
			
		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);
		Set_q(fq);

		/* setup workspace */
		dArrayT temp(1);

		for (int i = 0; i < NumElementNodes_Field(); i++)
		{
			fh1[i] = fh[i];
			fh2[i] = fh[i + NumElementNodes_Field()];

			fp1[i] = fp[i];
			fp2[i] = fp[i + NumElementNodes_Field()];

			fq1[i] = fq[i];
			fq2[i] = fq[i + NumElementNodes_Field()];
		}
		
		/* compute field using h shape functions */
		fh1.Multx(fLocField1, temp);
		fField_List[ip] = temp[0];
		fh2.Multx(fLocField2, temp);
		fField_List[ip] += temp[0];		

		fh1.Multx(fLocLastField1, temp);
		fField_last_List[ip] = temp[0];
		fh2.Multx(fLocLastField2, temp);
		fField_last_List[ip] += temp[0];

		fp1.Multx(fLocField1, temp);
		fGradField_List[ip] = temp[0];
		fp2.Multx(fLocField2, temp);
		fGradField_List[ip] += temp[0];		

		fp1.Multx(fLocLastField1, temp);
		fGradField_last_List[ip] = temp[0];
		fp2.Multx(fLocLastField2, temp);
		fGradField_last_List[ip] += temp[0];

		fq1.Multx(fLocField1, temp);
		fLapField_List[ip] = temp[0];
		fq2.Multx(fLocField2, temp);
		fLapField_List[ip] += temp[0];		

		fq1.Multx(fLocLastField1, temp);
		fLapField_last_List[ip] = temp[0];
		fq2.Multx(fLocLastField2, temp);
		fLapField_last_List[ip] += temp[0];
	}
	/********DEBUG*******/
	if (print_GlobalShape)
	{
		for (int nd_dof = 0; nd_dof < NumElementNodes_Field()*fNumDOF_Field1; nd_dof++)
			cout << "                 fLocField1[" << nd_dof << "]  : " << fLocField1[nd_dof] << endl;
		for (int nd_dof = 0; nd_dof < NumElementNodes_Field()*fNumDOF_Field2; nd_dof++)
			cout << "                 fLocField2[" << nd_dof << "]  : " << fLocField2[nd_dof] << endl;
		cout << endl;
		for (int ip = 0; ip < NumIP(); ip++)
			cout << "                 fField_List[" << ip << "]: " << fField_List[ip] << endl;
		for (int ip = 0; ip < NumIP(); ip++)
			cout << "                 fGradField_List[" << ip << "]: " << fGradField_List[ip] << endl;
		for (int ip = 0; ip < NumIP(); ip++)
			cout << "                 fLapField_List[" << ip << "]: " << fLapField_List[ip] << endl;
	}
	/*******************/
}

/* form the element stiffness matrix */
void GradSmallStrainT::FormStiffness(double constK)
{
	const char caller[] = "GradSmallStrainT::FormStiffness";
	int curr_group = ElementSupport().CurrentGroup();

	int fNumElementNodes_Disp = NumElementNodes();
	int fNumElementNodes_Field = NumElementNodes_Field();
	
	/* dmatrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole:
		dMatrixT::kUpperOnly;

	if(format != dMatrixT::kWhole)
		ExceptionT::BadInputValue(caller, "LHS must be NonSymmetric");

	/* integration rules (same for both fields)*/
	const double* Det = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	if (print_Stiffness)
	{
		cout << caller << endl;
		cout << "  element: " << CurrElementNumber() << endl;
	}
	/*******************/

	/* initialize */
	fK_bb = 0.0;
	fK_bh = 0.0;
	fK_hb = 0.0;
	fK_hh = 0.0;
	fK_hp = 0.0;
	fK_hq = 0.0;
	fK_ct = 0.0;

	fShapes->TopIP();

	/* integrate stiffness over the elements */
	while(fShapes->NextIP())
	{
		/* integration weight */
		double scale = constK*(*Det++)*(*Weight++);

		/* accumulate strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanGradient, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* accumulate shape functions */
		Set_h(fh);
		Set_p(fp);
		Set_q(fq);
		
		/* compute elastic stiffness matrix */
		fDM_bb.SetToScaled(scale, fCurrMaterial_Grad->dm_bb_ijkl());
		fK_bb.MultQTBQ(fB, fDM_bb, format, dMatrixT::kAccumulate);

		/* compute off-diagonal matrics */
		fOM_bh.SetToScaled(scale, fCurrMaterial_Grad->om_bh_ij());
		fK_bh.MultATBC(fB, fOM_bh, fh, format, dMatrixT::kAccumulate);

		fK_hb.Transpose(fK_bh);
		//		fOM_hb.SetToScaled(scale, fCurrMaterial_Grad->om_hb_ij());
		//		fK_hb.MultATBC(fh, fOM_hb, fB, format, dMatrixT::kAccumulate);

		/* compute non-symmetric, gradient dependent matrix */
		fGM_hh.SetToScaled(scale, fCurrMaterial_Grad->gm_hh());
		fK_hh.MultQTBQ(fh, fGM_hh, format, dMatrixT::kAccumulate);

		fGM_hp.SetToScaled(scale, fCurrMaterial_Grad->gm_hp());
		fK_hp.MultATBC(fh, fGM_hp, fp, format, dMatrixT::kAccumulate);

		fGM_hq.SetToScaled(scale, fCurrMaterial_Grad->gm_hq());
		fK_hq.MultATBC(fh, fGM_hq, fq, format, dMatrixT::kAccumulate);

		/* obtain yield condition residual */
		fYield_List[CurrIP()] = fCurrMaterial_Grad->yc();

		/********DEBUG*******/
		if (print_Stiffness)
		{
			cout<<"            ip: "             << CurrIP()                                           <<endl;
			cout<<"                 Field    : " << fField_List[CurrIP()]                              <<endl;
			cout<<"                 del_Field: " << fField_List[CurrIP()] - fField_last_List[CurrIP()] <<endl;
			cout<<"                 strain   : " << (fStrain_List[CurrIP()])[0]                        <<endl;
			cout<<"                 stress   : " << (fCurrMaterial_Grad->s_ij())[0]                    <<endl;
			cout<<"                 yc       : " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"                 om_bh_ij : " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"                 om_hb_ij : " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"                 gm_hh    : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp    : " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			cout<<"                 gm_hq    : " << fCurrMaterial_Grad->gm_hq()                        <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 h[" << i << "]     : " << fh[i]                                <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 p[" << i << "]     : " << fp[i]                                <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 q[" << i << "]     : " << fq[i]                                <<endl;
			for (int i = 0; i < NumElementNodes(); i++)
				cout<<"                 B[" << i << "]     : " << fB[i]                                <<endl;
			cout<<"                 scale    : " << scale                                              <<endl;
		}
		/*******************/
	}

	/* assemble into element stiffness matrix */
	fLHS.AddBlock(0,            0,            fK_bb);
	fLHS.AddBlock(0,            fK_bb.Cols(), fK_bh);
	fLHS.AddBlock(fK_bb.Rows(), 0,            fK_hb);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hh);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hp);
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_hq);

	dArrayT fLocYield(NumElementNodes());
	IP_ExtrapolateAll(fYield_List, fLocYield);
	
	/* add constraint to Krr if elastic */
	fK_ct = 0.;
	for (int nd = 0; nd < NumElementNodes_Field(); nd ++)
		if (fLocField1[nd] <= 0. || fLocField1[nd] <= fLocLastField1[nd]) // || fLocYield[nd] <= 0.)
			fK_ct(nd)[nd] = fNodalConstraint;
	fLHS.AddBlock(fK_bb.Rows(), fK_bb.Cols(), fK_ct);

	/********DEBUG*******/
	if (print_StiffnessMatrix)
	{
		cout << caller << " Stiffness Matrix" << endl;
		cout << "  element: " << CurrElementNumber() << endl;
		for (int i=0; i < fNumElementNodes_Disp * fNumDOF_Disp; i++)
		{
			cout<<"                 K("<< i << ",j): ";
			for (int j=0; j < fNumElementNodes_Disp * fNumDOF_Disp; j++)
				cout << fK_bb(j)[i] << "     ";
			for (int j=0; j < fNumElementNodes_Field * fNumDOF_Field; j++)
				cout << fK_bh(j)[i] << "         ";
			cout << endl;
		}
		for (int i=0; i < fNumElementNodes_Field * fNumDOF_Field; i++)
		{
			cout<<"                 K("<< i << ",j): ";
			for (int j=0; j < fNumElementNodes_Disp * fNumDOF_Disp; j++)
				cout << fK_hb(j)[i] << "     ";
			for (int j=0; j < fNumElementNodes_Field * fNumDOF_Field; j++)
			{
				cout << fK_hh(j)[i] << "+";
				if (fK_ct(j)[i] != 0)
					cout << "cst" << "     ";
				else
					cout << "0.0" << "     ";
			}
			cout << endl;
		}
	}
	/*******************/
}

/* calculate the internal force contribution ("-k*d") */
void GradSmallStrainT::FormKd(double constK)
{
	const char caller[] = "GradSmallStrainT::FormKd";
	int curr_group = ElementSupport().CurrentGroup();

	int fNumElementNodes_Disp = NumElementNodes();
	int fNumElementNodes_Field = NumElementNodes_Field();
	
	/* partition residual force vector */
	int neq_Disp = fNumElementNodes_Disp*fNumDOF_Disp;
	int neq_Field1 = fNumElementNodes_Field*fNumDOF_Field1;
	int neq_Field2 = fNumElementNodes_Field*fNumDOF_Field2;	

	dArrayT RHS_Disp(neq_Disp, fRHS.Pointer());
	dArrayT RHS_Field1(neq_Field1, fRHS.Pointer(neq_Disp));
	dArrayT RHS_Field2(neq_Field2, fRHS.Pointer(neq_Disp+neq_Field1));	

	/* integration rules (same for both fields)*/
	const double* Det = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/********DEBUG*******/
	if (print_Kd)
	{
		cout << caller << endl;
		cout<<"	 element: "<<CurrElementNumber()<<endl;
	}
	/*******************/

	/* initialize */
	fShapes->TopIP();

	/* integrate residuals over the elements */
	while(fShapes->NextIP())
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
		Set_q(fq);

		for (int i=0; i<2; i++)
		{
			fh1[i] = fh[i];
			fh2[i] = fh[i+2];
		}
		
		/* obtain yield condition residual */
		double yield = fCurrMaterial_Grad->yc();

		/* accumulate in rhs field equations */
		RHS_Field1.AddScaled(-scale*yield, fh1T.Transpose(fh1));
		RHS_Field2.AddScaled(-scale*yield, fh2T.Transpose(fh2));		

		/********DEBUG*******/
		if (print_Kd)
		{
			cout<<"            ip: "             << CurrIP()                                           <<endl;
			cout<<"                 Field    : " << fField_List[CurrIP()]                              <<endl;
			cout<<"                 del_Field: " << fField_List[CurrIP()] - fField_last_List[CurrIP()] <<endl;
			cout<<"                 strain   : " << (fStrain_List[CurrIP()])[0]                        <<endl;
			cout<<"                 stress   : " << (fCurrMaterial_Grad->s_ij())[0]                    <<endl;
			cout<<"                 yc       : " << fCurrMaterial_Grad->yc()                           <<endl;
			cout<<"                 om_bh_ij : " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"                 om_hb_ij : " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"                 gm_hh    : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp    : " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 h[" << i << "]     : " << fh[i]                                <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 p[" << i << "]     : " << fp[i]                                <<endl;
			for (int i = 0; i < NumElementNodes_Field()*NumDOF_Field(); i++)
				cout<<"                 q[" << i << "]     : " << fq[i]                                <<endl;
			for (int i = 0; i < NumElementNodes(); i++)
				cout<<"                 B[" << i << "]     : " << fB[i]                                <<endl;
			cout<<"                 scale    : " << scale                                              <<endl;
			for (int i=0; i < neq_Disp; i++)
				cout<<"                 RHS_Disp["<< i << "] : " << scale*fNEEvec[i] << endl;
			for (int i=0; i < neq_Field1; i++)
				cout<<"                 RHS_Field1["<< i << "]: " << -scale*yield*fh1T[i] << endl;
			for (int i=0; i < neq_Field2; i++)
				cout<<"                 RHS_Field2["<< i << "]: " << -scale*yield*fh2T[i] << endl;
		}
		/*******************/
	}
	/********DEBUG*******/
	if (print_KdMatrix)
	{
		cout << caller << " Kd Matrix" << endl;
		cout << "  element: " << CurrElementNumber() << endl;
		for (int i=0; i < neq_Disp; i++)
			cout<<"                 RHS_Disp["<< i << "]: " << RHS_Disp[i] << endl;
		for (int i=0; i < neq_Field1; i++)
			cout<<"                 RHS_Field1["<< i << "]: " << RHS_Field1[i] << endl;
		for (int i=0; i < neq_Field2; i++)
			cout<<"                 RHS_Field2["<< i << "]: " << RHS_Field2[i] << endl;
		cout << endl;
	}
	/*******************/
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

	out << " Associated field 1. . . . . . . . . . . . . . . . = \"" << fField1.Name() << "\"\n";
	out << " Associated field 2. . . . . . . . . . . . . . . . = \"" << fField2.Name() << "\"\n";	
	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode_Field << '\n';
	out << "	eq." << GeometryT::kLine	<< ", line\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP_Field	<< '\n';
}

/***********************************************************************
* Private
***********************************************************************/

/* accumulate shape function matrix */
void GradSmallStrainT::Set_h(dMatrixT& h) const
{
	const char caller[] = "GradSmallStrainT::Set_h";

	int  fNumElementNodes_Field = NumElementNodes_Field();
	const double* pShapes_Field = fShapes_Field->IPShapeU(fShapes->CurrIP());
	double* ph = h.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
			*ph++ = *pShapes_Field++;
	}
	else
		ExceptionT::GeneralFail(caller, "implemented for 1d only");
}

/* accumulate shape function matrix */
void GradSmallStrainT::Set_p(dMatrixT& p) const
{
	const char caller[] = "GradSmallStrainT::Set_p";

	int  fNumElementNodes_Field = NumElementNodes_Field();
	double* pp = p.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		const double* pNax = fShapes_Field->IPDShapeU(fShapes->CurrIP())(0);
		for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
			*pp++ = *pNax++;
	}
	else
		ExceptionT::GeneralFail(caller, "implemented for 1d only");
}

/* accumulate Laplacian C1 shape function matrix */
void GradSmallStrainT::Set_q(dMatrixT& q) const
{
	const char caller[] = "GradSmallStrainT::Set_q";

	int  fNumElementNodes_Field = NumElementNodes_Field();
	const double* pShapes_LapField = fShapes_Field->IPDDShapeU(fShapes->CurrIP());
	double* pq = q.Pointer();

	/* 1D */
	if (fNumSD == 1)
	{
		for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
			*pq++ = *pShapes_LapField++;
	}
	else
		ExceptionT::GeneralFail(caller, "implemented for 1d only");
}
