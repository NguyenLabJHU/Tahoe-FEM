/* $id: GradSmallStrainT.cpp,v 1.10 2004/06/17 20:45:13 rdorgan Exp $ */ 
#include "GradSmallStrainT.h"

/* shape functions */
#include "ShapeFunctionT.h"
#include "ShapeTools.h"
#include "C0ShapeTools.h"
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

/* parameters */
const double kYieldTol    = 1.0e-10;

/* constructor */
GradSmallStrainT::GradSmallStrainT(const ElementSupportT& support, 
									   const FieldT& disp, const FieldT& field):
	SmallStrainT(support), //pass the displacement field to the base class
	fDisplacement(disp),
	fField(field),
	
	fLocField(LocalArrayT::kDisp),
	fLocLastField(LocalArrayT::kLastDisp),

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

	fNumSD(NumSD()),
	fNumIP_Disp(NumIP()),
	fNumIP_Field(fNumIP_Disp),
	fNumDOF_Disp(fDisplacement.NumDOF()),
	fNumDOF_Field(fField.NumDOF()),
	
	fI(1)
{
ExceptionT::GeneralFail("GradSmallStrainT::GradSmallStrainT", "out of date");
#if 0
	const char caller[] = "GradSmallStrainT::GradSmallStrainT";

	ifstreamT& in = ElementSupport().Input();

	/* control parameters for field mesh */
	in >> fNumElementNodes_Field;
	in >> fNodalConstraint;
 
	/********DEBUG*******/
	in >> print_GlobalShape;
	in >> print_Kd;
	in >> print_KdMatrix;
	in >> print_Stiffness;
	in >> print_StiffnessMatrix;
	/*******************/

	/* plastic multiplier element type */
	if (fNumSD == 1) {
		if (fNumDOF_Field == 1)              // lambda
			fDegreeOfContinuity_Field = C0;
		else if (fNumDOF_Field == 2) {       // lambda, lambda_x
			if (fNumElementNodes_Field == 2) // vertex nodes only
				fDegreeOfContinuity_Field = C1;
			else
				ExceptionT::BadInputValue(caller, "fNumElementNodes_Field != 2 (fNumSD == 1, fNumDOF_Field == 2)");
		}
		else
			ExceptionT::BadInputValue(caller, "fNumDOF_Field != 1 && fNumDOF_Field != 2 (fNumSD == 1)");
	}	
	else
		ExceptionT::BadInputValue(caller, "fNumSD != 1");
#endif
}

/* destructor */
GradSmallStrainT::~GradSmallStrainT(void) 
{  
	delete fShapes_Field;
}


void GradSmallStrainT::Initialize(void)
{
ExceptionT::GeneralFail("GradSmallStrainT::Initialize", "out of date");
#if 0	
	const char caller[] = "GradSmallStrainT::Initialize";

	/* inherited */
	SmallStrainT::Initialize();

	/* element dimensions */
	fNumElementNodes_Disp = NumElementNodes();

	/* allocate lists */
	fField_List.Dimension(fNumIP_Field);           // field
	fField_last_List.Dimension(fNumIP_Field);      // "last" field

	fGradField_List.Dimension(fNumIP_Field);       // gradient field
	fGradField_last_List.Dimension(fNumIP_Field);  // "last" gradient field

	fLapField_List.Dimension(fNumIP_Field);        // Laplacian field
	fLapField_last_List.Dimension(fNumIP_Field);   // "last" Laplacian field

	fYield_List.Dimension(fNumIP_Disp);            // yield condition field
	
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

	/* allocate gradient of shape functions */
	fp.Dimension (fNumSD, fNumElementNodes_Field*fNumDOF_Field);

	/* allocate Laplacian shape functions */
	fq.Dimension (1, fNumElementNodes_Field*fNumDOF_Field);
#endif
}


void GradSmallStrainT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
								 AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)
	const char caller[] = "GradSmallStrainT::Equations";

	/* loop over connectivity blocks */
	for (int i = 0; i < fEqnos.Length(); i++)
	{
		/* connectivities */
		const iArray2DT& connects = *(fConnectivities[i]);
		const iArray2DT& connects_field = fConnectivities_Field[i];
		int fNumElements = connects.MajorDim();

		/* dimension */
		fEqnos[i].Dimension   (fNumElements, fNumEQ_Total);
		iArray2DT fEqnos_Disp (fNumElements, fNumElementNodes_Disp *fNumDOF_Disp );
		iArray2DT fEqnos_Field(fNumElements, fNumElementNodes_Field*fNumDOF_Field);

		/* get equation numbers */
		fDisplacement.SetLocalEqnos(connects, fEqnos_Disp);
		fField.SetLocalEqnos(connects_field, fEqnos_Field);

		/* write into one array */
		fEqnos[i].BlockColumnCopyAt(fEqnos_Disp, 0);
		fEqnos[i].BlockColumnCopyAt(fEqnos_Field, fEqnos_Disp.MinorDim());

		/* add to list of equation numbers */
		eq_1.Append(&fEqnos[i]);
	}

	/* reset pointers to element cards */
	SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
}

/***********************************************************************
* Protected
***********************************************************************/

/* echo element connectivity data. Calls ElementBaseT::ReadConnectivity */
void GradSmallStrainT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	const char caller[] = "GradSmallStrainT::EchoConnectivityData";

	/* initialize displacement connectivities */
	SmallStrainT::EchoConnectivityData(in, out);

	/* depending whether you need separate connectivities */
	bool need_separate_connectivities = true;
	//bool need_separate_connectivities = fNumElementNodes_Disp - fNumElementNodes_Field;
	
	fConnectivities_Field.Dimension(fConnectivities.Length());
	if (need_separate_connectivities) {
	
		//number of vertex nodes (fNumSD == 1)
		int num_vertex_nodes = 2;

		//list of the vertex nodes
		iArrayT vertex_nodes(num_vertex_nodes);
		vertex_nodes[0] = 0;
		vertex_nodes[1] = 1;

		fConnectivities_All.Dimension(NumElements(), num_vertex_nodes);
	
		/* translate blocks */
		int count = 0;
		for (int i = 0; i < fConnectivities.Length(); i++)
		{
			const iArray2DT& connects = *(fConnectivities[i]);
			iArray2DT& connects_field = fConnectivities_Field[i];
			connects_field.Alias(connects.MajorDim(), num_vertex_nodes, fConnectivities_All(count));
		
			/* extract */
			for (int j = 0; j < vertex_nodes.Length(); j++)
				connects_field.ColumnCopy(j, connects, vertex_nodes[j]);
				
			/* next block */
			count += connects.MajorDim();
		}
	}
	else /* same connectivities - make aliases */
	{
		for (int i = 0; i < fConnectivities.Length(); i++)
			fConnectivities_Field[i].Alias(*(fConnectivities[i]));
	}
}

/* construct a new material support and return a pointer */
MaterialSupportT* GradSmallStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	const char caller[] = "GradSmallStrainT::NewMaterialSupport";

	/* allocate new */
	if (!p) p = new GradSSMatSupportT(fNumSD, fNumDOF_Disp, fNumDOF_Field, fNumIP_Disp, fNumIP_Field);

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
	const char caller[] = "GradSmallStrainT::SetLocalArrays";

	/* inherited */
	SmallStrainT::SetLocalArrays();

	/* dimension local arrays */
	fLocField.Dimension(fNumElementNodes_Field, fNumDOF_Field);
	fLocLastField.Dimension(fNumElementNodes_Field, fNumDOF_Field);

	fLocFieldTranspose.Dimension(fLocField.Length());
	
	/* register local arrays */
	fField.RegisterLocal(fLocField);
	fField.RegisterLocal(fLocLastField);
}

/* initialization functions */
void GradSmallStrainT::SetShape(void)
{
	const char caller[] = "GradSmallStrainT::SetShape";

	/* inherited */
	SmallStrainT::SetShape();

	
	/* select shape function tools to use */
	switch(fDegreeOfContinuity_Field)
	{
		case C0:
		{
			fShapes_Field = new C0ShapeTools(GeometryCode(), fNumIP_Disp, fNumElementNodes_Field, fNumDOF_Field, fLocInitCoords);
			break;
		}
		case C1:
		{
			fShapes_Field = new C1ShapeTools(GeometryCode(), fNumIP_Disp, fNumElementNodes_Field, fNumDOF_Field, fLocInitCoords);
			break;
		}
		default:
			ExceptionT::GeneralFail(caller, "Bad fDegreeOfContinuity_Field");
	}

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
	fConnectivities_All.RowAlias(CurrElementNumber(), fNodesField);
	fLocField.SetLocal(fNodesField);
	fLocLastField.SetLocal(fNodesField);

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

		/* compute current fields using field shape functions */
		fLocField.ReturnTranspose(fLocFieldTranspose);

		fh.Multx(fLocFieldTranspose, temp);
		fField_List[ip] = temp[0];

		fp.Multx(fLocFieldTranspose, temp);
		fGradField_List[ip] = temp[0];

		fq.Multx(fLocFieldTranspose, temp);
		fLapField_List[ip] = temp[0];

		/* compute "last" fields using field shape functions */
		fLocLastField.ReturnTranspose(fLocFieldTranspose);

		fh.Multx(fLocFieldTranspose, temp);
		fField_last_List[ip] = temp[0];

		fp.Multx(fLocFieldTranspose, temp);
		fGradField_last_List[ip] = temp[0];

		fq.Multx(fLocFieldTranspose, temp);
		fLapField_last_List[ip] = temp[0];
	}

	/********DEBUG*******/
	if (print_GlobalShape)
	{
		for (int nd_dof = 0; nd_dof < fNumElementNodes_Field*fNumDOF_Field; nd_dof++)
			cout << "                 fLocField[" << nd_dof << "]  : " << fLocField[nd_dof] << " * " << fh[nd_dof] << endl;
		cout << endl;
		for (int ip = 0; ip < fNumIP_Disp; ip++)
			cout << "                 fField_List[" << ip << "]: " << fField_List[ip] << endl;
		for (int ip = 0; ip < fNumIP_Disp; ip++)
			cout << "                 fGradField_List[" << ip << "]: " << fGradField_List[ip] << endl;
		for (int ip = 0; ip < fNumIP_Disp; ip++)
			cout << "                 fLapField_List[" << ip << "]: " << fLapField_List[ip] << endl;
	}
	/*******************/
}

/* form the element stiffness matrix */
void GradSmallStrainT::FormStiffness(double constK)
{
	const char caller[] = "GradSmallStrainT::FormStiffness";
	int curr_group = ElementSupport().CurrentGroup();

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

		fOM_hb.SetToScaled(scale, fCurrMaterial_Grad->om_hb_ij());
		fK_hb.MultATBC(fh, fOM_hb, fB, format, dMatrixT::kAccumulate);

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
			cout<<"                 om_bh_ij : " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"                 om_hb_ij : " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"                 gm_hh    : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp    : " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			cout<<"                 gm_hq    : " << fCurrMaterial_Grad->gm_hq()                        <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 h[" << i << "]     : " << fh[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 p[" << i << "]     : " << fp[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 q[" << i << "]     : " << fq[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Disp; i++)
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

	dArrayT fLocYield(fNumElementNodes_Disp);
	IP_ExtrapolateAll(fYield_List, fLocYield);
	
	/* add constraint to Krr if elastic */
	fK_ct = 0.;

	for (int nd = 0; nd < fNumElementNodes_Field ; nd++)

		// BEST CONVERGENCE / BETTER TOLERANCE CONVERGENCE
		if (fLocField[nd] <= 0. || fLocField[nd] <= fLocLastField[nd])

		// 20ELEMENT KILLS EARLY / OTHER SOLUTIONS TAKE LONGER WITH WORSE TOLERANCE CONVERGENCE
		//		if (fLocField1[nd] <= 0. || fLocField1[nd] <= fLocLastField1[nd] || fabs(fLocYield[nd]) < kYieldTol) 

		// TAKES LONGER / WORSE TOLERANCE CONVERGENCE
		//		if (fabs(fLocYield[nd]) < kYieldTol)

			fK_ct(nd*fNumDOF_Field)[nd*fNumDOF_Field] = fNodalConstraint;

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

	/* partition residual force vector */
	int neq_Disp = fNumElementNodes_Disp*fNumDOF_Disp;
	int neq_Field = fNumElementNodes_Field*fNumDOF_Field;

	dArrayT RHS_Disp(neq_Disp, fRHS.Pointer());
	dArrayT RHS_Field(neq_Field, fRHS.Pointer(neq_Disp));

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

		/* obtain yield condition residual */
		double yield = fCurrMaterial_Grad->yc();

		/* accumulate in rhs field equations */
		RHS_Field.AddScaled(-scale*yield, fhT.Transpose(fh));

		/********DEBUG*******/
		if (print_Kd)
		{
			cout<<"            ip: "             << CurrIP()                                           <<endl;
			cout<<"                 om_bh_ij : " << (fCurrMaterial_Grad->om_hb_ij())[0]                <<endl;
			cout<<"                 om_hb_ij : " << (fCurrMaterial_Grad->om_bh_ij())[0]                <<endl;
			cout<<"                 gm_hh    : " << fCurrMaterial_Grad->gm_hh()                        <<endl;
			cout<<"                 gm_hp    : " << fCurrMaterial_Grad->gm_hp()                        <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 h[" << i << "]     : " << fh[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 p[" << i << "]     : " << fp[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Field*fNumDOF_Field; i++)
				cout<<"                 q[" << i << "]     : " << fq[i]                                <<endl;
			for (int i = 0; i < fNumElementNodes_Disp; i++)
				cout<<"                 B[" << i << "]     : " << fB[i]                                <<endl;
			cout<<"                 scale    : " << scale                                              <<endl;
			for (int i=0; i < neq_Disp; i++)
				cout<<"                 RHS_Disp["<< i << "] : " << scale*fNEEvec[i] << endl;
			for (int i=0; i < neq_Field; i++)
				cout<<"                 RHS_Field["<< i << "]: " << -scale*yield*fhT[i] << endl;
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
		for (int i=0; i < neq_Field; i++)
			cout<<"                 RHS_Field["<< i << "]: " << RHS_Field[i] << endl;
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

#if 0
/* print element group data */
void GradSmallStrainT::PrintControlData(ostream& out) const
{
	/* inherited */
	SmallStrainT::PrintControlData(out);

	out << " Associated field. . . . . . . . . . . . . . . . = \"" << fField.Name() << "\"\n";

	out << " Element degree of continuity. . . . . . . . . . = " << fDegreeOfContinuity_Field << '\n';
	out << " Number of element nodes . . . . . . . . . . . . = " << fNumElementNodes_Field    << '\n';
	out << " Number of integration points. . . . . . . . . . = " << fNumIP_Field              << '\n';
	out << " Number of degrees of freedom. . . . . . . . . . = " << fNumDOF_Field             << '\n';
	out << " Nodal constraint. . . . . . . . . . . . . . . . = " << fNodalConstraint          << '\n';
}
#endif

/***********************************************************************
* Private
***********************************************************************/

/* accumulate shape function matrix */
void GradSmallStrainT::Set_h(dMatrixT& h) const
{
	const char caller[] = "GradSmallStrainT::Set_h";

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
