/* $Id: TorsionKBCT.cpp,v 1.3 2003-08-18 03:45:17 paklein Exp $ */
#include "TorsionKBCT.h"
#include "NodeManagerT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* parameters */
const double Pi = acos(-1.0);

/* vector functions */
inline static void CrossProduct(const dArrayT& A, const dArrayT& B, dArrayT& AxB)
{
	AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

/* constructor */
TorsionKBCT::TorsionKBCT(NodeManagerT& node_manager, const double& time):
	KBC_ControllerT(node_manager),
	fTime(time),
	fStartTime(0.0),
	fw(0.0),
	fAxis(-1),
	fDummySchedule(1.0)
{
	SetName("torsion");
}

/* initialize data - called immediately after construction */
void TorsionKBCT::Initialize(ifstreamT& in)
{
	const char caller[] = "TorsionKBCT::Initialize";

	/* 3D only */
	if (fNodeManager.NumSD() != 3)
		ExceptionT::BadInputValue(caller, "3D only");

	/* rotation rate */
	in >> fw;
	
	/* direction of the axis of rotation */
	in >> fAxis;
	if (fAxis > 3 || fAxis < 1)
		ExceptionT::BadInputValue(caller, "axis is out of range: %d", fAxis);
	fAxis--;
	
	/* point on the axis of rotation */
	fPoint.Dimension(fNodeManager.NumSD());
	in >> fPoint;

	/* nodes */
	ReadNodes(in, fID_List, fNodes);

	/* constrained directions */
	int constrained_dirs[3][2] = {
		{1,2},
		{2,0},
		{0,1}};
	int* dir = constrained_dirs[fAxis];

	/* generate BC cards */
	int n_cards = 2;
	fKBC_Cards.Dimension(fNodes.Length()*n_cards);
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < fNodes.Length(); i++)
		for (int j = 0; j < n_cards; j++)
		{
			/* set values */
			pcard->SetValues(fNodes[i], dir[j], KBC_CardT::kDsp, 0, 0.0);
	
			/* dummy schedule */
			pcard->SetSchedule(&fDummySchedule);
			pcard++;
		}	
}

/* set to initial conditions */
void TorsionKBCT::InitialCondition(void)
{
	/* store start time */
	fStartTime = fTime;
}

void TorsionKBCT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

	out << "\n T o r s i o n   p a r a m e t e r s :\n\n";
	out << " Rotation rate . . . . . . . . . . . . . . . . . = " << fw << '\n';
	out << " Rotation axis direction (1:x, 2:y, 3:z) . . . . = " << fAxis+1 << '\n';
	out << " Point on the axis of rotation:\n";
	out << fPoint << '\n';
	out << " Number of group nodes . . . . . . . . . . . . . = " << fNodes.Length() << '\n';	
	iArrayT tmp;
	tmp.Alias(fNodes);
	tmp++;
	out << tmp.wrap(5) << '\n';
	tmp--;
}

/* initialize/finalize/reset step */
void TorsionKBCT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();

	/* rotation axes */
	double direction[3][3] = {
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0}};
	dArrayT axis(3, direction[fAxis]);
	
	/* work space */
	dArrayT v_op(3), R(3), xl(3), yl(3);
	dArrayT x(3), X(3), c(3);

	/* coordinates */
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();

	/* compute point by point */
	int dex = 0;
	double theta = (fTime-fStartTime)*fw;
	for (int i = 0; i < fNodes.Length(); i++)
	{
		/* node */
		int node = fNodes[i];

		/* coordinates */
		init_coords.RowAlias(node, X);

		/* from axis to point */
		v_op.DiffOf(X, fPoint);
		
		/* radial vector */
		R.SetToCombination(1.0, v_op, -dArrayT::Dot(v_op, axis), axis);
		double r = R.Magnitude();
		
		/* plane of rotation */
		c.DiffOf(X, R);
		xl.SetToScaled(1.0/r, R);
		CrossProduct(axis, xl, yl);
	
		/* new position */
		x.SetToCombination(1.0, c, r*cos(theta), xl, r*sin(theta), yl);
		
		KBC_CardT& card_1 = fKBC_Cards[dex++];
		int dof_1 = card_1.DOF();
		card_1.SetValues(node, dof_1, KBC_CardT::kDsp, 0, x[dof_1] - X[dof_1]);

		KBC_CardT& card_2 = fKBC_Cards[dex++];
		int dof_2 = card_2.DOF();
		card_2.SetValues(node, dof_2, KBC_CardT::kDsp, 0, x[dof_2] - X[dof_2]);
	}
}

/* describe the parameters needed by the interface */
void TorsionKBCT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	list.AddParameter(fw, "rotation_rate");
	list.AddParameter(fAxis, "rotation_axis");
}
