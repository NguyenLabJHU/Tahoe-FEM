/* $Id: TorsionKBCT.cpp,v 1.5 2004-07-22 21:07:49 paklein Exp $ */
#include "TorsionKBCT.h"
#include "NodeManagerT.h"

#include "ParameterUtils.h"

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
TorsionKBCT::TorsionKBCT(const BasicSupportT& support):
	KBC_ControllerT(support),
	fStartTime(0.0),
	fw(0.0),
	fAxis(-1),
	fDummySchedule(1.0)
{
	SetName("torsion");
}

/* set to initial conditions */
void TorsionKBCT::InitialCondition(void)
{
	/* store start time */
	fStartTime = fSupport.Time();
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
	const dArray2DT& init_coords = fSupport.InitialCoordinates();

	/* compute point by point */
	int dex = 0;
	double theta = (fSupport.Time() - fStartTime)*fw;
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
		
		/* compute new position */
		if (fabs(r) > kSmall) 
		{
			/* plane of rotation */
			c.DiffOf(X, R);
			xl.SetToScaled(1.0/r, R);
			CrossProduct(axis, xl, yl);
	
			/* new position */
			x.SetToCombination(1.0, c, r*cos(theta), xl, r*sin(theta), yl);
		}
		else /* motion */
			x = X;

		/* set prescribed motion */
		KBC_CardT& card_1 = fKBC_Cards[dex++];
		int dof_1 = card_1.DOF();
		card_1.SetValues(node, dof_1, KBC_CardT::kDsp, NULL, x[dof_1] - X[dof_1]);

		KBC_CardT& card_2 = fKBC_Cards[dex++];
		int dof_2 = card_2.DOF();
		card_2.SetValues(node, dof_2, KBC_CardT::kDsp, NULL, x[dof_2] - X[dof_2]);
	}
}

/* describe the parameters needed by the interface */
void TorsionKBCT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	list.AddParameter(fw, "rotation_rate");

	ParameterT axis(ParameterT::Integer, "rotation_axis");
	axis.AddLimit(1, LimitT::Only);
	axis.AddLimit(2, LimitT::Only);
	axis.AddLimit(3, LimitT::Only);
	list.AddParameter(axis);
}

/* information about subordinate parameter lists */
void TorsionKBCT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* point on the axis of rotation */
	sub_list.AddSub("point_on_axis");

	/* list of nodes to control */
	sub_list.AddSub("node_ID_list");
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* TorsionKBCT::NewSub(const StringT& name) const
{
	if (name == "point_on_axis")
		return new VectorParameterT(name, 3, 'x');
	else /* inherited */	
		return KBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void TorsionKBCT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "TorsionKBCT::TakeParameterList";

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* 3D only */
	if (fSupport.NumSD() != 3) ExceptionT::BadInputValue(caller, "3D only");

	/* rotation rate and axis */
	fw = list.GetParameter("rotation_rate");
	fAxis = list.GetParameter("rotation_axis");
	fAxis--;

	/* point on the axis of rotation */
	VectorParameterT vec("point_on_axis", 3, 'x');
	vec.TakeParameterList(list.GetList(vec.Name()));
	fPoint = vec;

	/* nodes */
	StringListT::Extract(list.GetList("node_ID_list"),  fID_List);	
	GetNodes(fID_List, fNodes);

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
			pcard->SetValues(fNodes[i], dir[j], KBC_CardT::kDsp, &fDummySchedule, 0.0);
			pcard++;
		}	

#if 0
//TEMP - write nodes???
	out << " Number of group nodes . . . . . . . . . . . . . = " << fNodes.Length() << '\n';	
	iArrayT tmp;
	tmp.Alias(fNodes);
	tmp++;
	out << tmp.wrap(5) << '\n';
	tmp--;
#endif
}
