/* $Id: SolidMaterialT.cpp,v 1.10.2.5 2004-03-04 06:45:37 paklein Exp $ */
/* created: paklein (11/20/1996) */
#include "SolidMaterialT.h"

#include <iostream.h>

#include "fstreamT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "LocalArrayT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
SolidMaterialT::SolidMaterialT(ifstreamT& in, const MaterialSupportT& support):
	ParameterInterfaceT("solid_material"),
	ContinuumMaterialT(support)
{
#pragma unused(in)
#if 0
	in >> fMassDamp;	if (fMassDamp  <  0.0) throw ExceptionT::kBadInputValue;
	in >> fStiffDamp;	if (fStiffDamp <  0.0) throw ExceptionT::kBadInputValue;
	in >> fDensity;		if (fDensity   <= 0.0) throw ExceptionT::kBadInputValue;
	fThermal = new ThermalDilatationT(in);
	if (!fThermal) throw ExceptionT::kOutOfMemory;

	SetName("solid_material");
#endif
}

SolidMaterialT::SolidMaterialT(void):
	ParameterInterfaceT("solid_material"),
	fThermal(NULL),
	fDensity(0.0),
	fConstraint(kNoConstraint),
	fMassDamp(0.0),
	fStiffDamp(0.0)
{

}

/* destructor */
SolidMaterialT::~SolidMaterialT(void) { delete fThermal; }

/* initialization */
void SolidMaterialT::Initialize(void)
{
	/* inherited */
	ContinuumMaterialT::Initialize();

	/* active multiplicative dilatation */
	if (fThermal->IsActive() && !SupportsThermalStrain())
	{
		cout << "\n SolidMaterialT::Initialize: material does not support\n"
		     <<   "     imposed thermal strain." << endl;
		throw ExceptionT::kBadInputValue;
	}
}

/* I/O functions */
void SolidMaterialT::Print(ostream& out) const
{
	/* inherited */
	ContinuumMaterialT::Print(out);
	
	out << " Mass damping coefficient. . . . . . . . . . . . = " << fMassDamp  << '\n';
	out << " Stiffness damping coefficient . . . . . . . . . = " << fStiffDamp << '\n';
	out << " Density . . . . . . . . . . . . . . . . . . . . = " << fDensity   << '\n';

//	fThermal->Print(out);
}

/* return the wave speeds */
void SolidMaterialT::WaveSpeeds(const dArrayT& normal, dArrayT& speeds)
{
#if __option(extended_errorcheck)
	if (normal.Length() != speeds.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* compute acoustical tensor */
	const dSymMatrixT& Q = AcousticalTensor(normal);

	/* get eigenvalues (sorted by magnitude) */
	Q.PrincipalValues(speeds);
	
	/* order results */
	if (speeds.Length() == 2)
	{
		/* order as {c_d, c_s} */
		double n_eig = Q.MultmBn(normal, normal);
		if (fabs(speeds[1] - n_eig) < fabs(speeds[0] - n_eig))
		{
			double temp = speeds[0];
			speeds[0] = speeds[1];
			speeds[1] = temp;
		}
		
		/* compute wave speeds */
		speeds[0] = (speeds[0] > 0.0) ? sqrt(speeds[0]/fDensity) : 0.0;
		speeds[1] = (speeds[1] > 0.0) ? sqrt(speeds[1]/fDensity) : 0.0;
	}
	else
	{
		/* order as {c_d, (c_s)_min, (c_s)_max} */
		double n_eig = Q.MultmBn(normal, normal);
		double d0 = fabs(speeds[0] - n_eig);
		double d1 = fabs(speeds[1] - n_eig);
		double d2 = fabs(speeds[2] - n_eig);

		int id, is1, is2;
		if (d0 < d1) {
			is1 = 1;			
			if (d0 < d2) {
				id  = 0;
				is2 = 2; }
			else {
				id  = 2;
				is2 = 0; } }
		else {
			is1 = 0;
			if (d1 < d2) {
				id  = 1;
				is2 = 2; }
			else {
				id  = 2;
				is2 = 1; } }
		
		/* sort */
		double temp[3];
		temp[0] = speeds[id];
		if (speeds[is1] < speeds[is2]) {
			temp[1] = speeds[is1];
			temp[2] = speeds[is2]; }
		else {
			temp[2] = speeds[is1];
			temp[1] = speeds[is2]; }
	
		/* compute wave speeds */
		speeds[0] = (temp[0] > 0.0) ? sqrt(temp[0]/fDensity) : 0.0;
		speeds[1] = (temp[1] > 0.0) ? sqrt(temp[1]/fDensity) : 0.0;
		speeds[2] = (temp[2] > 0.0) ? sqrt(temp[2]/fDensity) : 0.0;
	}
}

/* returns 1 if the strain localization conditions if satisfied,
* .ie if the acoustic tensor has zero (or negative eigenvalues),
* for the current conditions (current integration point and strain
* state). If localization is detected, the normal (current config)
* to the surface is returned in normal */
int SolidMaterialT::IsLocalized(dArrayT& normal)
{
#pragma unused(normal)

	/* by default, no localization */
	return 0;
}

/* describe the parameters needed by the interface */
void SolidMaterialT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ContinuumMaterialT::DefineParameters(list);

	/* density */
	ParameterT density(fDensity, "density");
	density.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(density);

	/* 2D constraint option */
	ParameterT constraint(ParameterT::Enumeration, "constraint_2D");
	constraint.AddEnumeration("none", kNoConstraint);
	constraint.AddEnumeration("plane_stress", kPlaneStress);
	constraint.AddEnumeration("plane_strain", kPlaneStrain);
	constraint.SetDefault(kNoConstraint);
	list.AddParameter(constraint);
}

/* information about subordinate parameter lists */
void SolidMaterialT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ContinuumMaterialT::DefineSubs(sub_list);
	
	/* thermal dilatation */
	sub_list.AddSub("thermal_dilatation", ParameterListT::ZeroOrOnce);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* SolidMaterialT::NewSub(const StringT& list_name) const
{
	if (list_name == "thermal_dilatation")
	{
		ParameterContainerT* thermal_dilatation = new ParameterContainerT(list_name);
		
		thermal_dilatation->AddParameter(ParameterT::Double, "percent_elongation");
		thermal_dilatation->AddParameter(ParameterT::Integer, "schedule_number");
		
		return thermal_dilatation;
	}
	else /* inherited */
		return ContinuumMaterialT::NewSub(list_name);
}

/* accept parameter list */
void SolidMaterialT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContinuumMaterialT::TakeParameterList(list);

	/* density */
	fDensity = list.GetParameter("density");

	/* 2D constraint - default to plane strain for 2D materials */
	list.GetParameter("constraint_2D", enum2int<ConstraintT>(fConstraint));
	if (NumSD() == 3)
		fConstraint = kNoConstraint;
	else if (NumSD() == 2 && fConstraint == kNoConstraint)
		fConstraint = kPlaneStrain;

	/* thermal dilatation */
	if (!fThermal) fThermal = new ThermalDilatationT;
	const ParameterListT* thermal = list.List("thermal_dilatation");
	if (thermal) {
		double elongation = thermal->GetParameter("percent_elongation");
		int schedule = thermal->GetParameter("schedule_number");
		schedule--;
		
		fThermal->SetPercentElongation(elongation);
		fThermal->SetScheduleNum(schedule);
	}
}
