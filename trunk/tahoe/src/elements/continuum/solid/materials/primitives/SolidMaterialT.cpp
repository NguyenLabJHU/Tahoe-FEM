/* $Id: SolidMaterialT.cpp,v 1.16 2005-01-28 23:22:04 raregue Exp $ */
/* created: paklein (11/20/1996) */
#include "SolidMaterialT.h"

#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "LocalArrayT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* dummy return values */
iArrayT ijunk;
dArrayT djunk;

/* constructor */
SolidMaterialT::SolidMaterialT(void):
	ParameterInterfaceT("solid_material"),
	fThermal(NULL),
	fDensity(0.0),
	fConstraint(kNoConstraint),
	fCTE(0.0)
{

}

/* destructor */
SolidMaterialT::~SolidMaterialT(void) { delete fThermal; }

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


/* assumes small strains. returns true if the strain localization condition is satisfied,
* .ie if the acoustic tensor has zero (or negative eigenvalues),
* for the current conditions (current integration point and strain
* state). If localization is detected, the normals (current config)
* to the surface and slip directions are returned */
bool SolidMaterialT::IsLocalized_SS(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs)
{
	/* stress tensor */
	const dSymMatrixT& stress = s_ij();
			
	/* consistent tangent modulus */
	const dMatrixT& modulus = c_ijkl();
	
	/* elastic modulus */
	const dMatrixT& modulus_e = ce_ijkl();

	/* localization condition checker */
	DetCheckT checker(stress, modulus, modulus_e);
	normals.Dimension(NumSD());
	slipdirs.Dimension(NumSD());
	return checker.IsLocalized_SS(normals,slipdirs);
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
	constraint.SetDefault(fConstraint);
	list.AddParameter(constraint);
	
	/* coefficient of thermal expansion */
	ParameterT CTE(fCTE, "CTE");
	CTE.SetDefault(fCTE);
	list.AddParameter(CTE);
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
ParameterInterfaceT* SolidMaterialT::NewSub(const StringT& name) const
{
	if (name == "thermal_dilatation")
	{
		ParameterContainerT* thermal_dilatation = new ParameterContainerT(name);
		
		thermal_dilatation->AddParameter(ParameterT::Double, "percent_elongation");
		thermal_dilatation->AddParameter(ParameterT::Integer, "schedule_number");
		
		return thermal_dilatation;
	}
	else /* inherited */
		return ContinuumMaterialT::NewSub(name);
}

/* accept parameter list */
void SolidMaterialT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "SolidMaterialT::TakeParameterList";

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

	/* coefficient of thermal expansion */
	fCTE = list.GetParameter("CTE");

	/* thermal dilatation */
	if (!fThermal) fThermal = new ThermalDilatationT;
	const ParameterListT* thermal = list.List("thermal_dilatation");
	if (thermal) {
		double elongation = thermal->GetParameter("percent_elongation");
		int schedule_num = thermal->GetParameter("schedule_number");
		schedule_num--;
		
		fThermal->SetPercentElongation(elongation);
		fThermal->SetScheduleNum(schedule_num);
		const ScheduleT* schedule = MaterialSupport().Schedule(schedule_num);
		if (!schedule) ExceptionT::GeneralFail(caller, "could not resolve schedule %d", schedule_num+1);		
		fThermal->SetSchedule(schedule);
	}

	/* active thermal dilatation */
	if (fThermal->IsActive() && !SupportsThermalStrain())
		ExceptionT::BadInputValue("SolidMaterialT::Initialize", 
			"material does not support imposed thermal strain");
}

const iArrayT& SolidMaterialT::InternalDOF(void) const {
	ExceptionT::GeneralFail("SolidMaterialT::InternalDOF", "not implemented");
	return ijunk;
}

const dArrayT& SolidMaterialT::InternalStressVars(void) {
	ExceptionT::GeneralFail("SolidMaterialT::InternalStressVars", "not implemented");
	return djunk;
};

const dArrayT& SolidMaterialT::InternalStrainVars(void) {
	ExceptionT::GeneralFail("SolidMaterialT::InternalStrainVars", "not implemented");
	return djunk;
};
