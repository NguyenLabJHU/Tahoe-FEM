/* $Id: SolidMaterialT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */

#include "SolidMaterialT.h"

#include <iostream.h>

#include "ThermalDilatationT.h"
#include "fstreamT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "ElasticT.h"
#include "LocalArrayT.h"

/* constructor */
SolidMaterialT::SolidMaterialT(ifstreamT& in,
	const ElasticT& element):
	ContinuumMaterialT(element)
//DEV
//	,fModuli(dSymMatrixT::NumValues(element.InitialCoordinates().MinorDim()))
{
	in >> fMassDamp;	if (fMassDamp  <  0.0) throw eBadInputValue;
	in >> fStiffDamp;	if (fStiffDamp <  0.0) throw eBadInputValue;
	in >> fDensity;		if (fDensity   <= 0.0) throw eBadInputValue;

	fThermal = new ThermalDilatationT(in);
	if (!fThermal) throw eOutOfMemory;
}

/* destructor */
SolidMaterialT::~SolidMaterialT(void) { delete fThermal; }

/* I/O functions */
void SolidMaterialT::Print(ostream& out) const
{
	/* inherited */
	ContinuumMaterialT::Print(out);
	
	out << " Mass damping coefficient. . . . . . . . . . . . = " << fMassDamp  << '\n';
	out << " Stiffness damping coefficient . . . . . . . . . = " << fStiffDamp << '\n';
	out << " Density . . . . . . . . . . . . . . . . . . . . = " << fDensity   << '\n';

	fThermal->Print(out);
}

/* return the wave speeds */
void SolidMaterialT::WaveSpeeds(const dArrayT& normal, dArrayT& speeds)
{
#if __option(extended_errorcheck)
	if (normal.Length() != speeds.Length()) throw eSizeMismatch;
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
			speeds[1] = speeds[0];
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

/* required parameter flags - (all FALSE by default) */
bool SolidMaterialT::NeedDisp(void) const     { return false; }
bool SolidMaterialT::NeedLastDisp(void) const { return false; }
bool SolidMaterialT::NeedVel(void) const      { return false; }

/* returns true if the material has internal forces in the unloaded
* configuration, ie thermal strains */
int SolidMaterialT::HasInternalStrain(void) const
{
	return fThermal->IsActive();
}

/* Thermal accessors */
int SolidMaterialT::ThermalLTfNumber(void) const
{
	return fThermal->LTfNumber();
}

void SolidMaterialT::SetThermalLTfPtr(const LoadTime* LTfPtr)
{
	fThermal->SetLTfPtr(LTfPtr);
}

//DEV - specific to 2D, small strain even, isotropic???
#if 0
double& SolidMaterialT::PlanarDilatationFactor(void)
{
	return fThermal->ScaleFactor();
}	
#endif

double SolidMaterialT::ThermalElongation(void) const //percentage
{
	return fThermal->PercentElongation();
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
