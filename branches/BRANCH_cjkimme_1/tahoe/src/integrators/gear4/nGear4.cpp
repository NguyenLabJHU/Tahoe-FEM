#include "nGear4.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "KBC_CardT.h"
#include "BasicFieldT.h"
#include "AndersenPressureT.h"

using namespace Tahoe;

/* constructor */
nGear4::nGear4(void):
	v_field(4)
{
	/* Gear constants */
	F02 = 1.0/12.0;
	F12 = 5.0/12.0;
	F22 = 1.0;
	F32 = 1.0;
}

/* consistent BC's - updates predictors and acceleration only */
void nGear4::ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC)
{
	const char caller[] = "nGear4::ConsistentKBC";

	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail(caller, "field must be order 3: %d", field.Order());

	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (field[0])(node, dof);
	double& v = (field[1])(node, dof);
	double& a = (field[2])(node, dof);

	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		{
			d = 0.0;
			v = 0.0; //correct?	
			a = 0.0; //correct?	
			break;
		}
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			d = KBC.Value();
		   	break;
			/* NOTE:  haven't figured out a correct way to
			   compute velocities and accelerations given a
			   prescribed displacement...*/
		}
		
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			v = KBC.Value();
			// NEED ACTUAL, NOT PREDICTED ACCELERATION TO DEFINE
			// THIS KBC!!! (update array?)
			break;
		}
		
		case KBC_CardT::kAcc: /* prescribed acceleration */
		{
#pragma message("Chris, fix KBC_CardT::kAcc")		
			double a_next  = KBC.Value();
			v -= F12 * (a - a_next);
			d -= F02 * (a - a_next);
			a = a_next;
			break;
		}

		case KBC_CardT::kNull: /* do nothing */
		{
			break;
		}

		default:
			ExceptionT::BadInputValue(caller, "unknown BC code: %d", KBC.Code() );
	}
}		

/* predictors - map ALL */
void nGear4::Predictor(BasicFieldT& field)
{
	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail("nGear4::Predictor", "field must be order 3: %d", field.Order());

	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();

	/* run through arrays */
	int len = field[0].Length();
	for (int i = 0; i < len; i++)
	{
		(*p0++) += fdt*(*p1) + fdt2*(*p2) + fdt3*(*p3);
		(*p1++) += fdt*(*p2) + fdt2*(*p3);
		(*p2++) += fdt*(*p3);
	}
	
	/* integrate the system volume */
	v_field[0] += fdt*v_field[1] + fdt2*v_field[2] + fdt3*v_field[3];
	v_field[1] += fdt*v_field[2] + fdt2*v_field[3];
	v_field[2] += fdt*v_field[3];
}		

/* corrector. Maps ALL degrees of freedom forward. */
void nGear4::Corrector(BasicFieldT& field, const dArray2DT& update)
{
	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail("nGear4::Corrector", "field must be order 3: %d", field.Order());
	
	/* get stuff from barostat */
	double dP(AndersenPressureT::MdP);
	double dVFactor = (v_field[2] - 2./3.*v_field[1]*v_field[1]/v_field[0])/v_field[0]/3.;
	if (!AndersenPressureT::fV_field)
		ExceptionT::GeneralFail("nGear4::Corrector", "cannot get field for volume");
	const dArrayT& v_field = *AndersenPressureT::fV_field;

	
	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();
	double* pu = update.Pointer();

	if (fabs(fdt) > kSmall)
	{
		/* run through arrays */
		int len = field[0].Length();
		for (int i = 0; i < len; i++)
		{	
			double error = *pu++ + *p0*dVFactor;	
			error = (*p2 - error)*fdt2;

			*p0 -= error*f02;
			*p1 -= error*f12;
			*p2 -= error*f22;
			*p3 -= error*f32;

			/* next */
			p0++; p1++; p2++; p3++;
		}
		
		double error = (v_field[2] - dP)*fdt2;
		v_field[0] -= error*f02;
		v_field[1] -= error*f12;
		v_field[2] -= error*f22;
		v_field[3] -= error*f32; 
	}
	else /* for dt -> 0.0 */
	{
		/* run through arrays */
		int len = field[0].Length();
		for (int i = 0; i < len; i++)
		{
			*p2 -= ((*p2) - (*pu++) - (*p0++)*dVFactor)*F22;
			p2++;
		}
	}
}

/* correctors - map ACTIVE */
void nGear4::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail("nGear4::Corrector", "field must be order 3: %d", field.Order());

	/* get stuff from barostat */
	double dP(AndersenPressureT::MdP);
	double dVFactor = (v_field[2] - 2./3.*v_field[1]*v_field[1]/v_field[0])/v_field[0]/3.;

	const iArray2DT& eqnos = field.Equations();

	/* add update - assumes that fEqnos maps directly into dva */
	int    *peq = eqnos.Pointer();
	
	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();

	if (fabs(fdt) > kSmall)
	{
		for (int i = 0; i < eqnos.Length(); i++)
		{
			int eq = *peq++ - eq_start;
		
			/* active dof */
			if (eq > -1 && eq < num_eq)
			{
				double a = update[eq];
				double error = a + *p0*dVFactor;	
				error = (*p2 - error)*fdt2;
				*p0 -= error*f02;
				*p1 -= error*f12;
				*p2 -= error*f22;
				*p3 -= error*f32; 
			}
		
			/* next */
			p0++; p1++; p2++; p3++;
		}
		
		double error = (v_field[2] - dP)*fdt2;
		v_field[0] -= error*f02;
		v_field[1] -= error*f12;
		v_field[2] -= error*f22;
		v_field[3] -= error*f32; 
	}
	else /* for dt -> 0.0 */
	{
		for (int i = 0; i < eqnos.Length(); i++)
		{
			int eq = *peq++ - eq_start;
		
			/* active dof */
			if (eq > -1 && eq < num_eq)
			{
				double a = update[eq];
				*p2 -= ((*p2) - a - (*p0++)*dVFactor)*F22;
			}
		
			/* next */
			p2++;
		}
	}
}

void nGear4::MappedCorrector(BasicFieldT& field, const iArrayT& map, 
	const iArray2DT& flags, const dArray2DT& update)
{
//TEMP
ExceptionT::Stop("nGear4::MappedCorrector", "not implemented");

	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail("nGear4::MappedCorrector", "field must be order 3: %d", field.Order());

	/* get stuff from barostat */
	double dP(AndersenPressureT::MdP);
	double dVFactor = (v_field[2] - 2./3.*v_field[1]*v_field[1]/v_field[0])/v_field[0]/3.;
	

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();

	/* fetch pointers */
	double* p3 = field[3].Pointer();

	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		int* pflags = flags(i);
		
		double* p0 = (field[0])(row);
		double* p1 = (field[1])(row);
		double* p2 = (field[2])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag > 0)
			{
				double a = *pupdate + *p0*dVFactor;
				double error = (*p2) - a;
				*p0 -= error * F02;
				*p1 -= error * F12;
				*p2 -= error;
				*p3 -= error * F32;
			}
			
			/* next */
			pflag++; pupdate++; 
			p0++; p1++; p2++; p3++;
		}
	}
	
	double error = (v_field[2] - dP)*fdt2;
	v_field[0] -= error*f02;
	v_field[1] -= error*f12;
	v_field[2] -= error*f22;
	v_field[3] -= error*f32; 
}

/* return the field array needed by nIntegratorT::MappedCorrector. */
const dArray2DT& nGear4::MappedCorrectorField(BasicFieldT& field) const
{
	/* check */
	if (field.Order() != 3)
		ExceptionT::GeneralFail("nGear4::MappedCorrectorField", "field must be order 3: %d", field.Order());

	return field[2];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nGear4::ExternalNodeCondition(void) const
{
	return KBC_CardT::kAcc;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nGear4::nComputeParameters(void)
{
	/* nothing implemented here */
	fdt2 =  fdt*fdt/2.0;
	fdt3 = fdt2*fdt/3.0;
	f02 = F02;
	f12 = F12/fdt;
	f22 = F22/fdt2;
	f32 = F32/fdt3;
}
