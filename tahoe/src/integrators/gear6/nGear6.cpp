/* $Id: nGear6.cpp,v 1.13 2004-12-26 21:09:05 d-farrell2 Exp $ */
#include "nGear6.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "KBC_CardT.h"
#include "BasicFieldT.h"

using namespace Tahoe;

/* constructor */
nGear6::nGear6(void)
{
	/* Gear constants */
	F02 = 3.0/16.0;
	F12 = 251.0/360.0;
	F22 = 1.0;
	F32 = 11.0/18.0;
	F42 = 1.0/6.0;
	F52 = 1.0/60.0;
}

/* consistent BC's - updates predictors and acceleration only */
void nGear6::ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC)
{
	const char caller[] = "nGear6::ConsistentKBC";

	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail(caller, "field must be order 6: %d", field.Order());

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
#pragma message("nGear6::Predictor, not implemented with limits yet, declaration changed to match others")
/* predictors - map ALL */
void nGear6::Predictor(BasicFieldT& field, int fieldstart /*= 0*/, int fieldend /*= -1*/)
{
	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail("nGear6::Predictor", "field must be order 6: %d", field.Order());

	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();
	double* p4 = field[4].Pointer();
	double* p5 = field[5].Pointer();

	/* run through arrays */
	int len = field[0].Length();
	for (int i = 0; i < len; i++)
	{
		(*p0++) += fdt*(*p1) + fdt2*(*p2) + fdt3*(*p3) + fdt4*(*p4) + fdt5*(*p5);
		(*p1++) += fdt*(*p2) + fdt2*(*p3) + fdt3*(*p4) + fdt4*(*p5);
		(*p2++) += fdt*(*p3) + fdt2*(*p4) + fdt3*(*p5);
		(*p3++) += fdt*(*p4) + fdt2*(*p5);
		(*p4++) += fdt*(*p5++);
	}
}		
#pragma message("nGear6::Corrector, not implemented with limits yet, declaration changed to match others")
/* corrector. Maps ALL degrees of freedom forward. */
void nGear6::Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart /*= 0*/, int fieldend /*= -1*/, int dummy /*= 0*/)
{
	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail("nGear6::Corrector", "field must be order 6: %d", field.Order());
	
	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();
	double* p4 = field[4].Pointer();
	double* p5 = field[5].Pointer();
	const double* pu = update.Pointer();

	if (fabs(fdt) > kSmall)
	{
		/* run through arrays */
		int len = field[0].Length();
		for (int i = 0; i < len; i++)
		{		
			double error = ((*p2) - (*pu++))*fdt2;

			*p0 -= (error*F02);
			*p1 -= (error*F12)/fdt;
			*p2 -= (error*F22)/fdt2;
			*p3 -= (error*F32)/fdt3; 
			*p4 -= (error*F42)/fdt4;
			*p5 -= (error*F52)/fdt5;

			/* next */
			p0++; p1++; p2++; p3++; p4++; p5++;
		}
	}
	else /* for dt -> 0.0 */
	{
		/* run through arrays */
		int len = field[0].Length();
		for (int i = 0; i < len; i++)
		{
			*p2 -= ((*p2) - (*pu++))*F22;
			p2++;
		}
	}
}

/* correctors - map ACTIVE */
void nGear6::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail("nGear6::Corrector", "field must be order 6: %d", field.Order());

	const iArray2DT& eqnos = field.Equations();

	/* add update - assumes that fEqnos maps directly into dva */
	const int *peq = eqnos.Pointer();
	
	/* fetch pointers */
	double* p0 = field[0].Pointer(); 
	double* p1 = field[1].Pointer();
	double* p2 = field[2].Pointer();
	double* p3 = field[3].Pointer();
	double* p4 = field[4].Pointer();
	double* p5 = field[5].Pointer();

	if (fabs(fdt) > kSmall)
	{
		for (int i = 0; i < eqnos.Length(); i++)
		{
			int eq = *peq++ - eq_start;
		
			/* active dof */
			if (eq > -1 && eq < num_eq)
			{
				double a = update[eq];
				double error = ((*p2) - a)*fdt2;
				*p0 -= (error*F02);
				*p1 -= (error*F12)/fdt;
				*p2 -= (error*F22)/fdt2;
				*p3 -= (error*F32)/fdt3; 
				*p4 -= (error*F42)/fdt4;
				*p5 -= (error*F52)/fdt5;
			}
		
			/* next */
			p0++; p1++; p2++; p3++; p4++; p5++;
		}
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
				*p2 -= ((*p2) - a)*F22;
			}
		
			/* next */
			p2++;
		}
	}
}

void nGear6::MappedCorrector(BasicFieldT& field, const iArrayT& map, 
	const iArray2DT& flags, const dArray2DT& update)
{
//TEMP
ExceptionT::Stop("nGear6::MappedCorrector", "not implemented");

	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail("nGear6::MappedCorrector", "field must be order 6: %d", field.Order());

	/* run through map */
	int minordim = flags.MinorDim();
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();

	/* fetch pointers */
	double* p3 = field[3].Pointer();
	double* p4 = field[4].Pointer();
	double* p5 = field[5].Pointer();

	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		const int* pflags = flags(i);
		
		double* p0 = (field[0])(row);
		double* p1 = (field[1])(row);
		double* p2 = (field[2])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag > 0)
			{
				double a = *pupdate;
				double error = (*p2) - a;
				*p0 -= error * F02;
				*p1 -= error * F12;
				*p2 -= error;
				*p3 -= error * F32;
				*p4 -= error * F42;
				*p5 -= error * F52;
			}
			
			/* next */
			pflag++; pupdate++; 
			p0++; p1++; p2++; p3++; p4++; p5++;
		}
	}
}

/* return the field array needed by nIntegratorT::MappedCorrector. */
const dArray2DT& nGear6::MappedCorrectorField(BasicFieldT& field) const
{
	/* check */
	if (field.Order() != 5)
		ExceptionT::GeneralFail("nGear6::MappedCorrectorField", "field must be order 6: %d", field.Order());

	return field[2];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nGear6::ExternalNodeCondition(void) const
{
	return KBC_CardT::kAcc;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nGear6::nComputeParameters(void)
{
	/* nothing implemented here */
	fdt2 =  fdt*fdt/2.0;
	fdt3 = fdt2*fdt/3.0;
	fdt4 = fdt3*fdt/4.0;
	fdt5 = fdt4*fdt/5.0;
}
