/* $Header: /home/cvs/t/ta/tahoe/tahoe/src/integrators/mixed/nMixed.cpp,v 1.4 2006-08-18 20:00:40 tdnguye Exp $ */
/* created: a-kopacz (08/08/2006) */

#include "nMixed.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "KBC_CardT.h"
#include "BasicFieldT.h"

using namespace Tahoe;

/* constructor */
nMixed::nMixed(void) 
{ }

/*All field arrays must be registered with nIntegratorT::Dimension before calling predictors, correctors, 
or boundary condition calculations. Note that the predictor is applied to all degrees of freedom, 
while the correctors are only applied to the active degrees of freedom
*/
void nMixed::Dimension(BasicFieldT& field)
{
	/*inherited*/
	nIntegratorT::Dimension(field);
	fNumDOF = field.NumDOF();
	dpred_v_.Dimension(fNumDOF);
	dcorr_v_.Dimension(fNumDOF);

}

/* consistent BC's */
void nMixed::ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC)
{
	/* destinations */
	int node = KBC.Node();
	int dof  = KBC.DOF();
	double& d = (field[0])(node, dof);
	double& v = (field[1])(node, dof);
		
	switch ( KBC.Code() )
	{
		case KBC_CardT::kFix: /* zero displacement */
		case KBC_CardT::kDsp: /* prescribed displacement */
		{
			double d_next = KBC.Value();

			if (fabs(dcorr_v_[dof]) > kSmall) /* for dt -> 0.0 */
				v = (d_next - d)/dcorr_v_[dof];
			else
				v = 0.0;
			d = d_next;
			break;
		}
		case KBC_CardT::kVel: /* prescribed velocity */
		{
			v  = KBC.Value();
			d += dcorr_v_[dof]*v;
			break;
		}
		case KBC_CardT::kNull:
		{
			break;
		}
		default:
			ExceptionT::BadInputValue("nMixed::ConsistentKBC",
				"unknown BC code: %d", KBC.Code());
	}
}		
#pragma message ("roll up redundancy after it works")
// predictors - map ALL, unless limit arguments are specified
void nMixed::Predictor(BasicFieldT& field, int fieldstart /*= 0*/, int fieldend /*= -1*/)
{
	int NumDOF = field.NumDOF();
	int NumNodes  = field.NumNodes();
	if (fieldend == -1) // operate on full arrays
	{
		for ( int d = 0; d< NumDOF; d++)
		{
			for ( int n = 0; n< NumNodes; n++)
			{	
				/* displacement predictor */	
				(field[0])(n,d) += dpred_v_[d]*(field[1])(n,d);

				/* velocity predictor */
				(field[1])(n,d) = 0.0;
			}	
		}	
	}
	else // operate on restricted contiguous block of the arrays
		ExceptionT::BadInputValue("nMixed::Predictor", "operate on full arrays ONLY");
}		

/* correctors - map ALL */
void nMixed::Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart /*= 0*/, int fieldend /*= -1*/, int dummy /*= 0*/)
{
	int NumDOF = field.NumDOF();
	int NumNodes  = field.NumNodes();
	int start = 0;
	int end = 0;
	if (fieldend == -1) // operate on full arrays
	{
/*		for ( int d = 0; d< NumDOF; d++)
		{
			for ( int n = 0; n< NumNodes; n++)
			{	
				(field[0])(n,d) += dcorr_v_[d]*update(n,d);

				(field[1])(n,d) += update(n,d);
			}
		}
*/
		start = 0;
		end = field[0].Length();
	}
	else // operate on restricted contiguous block of the arrays
	{
		start = fieldstart;
		end = fieldend + 1;		
	}
	
	double *pd  = field[0].Pointer();
	double *pv  = field[1].Pointer();
	for (int i = start; i < end; i++)
	{
		int dof = i % NumDOF;

		/* displacement corrector */	
		pd[i] += dcorr_v_[dof]*update[i];

		/* velocity corrector */	
		pv[i] += update[i];
	}
}

/* correctors - map ACTIVE */
void nMixed::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	const int NumDOF = field.NumDOF();
	
	const iArray2DT& eqnos = field.Equations();
	/*dimension of eqnos set in FieldT::InitEquations()*/
	/*eqnos dimensions are (num_nodes x num_dof)*/
	if (eqnos.Length() != field[0].Length())
		ExceptionT::GeneralFail("nMixed::Corrector",
				"length of eqnos array must be same in length as field arrays.");
	/* add update - assumes that fEqnos maps directly into FieldT */
	const int* peq = eqnos.Pointer();
	double *pd  = field[0].Pointer();
	double *pv  = field[1].Pointer();
	for (int i = 0; i < eqnos.Length(); i++)
	{
		int dof = i % NumDOF;
		int eq = *peq++ - eq_start;
		/* active dof */
		if (eq > -1 && eq < num_eq)
		{
			double v = update[eq];
			*pd += dcorr_v_[dof]*v;
			*pv += v;
		}
		pd++;
		pv++;
	}
}

void nMixed::MappedCorrector(BasicFieldT& field, const iArrayT& map, 
	const iArray2DT& flags, const dArray2DT& update)
{
	/* run through map */
	int minordim = flags.MinorDim();
	/*assume minordim = numdof*/
	if (minordim != field.NumDOF())
		ExceptionT::GeneralFail("nMixed::MappedCorrector",
				"minordim  of flags array equal number of DOFs.");
	
	const int* pflag = flags.Pointer();
	const double* pupdate = update.Pointer();
	for (int i = 0; i < map.Length(); i++)
	{
		int row = map[i];
		const int* pflags = flags(i);

		double* pd = (field[0])(row);
		double* pv = (field[1])(row);
		for (int j = 0; j < minordim; j++)
		{
			/* active */
			if (*pflag > 0)
			{
				double dv = *pupdate - *pv; 
/* NOTE: update is the total v_n+1, so we need to recover dv, the velocity
 *       increment. Due to truncation errors, this will not match the update
 *       applied in nTrapezoid::Corrector exactly. Therefore, ghosted nodes
 *       will not have exactly the same trajectory as their images. The solution
 *       would be to expand da to the full field, and send it. Otherwise, we
 *       could develop maps from the update vector to each communicated outgoing
 *       packet. Or {d,v} could be recalculated from t_n here and in Corrector,
 *       though this would require the data from t_n. */

				*pd += dcorr_v_[j]*dv;
				*pv  = *pupdate; /* a_n+1 matches exactly */
			}
			
			/* next */
			pflag++; pupdate++; pd++; pv++;
		}
	}
}

/* return the field array needed by nIntegratorT::MappedCorrector. */
const dArray2DT& nMixed::MappedCorrectorField(BasicFieldT& field) const
{
	return field[1];
}

/* pseudo-boundary conditions for external nodes */
KBC_CardT::CodeT nMixed::ExternalNodeCondition(void) const
{
	return KBC_CardT::kVel;
}

/***********************************************************************
* Protected
***********************************************************************/

/* recalculate time stepping constants */
void nMixed::nComputeParameters(void)
{
	/* predictor */
	
	/*sets different corrector and predictor for final degree of freedom*/
	for (int i = 0; i < fNumDOF-1; i++)
	{
		/* trapezoidal rule/midpoint rule/Crank-Nicolson */
		dpred_v_[i] = 0.5*fdt;
		dcorr_v_[i] = 0.5*fdt;
	}

	/* backward difference/backward Euler */
	dpred_v_[fNumDOF] = 1.0*fdt;
	dcorr_v_[fNumDOF] = 1.0*fdt;
}
