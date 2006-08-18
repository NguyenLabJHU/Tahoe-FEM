/* $Id: */
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
{ 
	/* fluid only for now :: {v_x, v_y, v_z, p} */
	int numDOF = 4; /* NEED TO GET THIS SOMEHOW */
	dpred_v_.Dimension(numDOF);
	dcorr_v_.Dimension(numDOF);
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
	if (fieldend == -1) // operate on full arrays
	{
		for ( int d = 0; d< NumDOF; d++)
		{
			for ( int n = 0; n< NumNodes; n++)
			{	
				/* displacement corrector */	
				(field[0])(n,d) += dcorr_v_[d]*update(n,d);

				/* velocity corrector */
				(field[1])(n,d) += update(n,d);
			}
		}
	}
	else // operate on restricted contiguous block of the arrays
		ExceptionT::BadInputValue("nMixed::Corrector", "operate on full arrays ONLY");
}

/* correctors - map ACTIVE */
void nMixed::Corrector(BasicFieldT& field, const dArrayT& update, 
	int eq_start, int num_eq)
{
	ExceptionT::BadInputValue("nMixed::Corrector", "bad call");
}

void nMixed::MappedCorrector(BasicFieldT& field, const iArrayT& map, 
	const iArray2DT& flags, const dArray2DT& update)
{
	ExceptionT::BadInputValue("nMixed::Corrector", "bad call");
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
	
	/* trapezoidal rule/midpoint rule/Crank-Nicolson */
	dpred_v_[0] = 0.5*fdt;
	dpred_v_[1] = 0.5*fdt;
	dpred_v_[2] = 0.5*fdt;
	
	/* backward difference/backward Euler */
	dpred_v_[3] = 1.0*fdt;

	/* corrector */
	
	/* trapezoidal rule/midpoint rule/Crank-Nicolson */
	dcorr_v_[0] = 0.5*fdt;
	dcorr_v_[1] = 0.5*fdt;
	dcorr_v_[2] = 0.5*fdt;
	
	/* backward difference/backward Euler */
	dcorr_v_[3] = 1.0*fdt;
}
