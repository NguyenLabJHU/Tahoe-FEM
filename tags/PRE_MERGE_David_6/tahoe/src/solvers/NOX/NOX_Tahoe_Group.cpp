/* $Id: NOX_Tahoe_Group.cpp,v 1.6 2003-12-28 08:24:18 paklein Exp $ */
#include "NOX_Tahoe_Group.h"

/* optional */
#ifdef __NOX__

#include "NOX_Tahoe_Vector.h"
#include "SolverInterfaceT.h"
#include "dArrayT.h"
#include "GlobalMatrixT.h"

using namespace Tahoe;
using namespace NOX;

/* constructor */
Group::Group(SolverInterfaceT& interface, dArrayT& X, GlobalMatrixT& J):
	fTahoeInterface(interface),
	fOwnJ(false),
	fJ(&J),
	fSolution(X),
	fRHS(fSolution, CopyShape),
	fGradient(fSolution, CopyShape),
	fNewton(fSolution, CopyShape),
	fIsRHS(false), fIsJacobian(false), fIsGrad(false), fIsNewton(false),
	fRHSNorm(0.0)
{

}

/* copy constructor */
Group::Group(const Group& source, NOX::CopyType type):
	fTahoeInterface(source.fTahoeInterface),
	fOwnJ(true),
	fJ(NULL)
{
#pragma unused(type)
	Group::operator=(source);
}

/* destructor */
Group::~Group(void)
{
	if (fOwnJ) delete fJ;
}

/* assignment operator */
NOX::Abstract::Group& Group::operator=(const Group& source)
{
	/* should not happen */
	if (!source.fJ) {
		cout << "\n Group::operator=: source has no Jacobian" << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* copy solution vector and Jacobian */
	if (fOwnJ) {
		if (fJ)
			*fJ = *(source.fJ);
		else
			fJ = source.fJ->Clone();
	}
	else {
	
		/* should not be */
		if (!fJ) {
			cout << "\n Group::operator=: shallow object has no Jacobian" << endl;
			throw ExceptionT::kGeneralFail;
		}
		
		/* copy */
		*fJ = *(source.fJ);
	}
	
	/* copy vectors */
	fRHS      = source.fRHS;
	fRHSNorm  = source.fRHSNorm;
	fSolution = source.fSolution;
	fGradient = source.fGradient;
	fNewton   = source.fNewton;

	/* copy flags */
	fIsRHS      = source.fIsRHS;
	fIsJacobian = source.fIsJacobian;
	fIsGrad     = source.fIsGrad;
	fIsNewton   = source.fIsNewton;

	return *this;
}

NOX::Abstract::Group& Group::operator=(const NOX::Abstract::Group& source)
{
#ifdef __NO_RTTI__
	cout << "\n Tahoe::Group::operator= : requires RTTI" << endl;
	throw ExceptionT::kGeneralFail;
#endif

	const Group* t_group = TB_DYNAMIC_CAST(const Group*, &source);
	if (!t_group) {
		cout << "\n Tahoe::Group::operator= : cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}
	return operator=(*t_group);
}

/* compute x where this.x = grp.x() + step * d */
bool Group::computeX(const Group& grp, const Vector& d, double step)
{
	/* dimens

	/* compute new solution */
	fSolution.get_dArrayT().SetToCombination(1.0, grp.fSolution, step, d);

	/* set flags */
	fIsRHS = fIsJacobian = fIsGrad = fIsNewton = false;
	return true;
}

bool Group::computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step)
{
#ifdef __NO_RTTI__
	cout << "\n Tahoe::Group::computeX: requires RTTI" << endl;
	throw ExceptionT::kGeneralFail;
#endif

	/* cast group */
	const Group* t_grp = TB_DYNAMIC_CAST(const Group*, &grp);
	if (!t_grp) {
		cout << "\n Tahoe::Group::computeX: Group cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* cast vector */
	const Vector* t_d = TB_DYNAMIC_CAST(const Vector*, &d);
	if (!t_d) {
		cout << "\n Tahoe::Group::computeX: Vector cast failed" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* do compute */
	return computeX(*t_grp, *t_d, step);
}

/* Compute and return RHS */
bool Group::computeRHS(void) 
{
	/* not up to date */
	if (!fIsRHS)
	{
		/* compute RHS */
		fIsRHS = fTahoeInterface.computeRHS(fSolution, fRHS);

		/* success */
		if (fIsRHS)
			fRHSNorm = fRHS.norm();
		else {
			fRHS = 0.0;
			fRHSNorm = 0.0;
		}
	}
	
	/* message */
	if (!fIsRHS) cout << "\n Tahoe::Group::computeRHS: failed" << endl;

	return fIsRHS;
}

/* Compute and return RHS */
bool Group::computeJacobian(void) 
{
	/* not up to date */
	if (!fIsJacobian)
	{
		/* Tahoe requires RHS be formed before asking for the Jacobian */
		fIsJacobian = computeRHS() && fTahoeInterface.computeJacobian(*fJ);
	}

	/* message */
	if (!fIsJacobian) cout << "\n Tahoe::Group::computeJacobian: failed" << endl;

	return fIsJacobian;
}

/* Compute and return gradient - J^T RHS*/
bool Group::computeGrad(void)
{
	/* not up to date */
	if (!fIsGrad)
	{
		/* Jacobian (and RHS) must be up to date */
		fIsGrad = computeJacobian() && fJ->MultTx(fRHS, fGradient);	
	}

	/* message */
	if (!fIsGrad) cout << "\n Tahoe::Group::computeGrad: failed" << endl;

	return fIsGrad;
}

/* compute the Newton direction */
bool Group::computeNewton(NOX::Parameter::List& params)
{
#pragma unused(params)

	/* not up to date */
	if (!fIsNewton)
	{
		/* Jacobian (and RHS) must be up to date */
		if (computeJacobian()) {

			/* initialize with RHS */
			fNewton = fRHS;
			
			/* copy Jacobian if needed */
			if (fJ->SolvePreservesData())
			{
				/* linear solve */
				if (!fJ->Solve(fNewton)) {
					cout << "\n Group::computeNewton: solve failed" << endl;
					throw ExceptionT::kBadJacobianDet;
				}
			}
			else /* need to make copy */
			{
				/* copy */
				GlobalMatrixT* A = fJ->Clone();

				/* linear solve */
				if (!A->Solve(fNewton)) {
					cout << "\n Group::computeNewton: solve failed" << endl;
					throw ExceptionT::kBadJacobianDet;
				}

				/* free */
				delete A;
			}
			
			/* success */
			fIsNewton = true;
		}
	}

	/* message */
	if (!fIsNewton) cout << "\n Tahoe::Group::computeNewton: failed" << endl;
		
	return fIsNewton;
}

/* compute result = Jacobian * input */
bool Group::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
	/* Jacobian must be up to date and cannot be calculated here */
	if (!fIsJacobian)
		return false;
	else
	{
		/* cast to Tahoe vectors */
		const Vector* t_input = TB_DYNAMIC_CAST(const Vector*, &input);
		if (!t_input) {
			cout << "\n Group::applyJacobian: cast of input failed" << endl;
			throw ExceptionT::kGeneralFail;
		}
		Vector* t_result = TB_DYNAMIC_CAST(Vector*, &result);
		if (!t_result) {
			cout << "\n Group::applyJacobian: cast of result failed" << endl;
			throw ExceptionT::kGeneralFail;
		}
	
		/* compute */
		return fJ->Multx(*t_input, *t_result);
	}
}

/* compute result = Jacobian^T * input */
bool Group::applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
	/* Jacobian must be up to date and cannot be calculated here */
	if (!fIsJacobian)
		return false;
	else
	{
		/* cast to Tahoe vectors */
		const Vector* t_input = TB_DYNAMIC_CAST(const Vector*, &input);
		if (!t_input) {
			cout << "\n Group::applyJacobianTranspose: cast of input failed" << endl;
			throw ExceptionT::kGeneralFail;
		}
		Vector* t_result = TB_DYNAMIC_CAST(Vector*, &result);
		if (!t_result) {
			cout << "\n Group::applyJacobianTranspose: cast of result failed" << endl;
			throw ExceptionT::kGeneralFail;
		}
	
		/* compute */
		return fJ->MultTx(*t_input, *t_result);
	}
}
  
 /* Applies the Jacobian Diagonal to the given input vector */
bool Group::applyJacobianDiagonalInverse(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
	/* Jacobian must be up to date and cannot be calculated here */
	if (!fIsJacobian)
		return false;
	else
	{
		/* cast to Tahoe vector */
		Vector* t_result = TB_DYNAMIC_CAST(Vector*, &result);
		if (!t_result) {
			cout << "\n Group::applyJacobianDiagonalInverse: cast of result failed" << endl;
			throw ExceptionT::kGeneralFail;
		}

		/* collect */
		dArrayT diags;
		if (fJ->CopyDiagonal(diags))
		{
			/* copy in */
			result = input;

			/* apply */
			t_result->get_dArrayT() /= diags;
			return true;	
		}
		else
			return false;
	}
}

/* create a new %Group of the same derived type as this one by cloning */
NOX::Abstract::Group* Group::clone(NOX::CopyType type) const
{
	Group* new_clone = new Group(*this, type);
	return new_clone;
}

#endif /* __NOX__ */
