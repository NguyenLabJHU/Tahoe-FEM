/* $Id: NOX_Tahoe_Group.cpp,v 1.1 2002-03-28 16:40:35 paklein Exp $ */
#include "NOX_Tahoe_Group.h"
#include "NOXInterfaceT.h"
#include "NOX_Tahoe_Vector.h"
#include "dArrayT.h"
#include "GlobalMatrixT.h"

using namespace NOX;
using namespace NOX::Tahoe;

/* constructor */
Group::Group(NOXInterfaceT& interface, dArrayT& X, GlobalMatrixT& J):
	fTahoeInterface(interface),
	fOwnJ(false),
	fSolution(X),
	fJ(&J),
	fJCopy(NULL),
	fIsRHS(false), fIsJacobian(false), fIsGrad(false), fIsNewton(false),
	fRHSNorm(0.0)
{

}

/* copy constructor */
Group::Group(const Group& source):
	fTahoeInterface(source.fTahoeInterface),
	fOwnJ(true),
	fJ(NULL),
	fJCopy(NULL)
{
	Group::operator=(source);
}

/* destructor */
Group::~Group(void)
{
	if (fOwnJ) delete fJ;
	delete fJCopy;
}

/* assignment operator */
Abstract::Group& Group::operator=(const Group& source)
{
	/* should not happen */
	if (!source.fJ) {
		cout << "\n Group::operator=: source has no Jacobian" << endl;
		throw eGeneralFail;
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
			throw eGeneralFail;
		}
		
		/* copy */
		*fJ = *(source.fJ);
	}

	/* copy of Jacobian */
	if (source.fJCopy) {
		if (fJCopy)
			*fJCopy = *(source.fJCopy);
		else
			fJCopy = source.fJCopy->Clone();
	}
	else {
		delete fJCopy;
		fJCopy = NULL;
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

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
#ifdef __NO_RTTI__
	cout << "\n NOX::Tahoe::Group::operator= : requires RTTI" << endl;
	throw eGeneralFail;
#endif

	const Group* t_group = dynamic_cast<const Group*>(&source);
	if (!t_group) {
		cout << "\n NOX::Tahoe::Group::operator= : cast failed" << endl;
		throw eGeneralFail;
	}
	return operator=(*t_group);
}

/* compute x where this.x = grp.x() + step * d */
bool Group::computeX(const Group& grp, const Vector& d, double step)
{
	/* compute new solution */
	fSolution.get_dArrayT().SetToCombination(1.0, grp.getX(), step, d);

	/* set flags */
	fIsRHS = fIsJacobian = fIsGrad = fIsNewton = false;
	return true;
}

bool Group::computeX(const Abstract::Group& grp, const Abstract::Vector& d, double step)
{
#ifdef __NO_RTTI__
	cout << "\n NOX::Tahoe::Group::computeX: requires RTTI" << endl;
	throw eGeneralFail;
#endif

	/* cast group */
	const Group* t_grp = dynamic_cast<const Group*>(&grp);
	if (!t_grp) {
		cout << "\n NOX::Tahoe::Group::computeX: Group cast failed" << endl;
		throw eGeneralFail;
	}

	/* cast vector */
	const Vector* t_d = dynamic_cast<const Vector*>(&d);
	if (!t_d) {
		cout << "\n NOX::Tahoe::Group::computeX: Vector cast failed" << endl;
		throw eGeneralFail;
	}

	/* do compute */
	return computeX(*t_grp, *t_d, step);
}

/** Compute and return RHS */
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
	if (!fIsRHS) cout << "NOX::Tahoe::Group::computeRHS: failed" << endl;

	return fIsRHS;
}

/** Compute and return RHS */
bool Group::computeJacobian(void) 
{
	/* not up to date */
	if (!fIsJacobian)
	{
		/* Tahoe requires RHS be formed before asking for the Jacobian */
		fIsJacobian = computeRHS() && fTahoeInterface.computeJacobian(*fJ);
	}

	/* message */
	if (!fIsRHS) cout << "NOX::Tahoe::Group::computeJacobian: failed" << endl;

	return fIsJacobian;
}

/* Compute and return gradient - J^T RHS*/
bool Group::computeGrad(void)
{
	//TEMP
	return false;

	/* not up to date */
	if (!fIsGrad) 
	{
	
		//fIsGrad = fIsRHS && fIsJacobian
	
	
	
	
	}
	return fIsGrad;
}

