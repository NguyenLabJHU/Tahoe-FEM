/* $Id: NLCSolverWrapper.cpp,v 1.2.10.1 2002-12-10 17:06:59 paklein Exp $ */
#include "NLCSolverWrapper.h"
#include "MaterialsConfig.h"
#include "ExceptionT.h"

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "PolyCrystalMatT.h"
#include "SlipHardening.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "EVPFDBaseT.h"
#endif

using namespace Tahoe;

/* base class: NLCSolverWrapper */
NLCSolverWrapper::~NLCSolverWrapper() { }

/* derived class: SolverWrapperPoly */
SolverWrapperPoly::SolverWrapperPoly(PolyCrystalMatT& poly):
  fpoly(poly) { }

SolverWrapperPoly::~SolverWrapperPoly() { }

void SolverWrapperPoly::FormRHS(dArrayT& x, dArrayT& rhs) 
{ 
#ifdef PLASTICITY_CRYSTAL_MATERIAL
	fpoly.FormRHS(x, rhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormRHS", "PLASTICITY_CRYSTAL_MATERIAL not enabled");
#endif
}

void SolverWrapperPoly::FormLHS(dArrayT& x, dMatrixT& lhs) 
{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
	fpoly.FormLHS(x, lhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormRHS", "PLASTICITY_CRYSTAL_MATERIAL not enabled");
#endif
}

/* derived class: SolverWrapperHard */
SolverWrapperHard::SolverWrapperHard(SlipHardening& hard):
  fhard(hard) { }

SolverWrapperHard::~SolverWrapperHard() { }

void SolverWrapperHard::FormRHS(dArrayT& x, dArrayT& rhs) 
{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
	fhard.FormRHS(x, rhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormRHS", "PLASTICITY_CRYSTAL_MATERIAL not enabled");
#endif
}

void SolverWrapperHard::FormLHS(dArrayT& x, dMatrixT& lhs)
{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
	fhard.FormLHS(x, lhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormLHS", "PLASTICITY_CRYSTAL_MATERIAL not enabled");
#endif
}

/* derived class: SolverWrapperEVPBase */
SolverWrapperEVPBase::SolverWrapperEVPBase(EVPFDBaseT& evp):
  fevp(evp) { }

SolverWrapperEVPBase::~SolverWrapperEVPBase() { }

void SolverWrapperEVPBase::FormRHS(dArrayT& x, dArrayT& rhs) 
{
#ifdef PLASTICITY_MACRO_MATERIAL
	fevp.FormRHS(x, rhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormRHS", "PLASTICITY_MACRO_MATERIAL not enabled");
#endif
}

void SolverWrapperEVPBase::FormLHS(dArrayT& x, dMatrixT& lhs)
{
#ifdef PLASTICITY_MACRO_MATERIAL
	fevp.FormLHS(x, lhs); 
#else
	ExceptionT::GeneralFail("SolverWrapperPoly::FormLHS", "PLASTICITY_MACRO_MATERIAL not enabled");
#endif
}
