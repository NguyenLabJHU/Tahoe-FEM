/* $Id: SSOptimize_MatT.cpp,v 1.1 2009-04-23 03:03:52 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "SSOptimize_MatT.h"
#include "SSMatSupportT.h"
//#include "SS_Optimize_Primal.h"

using namespace Tahoe;

/* constructor */
SSOptimize_MatT::SSOptimize_MatT(void):
	fSupport(NULL)
{

}

SSOptimize_MatT::~SSOptimize_MatT(void)
{

}

/* set the material support */
void SSOptimize_MatT::SetSupport(SSMatSupportT* support)
{
	fSupport = support;	
}

