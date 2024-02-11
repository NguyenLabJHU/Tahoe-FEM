/* $Id: FSFiberOptimize_MatT.cpp,v 1.1 2009-04-23 03:03:50 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberOptimize_MatT.h"
#include "FSFiberMatSupportT.h"
//#include "SS_Optimize_Primal.h"

using namespace Tahoe;

/* constructor */
FSFiberOptimize_MatT::FSFiberOptimize_MatT(void):
	fSupport(NULL)
{

}

FSFiberOptimize_MatT::~FSFiberOptimize_MatT(void)
{

}

/* set the material support */
void FSFiberOptimize_MatT::SetSupport(FSFiberMatSupportT* support)
{
	fSupport = support;	
}

