/* $Id: GlobalT.cpp,v 1.9.6.1 2005-05-18 18:30:46 paklein Exp $ */
/* created: paklein (04/01/2000) */
#include "GlobalT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* returns flag with precedence */
GlobalT::RelaxCodeT GlobalT::MaxPrecedence(GlobalT::RelaxCodeT code1,
	GlobalT::RelaxCodeT code2)
{
	RelaxCodeT result = kNoRelax;
	if (code1 == kReEQRelax || code2 == kReEQRelax)
		result = kReEQRelax;
	else if (code1 == kReEQ && code2 == kRelax)
		result = kReEQRelax;
	else if (code1 == kRelax && code2 == kReEQ)
		result = kReEQRelax;
	else if (code1 == kReEQ || code2 == kReEQ)
		result = kReEQ;
	else if (code1 == kRelax || code2 == kRelax)
		result = kRelax;
	else if (code1 == kNoRelax && code2 == kNoRelax)
		result = kNoRelax;
	else
		ExceptionT::GeneralFail("GlobalT::MaxPrecedence", "not expecting %d and %d",
			code1, code2);

	return result;
}

/* returns flag with precedence */
GlobalT::InitStatusT GlobalT::MaxPrecedence(
	GlobalT::InitStatusT status_1, GlobalT::InitStatusT status_2)
{
	return (status_1 > status_2) ? status_1 : status_2;
}

GlobalT::LoggingT GlobalT::int2LoggingT(int i)
{
	if (i == kVerbose)
		return kVerbose;
	else if (i == kModerate)
		return kModerate;
	else if (i == kSilent)
		return kSilent;
	else
		ExceptionT::GeneralFail("GlobalT::int2LoggingT", "could not translate %d", i);
	return kModerate;
}
