/* $Id: GlobalT.cpp,v 1.8 2003-06-09 07:01:17 paklein Exp $ */
/* created: paklein (04/01/2000) */

#include "GlobalT.h"
#include <iostream.h>
#include "ExceptionT.h"


using namespace Tahoe;

namespace Tahoe {

istream& operator>>(istream& in, GlobalT::AnalysisCodeT& code)
{
	int i_code;
	in >> i_code;
	switch (i_code)
	{
		case GlobalT::kNoAnalysis:
			code = GlobalT::kNoAnalysis;
			break;
		case GlobalT::kLinStatic:
			code = GlobalT::kLinStatic;
			break;
		case GlobalT::kLinDynamic:
			code = GlobalT::kLinDynamic;
			break;
		case GlobalT::kNLStatic:
			code = GlobalT::kNLStatic;
			break;
		case GlobalT::kNLDynamic:
			code = GlobalT::kNLDynamic;
			break;
		case GlobalT::kDR:
			code = GlobalT::kDR;
			break;
		case GlobalT::kLinExpDynamic:
			code = GlobalT::kLinExpDynamic;
			break;
		case GlobalT::kNLExpDynamic:
			code = GlobalT::kNLExpDynamic;
			break;
		case GlobalT::kCBStatic:
		{
			cout << "\n operator>>GlobalT::AnalysisCodeT: Cauchy-Born BC's converted to KBC controller" << endl;
			throw ExceptionT::kBadInputValue;
		}
		case GlobalT::kNLStaticKfield:
		{
			cout << "\n operator>>GlobalT::AnalysisCodeT: K-field converted to KBC controller" << endl;
			throw ExceptionT::kBadInputValue;
		}
		case GlobalT::kVarNodeNLStatic:
		case GlobalT::kVarNodeNLExpDyn:
		{
			cout << "\n operator>>GlobalT::AnalysisCodeT: analysis code is not longer\n"
			     <<   "     supported. Support for changing geometry is being re-\n"
			     <<   "     written: " << i_code << endl;
			throw ExceptionT::kBadInputValue;
		}
		case GlobalT::kAugLagStatic:
		{
			cout << "\n operator>>GlobalT::AnalysisCodeT: external degrees of freedom no longer\n" 
			     <<   "     require a specific analysis code: " << GlobalT::kAugLagStatic << endl;
			throw ExceptionT::kBadInputValue;
		}
		case GlobalT::kNLExpDynKfield:
		{
			cout << "\n operator>>GlobalT::AnalysisCodeT: K-field converted to KBC controller" << endl;
			throw ExceptionT::kBadInputValue;
		}
		case GlobalT::kLinStaticHeat:
			code = GlobalT::kLinStaticHeat;
			break;
		case GlobalT::kLinTransHeat:
			code = GlobalT::kLinTransHeat;
			break;
		case GlobalT::kNLStaticHeat:
			code = GlobalT::kNLStaticHeat;
			break;
		case GlobalT::kNLTransHeat:
			code = GlobalT::kNLTransHeat;
			break;
		case GlobalT::kPML:
			code = GlobalT::kPML;
			break;
		case GlobalT::kMultiField:
			code = GlobalT::kMultiField;
			break;
		default:
			cout << "\n operator>>GlobalT::AnalysisCodeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}
	return in;
}

}

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
	{
		cout << "\n GlobalT::MaxPrecedence: unexpected combination" << endl;	
		throw ExceptionT::kGeneralFail;
	}

	return result;
}
