/* $Id: ExceptionT.cpp,v 1.3 2002-11-13 08:28:14 paklein Exp $ */
#include "ExceptionT.h"
#include "ArrayT.h"
#include <iostream.h>
#include <iomanip.h>
#include <time.h>

/* initialize static data */
namespace Tahoe {
int ExceptionT::NumExceptions = 11;
const bool ArrayT<ExceptionT::CodeT>::fByteCopy = true;

/* exceptions strings */
const char* ExceptionT::fExceptionStrings[12] = 
{
/* 0 */ "no error",
/* 1 */ "general fail",
/* 2 */ "stop",
/* 3 */ "out of memory",
/* 4 */ "index out of range",
/* 5 */ "dimension mismatch",
/* 6 */ "invalid value read from input",
/* 7 */ "zero or negative jacobian",
/* 8 */ "MPI message passing error",
/* 9 */ "database read failure",
/*10 */ "bad MP heartbeat",
/*11 */ "unknown"};

/* write exception codes to output stream */
void ExceptionT::WriteExceptionCodes(ostream& out)
{
	out << "\nE x c e p t i o n   c o d e s :\n\n";	
	for (int i = 0; i < NumExceptions; i++)
	{
		out << setw(kIntWidth) << i << " : ";
		out << fExceptionStrings[i] << '\n';
	}	
		out << endl;
}

/* return exception string */
const char* ExceptionT::ToString(CodeT code)
{
	if (code >= 0 && code < NumExceptions)
		return fExceptionStrings[code];
	else
		return fExceptionStrings[NumExceptions];
}

void ExceptionT::Throw(ExceptionT::CodeT code, const char* caller, const char* message)
{
	/* write info */
	time_t t;
	time(&t);
	cout << "\n ExceptionT::Throw: " << ctime(&t) << '\n';
	cout << "      code: " << code << '\n';
	cout << " exception: " << ToString(code) << '\n';
	if (caller)  cout << "    caller: " << caller << '\n';
	if (message) cout << "   message: " << message << '\n';
	cout.flush();
	
	/* do the throw */
	throw code;
}

} /* namespace Tahoe */ 
