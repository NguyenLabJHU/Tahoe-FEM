/* $Id: ExceptionT.cpp,v 1.13 2004-02-05 18:42:35 paklein Exp $ */
#include "ExceptionT.h"
#include "ArrayT.h"
#include <iostream.h>
#include <iomanip.h>
#include <time.h>

#if defined(__SGI__) || defined(__DELMAR__) || defined(__PGI__) || defined(__JANUS__) || defined(__MWERKS__)
#include <stdio.h>
#include <stdarg.h>
#else
#include <cstdio>
#include <cstdarg> /* to handle the variable number of arguments in Throw() */
#endif

/* initialize static data */
namespace Tahoe {
int ExceptionT::NumExceptions = 12;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ExceptionT::CodeT>::fByteCopy = true;

/* exceptions strings */
const char* ExceptionT::fExceptionStrings[13] = 
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
/*11 */ "type mismatch",
/*12 */ "unknown"};

/* buffer for vsprintf */
static char message_buffer[255];

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

/* repeated here instead of passing a va_list to Throw so that the
 * cstdarg header is limited to this source file */
void ExceptionT::GeneralFail(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kGeneralFail, caller, message_buffer);
}

void ExceptionT::Stop(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kStop, caller, message_buffer);
}

void ExceptionT::OutOfMemory(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kOutOfMemory, caller, message_buffer);
}

void ExceptionT::OutOfRange(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kOutOfRange, caller, message_buffer);
}

void ExceptionT::SizeMismatch(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kSizeMismatch, caller, message_buffer);
}

void ExceptionT::BadInputValue(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kBadInputValue, caller, message_buffer);
}

void ExceptionT::BadJacobianDet(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kBadJacobianDet, caller, message_buffer);
}

void ExceptionT::MPIFail(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kMPIFail, caller, message_buffer);
}

void ExceptionT::DatabaseFail(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kDatabaseFail, caller, message_buffer);
}

void ExceptionT::BadHeartBeat(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kBadHeartBeat, caller, message_buffer);
}

void ExceptionT::TypeMismatch(const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(kTypeMismatch, caller, message_buffer);
}

void ExceptionT::Throw(ExceptionT::CodeT code, const char* caller, const char* fmt, ...)
{
	va_list argp;
	va_start(argp, fmt);
	vsprintf(message_buffer, fmt, argp);
	va_end(argp);
	Throw_(code, caller, message_buffer);
}

void ExceptionT::Throw_(ExceptionT::CodeT code, const char* caller, const char* message)
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* write info */
	time_t t;
	time(&t);
	cout << "\n ExceptionT::Throw: " << ctime(&t);
	cout <<   "        code: " << code << '\n';
	cout <<   "   exception: " << ToString(code) << '\n';
	if (caller)
	cout <<   "      caller: " << caller << '\n';
	if (message)
	cout <<   "     message: " << message << '\n';
#endif
	
	/* do the throw */
	throw code;
}

} /* namespace Tahoe */ 
