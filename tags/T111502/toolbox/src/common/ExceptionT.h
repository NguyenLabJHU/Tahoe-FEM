/* $Id: ExceptionT.h,v 1.4 2002-11-13 16:49:31 paklein Exp $ */
/* created: paklein (06/04/1996) */
#ifndef _EXCEPTION_T_H_
#define _EXCEPTION_T_H_

#include "ios_fwd_decl.h"
#include <stdlib.h>

namespace Tahoe {

/** enums and strings for exception codes */
class ExceptionT
{
  public:

	/** exception codes */
	enum CodeT {
		kNoError          = 0, /**< no error */
		kGeneralFail      = 1, /**< general unrecoverable error */
		kStop             = 2, /**< stop */ 
		kOutOfMemory      = 3, /**< out of memory */
		kOutOfRange       = 4, /**< index range error */
		kSizeMismatch     = 5, /**< (array) dimension mismatch */
		kBadInputValue    = 6, /**< bad input/construction parameter */
		kBadJacobianDet	  = 7, /**< ParentDomainT: bad jacobian determinant */
		kMPIFail          = 8, /**< general error on MPI call */
		kDatabaseFail     = 9, /**< general error reading/writing database */
		kBadHeartBeat     =10  /**< error detected on other processor */
	};

	/** convert ExceptionT::CodeT to a string */
	static const char* ToString(ExceptionT::CodeT code);

	/** write exception codes to output stream */
	static void WriteExceptionCodes(ostream& out);

	/** \name code to throw an exception, so you don't have to. This method writes
	 * some information about the exception to cout and then throws the appropriate
	 * ExceptionT::CodeT.
	 * \param caller the function that threw the exception
	 * \param message additional message to write */
	/*@{*/
	static void GeneralFail(const char* caller = NULL, const char* message = NULL);
	static void Stop(const char* caller = NULL, const char* message = NULL);
	static void OutOfMemory(const char* caller = NULL, const char* message = NULL);
	static void OutOfRange(const char* caller = NULL, const char* message = NULL);
	static void SizeMismatch(const char* caller = NULL, const char* message = NULL);
	static void BadInputValue(const char* caller = NULL, const char* message = NULL);
	static void BadJacobianDet(const char* caller = NULL, const char* message = NULL);
	static void MPIFail(const char* caller = NULL, const char* message = NULL);
	static void DatabaseFail(const char* caller = NULL, const char* message = NULL);
	static void BadHeartBeat(const char* caller = NULL, const char* message = NULL);

	/** general throw */
	static void Throw(ExceptionT::CodeT code, const char* caller, const char* message);
	/*@}*/

	/** number of exception codes */
	static int NumExceptions;

  private:
  
  	/** exception strings. One extra string to return "unknown". */
  	static const char* fExceptionStrings[12];
};

/* inline */
inline void ExceptionT::GeneralFail(const char* caller, const char* message)    { Throw(kGeneralFail, caller, message);    }
inline void ExceptionT::Stop(const char* caller, const char* message)           { Throw(kStop, caller, message);           }
inline void ExceptionT::OutOfMemory(const char* caller, const char* message)    { Throw(kOutOfMemory, caller, message);    }
inline void ExceptionT::OutOfRange(const char* caller, const char* message)     { Throw(kOutOfRange, caller, message);     }
inline void ExceptionT::SizeMismatch(const char* caller, const char* message)   { Throw(kSizeMismatch, caller, message);   }
inline void ExceptionT::BadInputValue(const char* caller, const char* message)  { Throw(kBadInputValue, caller, message);  }
inline void ExceptionT::BadJacobianDet(const char* caller, const char* message) { Throw(kBadJacobianDet, caller, message); }
inline void ExceptionT::MPIFail(const char* caller, const char* message)        { Throw(kMPIFail, caller, message);        }
inline void ExceptionT::DatabaseFail(const char* caller, const char* message)   { Throw(kDatabaseFail, caller, message);   }
inline void ExceptionT::BadHeartBeat(const char* caller, const char* message)   { Throw(kBadHeartBeat, caller, message);   }

} /* namespace Tahoe */
#endif /* _EXCEPTION_CODES_H_ */
