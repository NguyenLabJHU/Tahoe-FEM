/* $Id: ExceptionT.h,v 1.7 2004-01-21 17:02:33 paklein Exp $ */
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
		kBadHeartBeat     =10, /**< error detected on other processor */
		kTypeMismatch     =11  /**< type mismatch */
	};

	/** convert ExceptionT::CodeT to a string */
	static const char* ToString(ExceptionT::CodeT code);

	/** write exception codes to output stream */
	static void WriteExceptionCodes(ostream& out);

	/** \name code to throw an exception, so you don't have to. This method writes
	 * some information about the exception to cout and then throws the appropriate
	 * ExceptionT::CodeT.
	 * \param caller the method throwing the exception, or NULL to skip.
	 * \param fmt format additional message to write, as would be passed to printf, sprintf, etc,
	 *        or NULL to skip */
	/*@{*/
	static void GeneralFail(const char* caller = NULL);
	static void GeneralFail(const char* caller, const char* fmt, ...);

	static void Stop(const char* caller = NULL);
	static void Stop(const char* caller, const char* fmt, ...);

	static void OutOfMemory(const char* caller = NULL);
	static void OutOfMemory(const char* caller, const char* fmt, ...);

	static void OutOfRange(const char* caller = NULL);
	static void OutOfRange(const char* caller, const char* fmt, ...);

	static void SizeMismatch(const char* caller = NULL);
	static void SizeMismatch(const char* caller, const char* fmt, ...);

	static void BadInputValue(const char* caller = NULL);
	static void BadInputValue(const char* caller, const char* fmt, ...);

	static void BadJacobianDet(const char* caller = NULL);
	static void BadJacobianDet(const char* caller, const char* fmt, ...);

	static void MPIFail(const char* caller = NULL);
	static void MPIFail(const char* caller, const char* fmt, ...);

	static void DatabaseFail(const char* caller = NULL);
	static void DatabaseFail(const char* caller, const char* fmt, ...);

	static void BadHeartBeat(const char* caller = NULL);
	static void BadHeartBeat(const char* caller, const char* fmt, ...);

	static void TypeMismatch(const char* caller = NULL);
	static void TypeMismatch(const char* caller, const char* fmt, ...);

	static void Throw(ExceptionT::CodeT code, const char* caller = NULL);
	static void Throw(ExceptionT::CodeT code, const char* caller, const char* fmt, ...);
	/*@}*/

	/** number of exception codes */
	static int NumExceptions;

  private:

	/** general throw
	 * \param code the exception being thrown
	 * \param caller the method throwing the exception, or NULL to skip.
	 * \param message message string, or NULL to skip */
	static void Throw_(ExceptionT::CodeT code, const char* caller, const char* message);
  
  	/** exception strings. One extra string to return "unknown". */
  	static const char* fExceptionStrings[13];
};

/* inline */
inline void ExceptionT::Throw(ExceptionT::CodeT code, const char* caller) { Throw_(code, caller, NULL); }
inline void ExceptionT::GeneralFail(const char* caller)    { Throw_(kGeneralFail, caller, NULL);    }
inline void ExceptionT::Stop(const char* caller)           { Throw_(kStop, caller, NULL);           }
inline void ExceptionT::OutOfMemory(const char* caller)    { Throw_(kOutOfMemory, caller, NULL);    }
inline void ExceptionT::OutOfRange(const char* caller)     { Throw_(kOutOfRange, caller, NULL);     }
inline void ExceptionT::SizeMismatch(const char* caller)   { Throw_(kSizeMismatch, caller, NULL);   }
inline void ExceptionT::BadInputValue(const char* caller)  { Throw_(kBadInputValue, caller, NULL);  }
inline void ExceptionT::BadJacobianDet(const char* caller) { Throw_(kBadJacobianDet, caller, NULL); }
inline void ExceptionT::MPIFail(const char* caller)        { Throw_(kMPIFail, caller, NULL);        }
inline void ExceptionT::DatabaseFail(const char* caller)   { Throw_(kDatabaseFail, caller, NULL);   }
inline void ExceptionT::BadHeartBeat(const char* caller)   { Throw_(kBadHeartBeat, caller, NULL);   }
inline void ExceptionT::TypeMismatch(const char* caller)   { Throw_(kTypeMismatch, caller, NULL);   }

} /* namespace Tahoe */
#endif /* _EXCEPTION_CODES_H_ */
