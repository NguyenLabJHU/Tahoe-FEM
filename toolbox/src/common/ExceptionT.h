/* $Id: ExceptionT.h,v 1.1.2.1 2002-10-17 01:41:18 paklein Exp $ */
/* created: paklein (06/04/1996) */
#ifndef _EXCEPTION_T_H_
#define _EXCEPTION_T_H_

#include "ios_fwd_decl.h"

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

	/** number of exception codes */
	static int NumExceptions;

  private:
  
  	/** exception strings. One extra string to return "unknown". */
  	static const char* fExceptionStrings[12];
};
} /* namespace Tahoe */
#endif /* _EXCEPTION_CODES_H_ */
