/* $Id: ExceptionCodes.h,v 1.5 2002-08-21 07:15:12 paklein Exp $ */
/* created: paklein (06/04/1996) */

#ifndef _EXCEPTION_CODES_H_
#define _EXCEPTION_CODES_H_

namespace Tahoe {

/* number of exception codes */
const int eNumExceptions	= 11;

/* exception codes */
const int eNoError          = 0; // no error
const int eGeneralFail		= 1; // general unrecoverable error
const int eStop             = 2; // stop
const int eOutOfMemory      = 3; // out of memory
const int eOutOfRange       = 4; // index range error
const int eSizeMismatch     = 5; // (array) dimension mismatch
const int eBadInputValue    = 6; // bad input/construction parameter
const int eBadJacobianDet	= 7; // ParentDomainT:bad jacobian determinant
const int eMPIFail          = 8; // general error on MPI call
const int eDatabaseFail     = 9; // general error reading/writing database
const int eBadHeartBeat     =10; // error detected on other processor

} // namespace Tahoe 
#endif /* _EXCEPTION_CODES_H_ */
