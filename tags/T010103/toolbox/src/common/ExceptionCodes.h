/* $Id: ExceptionCodes.h,v 1.8 2002-11-04 21:34:44 paklein Exp $ */
/* created: paklein (06/04/1996) */

#ifndef _EXCEPTION_CODES_H_
#define _EXCEPTION_CODES_H_

/** \file  
 * Backward compatibility for exception codes.
 *
 * \deprecated This file provides backward compatibility for the "old" style
 * of Tahoe exception enums and output. See ExceptionT for revised definitions
 * and methods having to do with exceptions.
 */
 
 #include "ExceptionT.h"
 
using namespace Tahoe;

/* number of exception codes */
#define eNumExceptions   ExceptionT::NumExceptions

/* see ExceptionT for definitions */
#define eNoError         ExceptionT::kNoError         // no error
#define eGeneralFail     ExceptionT::kGeneralFail     // general unrecoverable error
#define eStop            ExceptionT::kStop            // stop
#define eOutOfMemory     ExceptionT::kOutOfMemory     // out of memory
#define eOutOfRange      ExceptionT::kOutOfRange      // index range error
#define eSizeMismatch    ExceptionT::kSizeMismatch    // (array) dimension mismatch
#define eBadInputValue   ExceptionT::kBadInputValue   // bad input/construction parameter
#define eBadJacobianDet  ExceptionT::kBadJacobianDet  // ParentDomainT:bad jacobian determinant
#define eMPIFail         ExceptionT::kMPIFail         // general error on MPI call
#define eDatabaseFail    ExceptionT::kDatabaseFail    // general error reading/writing database
#define eBadHeartBeat    ExceptionT::kBadHeartBeat    // error detected on other processor
 
#endif
