/* $Id: ExceptionCodes.h,v 1.9 2003-05-04 22:52:51 paklein Exp $ */
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
 
/* number of exception codes */
#define eNumExceptions   Tahoe::ExceptionT::NumExceptions

/* see ExceptionT for definitions */
#define eNoError         Tahoe::ExceptionT::kNoError         // no error
#define eGeneralFail     Tahoe::ExceptionT::kGeneralFail     // general unrecoverable error
#define eStop            Tahoe::ExceptionT::kStop            // stop
#define eOutOfMemory     Tahoe::ExceptionT::kOutOfMemory     // out of memory
#define eOutOfRange      Tahoe::ExceptionT::kOutOfRange      // index range error
#define eSizeMismatch    Tahoe::ExceptionT::kSizeMismatch    // (array) dimension mismatch
#define eBadInputValue   Tahoe::ExceptionT::kBadInputValue   // bad input/construction parameter
#define eBadJacobianDet  Tahoe::ExceptionT::kBadJacobianDet  // ParentDomainT:bad jacobian determinant
#define eMPIFail         Tahoe::ExceptionT::kMPIFail         // general error on MPI call
#define eDatabaseFail    Tahoe::ExceptionT::kDatabaseFail    // general error reading/writing database
#define eBadHeartBeat    Tahoe::ExceptionT::kBadHeartBeat    // error detected on other processor
#define eTypeMismatch    Tahoe::ExceptionT::kTypeMismatch    // type mismatch
 
#endif
