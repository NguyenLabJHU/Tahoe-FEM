C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   crayblas.f                                                              */
C/*                                                                           */
C/*   Written by Mahesh Joshi, U of MN.                                       */
C/*                                                                           */
C/*****************************************************************************/
C/*                                                                           */
C/* This code is meant to be used solely for educational, research, and       */
C/* benchmarking purposes by non-profit institutions and US government        */
C/* agencies only.  Use by any other organization requires prior written      */
C/* permission from both IBM Corporation and the University of Minnesota.     */
C/* The software may not be sold or redistributed.  One may make copies       */
C/* of the software or modify it for their use provided that the copies,      */
C/* modified or otherwise, are not sold or distributed, are used under the    */
C/* same terms and conditions, and this notice and any part of the source     */
C/* code that follows this notice are not separated.                          */
C/*                                                                           */
C/* As unestablished research software, this code is provided on an           */
C/* ``as is'' basis without warranty of any kind, either expressed or         */
C/* implied, including but not limited to implied warranties of               */
C/* merchantability and fitness for a particular purpose.  IBM does not       */
C/* warrant that the functions contained in this software will meet the       */
C/* user's requirements or that the operation of its routines will be         */
C/* uninterrupted or error-free.  Acceptance and use of this program          */
C/* constitutes the user's understanding that he/she will have no recourse    */
C/* to IBM for any actual or consequential damages, including, but not        */
C/* limited to, lost profits or savings, arising out of the use or inability  */
C/* to use these libraries.  Even if the user informs IBM of the possibility  */
C/* of such damages, IBM expects the user to accept the risk of any such      */
C/* harm, or the user shall not attempt to use these libraries for any        */
C/* purpose.                                                                  */
C/*                                                                           */
C/* The downloading, compiling, or executing any part of this software        */
C/* constitutes an implicit agreement to these terms.  These terms and        */
C/* conditions are subject to change at any time without prior notice.        */
C/*                                                                           */
C/*****************************************************************************/
C/* $Id: crayblas.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      double precision function ddot(n,x,incx,y,incy)
      double precision ddot
      integer n,incx,incy
      double precision x(*),y(*)

      ddot = sdot(n,x,incx,y,incy)
      end

      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      character*1 transa,transb
      integer m,n,k,lda,ldb,ldc
      double precision alpha,beta,a(*),b(*),c(*)

      call       sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      end

      subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      character*1 trans
      integer m,n,lda,incx,incy
      double precision alpha,beta,a(*),x(*),y(*)

      call       sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      end

      subroutine dpotrf(uplo,n,a,lda,info)
      character*1 uplo
      integer n,lda,info
      double precision a(*)

      call spotrf(uplo,n,a,lda,info)
      end

      subroutine dscal(n,alpha,x,incx)
      integer n,incx
      double precision alpha,x(*)

      call       sscal(n,alpha,x,incx)
      end

      subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
      character*1 uplo,trans
      integer n,k,lda,ldc
      double precision alpha,beta,a(*),c(*)

      call       ssyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
      end

      subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
      character*1 side,uplo,transa,diag
      integer m,n,lda,ldb
      double precision a(*),b(*),alpha

      call       strsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
      end

      subroutine dtrsv(uplo,trans,diag,n,a,lda,x,incx)
      character*1 uplo,trans,diag
      integer n,lda,incx
      double precision a(*),x(*)

      call       strsv(uplo,trans,diag,n,a,lda,x,incx)
      end

