
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   reordx.f                                                                */
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
C/* $Id: reordx.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine reordx(N,x,ldx,nrhs,rx,nown,dsc,dsd,drc,drd,
     +                  iown,rowdist,mynnodes,sloc,rloc,svals,
     +                  rvals,myid,comm)

      implicit none
      include 'mpif.h'

      integer N,nown,dsc(0:*),dsd(0:*),ldx
      integer drc(0:*),drd(0:*),mynnodes
      integer sloc(0:*),iown(0:*),rloc(0:*),rowdist(0:*)
      integer i,j,nrhs,comm
      double precision x(0:ldx-1,0:*),rx(0:N-1,0:*)
      double precision svals(0:*),rvals(0:*)
      integer myid

      do j=0,nrhs-1
        do i=0,nown-1
          svals(sloc(i)*nrhs+j) = rx(iown(2*i),j)
        end do
      end do

      call mpi_alltoallv(svals,dsc,dsd,MPI_DOUBLE_PRECISION,
     +                   rvals,drc,drd,MPI_DOUBLE_PRECISION,
     +                   comm,i)

      do j=0,nrhs-1
        do i=0,mynnodes-1
          x(rloc(i)-rowdist(myid),j) = rvals(i*nrhs+j)
        end do
      end do

      end

