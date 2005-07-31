C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   order.f                                                                 */
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
C/* $Id: order.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      subroutine porder(rowdist,aptrs,ainds,order,sizes,myid,pp,
     +      serialorder,dbgpp,comm)

      implicit none
      include 'mpif.h'

      integer rowdist(0:*),aptrs(2,0:*),ainds(*),order(0:*),sizes(0:*)
      integer, allocatable :: xadj(:),adjncy(:), sorder(:), counts(:)
      integer mynnodes,opts(5),i,j,k,l,m,ierr,ioasize,offdnz
      integer comm,myid,pp,serialorder
      integer dbgpp

      mynnodes = rowdist(myid+1)-rowdist(myid)

      ioasize = aptrs(2,mynnodes-1)+aptrs(1,mynnodes-1)-1

      allocate(xadj(0:mynnodes),stat=i)
      if(i.ne.0) then
        print *,myid,':memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(adjncy(0:ioasize-mynnodes-1),stat=i)
      if(i.ne.0) then
        print *,myid,':memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if

      ! form xadj and adjncy from aptrs, ainds

      m = rowdist(myid)
      l = 0
      do i=0,mynnodes-1
       xadj(i) = aptrs(1,i)-i-1
       do j=aptrs(1,i),aptrs(1,i)+aptrs(2,i)-1
         k = ainds(j)
         if(k.ne.m+i) then
           adjncy(l) = k
           l = l+1
         end if
       end do
      end do
      xadj(mynnodes) = xadj(mynnodes-1)+aptrs(2,mynnodes-1)-1

      call mpi_allreduce(xadj(mynnodes),offdnz,1,MPI_INTEGER,MPI_SUM,
     +                   comm,ierr);

      if(offdnz.eq.0) then
        
        do i=0,mynnodes-1
          order(i) = i + m
        end do

        do i=0,pp-2
          sizes(i) = rowdist(pp)/pp
        end do

        sizes(pp-1) = rowdist(pp) - (pp-1)*(rowdist(pp)/pp)

        do i=pp,2*pp-1
          sizes(i) = 0
        end do

      else

        opts(1) = 1
        opts(2) = 2
        opts(3) = dbgpp
        opts(4) = 0

        call parometisf(rowdist,xadj,adjncy,order,sizes,opts,
     +         serialorder,comm)

      end if

      deallocate(xadj)
      deallocate(adjncy)

      end
