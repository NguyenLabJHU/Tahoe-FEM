C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   alltoallu_hc.f                                                          */
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
C/* $Id: alltoallu_hc.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      subroutine all_to_all_union_hc ( src,tgt,tmp,mxsize,
     +                               srcsize,tgtsize,myid,dim,comm)
      
      implicit none

      include 'mpif.h'

      integer itag
      parameter(itag=1)

      integer srcsize,tgtsize,mxsize,dim,myid,ipartner,comm
      integer src(*), tgt(*), tmp(*)
      integer i,j,tmpsize,nbrecv
      integer ierr,mpistat(MPI_STATUS_SIZE)

      do i=1,srcsize
        tgt(i) = src(i)
      end do
      tgtsize = srcsize

      do i=0,dim-1

        ipartner = ieor(myid,ishft(1,i))

        if(ipartner.gt.myid) then
          call mpi_send(tgt,tgtsize,MPI_INTEGER,ipartner,itag,
     +                      comm,ierr)
          call mpi_recv(tmp,mxsize,MPI_INTEGER,ipartner,itag,
     +                  comm,mpistat,ierr)
        else
          call mpi_recv(tmp,mxsize,MPI_INTEGER,ipartner,itag,
     +                  comm,mpistat,ierr)
          call mpi_send(tgt,tgtsize,MPI_INTEGER,ipartner,itag,
     +                      comm,ierr)
        end if

        call mpi_get_count(mpistat,MPI_INTEGER,tmpsize,ierr)
        
        call mergelists(tmp,tmpsize,tgt,tgtsize,src,srcsize)

        do j=1,srcsize
          tgt(j) = src(j)
        end do
        tgtsize = srcsize

      end do
      

      end

      subroutine mergelists(as1,len1,as2,len2,at,lent)

      implicit none

      integer len1,len2,lent
      integer as1(*), as2(*), at(*)

      integer i,j,k

      i=1
      j=1
      k=1
      do while( (i.le.len1) .and. (j.le.len2) )
        if(as1(i) .lt. as2(j)) then
          at(k) = as1(i)
          i = i + 1
        else
          if(as1(i) .eq. as2(j)) then
            at(k) = as1(i)
            i = i + 1
            j = j + 1
          else
            at(k) = as2(j)
            j = j + 1
          end if
        end if
        k = k + 1
      end do
      do i=i,len1
        at(k) = as1(i)
        k = k + 1
      end do
      do j=j,len2
        at(k) = as2(j)
        k = k + 1
      end do
      lent = k-1

      end
