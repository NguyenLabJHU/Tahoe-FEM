C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   moveai.f                                                                */
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
C/* $Id: moveai.f,v 1.1 2004-12-10 20:26:45 paklein Exp $ */
C/*****************************************************************************/

      subroutine moveai(N,dd,pp,lgblk,myid,mynnodes,order,paptrs,
     +                  recvsizs,painds,aptrs,tainds,wrkint,ranmasks,
     +                  whichsnode,nsend,nrecvsizs,comm)

      implicit none

      include 'mpif.h'

      integer*8 loc

      integer order(0:*),aptrs(2,0:*),paptrs(2,0:*),painds(*)
      integer tainds(*),wrkint(0:*),whichsnode(0:*),ranmasks(5,0:*)
      integer N,dd,pp,lgblk,myid,mynnodes,nsend,comm

      integer, allocatable :: sendinds(:),sendsizs(:)
      integer recvsizs(0:*)

      integer proc,pgrsize,ierr,bmaskr,bmaskc,row,col
      integer i,j,k,l,m,ptr_r,fptr_r,ptr_c
      integer iwillsend_sizs,is1,ptr_sendsizs
      integer nrecvsizs,ptr_sendinds

      integer psci,psdi,prci,prdi,pscs,psds,prcs,prds,ppr,ppc,ppg

      psci = 0
      psdi = pp
      prci = 2*pp
      prdi = 3*pp
      pscs = 4*pp
      psds = 5*pp
      prcs = 6*pp
      prds = 7*pp
      ppr  = 8*pp
      ppc  = 9*pp
      ppg  = 10*pp

      pgrsize = ishft(1,ishft(dd,-1))

      wrkint(psdi) = 0
      wrkint(psds) = 0
      do proc=1,pp-1
        wrkint(psdi+proc) = wrkint(psdi+proc-1)+wrkint(psci+proc-1)
        wrkint(psds+proc) = wrkint(psds+proc-1)+wrkint(pscs+proc-1)
        wrkint(psci+proc-1) = 0
        wrkint(pscs+proc-1) = 0
      end do

      nsend = wrkint(psdi+pp-1)+wrkint(psci+pp-1)
      iwillsend_sizs = wrkint(psds+pp-1)+wrkint(pscs+pp-1)

      wrkint(psci+pp-1) = 0
      wrkint(pscs+pp-1) = 0

      allocate(sendinds(0:nsend-1),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if
      allocate(sendsizs(0:iwillsend_sizs-1),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if

      do i=0,mynnodes-1
        col = order(i)
        j = whichsnode(i)
        proc   = ranmasks(3,j)
        bmaskr = ranmasks(4,j)
        bmaskc = ranmasks(5,j)

        fptr_r = wrkint(ppr+proc)
        ptr_c  = wrkint(ppc+proc)+iand(ishft(col,-lgblk),bmaskc)

        do k=0,bmaskr
          proc = wrkint(ppg+ptr_c*pgrsize+fptr_r+k)
          ptr_sendsizs = wrkint(psds+proc)+wrkint(pscs+proc)
          sendsizs(ptr_sendsizs) = col
          wrkint(prcs+proc) = ptr_sendsizs+1
          sendsizs(ptr_sendsizs+1) = 0
          wrkint(pscs+proc) = wrkint(pscs+proc)+2
        end do

        do j=aptrs(1,i),aptrs(1,i)+aptrs(2,i)-1
          row = tainds(j)
          ptr_r = fptr_r+iand(ishft(row,-lgblk),bmaskr)
          proc = wrkint(ppg+ptr_c*pgrsize+ptr_r)
          ptr_sendinds = wrkint(psdi+proc)+wrkint(psci+proc)
          sendinds(ptr_sendinds) = row
          sendsizs(wrkint(prcs+proc)) = sendsizs(wrkint(prcs+proc))+1
          wrkint(psci+proc) = wrkint(psci+proc)+1
        end do
      end do 

      call mpi_alltoall(wrkint(pscs),1,MPI_INTEGER,wrkint(prcs),
     +                  1,MPI_INTEGER,comm,ierr)

      wrkint(prds) = 0
      do proc=1,pp-1
        wrkint(prds+proc) = wrkint(prds+proc-1)+wrkint(prcs+proc-1)
      end do
      nrecvsizs = wrkint(prds+pp-1)+wrkint(prcs+pp-1)

      call mpi_alltoallv(sendinds,wrkint(psci),
     +                   wrkint(psdi),MPI_INTEGER,painds,
     +                   wrkint(prci),wrkint(prdi),
     +                   MPI_INTEGER,comm,ierr)

      call mpi_alltoallv(sendsizs,wrkint(pscs),
     +                   wrkint(psds),MPI_INTEGER,recvsizs,
     +                   wrkint(prcs),wrkint(prds),
     +                   MPI_INTEGER,comm,ierr)

      do i=0,N-1
        paptrs(1,i) = 1
        paptrs(2,i) = 0
      end do

      i = 0
      j = 1
      do while (i.lt.nrecvsizs)
        col = recvsizs(i)
        paptrs(2,col) = recvsizs(i+1)
        i = i+2
        paptrs(1,col) = j
        j = j+paptrs(2,col)
      end do

      deallocate(sendinds)
      deallocate(sendsizs)

      end

