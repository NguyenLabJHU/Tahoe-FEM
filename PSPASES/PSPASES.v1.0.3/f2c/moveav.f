C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   moveav.f                                                                */
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
C/* $Id: moveav.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine moveav(N,dd,pp,lgblk,myid,rowdista,mynnodes,
     +                  order,aptrs,ainds,avals,pavals,wrkint,
     +                  maxnzpercol,ranmasks,comm)

      implicit none

      include 'mpif.h'

      integer*8 loc

      integer rowdista(0:*),order(0:*),aptrs(2,0:*),ainds(*)
      integer wrkint(0:*),ranmasks(5,0:*)
      integer N,dd,pp,lgblk,myid,mynnodes,maxnzpercol,comm
      double precision avals(*),pavals(*)

      integer, allocatable :: gorder(:),whichsnode(:),tainds(:)
      double precision, allocatable :: sendvals(:)

      integer proc,pgrsize,ierr,bmaskr,bmaskc,row,col
      integer i,j,k,l,m,ptr_r,fptr_r,ptr_c,itainds
      integer is1,nsend,ptr_sendvals
      integer beginrow

      integer pscv,psdv,prcv,prdv,ppr,ppc,ppg

      pscv = 0
      psdv = pp
      prcv = 2*pp
      prdv = 3*pp
      ppr  = 8*pp
      ppc  = 9*pp
      ppg  = 10*pp

      pgrsize = ishft(1,ishft(dd,-1))

      allocate(gorder(0:N-1),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if

      beginrow = rowdista(myid)

      do proc=0,pp-1
        wrkint(pscv+proc) = rowdista(proc+1)-rowdista(proc)
      end do

      call mpi_allgatherv(order,mynnodes,MPI_INTEGER,
     +                    gorder,wrkint(pscv),rowdista,MPI_INTEGER,
     +                    comm,ierr)

      do proc=0,pp-1
        wrkint(pscv+proc) = 0
      end do

      allocate(whichsnode(0:mynnodes-1),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if
      i = aptrs(1,mynnodes-1)+aptrs(2,mynnodes-1)-1
      allocate(tainds(2*i),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if

      itainds = 1
      do i=0,mynnodes-1

        col = gorder(beginrow+i)
        j = 0
        do while (col.lt.ranmasks(1,j) .or. col.gt.ranmasks(2,j))
          j = j+1
        end do
        proc   = ranmasks(3,j)
        bmaskr = ranmasks(4,j)
        bmaskc = ranmasks(5,j)

        fptr_r = wrkint(ppr+proc)
        ptr_c  = wrkint(ppc+proc)+iand(ishft(col,-lgblk),bmaskc)

        k = aptrs(1,i)
        l = aptrs(2,i)

        do j=0,l-1
          tainds(itainds+j) = gorder(ainds(k+j))
        end do

        call ikeysortf(l,tainds(itainds),tainds(itainds+l))

        m = 0
        do while(tainds(itainds+m).lt.col)
          m = m+1
        end do
        whichsnode(i) = m

        do j=m,l-1
          row = tainds(itainds+j)
          ptr_r = fptr_r+iand(ishft(row,-lgblk),bmaskr)
          proc = wrkint(ppg+ptr_c*pgrsize+ptr_r)
          wrkint(pscv+proc) = wrkint(pscv+proc)+1
          tainds(itainds+j) = proc
        end do

        itainds = itainds+2*l

      end do 

      deallocate(gorder)

      wrkint(psdv) = 0
      do proc=1,pp-1
        wrkint(psdv+proc) = wrkint(psdv+proc-1)+wrkint(pscv+proc-1)
        wrkint(pscv+proc-1) = 0
      end do
      nsend = wrkint(psdv+pp-1)+wrkint(pscv+pp-1)
      wrkint(pscv+pp-1) = 0

      allocate(sendvals(0:nsend-1),stat=is1)
      if(is1.ne.0) then
        print *,'Error in allocate'
        call mpi_abort(comm,1,ierr)
      end if

      itainds = 1
      do i=0,mynnodes-1

        k = aptrs(1,i)
        l = aptrs(2,i)

        m = whichsnode(i)
        do j=m,l-1
          proc = tainds(itainds+j)
          ptr_sendvals = wrkint(psdv+proc)+wrkint(pscv+proc)
          sendvals(ptr_sendvals) = avals(k+tainds(itainds+l+j))
          wrkint(pscv+proc) = wrkint(pscv+proc)+1
        end do

        itainds = itainds+2*l

      end do 

      deallocate(whichsnode)
      deallocate(tainds)

      call mpi_alltoall(wrkint(pscv),1,MPI_INTEGER,wrkint(prcv),
     +                  1,MPI_INTEGER,comm,ierr)

      wrkint(prdv) = 0
      do proc=1,pp-1
        wrkint(prdv+proc) = wrkint(prdv+proc-1)+wrkint(prcv+proc-1)
      end do

      call mpi_alltoallv(sendvals,wrkint(pscv),
     +                   wrkint(psdv),MPI_DOUBLE_PRECISION,pavals,
     +                   wrkint(prcv),wrkint(prdv),
     +                   MPI_DOUBLE_PRECISION,comm,ierr)

      deallocate(sendvals)

      end

