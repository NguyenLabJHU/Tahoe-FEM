C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   b_ax.f                                                                  */
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
C/* $Id: b_ax.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine db_ax(N,rowdista,rowdistb,nrhs,aptrs,ainds,avals,
     +                 b,ldb,x,ldx,myid,pp,emax,comm, tx,tb,bmax)
      implicit none
      include 'mpif.h'

      double precision zero
      parameter(zero=0.d0)

      integer N,rowdista(0:*),rowdistb(0:*),nrhs,aptrs(2,0:*),pp
      integer ldx,ldb
      integer ainds(*),comm,ierr,i,j,k,mynnodesb,mynnodesa,ofs,myid
      double precision avals(*),b(0:ldb-1,*),x(0:ldx-1,*),err,emax

C      double precision, allocatable :: tx(:,:),tb(:,:),bmax(:,:)
      double precision tx,tb,bmax
      dimension tx(0:N-1,nrhs)
      dimension tb(0:N-1,nrhs)
      dimension bmax(0:N-1,nrhs)

C      allocate(tx(0:N-1,nrhs),stat=k)
      if(k.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if
C      allocate(tb(0:N-1,nrhs),stat=k)
      if(k.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if
C     allocate(bmax(0:N-1,nrhs),stat=k)
      if(k.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if

      mynnodesb = rowdistb(myid+1)-rowdistb(myid)
      k = rowdistb(myid)

*      gather b.
      do j=1,nrhs
        do i=0,N-1
          bmax(i,j) = zero
        end do
      end do

      do j=1,nrhs
        do i=0,mynnodesb-1
          bmax(i+k,j) = b(i,j)
        end do
      end do

      call staged_mpirds(bmax,tb,N*nrhs,0,0,comm)

*      gather x.
      do j=1,nrhs
        do i=0,N-1
          bmax(i,j) = zero
        end do
      end do

      do j=1,nrhs
        do i=0,mynnodesb-1
          bmax(i+k,j) = x(i,j)
        end do
      end do

      call staged_mpirds(bmax,tx,N*nrhs,0,0,comm)

*      call b-Ax routine.
      do j=1,nrhs
        do i=0,N-1
          bmax(i,j) = zero
        end do
      end do

      mynnodesa = rowdista(myid+1)-rowdista(myid)
      ofs = rowdista(myid)

      call b_ax(aptrs,ainds,avals,tb,tx,N,mynnodesa,ofs,nrhs,bmax)

      call staged_mpirds(bmax,tx,N*nrhs,1,0,comm)

      if (myid.eq.0) then
        emax = 0.d0
        do j=1,nrhs
          do i=0,N-1
            err = dabs(tx(i,j))
            if(err.gt.emax) emax = err
          end do
        end do
c        print *,'max |B - AX| = ',emax
      end if

C      deallocate(tx)
C      deallocate(tb)
C      deallocate(bmax)

      end 

      subroutine b_ax(aptrs,ainds,avals,b,x,N,mynnodes,ofs,nrhs,ax)

      double precision avals(*),b(0:N-1,*),x(0:N-1,*),ax(0:N-1,*)
      integer aptrs(2,0:*),ainds(*),mynnodes,ofs

      do i=0,mynnodes-1
        ivptr = aptrs(1,i)
        ivsiz = aptrs(2,i)
        ig = i+ofs

        if(ivsiz.ne.0) then

          do while(ainds(ivptr).lt.ig)
            ivptr = ivptr+1
            ivsiz = ivsiz-1
          end do

          if(ainds(ivptr).eq.ig) then
            do k=1,nrhs
              ax(ig,k) = ax(ig,k) - avals(ivptr)*x(ig,k)
              ax(ig,k) = ax(ig,k) + b(ig,k)
            end do

            ivptr = ivptr+1
            ivsiz = ivsiz-1
          end if

          do k=1,nrhs
          do j=0,ivsiz-1
            node = ainds(ivptr+j)
            ax(node,k) = ax(node,k) - avals(ivptr+j)*x(ig,k)
            ax(ig,k) = ax(ig,k) - avals(ivptr+j)*x(node,k)
          end do
          end do

        end if

      end do
      end

      subroutine staged_mpirds(sbuf,rbuf,size,opt,root,comm)

      implicit none
      include 'mpif.h'

      integer eltlimit
      parameter(eltlimit=512*1024) ! corresponds to 4 MB 

      integer size,opt,root,comm
      double precision sbuf(*),rbuf(*)

      integer i,remains,ierr

      do i=1,size-eltlimit+1,eltlimit
        if(opt.eq.0) then
          call mpi_allreduce(sbuf(i),rbuf(i),eltlimit,
     +            MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        else
          call mpi_reduce(sbuf(i),rbuf(i),eltlimit,
     +            MPI_DOUBLE_PRECISION,MPI_SUM,root,comm,ierr)
        end if
      end do

      remains = size-i+1
      if(remains.gt.0) then
        if(opt.eq.0) then
          call mpi_allreduce(sbuf(i),rbuf(i),remains,
     +          MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
        else
          call mpi_reduce(sbuf(i),rbuf(i),remains,
     +            MPI_DOUBLE_PRECISION,MPI_SUM,root,comm,ierr)
        end if
      end if

      end
