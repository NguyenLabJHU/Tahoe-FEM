C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   ygentree.f                                                              */
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
C/* $Id: ygentree.f,v 1.1 2004-12-10 20:26:45 paklein Exp $ */
C/*****************************************************************************/

      subroutine ygentree(N,aptrs,ainds,beginleafnode,sizes,
     +                   parent,tptrs,tinds,dd,myid,pp,wrkint,
     +                   asize,recvsizs,nrsiz,comm)

      implicit none

      include 'mpif.h'

      integer N,beginleafnode,myid,dd,pp,asize,nrsiz,comm
      integer aptrs(2,0:*),ainds(*),sizes(0:*)
      integer parent(0:*),tptrs(3,0:*),tinds(*),recvsizs(0:*)
      integer wrkint(0:*),i,j,k,l,m,pn,ierr,endnode,endnodepar,col

      !assuming tinds size is >= N, it is used as a temporary array here

      endnode = beginleafnode+sizes(myid)-1
      endnodepar = parent(endnode)
      call genleaftree(beginleafnode,sizes(myid),parent,aptrs,ainds,
     +                 tinds)

      do i=beginleafnode,endnode-1
        if(parent(i).eq.-1) then
          parent(i) = endnode
        end if
      end do
      parent(endnode) = endnodepar

      call mpi_allgather(beginleafnode,1,MPI_INTEGER,
     +                   wrkint,1,MPI_INTEGER,comm,ierr)

c     +                   wrkint(N),1,MPI_INTEGER,comm,ierr)
      call mpi_allgatherv(parent(beginleafnode),sizes(myid),
     +                    MPI_INTEGER,parent,sizes,wrkint,
     +                    MPI_INTEGER,comm,ierr)

c     +                    MPI_INTEGER,parent,sizes,wrkint(N),
      do i=0,N-1
        tptrs(2,i) = 0
      end do

      !accumulate kid counts (tptrs(2,*))
      do i=0,N-2
        j = parent(i)
        tptrs(2,j) = tptrs(2,j)+1
      end do

      !form displacements array (tptrs(1,*))
      j = 1
      do i=0,N-1
        tptrs(1,i) = j
        j = j+tptrs(2,i)
        tptrs(2,i) = 0
      end do

      !fill tinds
      do i=0,N-2
        j = parent(i)
        tinds(tptrs(1,j)+tptrs(2,j)) = i
        tptrs(2,j) = tptrs(2,j)+1
      end do

      pn = 1

      i = 0
      do while (i.lt.nrsiz)
        col = recvsizs(i)
        j = aptrs(2,col)
        if(j.gt.0) then
          k = aptrs(1,col)
          l = 0
          do while (ainds(k+l).lt.col)
            l = l+1
          end do
          m = 0
          do l=l,j-1
            ainds(pn+m) = ainds(k+l)
            m = m+1
          end do
          aptrs(1,col) = pn
          aptrs(2,col) = m
          pn = pn+m
        end if
        i = i+2
      end do
      asize = pn-1

      end
