
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   preordbe.f                                                              */
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
C/* $Id: preordbe.f,v 1.1 2004-12-10 20:26:45 paklein Exp $ */
C/*****************************************************************************/

      subroutine preordbe(N,order,ranmasks,nown,isc,isd,irc,ird,
     +                    mynnodes,dd,myid,lgblk,pgrid,whichsnode,comm)
      implicit none
      include 'mpif.h'

      integer N,order(0:*),ranmasks(5,0:*),nown,isc(0:*),isd(0:*)
      integer irc(0:*),ird(0:*),mynnodes,dd,lgblk,comm
      integer pgrid(0:*),whichsnode(0:*)
      integer rbits,cbits,pgrsize,ppr,ppc,ppg,kr,kc,i,j,k,col,pp
      integer myid

      pp = ishft(1,dd)
      rbits = ishft(dd,-1)
      cbits = dd-rbits
      pgrsize = ishft(1,rbits)

      ppr = 0
      ppc = pp
      ppg = 2*pp

      do k=0,pp-1
        kr = 0
        i = ishft(k,-1)
        do j=0,rbits-1
          kr = ior(kr,ishft(iand(i,1),j))
          i = ishft(i,-2)
        end do
        kc = 0
        i = k
        do j=0,cbits-1
          kc = ior(kc,ishft(iand(i,1),j))
          i = ishft(i,-2)
        end do
        pgrid(ppr+k) = kr
        pgrid(ppc+k) = kc
        pgrid(ppg+kc*pgrsize+kr) = k
      end do

      do i=0,pp-1
        isc(i) = 0
      end do

      do i=0,mynnodes-1
        col = order(i)
        j = 0
        do while (col.lt.ranmasks(1,j) .or. col.gt.ranmasks(2,j))
          j = j+1
        end do
        whichsnode(i) = j
        k = ranmasks(3,j)
        kr = pgrid(ppr+k)+iand(ishft(col,-lgblk),ranmasks(4,j))
        kc = pgrid(ppc+k)+iand(ishft(col,-lgblk),ranmasks(5,j))
        k = pgrid(ppg+kc*pgrsize+kr)
        isc(k) = isc(k)+2
      end do

      isd(0) = 0
      do i=1,pp-1
        isd(i) = isd(i-1)+isc(i-1)
      end do

      call mpi_alltoall(isc,1,MPI_INTEGER,irc,1,MPI_INTEGER,comm,i)

      ird(0) = 0
      nown = 0
      do i=1,pp-1
        ird(i) = ird(i-1)+irc(i-1)
        nown = nown+irc(i-1)
      end do
      nown = nown+irc(pp-1)

      nown = ishft(nown,-1)

      end
