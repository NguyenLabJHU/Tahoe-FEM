
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   preordxc.f                                                              */
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
C/* $Id: preordxc.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine preordxc(N,iown,rowdist,mynnodes,isc,isd,
     +                    irc,ird,nown,slocx,rlocx,siprm,
     +                    dd,myid,whichproc,comm)
      implicit none
      include 'mpif.h'

      integer N,rowdist(0:*),nown,isc(0:*),isd(0:*)
      integer irc(0:*),ird(0:*),mynnodes,dd,comm
      integer whichproc(0:*),iown(0:*),slocx(0:*),rlocx(0:*)
      integer i,j,k,col,pp,siprm(0:*)
      integer myid

      pp = ishft(1,dd)

      do i=0,pp-1
        isc(i) = 0
      end do

      do i=0,nown-1
	col = iown(2*i+1)
	j = 1
        do while(col.ge.rowdist(j))
	  j = j+1
        end do
	k = j-1
	whichproc(i) = j-1
	isc(k) = isc(k)+1
      end do

      call mpi_alltoall(isc,1,MPI_INTEGER,irc,1,MPI_INTEGER,comm,i)

      isd(0) = 0
      ird(0) = 0
      do i=1,pp-1
	isd(i) = isd(i-1)+isc(i-1)
	ird(i) = ird(i-1)+irc(i-1)
        isc(i-1) = 0
      end do
      isc(pp-1) = 0

      do i=0,nown-1
	k = whichproc(i)
        j = isd(k)+isc(k)
	slocx(i) = isd(k)+isc(k)
	siprm(j) = iown(2*i+1)
	isc(k) = isc(k)+1
      end do

      call mpi_alltoallv(siprm,isc,isd,MPI_INTEGER,
     +                   rlocx,irc,ird,MPI_INTEGER,
     +                   comm,i)

      end
