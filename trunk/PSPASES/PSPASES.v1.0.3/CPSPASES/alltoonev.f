
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   alltoonev.f                                                             */
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
C/* $Id: alltoonev.f,v 1.1 2004-12-10 20:26:44 paklein Exp $ */
C/*****************************************************************************/

      subroutine alltoonev(tgt,tmp,size,psrc,lcsize,myidv,myid,comm)
      implicit none
      include 'mpif.h'

      integer lendp,itype
      parameter(lendp=8,itype=1)

      integer  size,psrc,lcsize,myid,myidv,partner,num,comm
      double precision tgt(size),tmp(size)
      integer diffbits,i,j,k,nbrecv,msglen,mids,midr
      integer mpistat(MPI_STATUS_SIZE),ierr

      diffbits = ieor(psrc,myidv)
      num = lcsize
      if(diffbits.ne.0) then
        ! the position of mostsignificant 1 in diffbits
        ! hints at how many stages of communication a
        ! processor will be involved in....
        num = lcsize+1
        do while(diffbits.ne.0)
          diffbits = ishft(diffbits,-1)
          num = num-1
        end do
      end if
      msglen = size*lendp

      k = ishft(lcsize,1)-1
      do i =1,num
        partner = ieor(myid,ishft(1,k))

        call mpi_sendrecv(tgt,msglen,MPI_BYTE,partner,itype,
     +                    tmp,msglen,MPI_BYTE,partner,itype,
     +                    comm,mpistat,ierr)

        do j=1,size
          tgt(j) = tgt(j) + tmp(j)
        end do
        k = k-2
      end do

      end
