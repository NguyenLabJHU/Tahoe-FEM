C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   eparfact1i.f                                                            */
C/*                                                                           */
C/*   Written by Anshul Gupta, IBM Corp.                                      */
C/*   Modified by Mahesh Joshi, U of MN.                                      */
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
C/* $Id: eparfact1i.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      subroutine EPARFACT1(N,aptrs,ainds,lptrs,linds,tptrs,tinds,
     1           sup,stak,nstak,root,dd,lgblk,myid,supinds,dimstak,
     2           wsize0,wsize1,ibuflen,dbuflen,iwspace,node,stakptr,
     4           nptr,lc,iptrs,lcsize,wsolvesize,wa1,comm) !Cmj

      integer KONSTANT
      parameter(KONSTANT=100000)

      integer N,root,dd,lgblk,blk,myid,nstak(*),supinds(*),comm
      integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*)
      integer tptrs(3,0:*),tinds(*),sup(*),stak(3,*)
      integer dimstak(*),wa1(*)

      integer rsuptr, csuptr, nptr, stakptr, rank, ibufptr, dbufptr
      integer hdim,vdim,dbuflen,wsize1,wsize0
      integer uptr, ibuflen
      integer lc(*),iptrs(*),lcsize,wsolvesize !Cmj

      include 'mpif.h'
      integer mpistat(MPI_STATUS_SIZE),ierr,ip,level

      vdim = ishft(dd,-1)
      hdim = dd - vdim

      idv = 0
      j = ishft(myid,-1)
      do 50 i = 0, vdim - 1
        idv = ior(idv,ishft(iand(j,1),i))
        j = ishft(j,-2)
50    continue

      idh = 0
      j = myid
      do 60 i = 0, hdim - 1
        idh = ior(idh,ishft(iand(j,1),i))
        j = ishft(j,-2)
60    continue

      wsize1 = 1
      wsize0 = 1
      mask1 = ishft(1,dd-1)
      node = root
      nptr = 1
      stakptr = 1
      
      do while ((hdim+vdim) .ne. 0)
        maskr = ishft(ishft(1,vdim)-1,lgblk)
        maskc = ishft(ishft(1,hdim)-1,lgblk)
        index_v = ishft(idv,lgblk)
        index_h = ishft(idh,lgblk)
        ncols_u = 0
        nrows_u = 0
        do while (tptrs(2,node) .eq. 1)
          nstak(nptr) = node
          nptr = nptr + 1
          k = tptrs(3,node)
          if (k .ne. 0) then
            if (sup(k) .eq. node) then
              csuptr = sup(k+2)
              l = sup(k+3)
              rsuptr = csuptr + l
              i = 0
              j = 0
              do 70 l = l-1, 0, -1
                if (iand(maskr,supinds(csuptr+l)).eq.index_v) goto 75
70            continue
75            do 30 m = 0, l 
                if (iand(maskc,supinds(csuptr+m)).eq.index_h) goto 35
30            continue
35            do 10 m = m, l
                inode = supinds(csuptr+m)
                if (iand(maskc,inode) .eq. index_h) then
                  supinds(csuptr+i) = inode
                  i = i + 1
                end if
                if (iand(maskr,inode) .eq. index_v) then
                  supinds(rsuptr+j) = inode
                  j = j + 1
                end if
10            continue
              stak(1,stakptr) = node
              stak(2,stakptr) = i
              stak(3,stakptr) = j
              stakptr = stakptr + 1
              if (i .gt. ncols_u) ncols_u = i
              if (j .gt. nrows_u) nrows_u = j
            end if
          end if
          node = tinds(tptrs(1,node))
        end do

        nstak(nptr) = node
        nptr = nptr + 1
        k = tptrs(3,node)
        csuptr = sup(k+2)
        l = sup(k+3)
        rsuptr = csuptr + l
        i = 0
        j = 0
        do 80 l = l-1, 0, -1
          if (iand(maskr,supinds(csuptr+l)).eq.index_v) goto 85
80      continue
85      do 40 m = 0, l 
          if (iand(maskc,supinds(csuptr+m)).eq.index_h) goto 45
40      continue
45      do 20 m = m, l
          inode = supinds(csuptr+m)
          if (iand(maskc,inode) .eq. index_h) then
            supinds(csuptr+i) = inode
            i = i + 1
          end if
          if (iand(maskr,inode) .eq. index_v) then
            supinds(rsuptr+j) = inode
            j = j + 1
          end if
20      continue
        stak(1,stakptr) = node
        stak(2,stakptr) = i
        stak(3,stakptr) = j
        stakptr = stakptr + 1
        if (i .gt. ncols_u) ncols_u = i
        if (j .gt. nrows_u) nrows_u = j
        dimstak(hdim+vdim) = nrows_u
        if (iand(myid,mask1) .eq. 0) then
          node = tinds(tptrs(1,node))
        else
          node = tinds(1+tptrs(1,node))
        end if
        mask1 = ishft(mask1,-1)
        if (hdim .ne. vdim) then
          hdim = hdim - 1
          idh = iand(idh,not(ishft(1,hdim)))
          wsize1 = max0(wsize1,ncols_u*nrows_u)
        else
          vdim = vdim - 1
          idv = iand(idv,not(ishft(1,vdim)))
          wsize0 = max0(wsize0,ncols_u*nrows_u)
        end if
      end do

      j = stak(1,stakptr-1)
      jtss = lptrs(2,j)
      do while (jtss .eq. 0 .and. j .lt. root)
      j = j + 1
        jtss = lptrs(2,j)
      end do
      itss = lptrs(2,sup(tptrs(3,node)))
      m = itss
      itss = itss * itss
      if (jtss .gt. 0) then
        jtss = jtss - ishft(jtss, -2)
        jtss = max(jtss*jtss,ishft(itss,-1))
      else
      jtss = itss
      end if
      dbuflen = jtss + 36000 * (dd + 4)
      ibuflen = max(N+ishft(N,-1)+1, ishft(m,3)+640*(dd + 1))

      is3 = 1
      iptrs(node+node+2) = 0 

      call gen_lc(node,lc,linds,lptrs,tinds,tptrs,sup,
     1            iptrs,lcsize,wa1,nstak(nptr),
     2            is3,is5,is4,iwstri)

      lcsize = is3-1
      wsolvesize = iwstri

      iwspace = is5 * is5 + is4
      wsize0 = max(wsize0,itss)
      iwspace = max(iwspace,wsize0+wsize1)

      dbuflen = max(min((dbuflen*3+1)/2,iwspace),KONSTANT)

      return
      end

