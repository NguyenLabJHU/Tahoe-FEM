C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   subpspaces.f                                                            */
C/*                                                                           */
C/*   Written by Anshul Gupta, IBM Corp.                                      */
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
C/* $Id: subpspaces.f,v 1.2 2004-12-15 01:14:19 paklein Exp $ */
C/*****************************************************************************/

C-------------------------------------------------------------------------------
      subroutine DPACK(sbuf,dbuf,blklen,offset,blknum)

      double precision sbuf(*), dbuf(*)
      integer blklen,offset,blknum

      jstart = 1
      istart = 1
      do 10 i = 1, blknum
        call mydc(blklen,sbuf(jstart),dbuf(istart))
        jstart = jstart + offset
        istart = istart + blklen
10    continue
      return
      end
C-------------------------------------------------------------------------------
      subroutine ASSIMILATE1(wmem,ibuf,dbuf,supinds,locinds,
     1           nrowlim,csuptr,rsuptr,ldf,ibuflen,dbuflen)

      integer ibuf(*),supinds(*),locinds(*)
      integer nrowlim,csuptr,rsuptr,ldf,ibuflen,dbuflen
      double precision wmem(*), dbuf(*)

      integer ibufptr, dbufptr, locf, locptr, jc, jr

      i = ibuf(1)
      ibufptr = i + 2
      if (i .eq. 0 .or. ibuf(ibufptr) .eq. -1) then
      return
      end if
      dbufptr = 1
      locf = 1
      locptr = -1
      jc = csuptr
      jr = rsuptr
      nextcol = ibuf(ibufptr)
      
      i = 2
      do while (ibuf(i) .lt. nextcol)
        i = i + 1
      end do
      k = rsuptr
      do while (i .lt. ibufptr .and. k .lt. nrowlim)
        locptr = locptr + 2
        do while (i .lt. ibufptr .and. ibuf(i) .eq. supinds(k))
          i = i + 1
          k = k + 1
        end do
        locinds(locptr) = k - 1
        if (i .lt. ibufptr) then
          do while (ibuf(i) .ne. supinds(k))
            k = k + 1
          end do
          locinds(locptr+1) = k - 1
        else
          locinds(locptr+1) = nrowlim - 1
        end if
      end do

      do while (nextcol .ge. 0)
        do while (supinds(jc) .ne. nextcol)
          jc = jc + 1
          locf = locf + ldf
        end do
        do while (supinds(jr) .lt. nextcol)
          jr = jr + 1
          locf = locf + 1
        end do
        k = jr
        locf = locf - jr
        do 10 l = 1, locptr, 2
          dbufptr = dbufptr - k
          do 25 k = k, locinds(l)
            wmem(locf+k) = wmem(locf+k) + dbuf(dbufptr+k)
25        continue
          dbufptr = dbufptr + k
          k = max0(k,locinds(l+1)+1)
10      continue
        locf = locf + jr
        ibufptr = ibufptr + 1
        nextcol = ibuf(ibufptr)
      end do

      return
      end
C-------------------------------------------------------------------------------
      SUBROUTINE MYFC (N,DX,DY)
      integer dx(*),dy(*)
      do i=1,n
        dy(i) = dx(i)
      end do
      return
      end
C-------------------------------------------------------------------------------
      SUBROUTINE MYDC (N,DX,DY)
      integer N
      double precision dx(*),dy(*)
      do i = 1, n
      dy(i) = dx(i)
      end do
      return
      end
C-------------------------------------------------------------------------------
      SUBROUTINE MYDDC (N,DX,DY,DZ)
      integer N
      double precision DX, DY, DZ
      call mydc(n,dx,dy)
      call mydc(n,dx,dz)
      return
      end
C------------------------------------------------------------------------------ 
      subroutine a2i(i,a,LEN)
      integer i
      character a(*)

c      parameter(LEN = 64)

      i = 0
      j = 1
      l = 1
      do while (j .lt. LEN .and. a(j) .ne. '0'
     1                     .and. a(j) .ne. '1'
     1                     .and. a(j) .ne. '2'
     1                     .and. a(j) .ne. '3'
     1                     .and. a(j) .ne. '4'
     1                     .and. a(j) .ne. '5'
     1                     .and. a(j) .ne. '6'
     1                     .and. a(j) .ne. '7'
     1                     .and. a(j) .ne. '8'
     1                     .and. a(j) .ne. '9')
        if (a(j) .eq. '-') l = -1
        j = j + 1
      end do
      do while (j .lt. LEN .and. (a(j) .eq. '0'
     1                     .or. a(j) .eq. '1'
     1                     .or. a(j) .eq. '2'
     1                     .or. a(j) .eq. '3'
     1                     .or. a(j) .eq. '4'
     1                     .or. a(j) .eq. '5'
     1                     .or. a(j) .eq. '6'
     1                     .or. a(j) .eq. '7'
     1                     .or. a(j) .eq. '8'
     1                     .or. a(j) .eq. '9'))
        i = i * 10
        if (a(j) .eq. '1') i = i + 1
        if (a(j) .eq. '2') i = i + 2
        if (a(j) .eq. '3') i = i + 3
        if (a(j) .eq. '4') i = i + 4
        if (a(j) .eq. '5') i = i + 5
        if (a(j) .eq. '6') i = i + 6
        if (a(j) .eq. '7') i = i + 7
        if (a(j) .eq. '8') i = i + 8
        if (a(j) .eq. '9') i = i + 9
        j = j + 1
      end do
      i = i * l
      return
      end

C------------------------------------------------------------------------------
      subroutine compmysan(sanity,n,lptrs,lvals,linds,cinfo)

c      include 'mpif.h'

      integer n, lptrs(3,0:*)
      double precision sanity, lvals(*)
      integer cinfo(0:n-1)

c      call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

c      write(10+myid,*)'cinfo:',(cinfo(i),i=0,n-1)

      sanity = 0.d0
      do 1130 is1 = 0, N - 1
          if (cinfo(is1) .eq. 1) then

c      write(10+myid,*)is1,':',(lvals(is2),
c     +			       linds(is2-lptrs(1,is1)+lptrs(3,is1)),
c     +			       is2 = lptrs(1,is1), 
c     +				     lptrs(1,is1)+lptrs(2,is1)-1)

          do 1170 is2 = lptrs(1,is1), lptrs(1,is1)+lptrs(2,is1)-1
            sanity = sanity + lvals(is2)
 1170     continue
        end if
1130  continue
      return
      end
C------------------------------------------------------------------------------
C    recursive
C    +
      subroutine symbolic6(aptrs,ainds,lptrs,
     1          linds,sup,myscrach,tptrs,tinds,
     2          nnz,root,scount,lspace,opcount,mystak,lctr,
     3          nnzxtra,opxtra,ssthresh,lsize,info) !Cmj added lsize and info

      integer scrsize,mystak(*),lspace
      integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*)
      integer sup(*),tptrs(3,0:*),tinds(*)
      integer jlctr,nnz,root,scount,lctr,nnzxtra
      double precision opcount,opxtra,xop,ssthresh
      integer myscrach(0:*)
      integer i,j,k,node,imyscr,itscr,supbot,suptop,nmiss
      integer lsize,info !Cmj

      if(info.eq.1) return !Cmj
      istk = 1
      node = root
      k = tptrs(2,node)
      do while (k .eq. 1)
        mystak(istk) = node
        node = tinds(tptrs(1,node))
        k = tptrs(2,node)
        istk = istk + 1
      end do
      mystak(istk) = node

      if (k .eq. 0) then
      jj = aptrs(1,node)
        do 10 itscr = 0, aptrs(2,node) - 1
          linds(lctr+itscr) = ainds(jj+itscr)
10      continue
      else
      do i = 0, k - 1
          call symbolic6_recursive(aptrs,ainds,lptrs,linds,
     1       sup,myscrach,tptrs,tinds,nnz,
     2       tinds(tptrs(1,node)+i),scount,lspace,opcount,
     3       mystak(1+istk),lctr,nnzxtra,opxtra,ssthresh,lsize,info) !Cmj
          if(info.eq.1) return !Cmj
      end do
      i = tinds(tptrs(1,node))
      j = tinds(tptrs(1,node)+1)
      call smerge(lptrs(2,i)-1,linds(lptrs(3,i)+1),
     1              lptrs(2,j)-1,linds(lptrs(3,j)+1),
     2              itscr,myscrach)
        do i = 2, k - 2, 2
        j = tinds(tptrs(1,node)+i)
        ii = itscr
        call smerge(ii,myscrach,lptrs(2,j)-1,linds(lptrs(3,j)+1),
     1                itscr,linds(lctr))
        j = tinds(tptrs(1,node)+i+1)
        ii = itscr
        call smerge(ii,linds(lctr),lptrs(2,j)-1,
     1                linds(lptrs(3,j)+1),itscr,myscrach)
      end do
        if (i .eq. k - 1) then
        j = tinds(tptrs(1,node)+i)
        ii = itscr
          call smerge(ii,myscrach,lptrs(2,j)-1,linds(lptrs(3,j)+1),
     1                itscr,linds(lctr))
        do j = 0, itscr - 1
          myscrach(j) = linds(lctr+j)
        end do
        end if

        imyscr = itscr
        itscr = 0
        ii = aptrs(1,node)
        kk = aptrs(2,node) + ii
        jj = 0
        do while (ii .lt. kk .and. jj .lt. imyscr)
          if (ainds(ii) .gt. myscrach(jj)) then
            linds(lctr+itscr) = myscrach(jj)
            jj = jj + 1
          else
            if (ainds(ii) .eq. myscrach(jj)) then
              linds(lctr+itscr) = myscrach(jj)
              jj = jj + 1
              ii = ii + 1
            else
              linds(lctr+itscr) = ainds(ii)
              ii = ii + 1
            end if
          end if
          itscr = itscr + 1
        end do
        do 444 ii = ii, kk - 1
          linds(lctr+itscr) = ainds(ii)
          itscr = itscr + 1
444     continue
        do 555 jj = jj, imyscr - 1
          linds(lctr+itscr) = myscrach(jj)
          itscr = itscr + 1
555     continue
      end if

      lptrs(1,node) = lspace
      lptrs(2,node) = itscr
      lptrs(3,node) = lctr
      jlctr = lctr + 1
      lctr = lctr + itscr
      if(lctr.gt.lsize) then
      info = 1
      return
      end if
      nnz = nnz + itscr
      opcount = opcount + dble(itscr * itscr)
      supbot = node
      isupbot = istk
      isupsize = 0
      lds = itscr
      lspace = lspace + lds

      do 30 j = istk-1, 1, -1
        suptop = node
        isupsize = isupsize + 1
        itscr = 0
        nmiss = 0
      jml = -1
        node = mystak(j)
        ii = jlctr
        jj = aptrs(1,node)
        k = aptrs(2,node) + jj
        do while (jj .lt. k .and. ii .lt. lctr)
          if (ainds(jj) .gt. linds(ii)) then
            linds(lctr+itscr) = linds(ii)
            ii = ii + 1
          else
            if (ainds(jj) .eq. linds(ii)) then
              linds(lctr+itscr) = linds(ii)
              ii = ii + 1
              jj = jj + 1
            else
            if (nmiss .eq. 0) jml = itscr
              linds(lctr+itscr) = ainds(jj)
              jj = jj + 1
              nmiss = nmiss + 1
            end if
          end if
          itscr = itscr + 1
        end do
      if (jml .eq. -1 .and. k .gt. jj) jml = itscr
        nmiss = nmiss + k - jj
        do 70 jj = jj, k - 1
            linds(lctr+itscr) = ainds(jj)
            itscr = itscr + 1
70      continue
        do 80 ii = ii, lctr - 1
          linds(lctr+itscr) = linds(ii)
          itscr = itscr + 1
80      continue
      if (jml .eq. -1) jml = itscr
        lptrs(2,node) = itscr
        nnz = nnz + itscr
        opcount = opcount + dble(itscr * itscr)

      kk = nmiss * isupsize
      xop = dble(kk * kk)
      if (nmiss .gt. 0 .and. xop*(opxtra/opcount) .lt. ssthresh) then
        lspace = lspace + kk
        jj = 0
        do i = isupbot, j+1, -1
          nnd = mystak(i)
          lptrs(1,nnd) = lptrs(1,nnd) + jj
          lptrs(2,nnd) = lptrs(2,nnd) + nmiss
          jj = jj + nmiss
        end do
        do i = jml, itscr - 1
          linds(jlctr+i) = linds(lctr+i)
        end do
        opxtra = opxtra + xop
        lds = lds + nmiss
        lctr = jlctr + itscr
        if(lctr.gt.lsize) then
          info = 1
          return
        end if
        nnzxtra = nnzxtra + kk
        nmiss = 0
      end if

        if (nmiss .ne. 0) then
          lptrs(1,node) = lspace
          lptrs(3,node) = lctr
          lds = itscr
          jlctr = lctr + 1
          lctr = lctr + itscr
        if(lctr.gt.lsize) then
          info = 1
          return
        end if
          tptrs(3,supbot) = scount
          tptrs(3,suptop) = scount
          sup(scount) = supbot
          sup(scount+1) = isupsize
          sup(scount+2) = node
          scount = scount + 4
          supbot = node
        isupbot = j
          isupsize = 0
        else
          lptrs(1,node) = lspace + isupsize
          lptrs(3,node) = jlctr
          jlctr = jlctr + 1
        end if
        lspace = lspace + lds
30    continue

      isupsize = isupsize + 1
      tptrs(3,supbot) = scount
      tptrs(3,node) = scount
      sup(scount) = supbot
      sup(scount+1) = isupsize
      sup(scount+2) = -1
      scount = scount + 4
      return
      end 
C-------------------------------------------------------------------------------
      subroutine smerge(i1,l1,i2,l2,i3,l3)

      integer i1,i2,i3,l1(*),l2(*),l3(0:*)

      ii = 1
      jj = 1
      i3 = 0
      do while (ii .le. i1 .and. jj .le. i2)
      if (l1(ii) .eq. l2(jj)) then
        l3(i3) = l1(ii)
        ii = ii + 1
        jj = jj + 1
        i3 = i3 + 1
      else if (l1(ii) .gt. l2(jj)) then
        l3(i3) = l2(jj)
        jj = jj + 1
        i3 = i3 + 1
      else
        l3(i3) = l1(ii)
        ii = ii + 1
        i3 = i3 + 1
      end if
      end do
      do i = ii, i1
      l3(i3) = l1(i)
      i3 = i3 + 1
      end do
      do i = jj, i2
      l3(i3) = l2(i)
      i3 = i3 + 1
      end do

      return
      end
C-------------------------------------------------------------------------------
      subroutine FSOLVE fsolve(N,lvals,linds,lptrs,tinds,tptrs,sup,rhs,
     1           nrhs,root,lc,iptrs,w)

      integer N,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*)
      integer sup(*),nrhs,root,lc(*),iptrs(2,0:*)
      double precision lvals(*),rhs(0:N-1,*),w(*)

      call rfsolve(N,root,lvals,linds,lptrs,tinds,tptrs,sup,rhs,
     1     nrhs,iptrs,lc,i,ldr,w)

      return
      end
C------------------------------------------------------------------------------
      subroutine BSOLVE(N,lvals,linds,lptrs,tinds,tptrs,sup,rhs,
     1           nrhs,root,lc,iptrs,w)

      integer N,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*)
      integer sup(*),nrhs,root,lc(*),iptrs(2,0:*)
      double precision lvals(*),rhs(0:N-1,*),w(*)

      call rbsolve(N,root,lvals,linds,lptrs,tinds,tptrs,sup,rhs,
     1     nrhs,iptrs,lc,w,0,w)
      return
      end
C------------------------------------------------------------------------------
C     recursive
C    +
      subroutine RBSOLVE(N,root,lvals,linds,lptrs,tinds,
     1          tptrs,sup,rhs,nrhs,iptrs,lc,w,ldw,wp)

      integer N,root,linds(*),lptrs(3,0:*),tinds(*),tptrs(3,0:*)
      integer sup(*),iptrs(2,0:*),lc(*),ldr,ldw
      double precision lvals(*),w(0:*),rhs(0:N-1,*),wp(0:*)

      integer rank

      i = tptrs(3,root)
      jbot = sup(i)
      rank = sup(i+1)
      ldr = lptrs(2,jbot)
      ldj = ldr - rank
      call compress(root,wp(rank),ldr,w,ldw,iptrs,lc,nrhs)
      call bgetrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,wp,rhs)
      if (nrhs .gt. 1) then
        call dgemm('T','N',rank,nrhs,ldj,-1.d0,lvals(lptrs(1,jbot)+rank)
     1     ,ldr,wp(rank),ldr,1.d0,wp,ldr)
        call dtrsm('L','L','T','N',rank,nrhs,1.d0,lvals(lptrs(1,jbot)),
     1     ldr,wp,ldr)
      else
        call dgemv('T',ldr-rank,rank,-1.d0,lvals(lptrs(1,jbot)+rank),
     1     ldr,wp(rank),1,1.d0,wp,1)
        call dtrsv('L','T','N',rank,lvals(lptrs(1,jbot)),ldr,wp,1)
      end if
      call putrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,wp,rhs)
      jloc = ldr * nrhs
      do i = tptrs(1, jbot), tptrs(1, jbot) + tptrs(2,jbot) - 1
        kid = tinds(i)
        call rbsolve_recursive(N,kid,lvals,linds,lptrs,tinds,tptrs,sup,rhs,nrhs,
     1       iptrs,lc,wp,ldr,wp(jloc))
      end do
      return
      end
C------------------------------------------------------------------------------
      subroutine BGETRHS(N,linds,rank,ldr,nrhs,w,rhs)
      integer N,linds(*),rank,ldr,nrhs
      double precision w(ldr,*),rhs(0:N-1,*)

      do 10 i = 1, rank
        j = linds(i)
        do 20 k = 1, nrhs
          w(i,k) = rhs(j,k)
20      continue
10    continue
      return
      end
C------------------------------------------------------------------------------
      subroutine BGETRHS1(N,linds,rank,w,rhs)
      integer N,linds(*),rank
      double precision w(*),rhs(0:*)

      do 10 i = 1, rank
        w(i) = rhs(linds(i))
10    continue
      return
      end
C------------------------------------------------------------------------------
      subroutine COMPRESS(node,dest,ldd,source,lds,iptrs,lc,nrhs)
      integer node,ldd,lds,nrhs
      integer iptrs(2,0:*),lc(*)
      double precision dest(ldd,*),source(lds,*)

      m = 1
      j = 1
      do i = iptrs(1,node), iptrs(1,node) + iptrs(2,node) - 1, 2
        do k = j, lc(i) + j - 1
          do nr = 1, nrhs
            dest(m,nr) = source(k,nr)
          end do
          m = m + 1
        end do
        j = k + lc(i+1)
      end do
      return
      end
C------------------------------------------------------------------------------
C     recursive
C    +
      subroutine RFSOLVE(N,root,lvals,linds,lptrs,tinds,
     1          tptrs,sup,rhs,nrhs,iptrs,lc,rank,ldr,w)

      integer N,root,linds(*),lptrs(3,0:*),tinds(N),tptrs(3,0:*)
      integer sup(*),iptrs(2,0:*),lc(*),rank,ldr
      double precision lvals(*),w(0:*),rhs(0:N-1,*)

      i = tptrs(3,root)
      jbot = sup(i)
      rank = sup(i+1)
      ldr = lptrs(2,jbot)
      jloc = ldr * nrhs
      call getrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,w,rhs)
      do i = tptrs(1,jbot), tptrs(1,jbot) + tptrs(2,jbot) - 1
        kid = tinds(i)
        call rfsolve_recursive(N,kid,lvals,linds,lptrs,tinds,tptrs,sup,
     1       rhs,nrhs,iptrs,lc,jrank,ldj,w(jloc))
        call mergerhs(ldr,w,ldj,w(jloc+jrank),nrhs,kid,iptrs,lc)
      end do
      if (nrhs .eq. 1) then
      call dtrsv('L','N','N',rank,lvals(lptrs(1,jbot)),ldr,w,1)
      call dgemv('N',ldr-rank,rank,-1.d0,
     1       lvals(lptrs(1,jbot)+rank),ldr,w,1,1.d0,w(rank),1)
      else
        call dtrsm('L','L','N','N',rank,nrhs,1.d0,lvals(lptrs(1,jbot)), 
     1     ldr,w,ldr)
        call dgemm('N','N',ldr-rank,nrhs,rank,-1.d0,
     1     lvals(lptrs(1,jbot)+rank),ldr,w,ldr,1.d0,w(rank),ldr)
      end if
      call putrhs(N,linds(lptrs(3,jbot)),rank,ldr,nrhs,w,rhs)
      return
      end
C------------------------------------------------------------------------------
      subroutine MERGERHS(ldr,w1,ldj,w2,nrhs,kid,iptrs,lc)
      integer ldr,ldj,nrhs,kid,iptrs(2,0:*),lc(*)
      double precision w1(ldr,*),w2(ldj,*)

      j = 1
      m = 1
      l = iptrs(1,kid)
      do i = l, l+iptrs(2,kid)-1, 2
        do k = j, lc(i) + j - 1
          do nr = 1, nrhs
            w1(k, nr) = w1(k,nr) + w2(m,nr)
          end do
          m = m + 1
        end do
        j = k + lc(i+1)
      end do
      return
      end
C------------------------------------------------------------------------------
      subroutine GETRHS(N,linds,rank,ldr,nrhs,w,rhs)
      integer N,linds(*),rank,ldr,nrhs
      double precision  w(ldr,*),rhs(0:N-1,*)

      do i = 1, rank
        j = linds(i)
        do k = 1, nrhs
          w(i,k) = rhs(j,k)
        end do
      end do
      do j = i, ldr
        do k = 1, nrhs
          w(j,k) = 0.d0
        end do
      end do
      return
      end
C------------------------------------------------------------------------------
      subroutine PUTRHS(N,linds,rank,ldr,nrhs,w,rhs)
      integer N,linds(*),rank,ldr,nrhs
      double precision  w(ldr,*),rhs(0:N-1,*)

      do i = 1, rank
        j = linds(i)
        do k = 1, nrhs
          rhs(j,k) = w(i,k)
        end do
      end do
      return
      end
C------------------------------------------------------------------------------
      subroutine DCHOL2(a,lda,b,ldb,ldim,lrank,lup,
     1           linds,dfopts,ifopts,info)
      implicit double precision (a-h,o-z)
      double precision a(0:lda-1,0:*),dfopts(7)
      double precision b(0:ldb-1,0:*)
      integer lda,ldim,lrank,lup,linds(0:*),ifopts(5),ldb
      parameter (one = 1.d0,zero = 0.d0)

      dsmall = dfopts(1)
      ddmax = dfopts(5)
      ddmin = dfopts(6)

      lr1 = mod(lrank,2)
      if (lr1 .eq. 1) then
      a11 = a(0,0)
      if (a11 .le. 0.d0) go to 99
      a11 = dsqrt(a11)
2111    if (a11 .gt. ddmax) ddmax = a11
      if (a11 .lt. ddmin) ddmin = a11
3111    rd1 = one/a11
        a(0,0) = a11
      do i = 1, ldim - 4, 4
        a(i,0) = a(i,0) * rd1
        a(i+1,0) = a(i+1,0) * rd1
        a(i+2,0) = a(i+2,0) * rd1
        a(i+3,0) = a(i+3,0) * rd1
      end do
      do j = i, ldim - 1
        a(j,0) = a(j,0) * rd1
      end do
      end if

      do i = lr1, lrank - 1, 2

      a11 = a(i,i)
      a21 = a(i+1,i)
      a22 = a(i+1,i+1)

      do j = 0, i - 2, 2
        t11 = a(i,j)
        t21 = a(i+1,j)
        t12 = a(i,j+1)
        t22 = a(i+1,j+1)

        a11 = a11 - t11 * t11
        a21 = a21 - t11 * t21
        a22 = a22 - t21 * t21
        a11 = a11 - t12 * t12
        a21 = a21 - t12 * t22
        a22 = a22 - t22 * t22
      end do

      if (j .eq. i-1) then
        t11 = a(i,j)
        t21 = a(i+1,j)
        a11 = a11 - t11 * t11
        a21 = a21 - t11 * t21
        a22 = a22 - t21 * t21
      end if

        if (a11 .le. 0.d0) go to 99
      a11 = dsqrt(a11)
2112    if (a11 .gt. ddmax) ddmax = a11
        if (a11 .lt. ddmin) ddmin = a11
3112    a(i,i) = a11
      rd1 = one/a11
      r21 = a21 * rd1
      a(i+1,i) = r21
      a22 = a22 - r21 * r21
      if (a22 .le. 0.d0) go to 99
      a22 = dsqrt(a22)
2113    if (a22 .gt. ddmax) ddmax = a22
        if (a22 .lt. ddmin) ddmin = a22
3113    a(i+1,i+1) = a22
      rd2 = one/a22

      do k = i+2, ldim - 4, 4

        a11 = a(k,i)
        a21 = a(k+1,i)
        a31 = a(k+2,i)
        a41 = a(k+3,i)
        a12 = a(k,i+1)
        a22 = a(k+1,i+1)
        a32 = a(k+2,i+1)
        a42 = a(k+3,i+1)
           
        do j = 0, i - 2, 2
          t11 = a(i,j)
          t21 = a(i+1,j)
          t12 = a(i,j+1)
          t22 = a(i+1,j+1)

          x11 = a(k,j)
          x21 = a(k+1,j)
          x12 = a(k,j+1)
          x22 = a(k+1,j+1)

          a11 = a11 - x11 * t11
          a21 = a21 - x21 * t11
          a12 = a12 - x11 * t21
          a22 = a22 - x21 * t21
          a11 = a11 - x12 * t12
          a21 = a21 - x22 * t12
          a12 = a12 - x12 * t22
          a22 = a22 - x22 * t22
 
          x31 = a(k+2,j)
          x41 = a(k+3,j)
          x32 = a(k+2,j+1)
          x42 = a(k+3,j+1)

          a31 = a31 - x31 * t11
          a41 = a41 - x41 * t11
          a32 = a32 - x31 * t21
          a42 = a42 - x41 * t21
          a31 = a31 - x32 * t12
          a41 = a41 - x42 * t12
          a32 = a32 - x32 * t22
          a42 = a42 - x42 * t22
        end do
        if (j .eq. i-1) then
          t11 = a(i,j)
          t21 = a(i+1,j)
          x11 = a(k,j)
          x21 = a(k+1,j) 
          x31 = a(k+2,j) 
          x41 = a(k+3,j) 
          a11 = a11 - x11 * t11
          a21 = a21 - x21 * t11
          a31 = a31 - x31 * t11
          a41 = a41 - x41 * t11
          a12 = a12 - x11 * t21
          a22 = a22 - x21 * t21
          a32 = a32 - x31 * t21
          a42 = a42 - x41 * t21
        end if

        a11 = a11 * rd1
        a21 = a21 * rd1
        a31 = a31 * rd1
        a41 = a41 * rd1
        a(k,i) = a11 
        a(k+1,i) = a21 
        a(k+2,i) = a31
        a(k+3,i) = a41

        a12 = a12 - a11 * r21
        a22 = a22 - a21 * r21
        a32 = a32 - a31 * r21
        a42 = a42 - a41 * r21
        a(k,i+1) = a12 * rd2
        a(k+1,i+1) = a22 * rd2
        a(k+2,i+1) = a32 * rd2
        a(k+3,i+1) = a42 * rd2
      end do

      do k1 = k, ldim - 1
        a11 = a(k1,i)
        a12 = a(k1,i+1)
        do j = 0, i - 1
          t11 = a(i,j)
          t21 = a(i+1,j)
          x11 = a(k1,j)
          a11 = a11 - t11 * x11
          a12 = a12 - t21 * x11
        end do
        a11 = a11 * rd1
        a(k1,i) = a11
        a12 = a12 - a11 * r21
        a(k1,i+1) = a12 * rd2
      end do      
      end do

      do i = lrank, lup - 2, 2

      a11 = zero
      a21 = zero
      a22 = zero

      do j = 0, lrank - 2, 2
        t11 = a(i,j)
        t21 = a(i+1,j)
        t12 = a(i,j+1)
        t22 = a(i+1,j+1)

        a11 = a11 - t11 * t11
        a21 = a21 - t11 * t21
        a22 = a22 - t21 * t21
        a11 = a11 - t12 * t12
        a21 = a21 - t12 * t22
        a22 = a22 - t22 * t22
      end do

      if (lr1 .eq. 1) then
        t11 = a(i,j)
        t21 = a(i+1,j)
        a11 = a11 - t11 * t11
        a21 = a21 - t11 * t21
        a22 = a22 - t21 * t21
      end if

      b(i,i) = a11
      b(i+1,i) = a21
      b(i+1,i+1) = a22

      do k = i+2, ldim - 4, 4

        a11 = zero
        a21 = zero
        a31 = zero
        a41 = zero
        a12 = zero
        a22 = zero
        a32 = zero
        a42 = zero
           
        do j = 0, lrank - 2, 2
          t11 = a(i,j)
          t21 = a(i+1,j)
          t12 = a(i,j+1)
          t22 = a(i+1,j+1)

          x11 = a(k,j)
          x21 = a(k+1,j)
          x12 = a(k,j+1)
          x22 = a(k+1,j+1)

          a11 = a11 - x11 * t11
          a21 = a21 - x21 * t11
          a12 = a12 - x11 * t21
          a22 = a22 - x21 * t21
          a11 = a11 - x12 * t12
          a21 = a21 - x22 * t12
          a12 = a12 - x12 * t22
          a22 = a22 - x22 * t22
 
          x31 = a(k+2,j)
          x41 = a(k+3,j)
          x32 = a(k+2,j+1)
          x42 = a(k+3,j+1)

          a31 = a31 - x31 * t11
          a41 = a41 - x41 * t11
          a32 = a32 - x31 * t21
          a42 = a42 - x41 * t21
          a31 = a31 - x32 * t12
          a41 = a41 - x42 * t12
          a32 = a32 - x32 * t22
          a42 = a42 - x42 * t22
        end do
        if (lr1 .eq. 1) then
          t11 = a(i,j)
          t21 = a(i+1,j)
          x11 = a(k,j)
          x21 = a(k+1,j) 
          x31 = a(k+2,j) 
          x41 = a(k+3,j) 
          a11 = a11 - x11 * t11
          a21 = a21 - x21 * t11
          a31 = a31 - x31 * t11
          a41 = a41 - x41 * t11
          a12 = a12 - x11 * t21
          a22 = a22 - x21 * t21
          a32 = a32 - x31 * t21
          a42 = a42 - x41 * t21
        end if

        b(k,i) = a11 
        b(k+1,i) = a21 
        b(k+2,i) = a31
        b(k+3,i) = a41

        b(k,i+1) = a12 
        b(k+1,i+1) = a22 
        b(k+2,i+1) = a32 
        b(k+3,i+1) = a42 
      end do

      do k1 = k, ldim - 1
        a11 = zero
        a12 = zero
        do j = 0, lrank - 1
          t11 = a(i,j)
          t21 = a(i+1,j)
          x11 = a(k1,j)
          a11 = a11 - t11 * x11
          a12 = a12 - t21 * x11
        end do
        b(k1,i) = a11
        b(k1,i+1) = a12
      end do      
      end do

      if (i .eq. lup-1) then
      do k = i, ldim - 1
        a11 = zero
        do j = 0, lrank - 1
          a11 = a11 - a(k,j) * a(i,j)
        end do
        b(k,i) = a11
      end do
      end if

      dfopts(5) = ddmax      
      dfopts(6) = ddmin      
      info=0 !Cmj
      return

99    dfopts(5) = ddmax
      dfopts(6) = ddmin
      info=1 !Cmj replaced "return 1"
      end
C-----------------------------------------------------------------------------
      subroutine PUTRHS1(N,linds,rank,w,rhs)
      integer N,linds(*),rank,ldr
      double precision  w(*),rhs(0:*)

      do i = 1, rank
        rhs(linds(i)) = w(i)
      end do
      return
      end
