C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   pparfact1i.f                                                            */
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
C/* $Id: pparfact1i.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine PPARFACT1(N,aptrs,ainds,avals,lptrs,linds,lvals,
     1           tptrs,tinds,sup,stak,nstak,root,dd,lgblk,blk,myid,
     2           cinfo,supinds,supindsize,dfopts,ifopts,dimstak,
     3           wsize0,wsize1,ibuflen,dbuflen,iwspace,node,stakptr,
     4           nptr,lc,iptrs,info,comm) 

      integer N,root,dd,lgblk,blk,myid,nstak(*),supinds(*)
      integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*)
      integer tptrs(3,0:*),tinds(*),sup(*),stak(3,*)
      integer supindsize,info,comm
      integer ifopts(*),dimstak(*)
      integer cinfo(0:*)
      double precision dfopts(*)
      double precision lvals(*),avals(*)
      integer AE_TYPE_I,AE_TYPE_D

      parameter (AE_TYPE_I=1,AE_TYPE_D=2)

      double precision, allocatable:: dbuf_s(:),wmem(:)      
      double precision wmem0(*),wmem1(*),dbuf_r(*)
      integer, allocatable:: ibuf_s(:)
      integer, allocatable:: locinds(:)
      integer ibuf_r(*)
      pointer (pdbufr,dbuf_r),(pibufr,ibuf_r),(pw0,wmem0),(pw1,wmem1)

      integer rsuptr, csuptr, nptr, stakptr, rank, ibufptr, dbufptr
      integer hdim,vdim,dbuflen,halfbuflen,udim,ldu,wsize1,wsize0
      integer fptr, uptr, rsuptr_u, csuptr_u, partner, ibuflen
      integer vcube(0:31),hcube(0:31)
      integer hcsize, vcsize, currdim, myleft, myright, myup, mydown
      integer diagproc, vcloc, hcloc
      integer lc(*),iptrs(*) !Cmj

      include 'mpif.h'
      integer mpistat(MPI_STATUS_SIZE),ierr,ip,level

      info = 0 !Cmj

      diagproc = 1
      vcsize = 1
      hcsize = 1
      vcube(0) = myid
      hcube(0) = myid
      myleft = myid
      myright = myid
      myup = myid
      mydown = myid

      halfbuflen = ishft(ishft(dbuflen-2,-2),1) + 1

      allocate(dbuf_s(dbuflen*2),stat=i)
      allocate(ibuf_s(ibuflen*2),stat=j)
      if (i+j .ne. 0) then
        print *,myid,': Unable to allocate working storage (dbuf+ibuf)'
        call mpi_abort(comm,0,ierr)
      end if
      pdbufr = loc(dbuf_s(dbuflen+1))
      pibufr = loc(ibuf_s(ibuflen+1))

      allocate(wmem(iwspace),stat=i) 
      if (i .ne. 0) then
        print *,myid,': Unable to allocate working storage (wmem)'
        call mpi_abort(comm,0,ierr)
      end if
      pw0 = loc(wmem)

      ldu = 0

      call factor6(wmem0,linds,lptrs,ainds,aptrs,avals, 
     1     lvals,tinds,tptrs,sup,rank,ldu,udim,kk,node,2*(iwspace+1),
     2     nstak(nptr),lc,iptrs,dfopts,ifopts,info)

      if (info.gt.0) return ! Cmj

      k = sup(tptrs(3,node))
      ncols_u = udim 
      nrows_u = udim
      csuptr_u = supindsize
      rsuptr_u = supindsize + udim
      call myfc(udim,linds(lptrs(3,k)),supinds(supindsize))
      call myfc(udim,supinds(supindsize),supinds(rsuptr_u))

      uptr = kk + rank * (ldu + 1)
      if (wsize1 .lt. uptr) then
      pw1 = loc(wmem)
      pw0 = loc(wmem(wsize1+1))
        uptr = uptr - wsize1
      else 
      j = kk + ldu*udim
      if (iwspace-j .ge. wsize1-1) then
        pw1 = loc(wmem(iwspace-wsize1+1))
      else
        pw1 = loc(wmem(wsize0+1))
        fptr = 1
        j = udim - rank 
        k = j - 1
        do i = 0, k
          do l = i, k
            wmem0(fptr+l) = wmem0(uptr+l)
          end do
          fptr = fptr + j
          uptr = uptr + ldu
        end do
        uptr = 1
        ldu = j
      end if
      end if
      nptr = nptr - 1
      hdim = 0
      vdim = 0
      currdim = 0
      mask1 = blk
      ncols_u = ncols_u - rank
      csuptr_u = csuptr_u + rank
      if (ncols_u .gt. 0) then
        k = supinds(csuptr_u)
        do 110 i = rsuptr_u, rsuptr_u + nrows_u - 1
          if (supinds(i) .ge. k) goto 120
110     continue
120     nrows_u = nrows_u - i + rsuptr_u
        rsuptr_u = i
      else
        nrows_u = 0
      end if

      allocate(locinds(0:N-1),stat=i)
      if(i.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      do while (currdim .ne. dd)

        node = nstak(nptr)
        partner = ieor(myid,ishft(1,currdim))
        stakptr = stakptr - 1
        ncols = stak(2,stakptr)
        nrows = stak(3,stakptr)
        i = tptrs(3,node) + 2
        csuptr = sup(i)
        rsuptr = csuptr + sup(i+1)

        currdim = currdim + 1
        ldf = dimstak(currdim)

        if (hdim .eq. vdim) then

          j = ishft(1,currdim-1)
          if (myid .lt. partner) then
            do 310 i = 0, hcsize - 1
              hcube(hcsize + i) = ieor(hcube(i), j)
310         continue
          else
            do 320 i = 0, hcsize - 1
              hcube(hcsize + i) = hcube(i)
              hcube(i) = ieor(hcube(i), j)
320         continue
          end if
          hcsize = hcsize + hcsize
          hcloc = 0
          do while (hcube(hcloc) .ne. myid)
            hcloc = hcloc + 1
          end do
          if (hcloc .eq. 0) then
            myleft = hcube(hcsize - 1)
          else
            myleft = hcube(hcloc - 1)
          end if
          if (hcloc .eq. hcsize-1) then
            myright = hcube(0)
          else 
            myright = hcube(hcloc + 1)
          end if

          ibuf_s(1) = nrows_u
          call myfc(nrows_u,supinds(rsuptr_u),ibuf_s(2))
          ibufptr = nrows_u + 2

          hdim = hdim + 1
          fptr = min0(wsize1,wsize1-ldf*(ncols-1)-nrows+1)
          i = rsuptr_u
          j = rsuptr_u + nrows_u
          k = rsuptr
          l = rsuptr + nrows

          locptr = -1
          do while (i .lt. j .and. supinds(i) .lt. supinds(k))
            i = i + 1
          end do
          do while (i .lt. j .and. k .lt. l)
            locptr = locptr + 2
            do while (i .lt. j .and. supinds(i) .eq. supinds(k))
              i = i + 1
              k = k + 1
            end do
            locinds(locptr) = k - 1
            if (i .lt. j) then
              do while (supinds(i) .ne. supinds(k))
                k = k + 1
              end do
              locinds(locptr+1) = k - 1
            else
              locinds(locptr+1) = l - 1
            end if
          end do
          if (locptr .eq. -1) then
            locinds(1) = rsuptr - 1
            locinds(2) = l - 1
            locptr = 1
          end if

          dbufptr = 1
          kk = csuptr
          ll = csuptr + ncols
          ncols_u = ncols_u + csuptr_u
          nrows_u = nrows_u + rsuptr_u
          locf = fptr
          ii = rsuptr
          do while (kk .lt. ll)

            do while (csuptr_u .lt. ncols_u .and. 
     1                supinds(csuptr_u) .lt. supinds(kk))
              do while (supinds(rsuptr_u) .lt. supinds(csuptr_u))
                rsuptr_u = rsuptr_u + 1
                uptr = uptr + 1
              end do
              m = nrows_u - rsuptr_u
              ibuf_s(ibufptr) = supinds(csuptr_u)
              ibufptr = ibufptr + 1
              call mydc(m,wmem0(uptr),dbuf_s(dbufptr))
              uptr = uptr + ldu
              dbufptr = dbufptr + m 
              csuptr_u = csuptr_u + 1
            end do

            if (csuptr_u .lt. ncols_u) then

              if (supinds(csuptr_u) .eq. supinds(kk)) then
                do while (supinds(rsuptr_u) .lt. supinds(csuptr_u))
                  rsuptr_u = rsuptr_u + 1
                  uptr = uptr + 1
                end do
                do while (supinds(ii) .lt. supinds(kk))
                  ii = ii + 1
                  locf = locf + 1
                end do
                i = ii
                locf = locf - ii
                j = uptr
                do 130 l = 1, locptr, 2
                  j = j - i
                  do 140 k = i, locinds(l)
                    wmem1(locf+k) = wmem0(j+k)
140               continue
                  j = j + k
                  do 150 i = k, locinds(l+1)
                    wmem1(locf+i) = 0.d0
150               continue
130             continue
                locf = locf + ii + ldf
                uptr = uptr + ldu
                csuptr_u = csuptr_u + 1
              else

                do while (supinds(ii) .lt. supinds(kk))
                  ii = ii + 1
                  locf = locf + 1
                end do
                do 160 i = locf, locinds(locptr+1)-ii+locf
                  wmem1(i) = 0.d0
160             continue
                locf = locf + ldf
              end if
              kk = kk + 1

            else
              do while (supinds(ii) .lt. supinds(kk))
                ii = ii + 1
                locf = locf + 1
              end do
              do 170 i = locf, locinds(locptr+1)-ii+locf
                wmem1(i) = 0.d0
170           continue
              locf = locf + ldf
              kk = kk + 1
            end if
          end do

          do while (csuptr_u .lt. ncols_u)
            do while (supinds(rsuptr_u) .lt. supinds(csuptr_u))
              rsuptr_u = rsuptr_u + 1
              uptr = uptr + 1
            end do
            m = nrows_u - rsuptr_u
            ibuf_s(ibufptr) = supinds(csuptr_u)
            ibufptr = ibufptr + 1
            call mydc(m,wmem0(uptr),dbuf_s(dbufptr))
            uptr = uptr + ldu
            dbufptr = dbufptr + m
            csuptr_u = csuptr_u + 1
          end do
          ibuf_s(ibufptr) = -1

          if (partner .gt. myid) then
            msgsize = ibufptr
            call mpi_send(ibuf_s,msgsize,MPI_INTEGER,partner,AE_TYPE_I,
     +        comm,ierr)
            msgsize = ishft(dbufptr-1,3)
            call mpi_send(dbuf_s,msgsize,MPI_BYTE,partner,AE_TYPE_D,
     +                      comm,ierr)
            msgsize = ibuflen
            call mpi_recv(ibuf_r,msgsize,MPI_INTEGER,partner,
     1                    AE_TYPE_I,comm,mpistat,ierr)
            call mpi_recv(dbuf_r,ishft(dbuflen,3),MPI_BYTE,partner,
     1                    AE_TYPE_D,comm,mpistat,ierr)
          else
            msgsize = ibuflen
            call mpi_recv(ibuf_r,msgsize,MPI_INTEGER,partner,
     1                    AE_TYPE_I,comm,mpistat,ierr)
            call mpi_recv(dbuf_r,ishft(dbuflen,3),MPI_BYTE,partner,
     1                    AE_TYPE_D,comm,mpistat,ierr)
            msgsize = ibufptr
            call mpi_send(ibuf_s,msgsize,MPI_INTEGER,partner,AE_TYPE_I,
     +                    comm,ierr)
            msgsize = ishft(dbufptr-1,3)
            call mpi_send(dbuf_s,msgsize,MPI_BYTE,partner,AE_TYPE_D,
     +                    comm,ierr)
          end if

          call assimilate1(wmem1(fptr),ibuf_r,dbuf_r,supinds,locinds,
     1         rsuptr+nrows,csuptr,rsuptr,ldf,ibuflen,dbuflen)

          call parelimv1(wmem1,dbuf_s,dbuf_s(halfbuflen),dbuf_r,
     1         dbuf_r(halfbuflen),supinds,tptrs,tinds,cinfo,
     2         stak,nstak,avals,aptrs,ainds,lvals,lptrs,linds,sup,
     3         stakptr,nptr,mask1,fptr,ldf,nrows,ncols,
     4         rsuptr,csuptr,myid, myleft,myright,myup,mydown,
     5         diagproc,halfbuflen-1,wsize1,locinds,N,
     6         ibuf_s,dfopts,ifopts,comm)

        else

          j = ishft(1,currdim-1)
          if (myid .lt. partner) then
            do 410 i = 0, vcsize - 1
              vcube(vcsize + i) = ieor(vcube(i), j)
410         continue
          else
            do 420 i = 0, vcsize - 1
              vcube(vcsize + i) = vcube(i)
              vcube(i) = ieor(vcube(i), j)
420         continue
          end if
          vcsize = vcsize + vcsize
          vcloc = 0
          do while (vcube(vcloc) .ne. myid)
            vcloc = vcloc + 1
          end do
          if (vcloc .eq. 0) then
            myup = vcube(vcsize - 1)
          else
            myup = vcube(vcloc - 1)
          end if
          if (vcloc .eq. vcsize-1) then
            mydown = vcube(0)
          else 
            mydown = vcube(vcloc + 1)
          end if
          if (hcloc .eq. vcloc) then
            diagproc = 1
          else
            diagproc = 0
          end if

          vdim = vdim + 1
          fptr = min0(wsize0,wsize0-ldf*(ncols-1)-nrows+1)
          nrows_u = rsuptr_u + nrows_u

          k = rsuptr
          l = rsuptr + nrows - 1
          ibufptr = 2
          if (nrows .gt. 0) then
            inode = supinds(k)
          else
            inode = 1000000000
          end if
          do 210 i = rsuptr_u, nrows_u - 1
            jnode = supinds(i)
            do while (k .lt. l)
              if (inode .lt. jnode) then
                k = k + 1
                inode = supinds(k)
              else
                goto 220
              end if
            end do
220         if (inode .eq. jnode) then
              locinds(i-rsuptr_u) = 1
            else
              ibuf_s(ibufptr) = jnode
              locinds(i-rsuptr_u) = 0
              ibufptr = ibufptr + 1
            end if
210       continue 
          ibuf_s(1) = ibufptr - 2
          dbufptr = 1
          locf = fptr

          k = csuptr 
          l = csuptr + ncols 
          kk = rsuptr
          m = rsuptr + nrows
          ii = rsuptr_u 
          if (ncols .gt. 0) then
            inode = supinds(k)
          else
            inode = 1000000000
          end if
          do 230 i = csuptr_u, csuptr_u + ncols_u - 1
            jnode = supinds(i)
            do while (supinds(ii) .lt. jnode)
              ii = ii + 1
              uptr = uptr + 1
            end do
            do while (k .lt. l)
              if (inode .lt. jnode) then
                locf = locf - kk
                do 260 j = kk, m - 1
                  wmem0(locf+j) = 0.d0
260             continue
                k = k + 1
                locf = locf + ldf + kk
                inode = supinds(k)
                do while (supinds(kk) .lt. inode)
                  locf = locf + 1
                  kk = kk + 1
                end do
              else
                goto 240
              end if
            end do
            kk = m
240         jj = 0
            uptr = uptr - ii
            do 250 j = ii, nrows_u - 1
              do 275 jj = jj, m-kk-1
                if (supinds(kk+jj) .lt. supinds(j)) then
                  wmem0(locf+jj) = 0.d0
                else
                  goto 270
                end if
275           continue
270           if (locinds(j-rsuptr_u) .eq. 0) then
                dbuf_s(dbufptr) = wmem1(uptr+j)
                dbufptr = dbufptr + 1
              else 
                wmem0(locf+jj) = wmem1(uptr+j)
                jj = jj + 1
              end if
250         continue
            do 255 jj = jj, m-kk-1
              wmem0(locf+jj) = 0.d0
255         continue
            if (inode .eq. jnode) then
              k = k + 1
              if (k .lt. l) then
                locf = locf + ldf
                inode = supinds(k)
                do while (supinds(kk) .lt. inode)
                  locf = locf + 1
                  kk = kk + 1
                end do
              else
                kk = m
              end if
            end if
            uptr = uptr + ldu + ii
            ibuf_s(ibufptr) = jnode
            ibufptr = ibufptr + 1
230       continue

          do 280 k = k, l-1, 1
            do while (supinds(kk) .lt. supinds(k))
              locf = locf + 1
              kk = kk + 1
            end do
            locf = locf - kk
            do 290 j = kk, m - 1
              wmem0(locf+j) = 0.d0
290         continue
            locf = locf + ldf + kk
280       continue

          i = ibuf_s(1)
          if (i .gt. 0) then
            do while (ibuf_s(ibufptr-1) .gt. ibuf_s(i+1))
              ibufptr = ibufptr - 1
            end do
          end if
          ibuf_s(ibufptr) = -1 
          if (partner .gt. myid) then
            msgsize = ibufptr
            call mpi_send(ibuf_s,msgsize,MPI_INTEGER,partner,AE_TYPE_I,
     +                    comm,ierr)
            msgsize = ishft(dbufptr-1,3)
            call mpi_send(dbuf_s,msgsize,MPI_BYTE,partner,AE_TYPE_D,
     +                    comm,ierr)
            msgsize = ibuflen
            call mpi_recv(ibuf_r,msgsize,MPI_INTEGER,partner,
     1                    AE_TYPE_I,comm,mpistat,ierr)
            call mpi_recv(dbuf_r,ishft(dbuflen,3),MPI_BYTE,partner,
     1                    AE_TYPE_D,comm,mpistat,ierr)
          else
            msgsize = ibuflen 
            call mpi_recv(ibuf_r,msgsize,MPI_INTEGER,partner,
     1                    AE_TYPE_I,comm,mpistat,ierr)
            call mpi_recv(dbuf_r,ishft(dbuflen,3),MPI_BYTE,partner,
     1                    AE_TYPE_D,comm,mpistat,ierr)
            msgsize = ibufptr
            call mpi_send(ibuf_s,msgsize,MPI_INTEGER,partner,AE_TYPE_I,
     +                    comm,ierr)
            msgsize = ishft(dbufptr-1,3)
            call mpi_send(dbuf_s,msgsize,MPI_BYTE,partner,AE_TYPE_D,
     +                    comm,ierr)
          end if

          call assimilate1(wmem0(fptr),ibuf_r,dbuf_r,supinds,locinds,
     1         rsuptr+nrows,csuptr,rsuptr,ldf,ibuflen,dbuflen)

          call parelimh1(wmem0,dbuf_s,dbuf_s(halfbuflen),dbuf_r,
     1         dbuf_r(halfbuflen),supinds,tptrs,tinds,cinfo,
     2         stak,nstak,avals,aptrs,ainds,lvals,lptrs,linds,sup,
     3         stakptr,nptr,mask1,fptr,ldf,nrows,ncols,
     4         rsuptr,csuptr,myid,myleft,myright,myup,mydown,
     5         diagproc,halfbuflen-1,wsize0,locinds,
     6         ibuf_s,dfopts,ifopts,comm)

          mask1 = ior(mask1,ishft(mask1,1))
        end if

        uptr = fptr
        ldu = ldf
        ncols_u = ncols
        nrows_u = nrows
        csuptr_u = csuptr
        rsuptr_u = rsuptr

      end do

      deallocate(dbuf_s)
      deallocate(ibuf_s)
      deallocate(wmem)
      deallocate(locinds)
      return

111   print *,'Bad news in serial factor!'
      write(*,*) (ifopts(i),i=1,5)
      do i = 1, 7
        print *, dfopts(i)
      end do
      return
      end
