C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   serialfactor.f                                                          */
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
C/* $Id: serialfactor.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      recursive
     +subroutine FACTOR6(wmem,linds,lptrs,ainds,
     1          aptrs,avals,lvals,tinds,tptrs,sup,rank,ldf,fdim,
     2          myfrontptr,root,wm,mystak,lc,iptrs,
     3          dfopts,ifopts,info)

      integer lptrs(3,0:*),aptrs(2,0:*),ainds(*),tptrs(3,0:*)
      integer linds(*),tinds(*),sup(*),mystak(*),lc(*)
      integer myfrontptr,rank,ldf,fdim,root,wm,iptrs(2,0:*)
      integer ifopts(5)
      double precision wmem(*),avals(*),lvals(*),dfopts(7)

      integer frontptr,mystakptr,ldu,udim,locptr

      parameter (LIMDENSE = 50, LIMDRK = 9)

      mystakptr = 1
      node = root
      k = sup(tptrs(3,node))

      ldf = max(ldf,sup(tptrs(3,root)+3))
      do while (tptrs(2,k) .eq. 1)
        mystak(mystakptr) = node
        node = tinds(tptrs(1,k))
        mystakptr = mystakptr + 2
        k = sup(tptrs(3,node))
      end do
      jfrontptr = ldf * ldf + 1
      fdim = lptrs(2,k)
      do 11 i = mystakptr-2, 1, -2
        j = mystak(i)
        j1 = lptrs(2,sup(tptrs(3,j)))
        mystak(i+1) = jfrontptr - j1*ldf + ldf - j1
11    continue
      myfrontptr = jfrontptr - fdim * ldf + ldf - fdim

      myinds = lptrs(3,k)
      mylim = myinds + fdim

      frontptr = 1
      ldu = ldf
      do 20 i = 0, tptrs(2,k) - 1
      if (i .gt. 0) then
        frontptr = jfrontptr
        ldu = 0
        end if
        kid = tinds(tptrs(1,k)+i)
      info = 0 
        call factor6(wmem(frontptr),linds,lptrs,ainds,aptrs,
     1       avals,lvals,tinds,tptrs,sup,rank,ldu,udim,kkk,kid,
     2       wm-frontptr+1,mystak(mystakptr),lc,iptrs,
     3       dfopts,ifopts,info)

      if (info.gt.0) goto 1111 

        locptr = iptrs(2,kid) - 1

        locu = frontptr + (ldu+1) * rank + kkk - 1
        locf = myfrontptr - 1
        kincr = 1
        kidincr = 1
        myj = 1
        locindj = 0

        if (i .eq. 0) then

          do 15 j1 = iptrs(1,kid), iptrs(1,kid) + locptr - 2, 2
          locindj = locindj + lc(j1)
            do 25 myj = myj, locindj
              myi = kincr
            locindi = locindj - lc(j1)
              do 35 i1 = j1, iptrs(1,kid) + locptr, 2
                locu = locu - myi
            locindi = locindi + lc(i1)
                do 46 myi = myi, locindi 
                  wmem(locf+myi) = wmem(locu+myi)
46              continue
                locu = locu + myi
            locindi = locindi + lc(i1+1)
                do 56 myi = myi, locindi
                  wmem(locf+myi) = 0.d0
56              continue
35            continue
              locu = locu + rank + kidincr + ldu - udim
              locf = locf + ldf
              kincr = kincr + 1
              kidincr = kidincr + 1
25          continue
          locindj = locindj + lc(j1+1)
            do 65 myj = myj, locindj
              do 76 i1 = kincr, fdim
                wmem(locf+i1) = 0.d0
76            continue
              locf = locf + ldf
              kincr = kincr + 1
65          continue
15        continue

        if (locf+kincr .ne. locu) then
          locindj = locindj + lc(j1)
            do 253 myj = myj, locindj
              myi = kincr
            locindi = locindj - lc(j1)
              do 353 i1 = j1, iptrs(1,kid) + locptr, 2
                locu = locu - myi
            locindi = locindi + lc(i1)
                do 463 myi = myi, locindi 
                  wmem(locf+myi) = wmem(locu+myi)
463             continue
                locu = locu + myi
              locindi = locindi + lc(i1+1)
                do 563 myi = myi, locindi
                  wmem(locf+myi) = 0.d0
563             continue
353           continue
              locu = locu + rank + kidincr + ldu - udim
              locf = locf + ldf
              kincr = kincr + 1
              kidincr = kidincr + 1
253         continue
          locindj = locindj + lc(j1+1)
            do 653 myj = myj, locindj
              do 763 i1 = kincr, fdim
                wmem(locf+i1) = 0.d0
763           continue
              locf = locf + ldf
              kincr = kincr + 1
653         continue
        end if
        else

          do 115 j1 = iptrs(1,kid), iptrs(1,kid) + locptr, 2
          locindj = locindj + lc(j1)
            do 125 myj = myj, locindj
              myi = kincr
            locindi = locindj - lc(j1)
              do 135 i1 = j1, iptrs(1,kid) + locptr, 2
                locu = locu - myi
            locindi = locindi + lc(i1)
                do 146 myi = myi, locindi
                  wmem(locf+myi)=wmem(locf+myi)+wmem(locu+myi)
146             continue
                locu = locu + myi
            locindi = locindi + lc(i1+1)
                myi = locindi + 1
135           continue
              locu = locu + rank + kidincr + ldu - udim
              locf = locf + ldf
              kincr = kincr + 1
              kidincr = kidincr + 1
125         continue
          locindj = locindj + lc(j1+1)
            j = locindj + 1 - myj
            locf = locf + ldf * j
            kincr = kincr + j
            myj = myj + j
115       continue
        end if
20    continue
      rank = sup(tptrs(3,k)+1)
      j = node

      if (tptrs(2,k) .eq. 0) then
        locf = lptrs(1,node)
        ldl = lptrs(2,k)
        do 80 j1 = 1, rank
          myi = lptrs(3,j)
          locf = locf - myi
          do 90 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1
            do while (linds(myi) .ne. ainds(i1))
              lvals(locf + myi) = 0.d0
              myi = myi + 1
            end do
            lvals(locf + myi) = avals(i1)
            myi = myi + 1
90        continue
          do 95 myi = myi, lptrs(3,j)+lptrs(2,j)-1
            lvals(locf + myi) = 0.d0
95        continue
          locf = locf + lptrs(3,j) - ldl - 1
          j = tinds(tptrs(1,j))
80      continue
        locf = lptrs(1,k)
      info = 0 
        call dchol2(lvals(locf),ldl,wmem(myfrontptr),ldf,fdim,rank,
     1       fdim,linds(lptrs(3,k)),dfopts,ifopts,info)
      if (info.gt.0) goto 1111 

      else
        locf = myfrontptr + (ldf+1) * (rank-1)
      newstakptr = mystakptr + rank
        do 130 j1 = rank, 1, -1
        mystak(newstakptr - j1) = j
          myi = lptrs(3,j)
          locf = locf - lptrs(3,j)
          do 140 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1
            do while (linds(myi) .ne. ainds(i1))
              myi = myi + 1
            end do
            wmem(locf + myi) = wmem(locf + myi) + avals(i1)
140       continue
          locf = locf + lptrs(3,j) - ldf - 1
          j = tinds(tptrs(1,j))
130     continue

        info = 0 
          call dpotrf('l',rank,wmem(myfrontptr),ldf,info)
        if (info.gt.0) goto 1111 
          call dtrsm('R','L','T','N',fdim-rank,rank,1.d0,
     1       wmem(myfrontptr),ldf,wmem(myfrontptr+rank),ldf)
      
          call dsyrk('L','N',fdim-rank,rank,-1.d0,
     1       wmem(myfrontptr+rank),ldf,1.d0,
     2       wmem(myfrontptr+(ldf+1)*rank),ldf)

      locf = myfrontptr
      do j1 = newstakptr - 1, mystakptr, -1
        j = mystak(j1)
        call mydc(lptrs(2,j),wmem(locf),lvals(lptrs(1,j)))
        locf = locf + ldf + 1
      end do
      end if

      do 10 i = mystakptr - 2, 1, -2
        udim = fdim
        kid = node
        frontptr = myfrontptr
        node = mystak(i)
        myfrontptr = mystak(i+1)
        k = sup(tptrs(3,node))
        myinds = lptrs(3,k)
        fdim = lptrs(2,k)
        myj = myinds
        mylim = myinds + fdim
        locptr = iptrs(2,kid) - 1

        locu = frontptr + (ldf+1) * rank
        locf = myfrontptr - 1
        kincr = 1
        kidincr = 1
        myj = 1
        locindj = 0
        do 215 j1 = iptrs(1,kid), iptrs(1,kid) + locptr - 2, 2
        locindj = locindj + lc(j1)
          do 225 myj = myj, locindj
            myi = kincr
          locindi = locindj - lc(j1)
            do 235 i1 = j1, iptrs(1,kid) + locptr, 2
              locu = locu - myi
            locindi = locindi + lc(i1)
              do 246 myi = myi, locindi
                wmem(locf+myi) = wmem(locu+myi)
246           continue
              locu = locu + myi
            locindi = locindi + lc(i1+1)
              do 256, myi = myi, locindi
                wmem(locf+myi) = 0.d0
256           continue
235         continue
            locu = locu + rank + kidincr + ldf - udim
            locf = locf + ldf
            kincr = kincr + 1
            kidincr = kidincr + 1
225       continue
        locindj = locindj + lc(j1+1)
          do 265 myj = myj, locindj
            do 276 i1 = kincr, fdim
              wmem(locf+i1) = 0.d0
276         continue
            locf = locf + ldf
            kincr = kincr + 1
265       continue
215     continue

        if ((locf+kincr) .ne. locu) then
        locindj = locindj + lc(j1)
          do 525 myj = myj, locindj
            myi = kincr
          locindi = locindj - lc(j1)
            do 535 i1 = j1, iptrs(1,kid) + locptr, 2
              locu = locu - myi
            locindi = locindi + lc(i1)
              do 546 myi = myi, locindi
                wmem(locf+myi) = wmem(locu+myi)
546           continue
              locu = locu + myi
            locindi = locindi + lc(i1+1)
              do 555 myi = myi, locindi
                wmem(locf+myi) = 0.d0
555           continue
535         continue
            locu = locu + rank + kidincr + ldf - udim
            kincr = kincr + 1
            locf = locf + ldf
            kidincr = kidincr + 1
525       continue
        locindj = locindj + lc(j1+1)
          do 565 myj = myj, locindj
            do 576 i1 = kincr, fdim
              wmem(locf+i1) = 0.d0
576         continue
            locf = locf + ldf
Cmj            kicr = kicr + 1
565       continue
        end if

        rank = sup(tptrs(3,k)+1)
        j = node
        locf = myfrontptr + (ldf+1) * (rank-1)
      newstakptr = mystakptr + rank
        do 30 j1 = rank, 1, -1
        mystak(newstakptr - j1) = j
          myi = lptrs(3,j)
          locf = locf - myi
          do 40 i1 = aptrs(1,j), aptrs(1,j) + aptrs(2,j) - 1
            do while (linds(myi) .ne. ainds(i1))
              myi = myi + 1
            end do
            wmem(locf + myi) = wmem(locf + myi) + avals(i1)
40        continue
          locf = locf + lptrs(3,j) - ldf - 1
          j = tinds(tptrs(1,j))
30      continue

        info = 0 
          call dpotrf('l',rank,wmem(myfrontptr),ldf,info)
          if (info.gt.0) goto 1111 
          call dtrsm('R','L','T','N',fdim-rank,rank,1.d0,
     1       wmem(myfrontptr),ldf,wmem(myfrontptr+rank),ldf)
          call dsyrk('L','N',fdim-rank,rank,-1.d0,
     1       wmem(myfrontptr+rank),ldf,1.d0,
     2       wmem(myfrontptr+(ldf+1)*rank),ldf)

      locf = myfrontptr
      do j1 = newstakptr - 1, mystakptr, -1
        j = mystak(j1)
        call mydc(lptrs(2,j),wmem(locf),lvals(lptrs(1,j)))
        locf = locf + ldf + 1
      end do

10    continue
      info=0 
      return
1111  continue
      info=1 
      end
