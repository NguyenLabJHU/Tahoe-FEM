C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   gen_lc.f                                                                */
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
C/* $Id: gen_lc.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      recursive
     +subroutine GEN_LC(root,lc,linds,lptrs,tinds,tptrs,sup,iptrs,
     1                 lcsize,wa1,wa2,lctr,wsmy,wstot,wstri)
      integer root,lcsize,lctr,lc(*),linds(*),lptrs(3,0:*),sup(*)
      integer tinds(*),tptrs(3,0:*),iptrs(2,0:*), wsmy, wstot, wstri
      integer wa1(0:*),wa2(0:*)

      itbot = tptrs(3,root)
      jbot = sup(itbot)
      wsmy = lptrs(2,jbot)
      wstri = wsmy
      do while (tptrs(2,jbot) .eq. 1)
        kid = tinds(tptrs(1,jbot))
        kbot = sup(tptrs(3,kid))
        iptrs(1,kid) = lctr
        iptrs(1,kbot) = lctr
        call gen_l(jbot,kid,ln,lc(lctr),linds,lptrs)
        lctr = lctr + ln
        iptrs(2,kid) = ln
        iptrs(2,kbot) = ln
        jbot = kbot
        j = lptrs(2,kbot)
        wsmy = max(wsmy,j)
        wstri = wstri + j

      end do

      wstot = 0
      itptr = tptrs(2,jbot)

      maxwstri = 0

      do i = 0, itptr - 1
        kid = tinds(tptrs(1,jbot)+i)
        kbot = sup(tptrs(3,kid))
        iptrs(1,kid) = lctr
        iptrs(1,kbot) = lctr

        call gen_l(jbot,kid,ln,lc(lctr),linds,lptrs)

        lctr = lctr + ln 
        iptrs(2,kid) = ln
        iptrs(2,kbot) = ln

        call gen_lc(kid,lc,linds,lptrs,tinds,tptrs,sup,iptrs,
     1              lcsize,wa1(itptr),wa2(itptr),lctr,wa1(i),wa2(i),
     2              kwstri)
        maxwstri = max(maxwstri,kwstri)

      end do

      wstri = wstri + maxwstri
      if (itptr .gt. 0) then
       kptr = tptrs(1,jbot)
       do i = 1, itptr - 1
        if (wa1(i) .gt. wa1(0)) then
          it = wa1(0)
          wa1(0) = wa1(i)
          wa1(i) = it
          it = wa2(0)
          wa2(0) = wa2(i)
          wa2(i) = it
          it = tinds(kptr)
          tinds(kptr) = tinds(kptr+i)
          tinds(kptr+i) = it
        end if
       end do

       wsmy = max(wsmy, wa1(0))
       wstot = wa2(0)
       do i = 1, itptr - 1
        wstot = max(wstot, wa1(i)*wa1(i)+wa2(i))
       end do
      end if

      sup(itbot+3) = wsmy       
      return
      end
C------------------------------------------------------------------------------
      subroutine GEN_L(par,kid,ln,lc,linds,lptrs)
      integer par,kid,ln,lc(0:*),linds(*),lptrs(3,0:*)

      i = lptrs(3,par)
      j = lptrs(3,kid) 
      ilim = i + lptrs(2,par)
      jlim = j + lptrs(2,kid)
      j = j + 1
      ln = 0
      if (j .eq. jlim) then
        ln = 2
        lc(0) = 0
        lc(1) = ilim - i
      else
        do while (j .lt. jlim .and. i .lt. ilim)
          k = i
          do while (j .lt. jlim .and. linds(j) .eq. linds(i))
            i = i + 1
            j = j + 1
          end do 
          lc(ln) = i - k
          if (j .lt. jlim) then
            k = i
            kl = linds(j)
            do while (linds(i) .ne. kl)
              i = i + 1
            end do
            lc(ln+1) = i - k
          else
            lc(ln+1) = ilim - i
          end if
          ln = ln + 2
        end do
      end if
      return
      end
