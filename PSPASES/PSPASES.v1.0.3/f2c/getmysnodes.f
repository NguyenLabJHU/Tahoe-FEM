C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   getmysnodes.f                                                           */
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
C/* $Id: getmysnodes.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine getmysnodes(root,sup,tinds,tptrs,N,supsize,mysnodes,
     +                 nsupnode,lrud,dd,maxhsize,maxvsize,ns,myid)

      implicit none
      integer root,N,supsize,nsupnode,dd,maxhsize,maxvsize,ns,myid
      integer sup(*),tinds(*),tptrs(3,0:*),mysnodes(*)
      integer lrud(*)

      integer mydown,myup,myright,myleft,level,level2,supptr,rowcol
      integer i,j,node,nbrp,jl,jr,ju,jd

      maxhsize = 0
      maxvsize = 0
      j = 1
      node = root
      mysnodes(j) = node
      j = j+1
      supptr = tptrs(3,node)
      maxhsize = max(sup(supptr+1),maxhsize)
      maxvsize = max(sup(supptr+3),maxvsize)
      node = sup(supptr) 
      do i=1,dd
        do while (tptrs(2,node).eq.1) 
          node = tinds(tptrs(1,node)) 
          mysnodes(j) = node
          j = j+1
          supptr = tptrs(3,node)
          maxhsize = max(sup(supptr+1),maxhsize)
          maxvsize = max(sup(supptr+3),maxvsize)
          node = sup(supptr) 
        end do
        if(iand(myid,ishft(1,dd-i)).eq.0) then       
          node = tinds(tptrs(1,node))
        else                                    
          node = tinds(tptrs(1,node)+1)
        end if
        mysnodes(j) = node
        j = j+1
        supptr = tptrs(3,node)
        maxhsize = max(sup(supptr+1),maxhsize)
        maxvsize = max(sup(supptr+3),maxvsize)
        node = sup(supptr)
      end do
      nsupnode = j-1 

      rowcol = 0
      myleft = myid
      myright = myid
      myup = myid
      mydown = myid
      jl = 0
      jr = 0
      ju = 1
      jd = 1

      nbrp = 1
      do level=dd-1,0,-1
        rowcol = 1-rowcol
        level2 = ishft(dd-level-1,-1)

        if(rowcol.eq.1) then      
          if(myleft.ge.myid) then
            myleft = ieor(myleft,ishft(1,jl))
            jl = jl+2
          endif
          if(myright.le.myid) then
            myright = ieor(myright,ishft(1,jr))
            jr = jr+2
          end if
        else                  
          if(myup.ge.myid) then
            myup = ieor(myup,ishft(1,ju))
            ju = ju+2
          end if
          if(mydown.le.myid) then
            mydown = ieor(mydown,ishft(1,jd))
            jd = jd+2
          end if
        endif
        lrud(nbrp)   = myleft
        lrud(nbrp+1) = myright
        lrud(nbrp+2) = myup
        lrud(nbrp+3) = mydown
        nbrp = nbrp + 4
      end do


      end
