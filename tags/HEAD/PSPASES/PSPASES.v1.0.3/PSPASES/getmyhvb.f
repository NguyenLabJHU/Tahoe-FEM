C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   getmyhvb.f                                                              */
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
C/* $Id: getmyhvb.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      subroutine getmyhvb(mysnodes,nsupnode,sup,ssize,supinds,
     +                    supindsize,tptrs,tinds,N,dd,lgblk,hvbndry,
     +                    hvbsize,myidh,myidv,myid)

      implicit none

      integer lendp
      parameter(lendp=8)

      integer nsupnode,ssize,supindsize,dd,lgblk
      integer hvbsize,myidv,myidh,N
      integer sup(*), supinds(*)
      integer tptrs(3,0:*),tinds(*)
      integer mysnodes(*),hvbndry(*)

      integer node,i,j,k,rowcol,level,level2,nb,nb_per_proc
      integer myleft,myright,myup,mydown,bmaskh,hsize,bmaskv,vsize
      integer hvbp,supptr,supbot,supsiz,lbotsiz,sptr,is,bnode
      integer bhs,bvs,bhp,bvp,pclean,currpos,pos,nvbg,pactual
      integer ph,pch,pv,pcv,nvb,nvbt,vsizer,bip,bis,pcleanh
      integer mymaskh,mymaskv,suptop,istflag
      integer myid,nbleft,nbleft1,nbleft2,poshdiff,posvdiff

      level = dd
      rowcol = 0
      bmaskh = 0
      bmaskv = 0
      hsize = 1
      vsize = 1
      mymaskh = 0
      mymaskv = 0

      hvbp = 1
      do is = nsupnode-1,1,-1

        suptop = mysnodes(is)
        supptr = tptrs(3,suptop)

        if (tptrs(2,sup(supptr)).ne.1) then
          level = level-1
          rowcol = 1-rowcol
          level2 = ishft(dd-level-1,-1)
          if(rowcol.eq.1) then      
            bmaskh = ior(bmaskh,ishft(1,level2))
            hsize = ishft(hsize,1)
            mymaskh = iand(myidh,bmaskh)
          else                  
            bmaskv = ior(bmaskv,ishft(1,level2))
            vsize = ishft(vsize,1)
            mymaskv = iand(myidv,bmaskv)
          end if
        end if
        
        supbot = sup(supptr)
        supsiz = sup(supptr+1)
        sptr   = sup(supptr+2)
        lbotsiz = sup(supptr+3)

        istflag = 0

        bhs = hvbp+6
        bhp = bhs
        i = 0
        bnode = ishft(supinds(sptr),-lgblk)
        pch = iand(bnode,bmaskh)
        pclean  = pch
        pactual = pch
        currpos = mod(mymaskh+hsize-pclean,hsize)
        pos = 0
        do while (i .lt. supsiz)
          j = i
          i = i+1
          do while (i.lt.supsiz)
            bnode = ishft(supinds(sptr+i),-lgblk)
            ph = iand(bnode,bmaskh)
            if(pch.ne.ph) goto 10
            i = i + 1
          end do

 10          if(pos.ge.currpos) then  
            hvbndry(bhp) = supinds(sptr+j)
            hvbndry(bhp+2) = supinds(sptr+i-1)
            if(pos.eq.currpos) then
              hvbndry(bhp+1) = i-j
            else
              hvbndry(bhp+1) = 0
            end if
            bhp = bhp+3
            currpos = currpos + hsize
          end if
          pclean = mod(pclean+1,hsize)
          pos = pos + 1 + mod(ph+hsize-pclean,hsize)
          pclean = ph
          pch = ph
        end do
        hvbndry(bhs-1) = (bhp-bhs)/3 

        bis = bhp+1
        bvs = bhp+4+lbotsiz
        bip = bis
        bvp = bvs

        node  = supinds(sptr)
        i = 0
        pos = 0
        bnode = ishft(node,-lgblk)
        nvb = 0
        nvbt = 0
        nvbg = 0
        pcv = iand(bnode,bmaskv)
        pch = iand(bnode,bmaskh)
        currpos = mod(mymaskv+vsize-pcv,vsize)
        pclean = pcv
        pcleanh = pch

        do while(i.lt.supsiz)
          j = i
          if(pos.eq.currpos) then
            hvbndry(bip) = node
            bip = bip+1
          end if
          i = i+1
          do while (i.lt.supsiz)
            node  = supinds(sptr+i)
            bnode = ishft(node,-lgblk)
            ph = iand(bnode,bmaskh)
            pv = iand(bnode,bmaskv)
            if((ph.ne.pch).or.(pv.ne.pcv)) goto 20

            if(pos.eq.currpos) then
              hvbndry(bip) = node
              bip = bip+1
            end if
            i = i + 1
          end do

 20          if(pos.ge.currpos) then  
            hvbndry(bvp) = supinds(sptr+j)
            hvbndry(bvp+2) = supinds(sptr+i-1)
            if(pos.eq.currpos) then
              hvbndry(bvp+1) = i-j
            else
              hvbndry(bvp+1) = 0
            end if
            bvp = bvp+3
            nvbt = nvbt+1
            currpos = currpos + vsize
          end if
          pcleanh = mod(pcleanh+1,hsize)
          if(pv.ne.pcv .and. i.lt.supsiz) then
            pclean = mod(pclean+1,vsize)
            posvdiff = mod(pv+vsize-pclean,vsize)
            poshdiff = mod(ph+hsize-pcleanh,hsize)
            pos = pos + 1 + posvdiff
            if(poshdiff.ne.posvdiff) then 
              hvbndry(bvp)   = node
              hvbndry(bvp+1) = 0
              hvbndry(bvp+2) = node
              bvp = bvp+3
              nvbt = nvbt+1
              nvbg = nvbg+vsize 
            end if
            pclean = pv
          else
            pos = pos + vsize
          end if
          pcleanh = ph
          nvbg = nvbg+1
          pch = ph
          pcv = pv
        end do
        nvb = nvbt

        vsizer = 0
        i = supsiz
        node  = supinds(sptr+i)
        bnode = ishft(node,-lgblk)
        pclean = iand(bnode,bmaskv)
        pactual = pclean
        currpos = mod(mymaskv+vsize-pclean,vsize)
        pos = 0
        do while(i.lt.lbotsiz)
          j = i
          if(pos.eq.currpos) then
            hvbndry(bip) = node
            bip = bip+1
          end if
          i = i+1 
          do while(i.lt.lbotsiz)
            node  = supinds(sptr+i)
            bnode = ishft(node,-lgblk)
            pv = iand(bnode,bmaskv)
            if(pv.ne.pactual) goto 30

            if(pos.eq.currpos) then
              hvbndry(bip) = node
              bip = bip+1
            end if
            i = i+1
          end do

 30          if(pos.ge.currpos) then
            hvbndry(bvp) = supinds(sptr+j)
            hvbndry(bvp+2) = supinds(sptr+i-1)
            if(pos.eq.currpos) then
              hvbndry(bvp+1) = i-j
            else
              hvbndry(bvp+1) = 0
            end if
            vsizer = vsizer+hvbndry(bvp+1)
            nvb = nvb+1
            bvp = bvp+3
            currpos = currpos + vsize
          end if
          pclean = mod(pclean+1,vsize)
          pactual = pv
          pos = pos + 1 + mod(pv+vsize-pclean,vsize)
          pclean = pv
        end do

        nb=1

        bnode = ishft(suptop,-lgblk)
        pv = iand(bnode,bmaskv)
        ph = iand(bnode,bmaskh)
        if(pv.eq.mymaskv) then
          istflag = 2
          if(ph.eq.mymaskh) istflag = 1
        end if

        hvbndry(bvs-1)  = nvb      
        hvbndry(bvs-2)  = nvb-nvbt 
        hvbndry(bvs-3)  = vsizer   
        hvbndry(hvbp)   = bvp+1
        hvbndry(hvbp+1) = bvs
        hvbndry(hvbp+2) = nb
        hvbndry(hvbp+3) = istflag
        hvbndry(hvbp+4) = lbotsiz
        hvbndry(bis-1)  = bip-bis  
        hvbndry(bvp) = hvbp
        hvbp = bvp+1

      end do
      hvbsize = hvbp-1

      end 
