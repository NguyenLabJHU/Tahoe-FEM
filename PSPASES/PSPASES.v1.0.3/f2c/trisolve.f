
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   trisolve.f                                                              */
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
C/* $Id: trisolve.f,v 1.2 2004-12-13 08:18:17 paklein Exp $ */
C/*****************************************************************************/

      subroutine trisolve(N,rowdist,order,lptrs,linds,lvals,
     +                    tptrs,tinds,sup,tsind,lc,iptrs,
     +                    ifopts,nrhs,options,rhso,ldo,rhsc,
     +                    ldc,ranmasks,comm,
     +                    hvbtemp, lrud,
     +                    wrkord0, wrkord1, wrkord2,
     +                    ty, dworkmj, ordvals,
     +                    mynnodes2, bnrhs2, ordvalsiz2, wsolvesize2, 
     +                    trhsize2, wrkord1siz2, wrkord2siz2, ns2, dd2, 
     +                    maxvsize2, hvbsize2)

      implicit none
      include 'mpif.h'

      double precision zero
      parameter(zero=0.d0)

      integer rowdist(0:*),order(0:*),lptrs(3,0:*),linds(*)
      integer tptrs(3,0:*),tinds(*),sup(*),tsind(*),lc(*),iptrs(2,0:*)
      integer ifopts(0:*),options(0:*),nrhs,ldo,ldc,ranmasks(5,0:*)
      double precision lvals(*),rhso(0:ldo-1,*),rhsc(0:ldc-1,*)

      integer N,dd,lgblk,wsolvesize,ns,supindsize,myid,myidh,myidv,pp
      integer comm,nsupnode,supsize,maxvsize,maxhsize,i,j,k,m,hvbsize
      integer uindptr,uvecptr,recvptr,rhsptr,trhsize,is1,nown,ierr
      integer mynnodes,uvlptr,bnrhs,sr,ir,rnr
      integer wrkord1siz,wrkord2siz,ordvalsiz
      integer pisc,pisd,pirc,pird,pdsc,pdsd,pdrc,pdrd
      integer piscx,pisdx,pircx,pirdx,pdscx,pdsdx,pdrcx,pdrdx,ppgr
      integer piown,psloc,pslocx
      integer psv,prv,psvx,prvx
C     integer, allocatable :: hvbtemp(:),lrud(:)
      integer hvbtemp, lrud
      integer ns2, dd2, maxvsize2, hvbsize2
      dimension hvbtemp(hvbsize2 + maxvsize2)
      dimension lrud(ns2+ 4*dd2)

C     integer, allocatable :: wrkord0(:), wrkord1(:), wrkord2(:)
      integer wrkord0, wrkord1, wrkord2
      integer mynnodes2, wrkord1siz2, wrkord2siz2
      dimension wrkord0(0:mynnodes2-1)
      dimension wrkord1(0:wrkord1siz2-1)
      dimension wrkord2(0:wrkord2siz2-1)

C     double precision, allocatable :: ty(:,:),dworkmj(:)
      double precision ty, dworkmj
      integer bnrhs2, wsolvesize2, trhsize2
      dimension ty(0:N-1,bnrhs2)
      dimension dworkmj(wsolvesize2 + 4*trhsize2)

C     double precision, allocatable :: ordvals(:)
      double precision ordvals
      integer ordvalsiz2
      dimension ordvals(0:ordvalsiz2-1)

      N = ifopts(0)
      dd = ifopts(1)
      lgblk = ifopts(2)
      wsolvesize = ifopts(3) 
      ns = ifopts(4)
      supindsize = ifopts(5)
      myid = ifopts(6)
      myidh = ifopts(7)
      myidv = ifopts(8)
      supsize = ifopts(9)

      bnrhs = min(options(0),nrhs)

      pp = ishft(1,dd)

C     allocate(lrud(ns+4*dd),stat=is1)
      if(is1.ne.0) then
        print *,'Allocate error'
        call mpi_abort(comm,112,ierr)
      end if

      call getmysnodes(N-1,sup,tinds,tptrs,N,supsize,lrud,nsupnode,
     +                 lrud(1+ns),dd,maxhsize,maxvsize,ns,myid)

      i = sup(tptrs(3,lrud(nsupnode))) 
      maxvsize = max(maxvsize,lptrs(2,i))

      hvbsize=nsupnode*(10+maxvsize+3*((maxhsize+1)/2+(maxvsize+1)/2))
C     allocate(hvbtemp(hvbsize+maxvsize),stat=is1)
      if(is1.ne.0) then
        print *,myid,': Cannot allocate memory for hvbtemp',hvbsize
        call mpi_abort(comm,113,ierr)
      end if
      uindptr = 1+hvbsize

      call getmyhvb(lrud,nsupnode,sup,supsize,tsind,
     +      supindsize,tptrs,tinds,N,dd,lgblk,hvbtemp,
     +      hvbsize,myidh,myidv,myid)

      wsolvesize = 4 *wsolvesize * nrhs
      trhsize = maxvsize * nrhs
C     allocate(dworkmj(wsolvesize+4*trhsize),stat=is1)
      if(is1.ne.0) then
        print *,myid,': Memory Allocation Failure'
        call mpi_abort(comm,111,ierr)
      end if
      uvecptr = 1+wsolvesize
      uvlptr  = uvecptr+trhsize
      recvptr = uvlptr+trhsize
      rhsptr  = recvptr+trhsize

C     allocate(ty(0:N-1,bnrhs),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      mynnodes = rowdist(myid+1)-rowdist(myid)
C     allocate(wrkord0(0:mynnodes-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      wrkord1siz = 19*pp
C     allocate(wrkord1(0:wrkord1siz-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      pisc  = 0
      pisd  = pp
      pirc  = 2*pp
      pird  = 2*pp+pp
      pdsc  = 4*pp
      pdsd  = 4*pp+pp
      pdrc  = 4*pp+2*pp
      pdrd  = 4*pp+2*pp+pp
      piscx = 8*pp
      pisdx = 8*pp+pp
      pircx = 8*pp+2*pp
      pirdx = 8*pp+2*pp+pp
      pdscx = 8*pp+4*pp
      pdsdx = 8*pp+4*pp+pp
      pdrcx = 8*pp+4*pp+2*pp
      pdrdx = 8*pp+4*pp+2*pp+pp
      ppgr  = 16*pp

      call preordbe(N,order,ranmasks,nown,wrkord1(pisc),wrkord1(pisd),
     +              wrkord1(pirc),wrkord1(pird),mynnodes,dd,myid,lgblk,
     +              wrkord1(ppgr),wrkord0,comm)
      
      wrkord2siz = 3*nown+mynnodes
C     allocate(wrkord2(0:wrkord2siz-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      piown  = 0
      psloc  = 2*nown
      pslocx = 2*nown+mynnodes

      ordvalsiz = 2*max(mynnodes,nown)*bnrhs
C     allocate(ordvals(0:ordvalsiz-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation error'
        call mpi_abort(comm,0,ierr)
      end if

      psv  = 0
      prv  = mynnodes*bnrhs
      psvx = 0
      prvx = nown*bnrhs

      call preordbc(N,order,ranmasks,nown,wrkord1(pisc),wrkord1(pisd),
     +              wrkord1(pirc),wrkord1(pird),wrkord2(piown),
     +              rowdist,mynnodes,wrkord2(psloc),ordvals(psv),
     +              dd,myid,lgblk,wrkord1(ppgr),
     +              wrkord0,comm)

      call preordxc(N,wrkord2(piown),rowdist,mynnodes,wrkord1(piscx),
     +              wrkord1(pisdx),wrkord1(pircx),wrkord1(pirdx),
     +              nown,wrkord2(pslocx),wrkord0,
     +              ordvals(psvx),dd,myid,ordvals(prvx),comm)

      do i=0,pp-1
        wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*bnrhs
        wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*bnrhs
        wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*bnrhs
        wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*bnrhs
        wrkord1(pdscx+i) = wrkord1(piscx+i)*bnrhs
        wrkord1(pdsdx+i) = wrkord1(pisdx+i)*bnrhs
        wrkord1(pdrcx+i) = wrkord1(pircx+i)*bnrhs
        wrkord1(pdrdx+i) = wrkord1(pirdx+i)*bnrhs
      end do

      do sr=1,nrhs-bnrhs+1,bnrhs  

      call reordb(N,rhso(0,sr),ldo,bnrhs,ty,nown,
     +            wrkord1(pdsc),wrkord1(pdsd),wrkord1(pdrc),
     +            wrkord1(pdrd),wrkord2(piown),mynnodes,
     +            wrkord2(psloc),ordvals(psv),ordvals(prv),
     +            comm)

      if(bnrhs.ne.1) then
        call pfsolvem(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,dd,lgblk,N,
     +            ty,bnrhs,dworkmj(rhsptr),dworkmj(uvecptr),
     +            dworkmj(uvlptr),dworkmj(recvptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      else
        call pfsolve1(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,dd,lgblk,N,
     +            ty,dworkmj(rhsptr),dworkmj(uvecptr),
     +            dworkmj(uvlptr),dworkmj(recvptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      end if

      if(bnrhs.ne.1) then
        call pbsolvem(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,myidv,dd,lgblk,N,
     +            ty,bnrhs,dworkmj(uvecptr),
     +            dworkmj(recvptr),dworkmj(rhsptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)

      else
        call pbsolve1(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,myidv,dd,lgblk,N,
     +            ty,dworkmj(uvecptr),
     +            dworkmj(recvptr),dworkmj(rhsptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      end if

      call reordx(N,rhsc(0,sr),ldc,bnrhs,ty,nown,wrkord1(pdscx),
     +            wrkord1(pdsdx),wrkord1(pdrcx),wrkord1(pdrdx),
     +            wrkord2(piown),rowdist,mynnodes,
     +            wrkord2(pslocx),wrkord0,ordvals(psvx),
     +            ordvals(prvx),myid,comm)

      end do 

      rnr = nrhs-sr+1

      if(rnr.ne.0) then

      do i=0,pp-1
        wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*rnr
        wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*rnr
        wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*rnr
        wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*rnr
        wrkord1(pdscx+i) = wrkord1(piscx+i)*rnr
        wrkord1(pdsdx+i) = wrkord1(pisdx+i)*rnr
        wrkord1(pdrcx+i) = wrkord1(pircx+i)*rnr
        wrkord1(pdrdx+i) = wrkord1(pirdx+i)*rnr
      end do

      call reordb(N,rhso(0,sr),ldo,rnr,ty,nown,
     +            wrkord1(pdsc),wrkord1(pdsd),wrkord1(pdrc),
     +            wrkord1(pdrd),wrkord2(piown),mynnodes,
     +            wrkord2(psloc),ordvals(psv),ordvals(prv),
     +            comm)

      if(rnr.ne.1) then
        call pfsolvem(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,dd,lgblk,N,
     +            ty,rnr,dworkmj(rhsptr),dworkmj(uvecptr),
     +            dworkmj(uvlptr),dworkmj(recvptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      else
        call pfsolve1(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,dd,lgblk,N,
     +            ty,dworkmj(rhsptr),dworkmj(uvecptr),
     +            dworkmj(uvlptr),dworkmj(recvptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      end if

      if(rnr.ne.1) then
        call pbsolvem(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,myidv,dd,lgblk,N,
     +            ty,rnr,dworkmj(uvecptr),
     +            dworkmj(recvptr),dworkmj(rhsptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)

      else
        call pbsolve1(lrud,nsupnode,sup,
     +            lptrs,linds,lvals,tptrs,tinds,
     +            myid,myidh,myidv,dd,lgblk,N,
     +            ty,dworkmj(uvecptr),
     +            dworkmj(recvptr),dworkmj(rhsptr),
     +            hvbtemp(uindptr),maxvsize,lrud(1+ns),
     +            hvbtemp,hvbsize,lc,dworkmj,iptrs,comm)
      end if

      call reordx(N,rhsc(0,sr),ldc,rnr,ty,nown,wrkord1(pdscx),
     +            wrkord1(pdsdx),wrkord1(pdrcx),wrkord1(pdrdx),
     +            wrkord2(piown),rowdist,mynnodes,
     +            wrkord2(pslocx),wrkord0,ordvals(psvx),
     +            ordvals(prvx),myid,comm)

      end if 

C      deallocate(lrud)
C      deallocate(hvbtemp)
C      deallocate(ty)
C      deallocate(dworkmj)
C      deallocate(wrkord0)
C      deallocate(wrkord1)
C      deallocate(wrkord2)
C      deallocate(ordvals)

      end
