                   
C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   pfsolve1.f                                                              */
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
C/* $Id: pfsolve1.f,v 1.1 2004-12-10 20:26:45 paklein Exp $ */
C/*****************************************************************************/

      subroutine pfsolve1(mysnodes,nsupnode,sup,
     +                  lptrs,linds,lvals,tptrs,tinds,
     +                  myid,myidh,dd,lgblk,N,
     +                  rhsc,rhs,uvec,uvl,recvec,
     +                  uinds,maxvsize,lrud,
     +                  hvbndry,hvbsize,lc,w,iptrs,comm)

      implicit none
      include 'mpif.h'

      integer N,dd,lgblk,myid,myidh,maxvsize
      integer nsupnode,hvbsize,comm

      integer lptrs(3,0:*), linds(*), iptrs(2,0:*)
      integer tptrs(3,0:*), tinds(*)
      integer sup(*), mysnodes(*)
      integer lrud(*), hvbndry(*)
      integer uinds(*), lc(*)

      double precision lvals(*)
      double precision rhsc(0:N-1),w(*)
      double precision uvec(*),uvl(*)
      double precision rhs(*),recvec(*)

      integer lendp,loglendp,itype
      parameter(lendp=8,loglendp=3,itype=1)
      double precision one,ngone,zero
      parameter(zero=0.d0,one=1.d0,ngone=-1.d0)

      integer root
      integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk
      integer myleft,myright,myup,mydown
      integer supptr,suptop,supbot,supsiz,lbotsiz
      integer bvs,nb,nhb,nvb,nbrp,hvbp,ih,bhp
      integer bhstrt,bhsize,bhend,lvalst,ldalb
      integer rhsst,uvecst,uvlst,itag
      integer iv,bvp,bvstrt,bvsize,bvend
      integer recvsize,cmpsize
      integer npending,npendings,npendings2
      integer mid,nbrecv
      integer lppar,kidtop,pparbot,kid
      integer diffbits,partner,ik,ldauk,lpkid,pkid
      integer i,j,k1,k2,is,usize
      integer mymaskh,mymaskhk,ldaup
      integer iamkid,iampar,kidr,kidl,pkido,bip,bip1,istflag
      integer vsizer,nbode,levelk,iaml,bnode,lvp,ivlim,nvinds,ntinds
      integer flagr,flags,indh,hsize
      integer cbsize,msizedp
      integer mpistat(MPI_STATUS_SIZE),req(5),ierr
      integer statall(MPI_STATUS_SIZE,5)

      do i=1,5
        req(i) = MPI_REQUEST_NULL
      end do

      root   = mysnodes(1)

      suptop = mysnodes(nsupnode)
      call fsolve(N,lvals,linds,lptrs,tinds,tptrs,sup,rhsc,1,
     +                suptop,lc,iptrs,w)

      supptr = tptrs(3,suptop)
      supbot = sup(supptr)
      bhsize = sup(supptr+1)
      ldalb  = lptrs(2,supbot)
      lpkid  = lptrs(3,supbot)

      ldauk  = ldalb-bhsize  
      lpkid = lpkid+bhsize   
                
      do i=1,ldauk
        uvl(i) = w(i+bhsize)
      end do

      level = dd  
      rowcol = 0  
      bmaskh = 0  
      bmaskv = 0  
      mymaskh = 0 
      hsize  = 1  

      hvbp = 1    
      nbrp = -3   

      do is = nsupnode-1,1,-1


        bnode = ishft(suptop,-lgblk)
        pkid = iand(bnode,bmaskh)

        suptop = mysnodes(is)        
        supptr = tptrs(3,suptop)      
        supbot = sup(supptr)            
        supsiz = sup(supptr+1)      
        bmaskhk = bmaskh            
        bmaskvk = bmaskv            
        mymaskhk = mymaskh            
        levelk  = level                  

        if (tptrs(2,supbot).ne.1) then 
          level = level-1
          rowcol = 1-rowcol
          level2 = ishft(dd-level-1,-1)
          if(rowcol.eq.1) then  
            bmaskh = ior(bmaskh,ishft(1,level2))
            mymaskh = iand(myidh,bmaskh)
            hsize = ishft(hsize,1)
          else                  
            bmaskv = ior(bmaskv,ishft(1,level2))
          end if
          nbrp = nbrp + 4 
        end if

        bnode = ishft(supbot,-lgblk)
        ldaup = lptrs(2,supbot)       
        lppar = lptrs(3,supbot)      
        pparbot = iand(bnode,bmaskh)

        iamkid = 0 
        iampar = 0 
        if(mymaskhk.eq.pkid) iamkid = 1
        if(mymaskh.eq.pparbot) iampar = 1
        kidl = tinds(tptrs(1,supbot)) 
        kidr = tinds(tptrs(1,supbot)+1) 
        iaml = 1 
        if(iand(myid,ishft(1,dd-level-1)).ne.0) iaml = 0

        cbsize  = ldauk*lendp 
        msizedp = maxvsize*lendp

        if(bmaskv.eq.bmaskvk) then
          partner  = myid
          if(iamkid.eq.1 .and. iampar.eq.1) then
            if(levelk.ne.level) then  
              kid = kidr
              if(iaml.eq.0) kid = kidl
              bnode = ishft(kid,-lgblk)
              pkido = iand(bnode,bmaskhk)
              diffbits = ieor(ieor(pkid,pkido),ishft(1,level2))
              do j = 0,level2
                partner = ieor(partner,ishft(iand(diffbits,1),j*2))
                diffbits = ishft(diffbits,-1)
              end do
              call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                        itype,comm,req(2),ierr)
              call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                        itype,comm,req(1),ierr)
              npending = 2
              do while(npending.gt.0)
                call mpi_waitany(2,req,mid,mpistat,ierr)
                if(mid.eq.1) then 
                  call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
                end if
                npending = npending-1
              end do
              call x10dad1(uvec,ldaup,uvl,ldauk,recvec,usize,
     +              linds(lppar),linds(lpkid),uinds)
            else
              call extend_op1(uvec,ldaup,uvl,ldauk,
     +              linds(lppar),linds(lpkid))
            end if
          end if
          if(iamkid.eq.0 .and. iampar.eq.1) then
            diffbits = ieor(pkid,mymaskhk)
            do j = 0,level2
              partner = ieor(partner,ishft(iand(diffbits,1),j*2))
              diffbits = ishft(diffbits,-1)
            end do
            call mpi_irecv(uvl,msizedp,MPI_BYTE,partner,
     +                     itype,comm,req(3),ierr)
            npending = 1

            if(levelk.ne.level) then
              kid = kidr
              if(iaml.eq.0) kid = kidl
              bnode = ishft(kid,-lgblk)
              pkido = iand(bnode,bmaskhk)
              diffbits = ieor(pparbot,pkido)
              if(iaml.eq.1) diffbits = ieor(diffbits,ishft(1,level2))
              partner = myid
              do j = 0,level2
                partner = ieor(partner,ishft(iand(diffbits,1),j*2))
                diffbits = ishft(diffbits,-1)
              end do
              call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                        itype,comm,req(2),ierr)
              call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                        itype,comm,req(1),ierr)
              npending = npending+2
              do while(npending.gt.0)
                call mpi_waitany(3,req,mid,mpistat,ierr)
                if(mid.eq.1) then
                  call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
                end if
                npending = npending-1
              end do
              call x10dad1(uvec,ldaup,uvl,ldauk,recvec,usize,
     +              linds(lppar),hvbndry(bip),uinds)
            else
              call mpi_wait(req(3),mpistat,ierr)
              call extend_op1(uvec,ldaup,uvl,ldauk,
     +              linds(lppar),hvbndry(bip))
            end if
          end if
          if(iamkid.eq.1 .and. iampar.eq.0) then
            diffbits = ieor(pparbot,mymaskh)
            i = diffbits
            do j = 0,level2
              partner = ieor(partner,ishft(iand(i,1),j*2))
              i = ishft(i,-1)
            end do
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                     itype,comm,req(4),ierr)
            npending = 1
            if(iand(diffbits,ishft(1,level2)).ne.0) then 
              call mpi_isend(linds(lpkid),ldauk,MPI_INTEGER,partner,
     +                          itype,comm,req(5),ierr)
              npending = 2
            end if
            do while(npending.gt.0)
              call mpi_waitany(2,req(4),mid,mpistat,ierr)
              npending = npending-1
            end do
          end if
        else
          if(iamkid.eq.1 .and. iampar.eq.1) then
            partner = ieor(myid,ishft(1,dd-level-1))
            call mpi_isend(linds(lpkid),ldauk,MPI_INTEGER,partner,
     +                     itype,comm,req(4),ierr)
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                     itype,comm,req(5),ierr)
            call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                     itype,comm,req(1),ierr)
            call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                     itype,comm,req(2),ierr)
            call mpi_wait(req(1),mpistat,ierr)
            call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
            call mpi_wait(req(2),mpistat,ierr)
            call x10dad1(uvec,ldaup,uvl,ldauk,recvec,usize,
     +              linds(lppar),linds(lpkid),uinds)
            call mpi_wait(req(4),mpistat,ierr)
            call mpi_wait(req(5),mpistat,ierr)
          end if
          if(iamkid.eq.0 .and. iampar.eq.1) then
            diffbits = ieor(pkid,mymaskhk)
            partner = myid
            do j=0,level2
              partner = ieor(partner,ishft(iand(diffbits,1),j*2))
              diffbits = ishft(diffbits,-1)
            end do
            call mpi_irecv(uvl,msizedp,MPI_BYTE,partner,
     +                     itype,comm,req(3),ierr)

            partner = ieor(myid,ishft(1,dd-level-1))
            call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                     itype,comm,req(1),ierr)
            call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                     itype,comm,req(2),ierr)
            call mpi_isend(hvbndry(bip),ldauk,MPI_INTEGER,partner,
     +                     itype,comm,req(5),ierr)
            call mpi_wait(req(3),mpistat,ierr)
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                     itype,comm,req(4),ierr)
            call mpi_wait(req(1),mpistat,ierr)
            call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
            call mpi_wait(req(2),mpistat,ierr)
            call x10dad1(uvec,ldaup,uvl,ldauk,recvec,usize,
     +              linds(lppar),hvbndry(bip),uinds)
            call mpi_wait(req(4),mpistat,ierr)
            call mpi_wait(req(5),mpistat,ierr)
          end if
          if(iamkid.eq.1 .and. iampar.eq.0) then
            diffbits = ieor(pparbot,mymaskh)
            partner = myid
            do j=0,level2
              partner = ieor(partner,ishft(iand(diffbits,1),j*2))
              diffbits = ishft(diffbits,-1)
            end do
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                     itype,comm,req(3),ierr)
            call mpi_wait(req(3),mpistat,ierr)
          end if
        end if
        
        bvs     = hvbndry(hvbp+1)
        nb      = hvbndry(hvbp+2)
        istflag = hvbndry(hvbp+3)
        lbotsiz = hvbndry(hvbp+4)
        nhb     = hvbndry(hvbp+5)
        vsizer  = hvbndry(bvs-3)
        bip     = hvbp+6+3*nhb+1
        bip1    = bip
        nvinds  = hvbndry(bip-1)
        ntinds  = nvinds-vsizer
        bip     = bip+ntinds

        myleft  = lrud(nbrp)            
        myright = lrud(nbrp+1)      
        myup    = lrud(nbrp+2)      
        mydown  = lrud(nbrp+3)      
        uvlst   = 1                  
        ldauk   = vsizer            
        npending = 0 

        if(iampar.eq.0) then 
          do i=1,nvinds
              uvec(i) = zero
            end do
        else
          do i=1,nvinds
            recvec(i) = zero
          end do
        end if
        bhp = hvbp+6 
        if(nhb.ne.0) then
        if(mod(dd-level-1,2).eq.0) then
          i = bhp
          k1 = bip1
          k2 = 1
          do ih=1,nhb
            bhstrt = hvbndry(i)
            bhsize = hvbndry(i+1)
            do while(hvbndry(k1).ne.bhstrt .and. k1-bip1.le.ntinds)
              k1 = k1+1
            end do
            if(bhstrt.eq.hvbndry(k1)) then
            do j=0,bhsize-1
              uinds(k2+j) = hvbndry(k1+j)
            end do
            k2 = k2+bhsize
            end if
            i = i+3
          end do
          ntinds = k2-1
          call bgetrhs1(N,uinds,ntinds,rhs,rhsc)
        else
          call bgetrhs1(N,hvbndry(bip1),ntinds,rhs,rhsc)
        end if
        rhsst = 1
        end if

        do ih = 1,nhb
          bhstrt = hvbndry(bhp)        
          bhsize = hvbndry(bhp+1)      
          bhend  = hvbndry(bhp+2)      

          lvalst = lptrs(1,bhstrt)  
          ldalb  = lptrs(2,bhstrt)  

          itag = MPI_ANY_TAG            
          npendings2 = 0            

          uvecst = 1                  
          iv = 1
          bvp = hvbndry(hvbp+1)
          nvb = hvbndry(bvp-1)

          do while (iv.le.nvb) 
            if(hvbndry(bvp).lt.bhstrt .or. hvbndry(bvp+1).eq.0) then
              bvsize = hvbndry(bvp+1)
              uvecst = uvecst+bvsize 
              iv = iv+1
              bvp = bvp+3
            else
              goto 10
            end if
          end do

 10          if(iv.le.nvb .and. bhsize.ne.0) then

          do while (iv.le.nvb)
            bvstrt = hvbndry(bvp)      
            bvsize = hvbndry(bvp+1)      
            bvend  = hvbndry(bvp+2)      

            flags = 0
            if(ih+1 .gt. nhb) then
              flags = 1
            else
              if(hvbndry(bhp+3).gt.bvstrt) flags = 1
            end if
            flagr = 0
            if(flags .eq. 1) then
              indh = bvend
              if(indh.gt.suptop) indh = suptop
              indh = iand(ishft(indh,-lgblk),bmaskh)
              if(mod(mymaskh+hsize-1,hsize).ne.indh .and. 
     +           bhstrt.ne.supbot) flagr=1
              if(bhend.eq.suptop) flags = 0
            end if 

            if(bhstrt.eq.bvstrt) then 

              if(bhstrt.ne.supbot) then 
                call mpi_irecv(recvec,msizedp,MPI_BYTE,myleft,
     +                         itype,comm,req(1),ierr)
                call mpi_wait(req(1),mpistat,ierr)
                call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr)
                do i=0,bhsize-1
                  uvec(uvecst+i) = uvec(uvecst+i) + recvec(i+1)
                  rhs(rhsst+i) = rhs(rhsst+i) + uvec(uvecst+i)
                end do
              else
                do i=0,bhsize-1
                  rhs(rhsst+i) = rhs(rhsst+i) + uvec(uvecst+i)
                end do
              end if

              call dtrsv('l','n','n',bhsize,
     +                    lvals(lvalst),ldalb,rhs(rhsst),1)

              itag = myid

              if(mydown.ne.myid .and. .not.(bvend.eq.suptop 
     +              .and. supsiz.eq.lbotsiz)) then
                call mpi_isend(rhs(rhsst),bhsize*lendp,MPI_BYTE,mydown,
     +                         itag,comm,req(5),ierr)
                npendings2 = 1
              end if

              call putrhs1(N,linds(lptrs(3,bhstrt)),bhsize,rhs(rhsst),
     +                     rhsc)
              lvalst = lvalst+bhsize
              uvecst = uvecst+bhsize

            else 

              npendings = 0
              if(itag.eq.MPI_ANY_TAG) then
                call mpi_irecv(rhs(rhsst),msizedp,MPI_BYTE,myup,
     +                       itag,comm,req(4),ierr)
                npendings = 1
              end if

              if(flagr.eq.1) then
                if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                call mpi_irecv(recvec,msizedp,MPI_BYTE,myleft,
     +                        itype,comm,req(1),ierr)
              end if
                
              if(bvstrt.gt.suptop) then 
                recvsize = 0 
                do while (iv.le.nvb)
                  recvsize = recvsize+hvbndry(bvp+1)
                  bvp = bvp+3
                  iv = iv+1
                end do
              else
                recvsize = bvsize
              end if

              if(npendings.eq.1) then
                call mpi_wait(req(4),mpistat,ierr)
                itag = mpistat(MPI_TAG)
                if(itag.ne.mydown .and. .not.(bvend.eq.suptop .and.
     +                lbotsiz.eq.supsiz)) then
                  call mpi_isend(rhs(rhsst),bhsize*lendp,MPI_BYTE,
     +                  mydown,itag,comm,req(5),ierr)
                  npendings2 = 1
                end if
                npendings = 0
              end if

              call dgemv('n',recvsize,bhsize,ngone,lvals(lvalst),
     +                     ldalb,rhs(rhsst),1,one,uvec(uvecst),1)

              if(flagr.eq.1) then
                call mpi_wait(req(1),mpistat,ierr)
                if(flags.eq.1) then
                  do i=0,recvsize-1
                    uvec(uvecst+i) = uvec(uvecst+i)+recvec(i+1)
                  end do
                  call mpi_isend(uvec(uvecst),recvsize*lendp,MPI_BYTE,
     +                       myright,itype,comm,req(1),ierr)
                  npending = 1
                else
                  do i = 0,recvsize-1
                    uvl(uvlst+i) = uvec(uvecst+i)+recvec(i+1)
                  end do
                  uvlst = uvlst+recvsize
                end if
              else
                if(flags.eq.1) then
                  if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                  call mpi_isend(uvec(uvecst),recvsize*lendp,MPI_BYTE,
     +                       myright,itype,comm,req(1),ierr)
                  npending = 1
                else
                  if(bhend.eq.suptop) then
                    do i = 0,recvsize-1
                      uvl(uvlst+i) = uvec(uvecst+i)
                    end do
                    uvlst = uvlst+recvsize
                  end if
                end if
              end if
              uvecst = uvecst+recvsize
              lvalst = lvalst+recvsize

            end if

            bvp = bvp+3
            iv = iv+1
            do while (iv.le.nvb)
              if(hvbndry(bvp+1).eq.0) then
                bvp = bvp+3
                iv = iv+1
              else
                goto 30
              end if
            end do

 30       end do 

          else 

            if(iv.gt.nvb .and. bhsize.ne.0) then

              if(itag.eq.MPI_ANY_TAG .and. lbotsiz.ne.supsiz) then
                call mpi_irecv(rhs(rhsst),msizedp,MPI_BYTE,myup,
     +                         itag,comm,req(4),ierr)
                call mpi_wait(req(4),mpistat,ierr)
                itag = mpistat(MPI_TAG)
                if(itag.ne.mydown) then
                  call mpi_isend(rhs(rhsst),bhsize*lendp,MPI_BYTE,
     +                     mydown,itag,comm,req(5),ierr)
                  npendings2 = 1
                end if
              end if

            else
              do while (iv.le.nvb)

                bvstrt = hvbndry(bvp)      
                bvsize = hvbndry(bvp+1)      
                bvend  = hvbndry(bvp+2)      

                flags = 0
                if(ih+1 .gt. nhb) then
                  flags = 1
                else
                  if(hvbndry(bhp+3).gt.bvstrt) flags = 1
                end if
                flagr = 0
                if(flags .eq. 1) then
                  indh = bvend
                  if(indh.gt.suptop) indh = suptop
                  indh = iand(ishft(indh,-lgblk),bmaskh)
                  if(mod(mymaskh+hsize-1,hsize).ne.indh) flagr=1
                end if 

                if(flagr.eq.1) then
                  if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                  call mpi_irecv(recvec,msizedp,MPI_BYTE,myleft,
     +                         itype,comm,req(1),ierr)
                  call mpi_wait(req(1),mpistat,ierr)
                  call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr)
                  recvsize = ishft(nbrecv,-loglendp)
                  do i=0,recvsize-1
                    uvec(uvecst+i) = uvec(uvecst+i)+recvec(i+1)
                  end do
                  npending = 0
                  if(bvstrt.gt.suptop) then 
                    cmpsize = 0
                    do while(cmpsize.lt.recvsize)
                      cmpsize = cmpsize+hvbndry(bvp+1)
                      bvp = bvp+3
                      iv = iv+1
                    end do
                  end if
                else
                  if(bvstrt.gt.suptop) then 
                    recvsize = 0 
                    do while (iv.le.nvb)
                      recvsize = recvsize+hvbndry(bvp+1)
                      bvp = bvp+3
                      iv = iv+1
                    end do
                  else
                    recvsize = bvsize
                  end if
                end if

                if(flags.eq.1) then
                  call mpi_isend(uvec(uvecst),recvsize*lendp,MPI_BYTE,
     +                       myright,itype,comm,req(1),ierr)
                  npending = 1
                end if

                uvecst = uvecst+recvsize
                bvp = bvp+3
                iv = iv+1
                do while (iv.le.nvb)
                  if(hvbndry(bvp+1).eq.0) then
                    bvp = bvp+3
                    iv = iv+1
                  else
                    goto 40
                  end if
                end do

 40           end do

            end if

          end if 

          if(npendings2.eq.1) then
            call mpi_wait(req(5),mpistat,ierr)
            npendings2 = 0
          end if
          if(npending.eq.1) then
            call mpi_wait(req(1),mpistat,ierr)
            npending = 0
          end if

          rhsst  = rhsst+bhsize
          bhp = bhp+3
        end do

        call mpi_waitall(5,req,statall,ierr)

        lpkid = lptrs(3,suptop)
        if(istflag.eq.1) lpkid = lpkid+1

        hvbp = hvbndry(hvbp)
      end do 

      end

      subroutine x10dad1(pvec,psize,kvec,ksize,rvec,rsize,
     +      indsp,indsk,indsr)
      integer psize,ksize,rsize,indsp(*),indsk(*),indsr(*)
      double precision pvec(*),kvec(*),rvec(*),zero
      parameter(zero=0.d0)

        ik = 1
        ir = 1
        do ip = 1,psize
          do while(indsr(ir).lt.indsp(ip) .and. ir.le.rsize)
            ir = ir+1
          end do
          do while(indsk(ik).lt.indsp(ip) .and. ik.le.ksize)
            ik = ik+1
          end do
          if(ir.le.rsize .and. indsr(ir).eq.indsp(ip)) then
            pvec(ip) = rvec(ir)
          else
            pvec(ip) = zero
          end if
          if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then
            pvec(ip) = pvec(ip)+kvec(ik)
          end if
        end do
      end 

      subroutine extend_op1(pvec,psize,kvec,ksize,indsp,indsk)
      integer psize,ksize,indsp(*),indsk(*)
      double precision pvec(*),kvec(*)

        ik = 1
        do ip = 1,psize
          do while(indsk(ik).lt.indsp(ip) .and. ik.le.ksize)
            ik = ik+1
          end do
          if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then
            pvec(ip) = kvec(ik)
          else
            pvec(ip) = 0.d0
          end if
        end do

      end
