C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   pfsolvem.f                                                              */
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
C/* $Id: pfsolvem.f,v 1.2 2004-12-15 01:14:19 paklein Exp $ */
C/*****************************************************************************/

      subroutine pfsolvem(mysnodes,nsupnode,sup,
     +                  lptrs,linds,lvals,tptrs,tinds,
     +                  myid,myidh,dd,lgblk,N,
     +                  rhsc,nrhs,rhs,uvec,uvl,recvec,
     +                  uinds,maxvsize,lrud,
     +                  hvbndry,hvbsize,lc,w,iptrs,comm)

      implicit none
      include 'mpif.h'

      integer N,dd,lgblk,myid,myidh,maxvsize,nrhs
      integer nsupnode,hvbsize,comm

      integer lptrs(3,0:*), linds(*), iptrs(2,0:*)
      integer tptrs(3,0:*), tinds(*)
      integer sup(*), mysnodes(*)
      integer lrud(*), hvbndry(*)
      integer uinds(*), lc(*)

      double precision lvals(*)
      double precision rhsc(0:N-1,nrhs),w(*)
      double precision uvec(*),uvl(*)
      double precision rhs(*),recvec(*)

      integer lendp,loglendp,itype
      parameter(lendp=8,loglendp=3,itype=1)
      double precision one,ngone,zero
      parameter(one=1.d0,ngone=-1.d0,zero=0.d0)

      integer root
      integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk
      integer myleft,myright,myup,mydown
      integer supptr,suptop,supbot,supsiz,lbotsiz
      integer bvs,nb,nhb,nvb,nbrp,hvbp,ih,bhp
      integer bhstrt,bhsize,bhend,lvalst,ldalb
      integer uvecst,uvlst,itag
      integer iv,bvp,bvstrt,bvsize,bvend
      integer recvsize,cmpsize
      integer npending,npendings,npendings2
      integer mid,nbrecv
      integer lppar,kidtop,pparbot,kid
      integer diffbits,partner,ik,ldauk,lpkid,pkid
      integer i,j,ks,kd,is,usize
      integer mymaskh,mymaskhk,ldaup
      integer iamkid,iampar,kidr,kidl,pkido,bip,istflag
      integer vsizer,nbode,levelk,iaml,bnode,lvp,ivlim,nvinds
      integer brhsr,msizedp,vcsize,cbsize
      integer flagr,flags,indh,hsize,ldauvec,ks1
      integer mpistat(MPI_STATUS_SIZE),req(5),ierr
      integer statall(MPI_STATUS_SIZE,5)

      do i=1,5
        req(i) = MPI_REQUEST_NULL
      end do
      root   = mysnodes(1)

      suptop = mysnodes(nsupnode)
      call fsolve(N,lvals,linds,lptrs,tinds,tptrs,sup,rhsc,nrhs,
     +                suptop,lc,iptrs,w)

      supptr = tptrs(3,suptop)
      supbot = sup(supptr)
      bhsize = sup(supptr+1)
      ldalb  = lptrs(2,supbot)
      lpkid  = lptrs(3,supbot)

      ldauk  = ldalb-bhsize
      lpkid = lpkid+bhsize
                
      do j=0,nrhs-1
        kd = j*ldauk
        ks = j*ldalb+bhsize
        do i=1,ldauk
          uvl(kd+i) = w(ks+i)
        end do
      end do

      level = dd
      rowcol = 0
      bmaskh = 0
      bmaskv = 0
      mymaskh = 0
      hsize  = 1

      brhsr = nrhs*lendp
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
            hsize  = ishft(hsize,1)
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

        cbsize  = ldauk*brhsr
        msizedp = maxvsize*brhsr

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

              call x10dad(uvec,ldaup,uvl,ldauk,recvec,usize,nrhs,
     +              linds(lppar),linds(lpkid),uinds)
            else
              call extend_op(uvec,ldaup,uvl,ldauk,nrhs,
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
     +                  itype,comm,req(3),ierr)
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
              call x10dad(uvec,ldaup,uvl,ldauk,recvec,usize,nrhs,
     +              linds(lppar),hvbndry(bip),uinds)
            else
              call mpi_wait(req(3),mpistat,ierr)
              call extend_op(uvec,ldaup,uvl,ldauk,nrhs,
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
     +                  itype,comm,req(4),ierr)
            npending = 1
            if(iand(diffbits,ishft(1,level2)).ne.0) then 
              call mpi_isend(linds(lpkid),ldauk,MPI_INTEGER,
     +                  partner,itype,comm,req(5),ierr)
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
     +                  itype,comm,req(4),ierr)
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                  itype,comm,req(5),ierr)
            call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                  itype,comm,req(1),ierr)
            call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                  itype,comm,req(2),ierr)
            call mpi_wait(req(1),mpistat,ierr)
            call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
            call mpi_wait(req(2),mpistat,ierr)
            call x10dad(uvec,ldaup,uvl,ldauk,recvec,usize,nrhs,
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
     +                  itype,comm,req(3),ierr)

            partner = ieor(myid,ishft(1,dd-level-1))
            call mpi_irecv(uinds,maxvsize,MPI_INTEGER,partner,
     +                  itype,comm,req(1),ierr)
            call mpi_irecv(recvec,msizedp,MPI_BYTE,partner,
     +                  itype,comm,req(2),ierr)
            call mpi_isend(hvbndry(bip),ldauk,MPI_INTEGER,partner,
     +                  itype,comm,req(5),ierr)
            call mpi_wait(req(3),mpistat,ierr)
            call mpi_isend(uvl,cbsize,MPI_BYTE,partner,
     +                  itype,comm,req(4),ierr)
            call mpi_wait(req(1),mpistat,ierr)
            call mpi_get_count(mpistat,MPI_INTEGER,usize,ierr)
            call mpi_wait(req(2),mpistat,ierr)
            call x10dad(uvec,ldaup,uvl,ldauk,recvec,usize,nrhs,
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
     +                  itype,comm,req(3),ierr)
            call mpi_wait(req(3),mpistat,ierr)
          end if
        end if
        
        bvs     = hvbndry(hvbp+1)
        nb      = hvbndry(hvbp+2)
        istflag = hvbndry(hvbp+3)
        lbotsiz = hvbndry(hvbp+4)
        nhb     = hvbndry(hvbp+5)
        bip     = hvbp+6+3*nhb+1
        nvinds  = hvbndry(bip-1)
        vsizer  = hvbndry(bvs-3)
        bip     = bip+nvinds-vsizer

        myleft  = lrud(nbrp)
        myright = lrud(nbrp+1)
        myup    = lrud(nbrp+2)
        mydown  = lrud(nbrp+3)
        uvlst   = 1
        ldauk   = vsizer
        ldauvec = nvinds
        npending = 0 

        if(iampar.eq.0) then
          do i=1,nvinds*nrhs
            uvec(i) = zero
          end do
        else
          do i=1,nvinds*nrhs
            recvec(i) = zero
          end do
        end if

        bhp = hvbp+6 
        do ih = 1,nhb
          bhstrt = hvbndry(bhp)
          bhsize = hvbndry(bhp+1)
          bhend  = hvbndry(bhp+2)

          lvalst = lptrs(1,bhstrt)
          ldalb  = lptrs(2,bhstrt)

          itag = MPI_ANY_TAG
          npendings2 = 0
          vcsize = bhsize*brhsr

          uvecst = 1
          iv = 1
          bvp = hvbndry(hvbp+1)
          nvb = hvbndry(bvp-1)

          do while (iv.le.nvb) 
            if(hvbndry(bvp).lt.bhstrt .or. hvbndry(bvp+1).eq.0) then
              uvecst = uvecst+hvbndry(bvp+1)
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

              call bgetrhs(N,linds(lptrs(3,bhstrt)),bhsize,
     +                     bhsize,nrhs,rhs,rhsc)
              if(bhstrt.ne.supbot) then 
                call mpi_irecv(recvec,msizedp,MPI_BYTE,myleft,
     +                  itype,comm,req(1),ierr)
                call mpi_wait(req(1),mpistat,ierr)
                do j=0,nrhs-1
                  ks = j*bhsize+1
                  ks1 = uvecst+j*ldauvec
                  do i=0,bhsize-1
                    uvec(ks1+i) = uvec(ks1+i) + recvec(ks+i)
                    rhs(ks+i) = rhs(ks+i) + uvec(ks1+i)
                  end do
                end do
              else
                do j=0,nrhs-1
                  kd = j*bhsize+1
                  ks = uvecst+j*ldauvec
                  do i=0,bhsize-1
                    rhs(kd+i) = rhs(kd+i) + uvec(ks+i)
                  end do
                end do
              end if

              call dtrsm('l','l','n','n',bhsize,nrhs,one,
     +                    lvals(lvalst),ldalb,rhs,bhsize)

              itag = myid
              if(mydown.ne.myid .and. .not.(bvend.eq.suptop 
     +              .and. supsiz.eq.lbotsiz)) then
                call mpi_isend(rhs,vcsize,MPI_BYTE,mydown,
     +                  itag,comm,req(5),ierr)
                npendings2 = 1
              end if
              call putrhs(N,linds(lptrs(3,bhstrt)),bhsize,
     +                      bhsize,nrhs,rhs,rhsc)
              uvecst = uvecst+bhsize
              lvalst = lvalst+bhsize

            else 

              npendings = 0
              if(itag.eq.MPI_ANY_TAG) then
                call mpi_irecv(rhs,msizedp,MPI_BYTE,myup,
     +                  itag,comm,req(4),ierr)
                npendings = 1
              end if

              if(flagr.eq.1) then
                if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                call mpi_irecv(recvec,msizedp,MPI_BYTE,myleft,
     +                  itype,comm,req(1),ierr)
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
                  call mpi_isend(rhs,vcsize,MPI_BYTE,mydown,
     +                        itag,comm,req(5),ierr)
                  npendings2 = 1
                end if
                npendings = 0
              end if

              call dgemm('n','n',recvsize,nrhs,bhsize,ngone,
     +                   lvals(lvalst),ldalb,rhs,bhsize,one,
     +                   uvec(uvecst),ldauvec)

              if(flagr.eq.1) then
                call mpi_wait(req(1),mpistat,ierr)
                if(flags.eq.1) then
                  do j=0,nrhs-1
                    kd = uvecst+j*ldauvec
                    ks = j*recvsize+1
                    do i=0,recvsize-1
                      recvec(ks+i) = recvec(ks+i) + uvec(kd+i)
                    end do
                  end do
                  call mpi_isend(recvec,recvsize*brhsr,MPI_BYTE,myright,
     +                         itype,comm,req(1),ierr)
                  npending = 1
                else
                  do j=0,nrhs-1
                    kd = j*ldauk+uvlst
                    ks = j*ldauvec+uvecst
                    ks1 = j*recvsize+1
                    do i = 0,recvsize-1
                      uvl(kd+i) = uvec(ks+i) + recvec(ks1+i)
                    end do
                  end do
                  uvlst = uvlst+recvsize
                end if
              else
                if(flags.eq.1) then
                  if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                  do j=0,nrhs-1
                    kd = 1+recvsize*j
                    ks = uvecst+j*ldauvec
                    do i=0,recvsize-1
                      recvec(kd+i) = uvec(ks+i)
                    end do
                  end do

                  call mpi_isend(recvec,recvsize*brhsr,MPI_BYTE,myright,
     +                         itype,comm,req(1),ierr)
                  npending = 1
                else
                  if(bhend.eq.suptop) then
                    do j=0,nrhs-1
                      kd = j*ldauk+uvlst
                      ks = j*ldauvec+uvecst
                      do i = 0,recvsize-1
                        uvl(kd+i) = uvec(ks+i)
                      end do
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
                call mpi_irecv(rhs,msizedp,MPI_BYTE,myup,
     +                        itag,comm,req(4),ierr)
                call mpi_wait(req(4),mpistat,ierr)
                itag = mpistat(MPI_TAG)
                call mpi_get_count(mpistat,MPI_BYTE,mpistat,ierr)
                itag = mpistat(MPI_TAG)
                call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr)
                if(itag.ne.mydown) then
                  call mpi_isend(rhs,nbrecv,MPI_BYTE,mydown,
     +                          itag,comm,req(5),ierr)
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
     +                    itype,comm,req(1),ierr)
                  call mpi_wait(req(1),mpistat,ierr)
                  call mpi_get_count(mpistat,MPI_BYTE,nbrecv,ierr)
                  recvsize = ishft(nbrecv,-loglendp)/nrhs
                  do j=0,nrhs-1
                    kd = uvecst+j*ldauvec
                    ks = j*recvsize+1
                    do i=0,recvsize-1
                      uvec(kd+i) = uvec(kd+i) + recvec(ks+i)
                    end do
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
                  if(npending.eq.1) call mpi_wait(req(1),mpistat,ierr)
                  do j=0,nrhs-1
                    kd = 1+recvsize*j
                    ks = uvecst+j*ldauvec
                    do i=0,recvsize-1
                      recvec(kd+i) = uvec(ks+i)
                    end do
                  end do

                  call mpi_isend(recvec,recvsize*brhsr,MPI_BYTE,myright,
     +                  itype,comm,req(1),ierr)
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

          bhp = bhp+3
        end do

        call mpi_waitall(5,req,statall,ierr)

        lpkid = lptrs(3,suptop)
        if(istflag.eq.1) lpkid = lpkid+1

        hvbp = hvbndry(hvbp)

      end do 

      end

      subroutine extend_op(pvec,psize,kvec,ksize,nrhs,indsp,indsk)
      integer psize,ksize,indsp(*),indsk(*),nrhs
      double precision pvec(*),kvec(*)

        ik = 1
        do ip = 1,psize
          do while(indsk(ik).lt.indsp(ip) .and. ik.le.ksize)
            ik = ik+1
          end do
          if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then
            do k=0,nrhs-1
              pvec(psize*k+ip) = kvec(ksize*k+ik)
            end do
          else
            do k=0,nrhs-1
              pvec(psize*k+ip) = 0.d0
            end do
          end if
        end do

      end

c@process opt(3) strict debug(inline)
      subroutine x10dad(pvec,psize,kvec,ksize,rvec,rsize,nrhs,
     +      indsp,indsk,indsr)
      integer psize,ksize,rsize,nrhs,indsp(*),indsk(*),indsr(*)
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
            do k=0,nrhs-1
              pvec(psize*k+ip) = rvec(rsize*k+ir)
            end do
          else
            do k=0,nrhs-1
              pvec(psize*k+ip) = zero
            end do
          end if
          if(ik.le.ksize .and. indsk(ik).eq.indsp(ip)) then
            do k=0,nrhs-1
              idst = psize*k+ip
              pvec(idst) = pvec(idst)+kvec(ksize*k+ik)
            end do
          end if
        end do
      end

c@process nosave
C     recursive 
      subroutine comp_sty(root,sanity,tptrs,tinds,rhsc,N,nrhs)
      implicit none

      integer N,nrhs,i,j,kid,root,size
      integer tptrs(3,0:N-1),tinds(N)
      double precision rhsc(0:N-1,nrhs),sanity

      size = tptrs(2,root)
      do i = 0,size-1
        kid = tinds(tptrs(1,root)+i)
        call comp_sty_recursive(kid,sanity,tptrs,tinds,rhsc,N,nrhs)
        do j=1,nrhs
         sanity = sanity + rhsc(kid,j)
        end do
      end do
      end
