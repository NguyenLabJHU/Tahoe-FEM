C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   pbsolve1.f                                                              */
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
C/* $Id: pbsolve1.f,v 1.1 2004-12-10 20:28:27 paklein Exp $ */
C/*****************************************************************************/

      subroutine pbsolve1(mysnodes,nsupnode,sup,
     +                  lptrs,linds,lvals,tptrs,tinds,
     +                  myid,myidh,myidv,dd,lgblk,N,
     +                  rhsc,uvec,recvec,rhs,uinds,maxvsize,
     +                  lrud,hvbndry,hvbsize,lc,w,iptrs,comm)

      implicit none
      include 'mpif.h'

      integer N,dd,lgblk,myid,myidh,myidv,maxvsize
      integer nsupnode,hvbsize,comm

      integer lptrs(3,0:*), linds(*), iptrs(2,0:*)
      integer tptrs(3,0:*), tinds(*)
      integer sup(*), mysnodes(*)
      integer lrud(*), hvbndry(*)
      integer uinds(*), lc(*)

      double precision lvals(*), w(*)
      double precision rhsc(0:*),uvec(*),rhs(*)
      double precision recvec(*)

      integer lendp,loglendp,itype
      parameter(lendp=8,loglendp=3,itype=1)
      double precision one,ngone,zero
      parameter(zero=0.d0,one=1.d0,ngone=-1.d0)

      integer level,level2,rowcol,bmaskh,bmaskv,bmaskhk,bmaskvk
      integer myleft,myright,myup,mydown
      integer supptr,supbot
      integer nhb,nvb,nbrp,hvbp,ih,iv,ivlim,bhp,bvp,bvs
      integer bhstrt,bhsize,bvstrt,bvsize
      integer uvecst,i,j,is,nbrecv,k1,k2
      integer lhsize,lvsize,hsize,vsize,vfac
      integer psrcv,mymaskv,bnode

      integer vsizek,nuinds,partner,lvalst,ldalb,itag
      integer ir,ip,ik,indk,krhsst,kbipi,mid,npending,npendings
      integer ksupbot,klbotsiz,kbip,kvsizer,kbvs
      integer trhs,ksupptr,khvbp,kid,knvinds,nvbr,nvbt
      integer istflag,rhsst,lbotsiz,suptop,nvinds
      integer vsizer,bip,tstf,ksupsiz,supsiz,szflag,rrec
      integer mpistat(MPI_STATUS_SIZE),req(2),ierr

      req(1) = MPI_REQUEST_NULL
      req(2) = MPI_REQUEST_NULL

      level  = 0
      rowcol = mod(dd,2)
      bmaskh = (dd+1)/2
      bmaskv = dd-bmaskh
      hsize  = ishft(1,bmaskh)
      vsize  = ishft(1,bmaskv)
      lhsize = bmaskh
      lvsize = bmaskv
      bmaskh = hsize-1
      bmaskv = vsize-1

      nbrp = 1+(dd-1)*4
      khvbp = hvbsize+1

      kid      = mysnodes(1)

      ksupptr  = tptrs(3,kid)
      ksupbot  = sup(ksupptr)
      ksupsiz  = sup(ksupptr+1)
      khvbp    = hvbndry(khvbp-1)
      kbvs     = hvbndry(khvbp+1)
      klbotsiz = hvbndry(khvbp+4)
      kbip     = kbvs-3-klbotsiz
      kvsizer  = hvbndry(kbvs-3)
      knvinds  = hvbndry(kbip-1)
      kbipi    = kbip+knvinds-kvsizer

      do is=1,nsupnode-1 

        hvbp = khvbp

        myleft = lrud(nbrp) 
        myright= lrud(nbrp+1) 
        myup   = lrud(nbrp+2) 
        mydown = lrud(nbrp+3) 

        suptop  = kid
        supbot  = ksupbot
        lbotsiz = klbotsiz
        supsiz  = ksupsiz
        bvs     = kbvs
        vsizer  = kvsizer
        nvinds  = knvinds
        bip     = kbipi
        bnode   = ishft(suptop,-lgblk)
        mymaskv = iand(myidv,bmaskv)

        istflag = hvbndry(hvbp+3)
        nhb     = hvbndry(hvbp+5)
        nvb     = hvbndry(bvs-1)
        nvbr    = hvbndry(bvs-2)
        bhp     = hvbp+6+(nhb-1)*3      
        nvbt    = nvb-nvbr
        bvp     = bvs+3*(nvbt-1) 
        psrcv   = iand(bnode,bmaskv)
        vfac    = hsize/vsize
        iv      = nvbt
        szflag  = supsiz-lbotsiz

        if((nvb.ne.0).and.(nhb.ne.0)) then
        

        do ih = nhb,1,-1
          bhstrt = hvbndry(bhp)
          bhsize = hvbndry(bhp+1)

          rhsst = nvinds-vsizer+1
          tstf  = istflag
          ivlim = max(iv-vfac+1,1)

          if(bhsize.ne.0) then

          ldalb  = lptrs(2,bhstrt)
          lvalst = lptrs(1,bhstrt)+ldalb-vsizer

          do i=1,bhsize
            recvec(i) = zero
          end do

          if(vsizer.ne.0) then
            
            call dgemv('t',vsizer,bhsize,one,lvals(lvalst),
     +        ldalb,rhs(rhsst),1,one,recvec,1)
          end if

          if (szflag.ne.0) then
            call alltoonev(recvec,uvec,bhsize,psrcv,lvsize,
     +          mymaskv,myid,comm)
          else
            szflag = 1
          end if
        
          do iv = iv,ivlim,-1
            bvstrt = hvbndry(bvp)
            bvsize = hvbndry(bvp+1)

            lvalst = lvalst-bvsize
            rhsst  = rhsst-bvsize
            bip    = bip-bvsize
            itag   = MPI_ANY_TAG

            if(bvsize.ne.0) then

            if(tstf.eq.1 .or. tstf.eq.2) then
            
              if(bhstrt.eq.bvstrt) then 

              call bgetrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc)
              do i=0,bvsize-1
                rhs(rhsst+i) = rhs(rhsst+i) - recvec(i+1)
              end do
              call dtrsv('l','t','n',bvsize,lvals(lvalst),ldalb,
     +          rhs(rhsst),1)
              
              
              itag = myid
              call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,myleft,
     +                       itag,comm,req(1),ierr)
              call putrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc)
              call mpi_wait(req(1),mpistat,ierr)

              else
              
                npendings = 0
                call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,myright,
     +            itag,comm,req(1),ierr)
                call mpi_wait(req(1),mpistat,ierr)
                itag = mpistat(MPI_TAG)
                if(itag.ne.myleft) then
                  call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,
     +              myleft,itag,comm,req(1),ierr)
                  npendings = 1
                end if
                call dgemv('t',bvsize,bhsize,one,lvals(lvalst),
     +            ldalb,rhs(rhsst),1,one,recvec,1)
                if (myup.ne.myid) then
                  call mpi_isend(recvec,bhsize*lendp,MPI_BYTE,myup,
     +                     itype,comm,req(2),ierr)
                           npendings = npendings + 1
                end if
                do while (npendings.gt.0)
                  call mpi_waitany(2,req,mid,mpistat,ierr)
                  npendings = npendings-1
                end do

              end if
              tstf = 0

            else  

              if(bhstrt.eq.bvstrt) then 

                if(mydown.ne.myid) then
                  call mpi_irecv(recvec,bhsize*lendp,MPI_BYTE,mydown,
     +                     itype,comm,req(2),ierr)
                end if
                call bgetrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc)
                if(mydown.ne.myid) then
                  call mpi_wait(req(2),mpistat,ierr)
                end if
                do i=0,bvsize-1
                  rhs(rhsst+i) = rhs(rhsst+i) - recvec(i+1)
                end do
                call dtrsv('l','t','n',bvsize,lvals(lvalst),ldalb,
     +            rhs(rhsst),1)
                itag = myid
                
                call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,myleft,
     +                         itag,comm,req(1),ierr)
                call putrhs1(N,hvbndry(bip),bvsize,rhs(rhsst),rhsc)
                call mpi_wait(req(1),mpistat,ierr)

              else   

                if(bvstrt.gt.bhstrt) then

                  npendings = 0
                  rrec = 0
                  call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,
     +                      myright,itag,comm,req(1),ierr)
                  npending = 1
                  if(mydown.ne.myid) then
                    call mpi_irecv(recvec,bhsize*lendp,MPI_BYTE,
     +                 mydown,itype,comm,req(2),ierr)
                    npending = 2
                  end if
                  do while(npending.gt.0)
                    call mpi_waitany(2,req,mid,mpistat,ierr)
                    if(mid.eq.1 .and. rrec.eq.0) then
                      itag = mpistat(MPI_TAG)
                      if(itag.ne.myleft) then
                        call mpi_isend(rhs(rhsst),bvsize*lendp,
     +                          MPI_BYTE,myleft,itag,comm,
     +                     req(1),ierr)
                        npending = npending+1
                      end if
                      rrec = 1
                    end if
                    npending = npending-1
                  end do
                  call dgemv('t',bvsize,bhsize,one,lvals(lvalst),
     +              ldalb,rhs(rhsst),1,one,recvec,1)
                  if(myup.ne.myid) then
                         call mpi_isend(recvec,bhsize*lendp,MPI_BYTE,
     +                  myup,itype,comm,req(2),ierr)
                    call mpi_wait(req(2),mpistat,ierr)
                  end if

                else 
                
                  call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,
     +               myright,itag,comm,req(1),ierr)
                  call mpi_wait(req(1),mpistat,ierr)
                  itag = mpistat(MPI_TAG)
                  if(itag.ne.myleft) then
                    call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,
     +                  myleft,itag,comm,req(1),ierr)
                    call mpi_wait(req(1),mpistat,ierr)
                  end if

                end if 

              end if 

            end if 

            else 
            
              if(myup.ne.myid .and. bhstrt.lt.bvstrt) then
                if(tstf.eq.2) then
                  tstf = 0
                else
                  call mpi_irecv(recvec,bhsize*lendp,MPI_BYTE,mydown,
     +                     itype,comm,req(2),ierr)
                  call mpi_wait(req(2),mpistat,ierr)
                end if

                call mpi_isend(recvec,bhsize*lendp,MPI_BYTE,myup,
     +                         itype,comm,req(2),ierr)
                call mpi_wait(req(2),mpistat,ierr)

              end if

            end if

            vsizer = vsizer+bvsize
            bvp = bvp-3
          end do  

          else

            do iv = iv,ivlim,-1
              bvsize = hvbndry(bvp+1)

              lvalst = lvalst-bvsize
              rhsst  = rhsst-bvsize
              bip    = bip-bvsize
              itag   = MPI_ANY_TAG

              if(bvsize.ne.0) then
                call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,myright,
     +                         itag,comm,req(1),ierr)
                call mpi_wait(req(1),mpistat,ierr)
                itag = mpistat(MPI_TAG)
                if(itag.ne.myleft) then
                  call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,
     +              myleft,itag,comm,req(1),ierr)
                  call mpi_wait(req(1),mpistat,ierr)
                end if
              end if

              vsizer = vsizer+bvsize
              bvp = bvp-3
            end do
          end if

          bhp = bhp-3
        end do 

        do iv = ivlim-1,1,-1
          bvsize = hvbndry(bvp+1)
          if (bvsize.ne.0) then
            rhsst  = rhsst-bvsize
            bip    = bip-bvsize
            itag   = MPI_ANY_TAG
            call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,myright,
     +                       itag,comm,req(1),ierr)
            call mpi_wait(req(1),mpistat,ierr)
            itag = mpistat(MPI_TAG)
            if(itag.ne.myleft) then
              call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,myleft,
     +                       itag,comm,req(1),ierr)
              call mpi_wait(req(1),mpistat,ierr)
            end if
          end if
          bvp = bvp-3
        end do

        else if((nvb.eq.0).and.(nhb.ne.0)) then
         
          if(szflag.ne.0) then 
            bhsize = hvbndry(bhp+1)
            if(bhsize.ne.0) then
            do i=1,bhsize
              uvec(i) = zero
            end do
            call alltoonev(uvec,recvec,bhsize,psrcv,lvsize,
     +        mymaskv,myid,comm)
            end if
            bhp = bhp-3
          end if

          do ih=nhb-1,1,-1
            bhsize = hvbndry(bhp+1)
            if(bhsize.ne.0) then
            do i=1,bhsize
              uvec(i) = zero
            end do
            call alltoonev(uvec,recvec,bhsize,psrcv,lvsize,
     +        mymaskv,myid,comm)
            end if
            bhp = bhp-3
          end do

        else if((nhb.eq.0).and.(nvb.ne.nvbr)) then
        
          rhsst = nvinds-vsizer+1
          do iv = nvbt,1,-1
            bvsize = hvbndry(bvp+1)
            if(bvsize.ne.0) then
              rhsst  = rhsst-bvsize
              bip    = bip-bvsize
              itag   = MPI_ANY_TAG
              call mpi_irecv(rhs(rhsst),bvsize*lendp,MPI_BYTE,myright,
     +                           itag,comm,req(1),ierr)
              call mpi_wait(req(1),mpistat,ierr)
              itag = mpistat(MPI_TAG)
              if(itag.ne.myleft) then
                call mpi_isend(rhs(rhsst),bvsize*lendp,MPI_BYTE,myleft,
     +                              itag,comm,req(1),ierr)
                call mpi_wait(req(1),mpistat,ierr)
              end if
            end if
            bvp = bvp-3
          end do

        end if

        vsizek = vsize

        if(tptrs(2,supbot).ne.1) then 
          level2 = ishft(dd-level-1,-1)
          if(rowcol.eq.1) then
            bmaskh = ieor(bmaskh,ishft(1,level2))
            hsize  = ishft(hsize,-1)
            lhsize = lhsize-1
          else
            bmaskv = ieor(bmaskv,ishft(1,level2))
            vsize  = ishft(vsize,-1)
            lvsize = lvsize-1
          end if
          level = level+1
          rowcol = 1-rowcol
          nbrp = nbrp-4  
        end if

        bip = kbip
        kid = mysnodes(is+1)

        if(is.ne.nsupnode-1) then

          khvbp    = hvbndry(khvbp-1)
          ksupptr  = tptrs(3,kid)
          ksupbot  = sup(ksupptr)
          ksupsiz  = sup(ksupptr+1)
          kbvs     = hvbndry(khvbp+1)
          klbotsiz = hvbndry(khvbp+4)
          kbip     = kbvs-3-klbotsiz
          kvsizer  = hvbndry(kbvs-3)
          knvinds  = hvbndry(kbip-1)
          kbipi    = kbip+knvinds-kvsizer

          if(vsize.ne.vsizek) then  

            partner = ieor(myid,ishft(1,dd-level))
            if(partner.gt.myid) then
              call mpi_send(hvbndry(bip),nvinds,MPI_INTEGER,partner,
     +                      itype,comm,ierr)
              call mpi_recv(uinds,maxvsize,MPI_INTEGER,partner,
     +                      itype,comm,mpistat,ierr)
            else
              call mpi_recv(uinds,maxvsize,MPI_INTEGER,partner,
     +                      itype,comm,mpistat,ierr)
              call mpi_send(hvbndry(bip),nvinds,MPI_INTEGER,partner,
     +                      itype,comm,ierr)
            end if
            call mpi_get_count(mpistat,MPI_INTEGER,nuinds,ierr)
            if(partner.gt.myid) then
              call mpi_send(rhs,nvinds*lendp,MPI_BYTE,partner,
     +                      itype,comm,ierr)
              call mpi_recv(recvec,maxvsize*lendp,MPI_BYTE,partner,
     +                      itype,comm,mpistat,ierr)
            else
              call mpi_recv(recvec,maxvsize*lendp,MPI_BYTE,partner,
     +                      itype,comm,mpistat,ierr)
              call mpi_send(rhs,nvinds*lendp,MPI_BYTE,partner,
     +                      itype,comm,ierr)
            end if

            ir = 1
            ip = 1
            do ik=1,kvsizer
              indk = hvbndry(kbipi+ik-1)
              do while(uinds(ir).lt.indk .and. ir.le.nuinds)
                ir = ir+1
              end do
              do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds)
                ip = ip+1
              end do
              if(ip.le.nvinds .and. indk.eq.hvbndry(bip+ip-1)) then
              
                uvec(ik) = rhs(ip)
              end if
              if(ir.le.nuinds .and. indk.eq.uinds(ir)) then
              
                uvec(ik) = recvec(ir)
              end if
            end do

          else 

            ip = 1
            do ik=1,kvsizer
              indk = hvbndry(kbipi+ik-1)
              do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds)
                ip = ip+1
              end do
              if(indk.eq.hvbndry(bip+ip-1)) then
                uvec(ik) = rhs(ip)
              end if
            end do

          end if

          krhsst = knvinds-kvsizer
          do i=1,krhsst
            rhs(i) = zero
          end do
          do i=1,kvsizer
            rhs(krhsst+i) = uvec(i)
          end do

        else 

          i = tptrs(3,kid)
          knvinds = lptrs(2,sup(i))
          kvsizer = lptrs(2,kid)-1
          kbipi   = lptrs(3,kid)+1
          krhsst = knvinds-kvsizer

          ip = 1
          do ik=1,kvsizer
            indk = linds(kbipi+ik-1)
            do while(hvbndry(bip+ip-1).lt.indk .and. ip.le.nvinds)
              ip = ip+1
            end do

            if(ip.le.nvinds .and. indk.eq.hvbndry(bip+ip-1)) then
              w(krhsst+ik) = rhs(ip)
            end if
          end do

        end if

      end do  

      call bsolve(N,lvals,linds,lptrs,tinds,tptrs,sup,rhsc,
     +            1,mysnodes(nsupnode),lc,iptrs,w)

      end
