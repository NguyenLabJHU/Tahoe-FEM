C/*****************************************************************************/
C/*                                                                           */
C/*   (C) Copyright IBM Corporation, 1997                                     */
C/*   (C) Copyright Regents of the University of Minnesota, 1997              */
C/*                                                                           */
C/*   parsymb.f                                                               */
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
C/* $Id: parsymb.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
C/*****************************************************************************/

      subroutine parsymb(root,aptrs,ainds,tptrs,tinds,lptrs,
     +                   linds,lindsptr,lvalsptr,nnz,lbal,sizes,
     +                   sup,supinds,supptr,supindsptr,iu,lsize,
     +                   dd,N,myid,lgblk,myidr,myidc,mystak,
     +                   resdcol,info,cinfo,comm)
      
      implicit none

      include 'mpif.h'

      integer itag
      parameter(itag=1)
      double precision SSTDEF
      parameter(SSTDEF=2.0d2)

      integer lsize,info
      integer node,level,rowcol,dd,N,myid,lgblk,root
      integer supptr,supindsptr,iu,lindsptr,lvalsptr
      integer myidr,myidc,bmaskr,bmaskc,mymaskr,mymaskc

      integer stakptr,kidnode,scptr,comm
      integer ns,is1,is4,nnz,nnzxtra,nnzl
      double precision opcount,opcountl,opxtra,ssthresh,lbal(0:*)
      double precision osnode
      integer i,j,k,l,ipartner,nbrecv,level2,leafroot
      integer lbotsiz , resdcolsiz, recvcolsiz
      integer stakst,supbot,suptop,supsiz

      integer aptrs(2,0:*),ainds(*)
      integer tptrs(3,0:*),tinds(*)
      integer mystak(*) , resdcol(*)
      integer lptrs(3,0:*),linds(*)
      integer sup(*),supinds(*)
      integer cinfo(0:*)
      logical ihadlast
      integer firstindex
      integer nuptopl

      integer, allocatable :: temp1(:),temp2(:),sinds(:),sptrs(:)
      double precision, allocatable :: temp3(:),temp4(:),svals(:)
      integer sizes(0:*),nrecv,maxisizel,maxisizeg,sisize,isize
      integer tsize,siptr,n1,n2,m,ti

      integer mpistat(MPI_STATUS_SIZE),ierr

      node = root
      j=1
      do i=1,dd
        k = j
        do while(tptrs(2,node).eq.1)
          mystak(j) = node
          node = tinds(tptrs(1,node))
          j = j+1
        end do
        mystak(j) = node
        j = j+1
        mystak(j) = j-k
        j = j+1
        if(tptrs(2,node).gt.2) then
          print *,'error - tree not binarized at level',i,
     +            ' for node ',node
          call mpi_abort(comm,0,ierr)
        end if
        if(iand(myid,ishft(1,dd-i)).eq.0) then 
          node = tinds(tptrs(1,node)) 
        else 
          node = tinds(tptrs(1,node)+1) 
        end if
      end do
      stakptr = j-1 

      do i=0,N-1
        lptrs(2,i) = 0
      end do

      lvalsptr = 1
      nnzxtra = 0
      nnzl = 1
      opxtra = 1.d0
      opcountl = 0.d0
      supptr = 1
      lindsptr = 1
      ssthresh = SSTDEF

      info = 0
      leafroot = node
      call symbolic6(aptrs,ainds,lptrs,linds,sup,resdcol,
     +               tptrs,tinds,nnzl,node,supptr,lvalsptr,
     +               opcountl,mystak(j),lindsptr,nnzxtra,opxtra,
     +               ssthresh,lsize,info)

      if(info.eq.1) return

      j = sup(supptr-4)  
      iu = ishft(lptrs(2,j),1)
      supindsptr = 1

      lbal(myid) = opcountl

      if (dd.ne.0) then 

      rowcol = 0 
      bmaskc = 0
      bmaskr = 0
      mymaskc = 0
      mymaskr = 0
      ipartner = ieor(myid,1)

      scptr = lptrs(3,node)+1
      lbotsiz = lptrs(2,node)-1

      call mpi_sendrecv(linds(scptr),lbotsiz,MPI_INTEGER,ipartner,
     +                 itag,supinds(supindsptr),N,MPI_INTEGER, 
     +                 ipartner,itag,comm, mpistat,ierr)
      call mpi_get_count(mpistat,MPI_INTEGER,recvcolsiz,ierr)

      call mergelists(linds(scptr),lbotsiz,
     +                  supinds(supindsptr),recvcolsiz,
     +                  resdcol,resdcolsiz)

      nuptopl = ishft(1,dd)

      do level=dd-1,0,-1

        stakst = stakptr
        supsiz = mystak(stakptr)

        if(mod(supsiz,2).eq.0) then 
          do stakptr = stakst-1,stakst-supsiz,-2
            kidnode = mystak(stakptr)
            call mergelists(resdcol,resdcolsiz,
     +                      ainds(aptrs(1,kidnode)),aptrs(2,kidnode),
     +                      linds(lindsptr), lbotsiz)
            kidnode = mystak(stakptr-1)
            call mergelists(linds(lindsptr),lbotsiz,
     +                      ainds(aptrs(1,kidnode)),aptrs(2,kidnode),
     +                      resdcol, resdcolsiz)
          end do
          call all_to_all_union_hc(resdcol,supinds(supindsptr),
     +               linds(lindsptr),N,resdcolsiz,lbotsiz,myid,dd-level,
     +               comm)
        else 
          do stakptr = stakst-1,stakst-supsiz+1,-2
            kidnode = mystak(stakptr)
            call mergelists(resdcol,resdcolsiz,
     +                      ainds(aptrs(1,kidnode)),aptrs(2,kidnode),
     +                      linds(lindsptr), lbotsiz)
            kidnode = mystak(stakptr-1)
            call mergelists(linds(lindsptr),lbotsiz,
     +                      ainds(aptrs(1,kidnode)),aptrs(2,kidnode),
     +                      resdcol, resdcolsiz)
          end do
          kidnode = mystak(stakptr)
          stakptr = stakptr-1
          call mergelists(resdcol,resdcolsiz,
     +                    ainds(aptrs(1,kidnode)),aptrs(2,kidnode),
     +                    linds(lindsptr), lbotsiz)
          call all_to_all_union_hc(linds(lindsptr),supinds(supindsptr),
     +               resdcol,N,lbotsiz,lbotsiz,myid,dd-level,
     +               comm)
        end if

        rowcol = 1-rowcol 
        level2 = ishft(dd-level-1,-1)
        if(rowcol.eq.1) then 
          bmaskc = ior(bmaskc,ishft(1,level2))
        else
          bmaskr = ior(bmaskr,ishft(1,level2))
        end if

        l = supindsptr+lbotsiz
        k=1
        do i=0,lbotsiz-1
          j = supinds(supindsptr+i)
          if(iand(myidr,bmaskr).eq.iand(ishft(j,-lgblk),bmaskr)) then
            supinds(l+k-1) = j
            k = k+1
          end if
        end do
        recvcolsiz = k-1

        if(iand(myidr,bmaskr).eq.0 .and. iand(myidc,bmaskc).eq.0) then
          osnode = 0.d0
          j = lbotsiz
          do i=0,supsiz-1
            osnode = osnode + dble(j*j)
            j = j-1
          end do
          j = ishft(myid,level-dd)+nuptopl
          nuptopl = nuptopl+ishft(1,level)
          lbal(j) = osnode
        end if

        k=1
        j=0
        ihadlast = .false.
        do i = stakst-1,stakst-supsiz,-1
          kidnode = mystak(i)
          if(iand(myidc,bmaskc) .eq.
     +       iand(ishft(kidnode,-lgblk),bmaskc)) then
            do while((supinds(l+k-1).lt.kidnode).and.
     +                   (k.le.recvcolsiz))
              k = k+1
            end do
            if(j.eq.0) j=k
            lptrs(2,kidnode) = recvcolsiz-k+1
            lptrs(3,kidnode) = lindsptr+k-j
            cinfo(kidnode) = 1 
            if(.not.ihadlast) then
              lptrs(1,kidnode) = lvalsptr
              lvalsptr = lvalsptr+lptrs(2,kidnode)
              ihadlast = .true.
              firstindex = k
            else
              lvalsptr = lvalsptr+k-firstindex
              lptrs(1,kidnode) = lvalsptr
              lvalsptr = lvalsptr+lptrs(2,kidnode)
            end if
          else
            lptrs(2,kidnode) = 0
            lptrs(3,kidnode) = lindsptr
            cinfo(kidnode) = 0 
            ihadlast = .false.
          end if
          nnzl = nnzl+lptrs(2,kidnode)
        end do

        if(j.ne.0) then
          k = recvcolsiz-j+1
          do i=0,k-1
            linds(lindsptr+i) = supinds(l+j-1+i)
          end do
          lindsptr = lindsptr+k
        end if

        supbot = mystak(stakst-1)
        suptop = mystak(stakst-supsiz)
        tptrs(3,suptop) = supptr
        tptrs(3,supbot) = supptr

        sup(supptr) = supbot           
        sup(supptr+1) = supsiz         
        sup(supptr+2) = supindsptr         
        sup(supptr+3) = lbotsiz         
        supptr = supptr + 4

        scptr = supindsptr+supsiz
        supindsptr = l+lbotsiz
        lbotsiz = lbotsiz-supsiz

        if(level.ne.0) then
          ipartner = ieor(myid,ishft(1,dd-level))

          call mpi_sendrecv(supinds(scptr),lbotsiz,MPI_INTEGER,
     +                    ipartner,itag,
     +                    supinds(supindsptr),N,MPI_INTEGER,
     +                    ipartner,itag,comm,mpistat,ierr)
          call mpi_get_count(mpistat,MPI_INTEGER,recvcolsiz,ierr)
          call mergelists(supinds(scptr),lbotsiz,
     +                    supinds(supindsptr),recvcolsiz,
     +                    resdcol,resdcolsiz)

        end if

      end do

      end if 

      call mpi_reduce(nnzl,nnz,1,MPI_INTEGER,MPI_SUM,0,comm,ierr)

      end

