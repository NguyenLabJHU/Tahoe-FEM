      program pspases_test

***** This code reads in a matrix file in .rsa format (either HB or RB), and 
***** tests PSPASES functionality for various parameters and solution paths.
***** 					        written by: Mahesh Joshi
*****  $Id: pspases_testfr.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $

      implicit none

      include 'mpif.h'

      integer lendp,itag
      parameter(lendp=8,itag=1)
      double precision zero
      parameter(zero=0.d0)

      integer DPSPACEN_REPS
      parameter(DPSPACEN_REPS=3)

      integer, allocatable :: aptrs(:,:),ainds(:),rowdist(:)
      integer, allocatable :: rowdistbx(:),order(:),sizes(:)

      double precision, allocatable :: avals(:),b(:,:)
      double precision, allocatable :: x(:,:)

      integer, allocatable :: taptrs(:),tainds(:),aisizes(:)
      integer, allocatable :: ptrex(:),index(:),painds(:)
      integer, allocatable :: wrksb(:),wrkrb(:),wrkint(:),paptrs(:)
      double precision, allocatable :: tavals(:),tb(:),wrksbv(:)
      double precision, allocatable :: valex(:),pavals(:)
      double precision, allocatable :: wrkrbv(:)

      integer iargc

      integer pp,myid,blk,nrhs,N,nbuf(10),msglen,br,mynb,ldb,ldx
      integer Npp,Npp1,asize,is1,i,j,k,l,j1,j2,j3,size,piv,proc,nown
      integer mpistat(MPI_STATUS_SIZE),ierr,ioasize,m,maxm
      integer options(16),checksymm,sortinds,serialorder,runopt
      integer beginOpt,endOpt,clean_option,mynnodes,comm
      integer(8) pspcomm,pspcommy,pspcommf,pspcommn(DPSPACEN_REPS)
      integer fmt,mxasize,nextcolptr,eptr,eind,eval,iex,lenex
      character*64 inpfile
      character*16 inpfile1
      double precision doptions(16),emax
      double precision time0,Ftime,Otime,Ytime,Ntime,Ttime,Ltime
      double precision, allocatable :: memPM(:),memFM(:)
      double precision minPM,maxPM,sumPM,minFM,maxFM,sumFM,minTM,maxTM
      integer psc,pscv,psd,psdv,prc,prcv,prd,prdv,pcl
      integer sendsiz,sendsizv,recvsiz,recvsizv

C      *----------- PART FROM hbcode1.f at HBSMC
      CHARACTER      TITLE*72, MTRXID*8,  MXTYPE*3,  RHSTYP*3, 
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD,  VALCRD,  RHSCRD,
     1               NROW,  NCOL,  NNZERO, NELTVL, NRHSIX
C      *----------- 

      call mpi_init(ierr)

      comm = MPI_COMM_WORLD
      call mpi_comm_size(comm,pp,ierr)
      call mpi_comm_rank(comm,myid,ierr)

      if(myid.eq.0) then

        call getarg(0,inpfile1)
        i = iargc()
        if(i.lt.5) then
          write(*,*)'Usage: ',inpfile1,'<file> <fmt> <ptr_npl>'
          write(*,51)'<ind_npl> <val_npl> [<blk> [<nrhs> [<br>        '
          write(*,51)'[<symm> [<sort> [<oopt> [<ropt>] ]]]]]]         '
          print *,'<file>    : input file containing A in .rsa format'
          print *,'<fmt>     : if 0, Harwell-Boeing    (.hb) Format'
          print *,'            if 1, Rutherford-Boeing (.rb) Format'
          print *,'<ptr_npl> : number of items per line in PTRFMT'
          print *,'<ind_npl> : number of items per line in INDFMT'
          print *,'<val_npl> : number of items per line in VALFMT'
          print *
          print *,' NOTE: The above three values should be obtained'
          print *,'       by examining the first three format specs'
          print *,'       on the fourth line of the input file.' 
          print *
          print *,'Optional Parameters --'
          print *,'<blk>    : block size for distribution of A'
          print *,'<nrhs>   : number of right hand sides'
          print *,'<br>     : blocking for nrhs'
          print *,'<symm>   : if 1, check symmetry of input matrix'
          print *,'<sort>   : if 1, the indices in ainds need sorting'
          print *,'<oopt>   : if 1, serial METIS is used for ordering'
          print *,'           if 0, parallel METIS is used for ordering'
          print *,'<ropt>   : if 0, DPSPACEF+DPSPACET'
          print *,'           if 1, PSPACEO+PSPACEY+DPSPACEN+DPSPACET'
          print *,'           if 2, PSPACEO+PSPACEY'
          print *,'Defaults: blk=64, nrhs=1, br=nrhs'
          print *,'          symm,sort,oopt=0, ropt = 2 -> 1 -> 0'
          call mpi_abort(comm,0,ierr)
        end if
 51   FORMAT(24X,A48)

        call getarg(1,inpfile)
        open (unit = 90, file = inpfile, status = 'old')

        print *
        print *,'* Testing ',inpfile

        call getarg(2,inpfile)
        call a2i(fmt,inpfile,64)
        if(fmt.ne.1 .and. fmt.ne.0) then
          print *,'<fmt> has to be 0 (HB) or 1 (RB)'
          call mpi_abort(comm,0,ierr)
        end if

C      *----------- PART FROM hbcode1.f at HBSMC
        if( fmt.eq.0 ) then ! Harwell-Boeing (HB) format header lines
          READ ( 90, 1000 ) TITLE , MTRXID,
     1                      TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, 
     2                      MXTYPE, NROW , NCOL , NNZERO, NELTVL,
     3                      PTRFMT, INDFMT, VALFMT, RHSFMT
          if(RHSCRD.ne.0) then
            READ ( 90, 1001 ) RHSTYP , NRHS , NRHSIX
            print *,'WARNING: Although the matrix file has RHS '
            print *,'         information, we will ignore it and'
            print *,'         generate our own random RHS.'
          end if
        end if

        if( fmt.eq.1 ) then ! Rutherford-Boeing (RB) format header lines
          READ ( 90, 2000 ) TITLE , MTRXID,
     1                      TOTCRD, PTRCRD, INDCRD, VALCRD,
     2                      MXTYPE, NROW , NCOL , NNZERO, NELTVL,
     3                      PTRFMT, INDFMT, VALFMT
        end if

 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
 1001 FORMAT ( A3, 11X, 2I14 )

 2000 FORMAT ( A72,A8 / I14,3(1X,I13) / A3,11X,4(1X,I13) / 2A16,A20 )
C      *----------- 

        if(NROW.ne.NCOL) then
          print *,'Matrix is not symmetric in dimensions'
          call mpi_abort(comm,0,ierr)
        end if
        if( (MXTYPE(1:1).ne.'r' .and. MXTYPE(1:1).ne.'R') .or.
     +      (MXTYPE(2:2).ne.'s' .and. MXTYPE(2:2).ne.'S') .or.
     +      (MXTYPE(3:3).ne.'a' .and. MXTYPE(3:3).ne.'A') ) then
          print *,'The matrix must be in RSA format, where'
          print *,'  RSA = Real Symmetric Assembled'
          call mpi_abort(comm,0,ierr)
        end if
        if(VALCRD.le.0) then
          print *,'Values must be specified in the file'
          call mpi_abort(comm,0,ierr)
        end if

        N = NROW
 
        call getarg(3,inpfile)
        call a2i(eptr,inpfile,64)
        call getarg(4,inpfile)
        call a2i(eind,inpfile,64)
        call getarg(5,inpfile)
        call a2i(eval,inpfile,64)

        nrhs = 0
        blk = 0
        br = 0
        checksymm = 0
        sortinds = 0
        serialorder = 0
        runopt = -1
        if(i.gt.5) then
         call getarg(6,inpfile)
         call a2i(blk,inpfile,64)
         if(i.gt.6) then
          call getarg(7,inpfile)
          call a2i(nrhs,inpfile,64)
          if(i.gt.7) then
           call getarg(8,inpfile)
           call a2i(br,inpfile,64)
           if(i.gt.8) then
            call getarg(9,inpfile)
            call a2i(checksymm,inpfile,64)
            if(i.gt.9) then
             call getarg(10,inpfile)
             call a2i(sortinds,inpfile,64)
             if(i.gt.10) then
              call getarg(11,inpfile)
              call a2i(serialorder,inpfile,64)
              if(i.gt.11) then
               call getarg(12,inpfile)
               call a2i(runopt,inpfile,64)
              end if
             end if
            end if
           end if
          end if
         end if
        end if
        if(blk.le.0) blk = min(64,N)
        if(nrhs.le.0) nrhs = 1
        if(br.le.0) br = nrhs
        if(runopt.lt.0 .or. runopt.gt.2) runopt = -1

        nbuf(1) = N
        nbuf(2) = blk
        nbuf(3) = nrhs
        nbuf(4) = br
        nbuf(5) = checksymm
        nbuf(6) = sortinds
        nbuf(7) = serialorder
        nbuf(8) = runopt

        print *,'Dimension = ',N,' #processors = ',pp
        print *,'Parameters: blk  = ',blk,' nrhs = ',nrhs,' br = ',br
        print *,'            symm = ',checksymm,' sort = ',sortinds,
     +          ' serialorder =',serialorder
        print *
      end if

      call mpi_bcast(nbuf,10,MPI_INTEGER,0,comm,ierr)

      if(myid.ne.0) then
        N = nbuf(1)
        blk = nbuf(2)
        nrhs = nbuf(3)
        br = nbuf(4)
        checksymm = nbuf(5)
        sortinds = nbuf(6)
        serialorder = nbuf(7)
        runopt = nbuf(8)
      end if

      Npp = N/pp
      Npp1 = (N+pp-1)/pp
      if(N-Npp1*(pp-1).gt.0 .and. (Npp1*pp-N .lt. N-Npp*pp)) Npp=Npp1

      allocate(rowdist(0:pp),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(1)'
        call mpi_abort(comm,0,ierr)
      end if

      do i=0,pp-1
        rowdist(i) = i*Npp
      end do
      rowdist(pp) = N

*      read data      
      
      m = rowdist(myid+1)-rowdist(myid)
      maxm = max(Npp,N-(pp-1)*Npp)

      allocate(paptrs(0:m),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(2)'
        call mpi_abort(comm,0,ierr)
      end if

      allocate(aisizes(0:pp-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(3)'
        call mpi_abort(comm,0,ierr)
      end if

      if(myid.eq.0) then

       allocate(ptrex(0:eptr),stat=is1)
       if(is1.ne.0) then
         print *,myid,': memory allocation failure(4)'
         call mpi_abort(comm,0,ierr)
       end if

       m = rowdist(1)-rowdist(0)
       lenex = eptr-mod(m+1,eptr)+1

       READ ( 90, PTRFMT ) ( paptrs (i), i = 0, m-1 ), 
     +                     ( ptrex (i), i = 0, lenex-1 )

       paptrs(m) = ptrex(0)

       aisizes(0) = paptrs(m)-paptrs(0)
       mxasize = aisizes(0)

       allocate(taptrs(0:maxm),stat=is1)
       if(is1.ne.0) then
        print *,myid,': memory allocation failure(5)'
        call mpi_abort(comm,0,ierr)
       end if

       do proc=1,pp-1

        m = rowdist(proc+1)-rowdist(proc)
        j1 = ptrex(0)-1

        if(lenex.gt.m) then

          do iex=0,m
            taptrs(iex) = ptrex(iex)-j1
          end do

          j = 0
          do iex=m,lenex-1
            ptrex(j) = ptrex(iex)
            j = j+1
          end do

          lenex = lenex-m

        else

          do iex=0,lenex-1
            taptrs(iex) = ptrex(iex)-j1
          end do

          if(proc.eq.pp-1) then
            lenex = 1
          else
            lenex = eptr-mod(m-iex+1,eptr)+1
          end if

          READ ( 90, PTRFMT ) ( taptrs (i), i = iex, m-1 ),
     +                      ( ptrex (i), i = 0,lenex-1)

          do i=iex,m-1
            taptrs(i) = taptrs(i)-j1
          end do
          taptrs(m) = ptrex(0)-j1

        end if

        aisizes(proc) = taptrs(m)-taptrs(0)
        mxasize = max(aisizes(proc),mxasize)

        call mpi_send(taptrs,m+1,MPI_INTEGER,proc,itag,
     +                comm,ierr)
       end do

       deallocate(taptrs)
       deallocate(ptrex)

      else

        m = rowdist(myid+1)-rowdist(myid)
        call mpi_recv(paptrs,m+1,MPI_INTEGER,0,itag,
     +                comm,mpistat,ierr)

      end if

      call mpi_bcast(aisizes,pp,MPI_INTEGER,0,comm,ierr)
      asize = aisizes(myid)

      allocate(painds(asize),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(6)'
        call mpi_abort(comm,0,ierr)
      end if

      if(myid.eq.0) then

       allocate(index(eind),stat=is1)
       if(is1.ne.0) then
         print *,myid,': memory allocation failure(7)'
         call mpi_abort(comm,0,ierr)
       end if

       lenex = eind-mod(asize,eind)

       READ ( 90, INDFMT ) ( painds (i), i = 1, asize ), 
     +                     ( index (i), i = 1, lenex)

       allocate(tainds(mxasize),stat=is1)
       if(is1.ne.0) then
        print *,myid,': memory allocation failure(8)'
        call mpi_abort(comm,0,ierr)
       end if

       do proc=1,pp-1

        j = aisizes(proc)

        if(lenex.ge.j) then

          do iex=1,j
            tainds(iex) = index(iex)
          end do

          i = 1
          do iex=j+1,lenex
            index(i) = index(iex)
            i = i+1
          end do

          lenex = lenex-j

        else

          do iex=1,lenex
            tainds(iex) = index(iex)
          end do

          if(proc.eq.pp-1) then
            lenex = 0
          else
            lenex = eind-mod(j-lenex,eind)
          end if

          READ ( 90, INDFMT ) ( tainds (i), i = iex, j ),
     +                      ( index (i), i = 1,lenex )

        end if

        call mpi_send(tainds,j,MPI_INTEGER,proc,itag,
     +                comm,ierr)

       end do

       deallocate(tainds)
       deallocate(index)

      else

        call mpi_recv(painds,asize,MPI_INTEGER,0,itag,
     +                comm,mpistat,ierr)
      end if

      allocate(pavals(asize),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(9)'
        call mpi_abort(comm,0,ierr)
      end if

      if(myid.eq.0) then

       allocate(valex(eval),stat=is1)
       if(is1.ne.0) then
         print *,myid,': memory allocation failure(10)'
         call mpi_abort(comm,0,ierr)
       end if

       lenex = eval-mod(asize,eval)

       READ ( 90, VALFMT ) ( pavals (i), i = 1, asize ), 
     +                     ( valex (i), i = 1, lenex)

       allocate(tavals(mxasize),stat=is1)
       if(is1.ne.0) then
        print *,myid,': memory allocation failure(11)'
        call mpi_abort(comm,0,ierr)
       end if

       do proc=1,pp-1

        j = aisizes(proc)

        if(lenex.ge.j) then

          do iex=1,j
            tavals(iex) = valex(iex)
          end do

          i = 1
          do iex=j+1,lenex
            valex(i) = valex(iex)
            i = i+1
          end do

          lenex = lenex-j

        else

          do iex=1,lenex
            tavals(iex) = valex(iex)
          end do

          if(proc.eq.pp-1) then
            lenex = 0
          else
            lenex = eval-mod(j-lenex,eval)
          end if

          READ ( 90, VALFMT ) ( tavals (i), i = iex, j ),
     +                        ( valex (i), i = 1,lenex )

        end if

        call mpi_send(tavals,lendp*j,MPI_BYTE,proc,itag,
     +                comm,ierr)

       end do

       deallocate(tavals)
       deallocate(valex)
       close(90)

      else

        call mpi_recv(pavals,asize*lendp,MPI_BYTE,0,itag,
     +                comm,mpistat,ierr)
      end if

      deallocate(aisizes)

*     convert from symmetric storage to full storage as required by PSPASES.
      
      allocate(wrkint(0:9*pp-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(12)'
        call mpi_abort(comm,0,ierr)
      end if
      psc = 0
      pscv = pp
      psd = 2*pp
      psdv = 3*pp
      prc = 4*pp
      prcv = 5*pp
      prd = 6*pp
      prdv = 7*pp
      pcl = 8*pp

      m = rowdist(myid+1)-rowdist(myid)

      do proc=0,pp-1
        wrkint(psc+proc) = 0
        wrkint(pscv+proc) = 0
        wrkint(psd+proc) = 0
      end do

      do i=0,m-1
        k = paptrs(i)
        l = paptrs(i+1) - k
        do proc=myid,pp-1
          wrkint(psd+proc) = 0
        end do
        proc = myid
        do j=1,l-1
          j1 = painds(k+j)-1
          do while(j1.ge.rowdist(proc+1))
            proc = proc+1
          end do
          if(wrkint(psd+proc).eq.0) then
            wrkint(psc+proc) = wrkint(psc+proc)+2
            wrkint(psd+proc) = 1
          end if
          wrkint(psc+proc) = wrkint(psc+proc)+1
          wrkint(pscv+proc) = wrkint(pscv+proc)+1
        end do
      end do

      call mpi_alltoall(wrkint(psc),1,MPI_INTEGER,wrkint(prc),
     +                  1,MPI_INTEGER,comm,ierr)
      call mpi_alltoall(wrkint(pscv),1,MPI_INTEGER,wrkint(prcv),
     +                  1,MPI_INTEGER,comm,ierr)

      wrkint(psd) = 0
      wrkint(psdv) = 0
      sendsiz = 0
      sendsizv = 0
      do proc=1,pp-1
        wrkint(psd+proc) = wrkint(psd+proc-1)+wrkint(psc+proc-1)
        wrkint(psdv+proc) = wrkint(psdv+proc-1)+wrkint(pscv+proc-1)
        sendsiz = sendsiz+wrkint(psc+proc-1)
        sendsizv = sendsizv+wrkint(pscv+proc-1)
        wrkint(psc+proc-1) = 0
        wrkint(pscv+proc-1) = 0
      end do
      sendsiz = sendsiz+wrkint(psc+pp-1)
      sendsizv = sendsizv+wrkint(pscv+pp-1)
      wrkint(psc+pp-1) = 0
      wrkint(pscv+pp-1) = 0

      allocate(wrksb(0:sendsiz-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(13)'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(wrksbv(0:sendsizv-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(14)'
        call mpi_abort(comm,0,ierr)
      end if

      do i=0,m-1
        k = paptrs(i)
        l = paptrs(i+1) - k
        do proc=myid,pp-1
          wrkint(prd+proc) = 0
        end do
        proc = myid
        do j=1,l-1
          j1 = painds(k+j)-1
          do while(j1.ge.rowdist(proc+1))
            proc = proc+1
          end do
          j2 = wrkint(psd+proc)+wrkint(psc+proc)
          if(wrkint(prd+proc).eq.0) then
            wrksb(j2)   = i+rowdist(myid)
            wrksb(j2+1) = 0
            wrkint(pcl+proc) = j2+1
            wrkint(psc+proc) = wrkint(psc+proc)+2
            j2 = j2+2
            wrkint(prd+proc) = 1
          end if
          wrksb(j2) = j1
          wrkint(psc+proc) = wrkint(psc+proc)+1
          wrksbv(wrkint(psdv+proc)+wrkint(pscv+proc)) = pavals(k+j)
          wrkint(pscv+proc) = wrkint(pscv+proc)+1
          j2 = wrkint(pcl+proc)
          wrksb(j2) = wrksb(j2)+1
        end do
      end do

      wrkint(prd) = 0
      wrkint(prdv) = 0
      recvsiz = 0
      recvsizv = 0
      do proc=1,pp-1
        wrkint(prd+proc) = wrkint(prd+proc-1)+wrkint(prc+proc-1)
        wrkint(prdv+proc) = wrkint(prdv+proc-1)+wrkint(prcv+proc-1)
        recvsiz = recvsiz+wrkint(prc+proc-1)
        recvsizv = recvsizv+wrkint(prcv+proc-1)
      end do
      recvsiz = recvsiz+wrkint(prc+pp-1)
      recvsizv = recvsizv+wrkint(prcv+pp-1)

      allocate(wrkrb(0:recvsiz-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(15)'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(wrkrbv(0:recvsizv-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(16)'
        call mpi_abort(comm,0,ierr)
      end if

      call mpi_alltoallv(wrksb,wrkint(psc),wrkint(psd),MPI_INTEGER,
     +                   wrkrb,wrkint(prc),wrkint(prd),MPI_INTEGER,
     +                   comm,ierr)

      call mpi_alltoallv(wrksbv,wrkint(pscv),wrkint(psdv),
     +                   MPI_DOUBLE_PRECISION,wrkrbv,
     +                   wrkint(prcv),wrkint(prdv),
     +                   MPI_DOUBLE_PRECISION,comm,ierr)

      deallocate(wrkint)
      deallocate(wrksb)
      deallocate(wrksbv)

      allocate(aptrs(2,0:m-1),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(17)'
        call mpi_abort(comm,0,ierr)
      end if

      do k=0,m-1
        aptrs(2,k) = paptrs(k+1)-paptrs(k)
      end do

      j1 = 0
      do while(j1 .lt. recvsiz)
        l = wrkrb(j1+1)
        j1 = j1+2
        do j=0,l-1
          k = wrkrb(j1+j)-rowdist(myid)
          aptrs(2,k) = aptrs(2,k)+1
        end do
        j1 = j1+l
      end do

      aptrs(1,0) = 1
      do i=1,m-1
        aptrs(1,i) = aptrs(1,i-1)+aptrs(2,i-1)
        aptrs(2,i-1) = 0
      end do
      aptrs(2,m-1) = 0

      allocate(ainds(asize+recvsizv),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(18)'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(avals(asize+recvsizv),stat=is1)
      if(is1.ne.0) then
        print *,myid,': memory allocation failure(19)'
        call mpi_abort(comm,0,ierr)
      end if

      j1 = 0
      j2 = 0
      do while(j1 .lt. recvsiz)
        i = wrkrb(j1)
        l = wrkrb(j1+1)
        j1 = j1+2
        do j=0,l-1
          k = wrkrb(j1+j)-rowdist(myid)
          j3 = aptrs(1,k)+aptrs(2,k)
          ainds(j3) = i
          avals(j3) = wrkrbv(j2+j)
          aptrs(2,k) = aptrs(2,k)+1
        end do
        j1 = j1+l
        j2 = j2+l
      end do

      do i=0,m-1
        j2 = paptrs(i)
        j3 = aptrs(1,i)+aptrs(2,i)
        l = paptrs(i+1)-j2
        do j=0,l-1
          ainds(j3+j) = painds(j2+j)-1
          avals(j3+j) = pavals(j2+j)
        end do
        aptrs(2,i) = aptrs(2,i)+l
      end do

      deallocate(paptrs)
      deallocate(painds)
      deallocate(pavals)
      deallocate(wrkrb)
      deallocate(wrkrbv)

* test PSPASES library functions.

      options(1) = blk        ! block size
      options(2) = checksymm  ! check symmetry (1: yes,0: no, Default:0)
      options(3) = sortinds   ! sort indices (1: yes,0: no, Default:0)
      options(4) = serialorder! use serial METIS (1: yes, 0: no, Default:0)

      if(runopt.eq.-1) then
        beginOpt = 2
        endOpt = 0
      else
        beginOpt = runopt
        endOpt = runopt
      end if

      do runopt=beginOpt,endOpt,-1

      call mpi_barrier(comm,ierr)
      if(myid.eq.0) then
      print *
      if(runopt.eq.0) 
     +  print *,'-------> Testing DPSPACEF+DPSPACET'
      if(runopt.eq.1) 
     +  print *,'-------> Testing PSPACEO+PSPACEY+DPSPACEN+DPSPACET'
      if(runopt.eq.2) 
     +  print *,'-------> Testing PSPACEO+PSPACEY'
      print *
      end if

      if(runopt.eq.0) then
        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +  print *,'calling DPSPACEF (Ordering+sYmbolic+Numerical)..'

        call mpi_barrier(comm,ierr)
        time0 = mpi_wtime()

        call DPSPACEF(rowdist,aptrs,ainds,avals,options,
     +                doptions,pspcommf,comm)

        call mpi_barrier(comm,ierr)
        Ftime = mpi_wtime()-time0

        clean_option = 1

      else

        mynnodes = rowdist(myid+1)-rowdist(myid)
        allocate(order(0:mynnodes-1),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(20)'
          call mpi_abort(comm,0,ierr)
        end if
        allocate(sizes(0:2*pp-1),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(21)'
          call mpi_abort(comm,0,ierr)
        end if

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +  print *,'calling  PSPACEO (Compute fill-reducing Ordering)..'

        call mpi_barrier(comm,ierr)
        time0 = mpi_wtime()

        call PSPACEO(rowdist,aptrs,ainds,order,sizes,options,
     +               comm)

        call mpi_barrier(comm,ierr)
        Otime = mpi_wtime()-time0

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +  print *,'calling  PSPACEY (sYmbolic Factorization)..'

        call mpi_barrier(comm,ierr)
        time0 = mpi_wtime()

        call PSPACEY(rowdist,aptrs,ainds,order,sizes,options,
     +               doptions,pspcommy,comm)

        call mpi_barrier(comm,ierr)
        Ytime = mpi_wtime()-time0

        deallocate(order)
        deallocate(sizes)

        clean_option = 0

      end if

      if(myid.eq.0) then
        allocate(memPM(0:pp-1),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(22)'
          call mpi_abort(comm,0,ierr)
        end if
        allocate(memFM(0:pp-1),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(23)'
          call mpi_abort(comm,0,ierr)
        end if
      end if

      doptions(1) = doptions(1)/dble(1024*1024)
      doptions(2) = doptions(2)/dble(1024*1024)
      call mpi_gather(doptions(1),1,MPI_DOUBLE_PRECISION,memPM,1,
     +                MPI_DOUBLE_PRECISION,0,comm,ierr)
      call mpi_gather(doptions(2),1,MPI_DOUBLE_PRECISION,memFM,1,
     +                MPI_DOUBLE_PRECISION,0,comm,ierr)

      call mpi_barrier(comm,ierr)
      if(myid.eq.0) then
        minPM = memPM(0)
        maxPM = memPM(0)
        sumPM = memPM(0)
        do i=1,pp-1
          minPM = min(minPM,memPM(i))
          maxPM = max(maxPM,memPM(i))
          sumPM = sumPM + memPM(i)
        end do

        write(*,*)

        write(*,100)'                               ','Max','Min','Sum'
        write(*,101)'  Memory Consumed by PSPCOMM = ',maxPM,minPM,sumPM

        if(runopt.ne.0) then
          minFM = memFM(0)
          maxFM = memFM(0)
          sumFM = memFM(0)
          minTM = memFM(0)+memPM(0)
          maxTM = memFM(0)+memPM(0)
          do i=1,pp-1
            minFM = min(minFM,memFM(i))
            maxFM = max(maxFM,memFM(i))
            sumFM = sumFM + memFM(i)
            minTM = min(minTM,memFM(i)+memPM(i))
            maxTM = max(maxTM,memFM(i)+memPM(i))
          end do

          write(*,101)'  More Factor Memory Needed  = ',maxFM,minFM,
     +                sumFM
          write(*,101)'  Total Factor Memory Needed = ',maxTM,minTM,
     +                sumPM+sumFM
        end if

        deallocate(memPM)
        deallocate(memFM)

        write(*,*)

        write(*,102)'  NNZ_Lower_Triangle_of_A = ',doptions(3),
     +              '  NNZ_L          = ',doptions(4)
        write(*,102)'  Tree_Opcount_Imbalance  = ',doptions(6),
     +              '  Factor_Opcount = ',doptions(5)

        write(*,*)

        if(runopt.eq.0) then
          Ltime = Ftime
          write(*,103)'  DPSPACEF Time (Ftime)      = ',Ftime
          write(*,104)'  Factor_Opcount / Ftime     = ',
     +            doptions(5)/Ftime*1.d-6
        else
          Ltime = Otime+Ytime
          write(*,103)'  Order Time                 = ',Otime
          write(*,103)'  sYmbolic Time              = ',Ytime
        end if

        write(*,*)

      end if

      if(runopt.eq.1) then

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +    print *,'calling DPSPACEN (Numerical Factorization)..'

        call mpi_barrier(comm,ierr)
        time0 = mpi_wtime()

        do i=1,DPSPACEN_REPS
          call DPSPACEN(rowdist,aptrs,ainds,avals,
     +                  pspcommy,pspcommn(i),comm)
        end do

        call mpi_barrier(comm,ierr)
        Ntime = (mpi_wtime()-time0)/dble(DPSPACEN_REPS)

        clean_option = 1

        if(myid.eq.0) then
          Ltime = Ltime + Ntime
          write(*,*)
          write(*,103)'  DPSPACEN Time              = ',Ntime
          write(*,104)'  Numerical Factor Perf      = ',
     +                doptions(5)/Ntime*1.e-6
          write(*,*)
        end if

      end if

      if(runopt.eq.1) then

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) print *,'calling PSPACEC (option=1) on pspcommY..'
	clean_option = 1;
        call PSPACEC(pspcommy,clean_option)

	do i=2,DPSPACEN_REPS
          call mpi_barrier(comm,ierr)
          if(myid.eq.0) 
     +	    print *,'calling PSPACEC (option=0) on pspcommN(',i,')..'
	  clean_option = 0;
          call PSPACEC(pspcommn(i),clean_option)
	end do

	if(myid.eq.0) print *

      else if(runopt.eq.2) then

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) print *,'calling PSPACEC (option=0) on pspcommY..'
	clean_option = 0;
        call PSPACEC(pspcommy,clean_option)

      end if

      if(runopt.ne.2) then

        allocate(rowdistbx(0:pp),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(24)'
          call mpi_abort(comm,0,ierr)
        end if

        j = 4
        do i=0,pp-1
          rowdistbx(i) = i*(Npp-Npp/j)
        end do
        rowdistbx(pp) = N
        maxm = max(Npp-Npp/j,N-(pp-1)*(Npp-Npp/j))

        mynb = rowdistbx(myid+1)-rowdistbx(myid)

        ldb = mynb
        allocate(b(0:ldb-1,nrhs),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(25)'
          call mpi_abort(comm,0,ierr)
        end if

        if(myid.eq.0) then
         i = 13
         call initrnd(i)

         m = rowdistbx(1)-rowdistbx(0)
         do j=1,nrhs
          do i=0,m-1
           call mydrand48n(b(i,j))
           b(i,j) = 2.d0*(b(i,j)-0.5d0)
          end do
         end do

         allocate(tb(0:maxm-1),stat=is1)
         if(is1.ne.0) then
          print *,myid,': memory allocation failure(26)'
          call mpi_abort(comm,0,ierr)
         end if

         do proc=1,pp-1
          m = rowdistbx(proc+1)-rowdistbx(proc)
          do j=1,nrhs
           do i=0,m-1
            call mydrand48n(tb(i))
            tb(i) = 2.d0*(tb(i)-0.5d0)
           end do
           call mpi_send(tb,lendp*m,MPI_BYTE,proc,itag,
     +                   comm,ierr)
          end do
         end do
         deallocate(tb)

        else

          do i=1,nrhs
            call mpi_recv(b(0,i),lendp*mynb,MPI_BYTE,0,itag,
     +                    comm,mpistat,ierr)
          end do

        end if

        ldx = mynb
        allocate(x(0:ldx-1,nrhs),stat=is1)
        if(is1.ne.0) then
          print *,myid,': memory allocation failure(27)'
          call mpi_abort(comm,0,ierr)
        end if

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +    print *,'calling DPSPACET (Triangular Systems Solution)..'

	if(runopt.eq.0) then
	  pspcomm = pspcommf;
	else
	  pspcomm = pspcommn(1);
	end if

        options(1) = br ! blocking factor on nrhs.

        call mpi_barrier(comm,ierr)
        time0 = mpi_wtime()

        call DPSPACET(rowdistbx,nrhs,b,ldb,x,ldx,
     +                options,pspcomm,comm)

        call mpi_barrier(comm,ierr)
        Ttime = mpi_wtime()-time0

        if(myid.eq.0) then
          write(*,*)
          write(*,103)'  DPSPACET Time              = ',Ttime
          write(*,104)'  Triangular Solution Perf   = ',
     +                  doptions(4)*dble(nrhs)*4.d-6/Ttime
          write(*,*)
          Ltime = Ltime+Ttime
          write(*,103)'  Total Solver Time          = ',Ltime
          write(*,*)
        end if

        clean_option = 0
	if(runopt.eq.0) then
          call mpi_barrier(comm,ierr)
          if(myid.eq.0) 
     +       print *,'calling PSPACEC (option=0) on pspcommF..'
	  call PSPACEC(pspcommf,clean_option);  
        else
          call mpi_barrier(comm,ierr)
          if(myid.eq.0) 
     +       print *,'calling PSPACEC (option=0) on pspcommN( 1 )..'
	  call PSPACEC(pspcommn(1),clean_option);  

          call mpi_barrier(comm,ierr)
          if(myid.eq.0) 
     +       print *,'calling PSPACEC (option=0) on pspcommY..'
	  call PSPACEC(pspcommy,clean_option);  
	end if

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) 
     +    print *,'all stages of PSPASES are done! checking B-AX ..'

        call CHECKB_AX(rowdist,aptrs,ainds,avals,
     +                 rowdistbx,nrhs,b,ldb,x,ldx,
     +                 emax,comm)

        call mpi_barrier(comm,ierr)
        if(myid.eq.0) print *,'max |B - AX| = ',emax

        deallocate(rowdistbx)
        deallocate(b)
        deallocate(x)

      end if !if runopt.ne.2

      end do

      deallocate(rowdist)
      deallocate(aptrs)
      deallocate(ainds)
      deallocate(avals)

      call mpi_finalize(ierr)

 100  format(A31,3A10)
 101  format(A31,3F10.3,' MB')
 102  format(A28,G10.3,A19,G12.4)
 103  format(A31,F10.3,' sec')
 104  format(A31,F10.3,' MFLOPS')

      end
