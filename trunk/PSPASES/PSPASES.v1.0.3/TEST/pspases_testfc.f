      program pspases_test

***** This code reads in a matrix file in .fcc format, and tests PSPASES
***** functionality for various parameters and solution paths.
***** 					        written by: Mahesh Joshi
*****  $Id: pspases_testfc.f,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $

      implicit none

      include 'mpif.h'

      integer lendp,itag
      parameter(lendp=8,itag=1)
      double precision zero
      parameter(zero=0.d0)

      integer NZPERCOL,DPSPACEN_REPS
      parameter(NZPERCOL=200,DPSPACEN_REPS=3)

      integer, allocatable :: aptrs(:,:),ainds(:),rowdist(:)
      integer, allocatable :: rowdistbx(:),order(:),sizes(:)

      double precision, allocatable :: avals(:),b(:,:)
      double precision, allocatable :: x(:,:)

      integer, allocatable :: taptrs(:,:),tainds(:)
      double precision, allocatable :: tavals(:),tb(:)

      integer iargc

      integer pp,myid,blk,nrhs,N,nbuf(10),msglen,br,mynb,ldb,ldx
      integer Npp,Npp1,asize,is1,i,j,k,size,piv,proc,nown
      integer mpistat(MPI_STATUS_SIZE),ierr,ioasize,m,maxm,maxasize
      integer options(16),checksymm,sortinds,serialorder,runopt
      integer beginOpt,endOpt,clean_option,mynnodes,comm
      integer(8) pspcomm,pspcommy,pspcommf,pspcommn(DPSPACEN_REPS)
      character*64 inpfile
      character*16 inpfile1
      double precision doptions(16),emax
      double precision time0,Ftime,Otime,Ytime,Ntime,Ttime,Ltime
      double precision, allocatable :: memPM(:),memFM(:)
      double precision minPM,maxPM,sumPM,minFM,maxFM,sumFM,minTM,maxTM

      call mpi_init(ierr)

      comm = MPI_COMM_WORLD
      call mpi_comm_size(comm,pp,ierr)
      call mpi_comm_rank(comm,myid,ierr)

      if(myid.eq.0) then

        call getarg(0,inpfile1)
        i = iargc()
        if(i.lt.1) then
          print *,'Usage: ',inpfile1,' <file> [<blk> [<nrhs> [<br>'
          write(*,51)'[<symm> [<sort> [<oopt> [<ropt>] ]]]]]]         '
          print *,'<file> : input file containing A in .fcc format'
          print *
          print *,'Optional Parameters --'
          print *,'<blk>  : block size for distribution of A'
          print *,'<nrhs> : number of right hand sides'
          print *,'<br>   : blocking for nrhs'
          print *,'<symm> : if 1, checks symmetry of input matrix'
          print *,'<sort> : if 1, sorts the indices in ainds'
          print *,'<oopt> : if 1, serial METIS is used for ordering'
          print *,'         if 0, parallel METIS is used for ordering'
          print *,'<ropt> : if 0, DPSPACEF+DPSPACET'
          print *,'         if 1, PSPACEO+PSPACEY+DPSPACEN+DPSPACET'
          print *,'         if 2, PSPACEO+PSPACEY'
          print *,'Defaults: blk=64, nrhs=1, br=nrhs'
          print *,'          symm,sort,oopt=0, ropt = 2 -> 1 -> 0'
          call mpi_abort(comm,0,ierr)
        end if
 51   FORMAT(24X,A48)

        call getarg(1,inpfile)
        open (unit = 90, file = inpfile, status = 'old')

        print *
        print *,'* Testing ',inpfile

        read(90,*)N

        nrhs = 0
        blk = 0
        br = 0
        checksymm = 0
        sortinds = 0
        serialorder = 0
        runopt = -1
        if(i.gt.1) then
         call getarg(2,inpfile)
         call a2i(blk,inpfile,64)
         if (i.gt.2) then
          call getarg(3,inpfile)
          call a2i(nrhs,inpfile,64)
          if(i.gt.3) then
           call getarg(4,inpfile)
           call a2i(br,inpfile,64)
           if(i.gt.4) then
            call getarg(5,inpfile)
            call a2i(checksymm,inpfile,64)
            if(i.gt.5) then
             call getarg(6,inpfile)
             call a2i(sortinds,inpfile,64)
             if(i.gt.6) then
              call getarg(7,inpfile)
              call a2i(serialorder,inpfile,64)
              if(i.gt.7) then
               call getarg(8,inpfile)
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
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if

      do i=0,pp-1
        rowdist(i) = i*Npp
      end do
      rowdist(pp) = N

      m = rowdist(myid+1)-rowdist(myid)
      asize = m*NZPERCOL

      maxm = max(Npp,N-(pp-1)*Npp)
      maxasize = maxm*NZPERCOL

      allocate(ainds(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(avals(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if
      allocate(aptrs(2,0:m-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
      end if

*      read data      
      
      if(myid.eq.0) then
       piv = 1
       do i=0,m-1
        read(90,*) size,(avals(k),ainds(k),k=piv,piv+size-1)
        if(piv+size-1 .gt. asize) then
         print *,'Insufficient memory allocated for ainds.'
         print *,'Increase the value of NZPERCOL variable.'
         call mpi_abort(comm,0,ierr)
        end if
        aptrs(1,i) = piv
        aptrs(2,i) = size
        piv = piv+size
       end do
       ioasize = piv-1

       allocate(tainds(maxasize),stat=is1)
       if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
       end if
       allocate(tavals(maxasize),stat=is1)
       if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
       end if
       allocate(taptrs(2,0:maxm-1),stat=is1)
       if(is1.ne.0) then
        print *,'memory allocation failure'
        call mpi_abort(comm,0,ierr)
       end if

       do proc=1,pp-1

        m = rowdist(proc+1)-rowdist(proc)
        piv = 1
        do i=0,m-1
         read(90,*) size,(tavals(k),tainds(k),k=piv,piv+size-1)
         if(piv+size-1 .gt. maxasize) then
          print *,'insufficient memory allocated for asize.'
          print *,'Increase the value of NZPERCOL variable.'
          call mpi_abort(comm,0,ierr)
         end if
         taptrs(1,i) = piv
         taptrs(2,i) = size
         piv = piv+size
        end do

        call mpi_send(tavals,lendp*(piv-1),MPI_BYTE,proc,itag,
     +                comm,ierr)
        call mpi_send(tainds,(piv-1),MPI_INTEGER,proc,itag,
     +                comm,ierr)
        call mpi_send(taptrs,2*m,MPI_INTEGER,proc,itag,
     +                comm,ierr)
       end do

       close(90)
       deallocate(tavals)
       deallocate(tainds)
       deallocate(taptrs)

      else

        call mpi_recv(avals,lendp*asize,MPI_BYTE,0,itag,
     +                comm,mpistat,ierr)
        call mpi_recv(ainds,asize,MPI_INTEGER,0,itag,
     +                comm,mpistat,ierr)
        call mpi_get_count(mpistat,MPI_INTEGER,ioasize,ierr)
        call mpi_recv(aptrs,2*m,MPI_INTEGER,0,itag,
     +                comm,mpistat,ierr)

      end if

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

      else

        mynnodes = rowdist(myid+1)-rowdist(myid)
        allocate(order(0:mynnodes-1),stat=is1)
        if(is1.ne.0) then
          print *,'memory allocation failure'
          call mpi_abort(comm,0,ierr)
        end if
        allocate(sizes(0:2*pp-1),stat=is1)
        if(is1.ne.0) then
          print *,'memory allocation failure'
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

      end if

      if(myid.eq.0) then
        allocate(memPM(0:pp-1),stat=is1)
        if(is1.ne.0) then
          print *,'memory allocation failure'
          call mpi_abort(comm,0,ierr)
        end if
        allocate(memFM(0:pp-1),stat=is1)
        if(is1.ne.0) then
          print *,'memory allocation failure'
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

        if(myid.eq.0) then
          Ltime = Ltime+Ntime
          write(*,*)
          write(*,103)'  DPSPACEN Time              = ',Ntime
          write(*,104)'  Numerical Factor Perf      = ',
     +            doptions(5)/Ntime*1.e-6
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
          print *,'memory allocation failure'
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
          print *,'memory allocation failure'
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
            print *,'memory allocation failure'
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
     +                      comm,ierr)
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
          print *,'memory allocation failure'
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
