      program rsa2fcc

**** This code converts a matrix file  in .rsa (either HB or RB) format to
**** .fcc format file.
****    $Id: rsa2fcc.f,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $

      CHARACTER      TITLE*72, MTRXID*8,  MXTYPE*3,  RHSTYP*3, 
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD,  VALCRD,  RHSCRD,
     1               NROW,  NCOL,  NNZERO, NELTVL, NRHSIX

      INTEGER, allocatable :: COLPTR (:), ROWIND (:)
      DOUBLE PRECISION, allocatable :: VALUES (:)

      integer, allocatable :: aptrs (:,:), ainds (:)
      double precision, allocatable :: avals(:)

      integer iargc,asize,fmt
      character*64 inpfile
      character*10 inpfile1

      call getarg(0,inpfile1)
      i = iargc()
      if(i.lt.3) then
        print *,'Usage: ',inpfile1,' <rsa_file> <fmt> <fcc_file>'
        print *,'  <fmt>     : if 0, Harwell-Boeing    (.hb) Format'
        print *,'              if 1, Rutherford-Boeing (.rb) Format'
        stop
      end if

      call getarg(1,inpfile)
      IDIN = 90
      open (unit = IDIN, file = inpfile, status = 'old')

      call getarg(2,inpfile)
      call a2i(fmt,inpfile,64)
      if(fmt.ne.1 .and. fmt.ne.0) then
        print *,'<fmt> has to be 0 (HB) or 1 (RB)'
        stop
      end if

      call getarg(3,inpfile)
      IDOUT = 80
      open (unit = IDOUT, file = inpfile, status = 'unknown')

C    ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------

      if( fmt.eq.0 ) then ! Harwell-Boeing (HB) format header lines
        READ ( IDIN, 1000 ) TITLE , MTRXID,
     1                      TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, 
     2                      MXTYPE, NROW , NCOL , NNZERO, NELTVL,
     3                      PTRFMT, INDFMT, VALFMT, RHSFMT
        if(RHSCRD.ne.0) then
          READ ( IDIN, 1001 ) RHSTYP , NRHS , NRHSIX
        end if
      end if

      if( fmt.eq.1 ) then ! Rutherford-Boeing (RB) format header lines
        READ ( IDIN, 2000 ) TITLE , MTRXID,
     1                      TOTCRD, PTRCRD, INDCRD, VALCRD,
     2                      MXTYPE, NROW , NCOL , NNZERO, NELTVL,
     3                      PTRFMT, INDFMT, VALFMT
      end if

 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
 1001 FORMAT ( A3, 11X, 2I14 )

 2000 FORMAT ( A72,A8 / I14,3(1X,I13) / A3,11X,4(1X,I13) / 2A16,A20 )

      if(NROW.ne.NCOL) then
        print *,'Matrix is not symmetric in dimensions'
        stop
      end if
      if( (MXTYPE(1:1).ne.'r' .and. MXTYPE(1:1).ne.'R') .or.
     +    (MXTYPE(2:2).ne.'s' .and. MXTYPE(2:2).ne.'S') .or.
     +    (MXTYPE(3:3).ne.'a' .and. MXTYPE(3:3).ne.'A') ) then
        print *,'The matrix must be in RSA format, where'
        print *,'  RSA = Real Symmetric Assembled'
        stop
      end if
      if(VALCRD.le.0) then
        print *,'Values must be specified in the file'
        stop
      end if

      allocate(colptr(0:ncol),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if

C     -------------------------
C     ... READ MATRIX STRUCTURE
C     -------------------------

      READ ( IDIN, PTRFMT ) ( COLPTR (I), I = 0, NCOL )

      asize = colptr(ncol)-colptr(0)
      if(asize.ne.nnzero) then
        print *,'error in colptr'
        stop
      end if

      allocate(rowind(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if

      READ ( IDIN, INDFMT ) ( ROWIND (I), I = 1, NNZERO )

      IF  ( VALCRD .LE. 0 )  THEN

C         ----------------------
C         ... READ MATRIX VALUES
C         ----------------------

        print *,'values are not provided!'
        stop

      ENDIF

      allocate(values(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if

      READ ( IDIN, VALFMT ) ( VALUES (I), I = 1, NNZERO )

      close(IDIN)

      allocate(aptrs(2,0:nrow-1),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if
      asize = (nnzero-nrow)*2+nrow
      allocate(ainds(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if
      allocate(avals(asize),stat=is1)
      if(is1.ne.0) then
        print *,'memory allocation failure'
        stop 
      end if

      do i=0,nrow-1
        aptrs(2,i) = colptr(i+1)-colptr(i)
      end do

      do i=0,nrow-1
        do j = colptr(i)+1,colptr(i+1)-1
          m = rowind(j)-1
          aptrs(2,m) = aptrs(2,m)+1
        end do
      end do

      aptrs(1,0) = 1
      do i=1,nrow-1
        aptrs(1,i) = aptrs(1,i-1)+aptrs(2,i-1)
        aptrs(2,i-1) = 0
      end do
      aptrs(2,nrow-1) = 0

      do i=0,nrow-1
        do j = colptr(i)+1,colptr(i+1)-1
          m = rowind(j)-1
          n = aptrs(1,m)+aptrs(2,m)
          ainds(n) = i
          avals(n) = values(j)
          aptrs(2,m) = aptrs(2,m)+1
        end do
      end do

      do i=0,nrow-1
        do j=colptr(i),colptr(i+1)-1
          n = aptrs(1,i)+aptrs(2,i)
          ainds(n) = rowind(j)-1
          avals(n) = values(j)
          aptrs(2,i) = aptrs(2,i)+1
        end do
      end do

      write(IDOUT,*)nrow
      do i=0,nrow-1
        write(IDOUT,*)aptrs(2,i),(avals(j),ainds(j),
     +                j=aptrs(1,i),aptrs(1,i)+aptrs(2,i)-1)
      end do

      close(IDOUT)

      end

      subroutine a2i(i,a,LEN)
      integer i
      character a(*)

c      parameter(LEN = 64)

      i = 0
      j = 1
      l = 1
      do while (j .lt. LEN .and. a(j) .ne. '0'
     1                     .and. a(j) .ne. '1'
     1                     .and. a(j) .ne. '2'
     1                     .and. a(j) .ne. '3'
     1                     .and. a(j) .ne. '4'
     1                     .and. a(j) .ne. '5'
     1                     .and. a(j) .ne. '6'
     1                     .and. a(j) .ne. '7'
     1                     .and. a(j) .ne. '8'
     1                     .and. a(j) .ne. '9')
        if (a(j) .eq. '-') l = -1
        j = j + 1
      end do
      do while (j .lt. LEN .and. (a(j) .eq. '0'
     1                     .or. a(j) .eq. '1'
     1                     .or. a(j) .eq. '2'
     1                     .or. a(j) .eq. '3'
     1                     .or. a(j) .eq. '4'
     1                     .or. a(j) .eq. '5'
     1                     .or. a(j) .eq. '6'
     1                     .or. a(j) .eq. '7'
     1                     .or. a(j) .eq. '8'
     1                     .or. a(j) .eq. '9'))
        i = i * 10
        if (a(j) .eq. '1') i = i + 1
        if (a(j) .eq. '2') i = i + 2
        if (a(j) .eq. '3') i = i + 3
        if (a(j) .eq. '4') i = i + 4
        if (a(j) .eq. '5') i = i + 5
        if (a(j) .eq. '6') i = i + 6
        if (a(j) .eq. '7') i = i + 7
        if (a(j) .eq. '8') i = i + 8
        if (a(j) .eq. '9') i = i + 9
        j = j + 1
      end do
      i = i * l
      return
      end

