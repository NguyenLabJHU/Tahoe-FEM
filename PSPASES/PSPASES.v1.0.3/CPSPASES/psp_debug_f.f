C
C debugging wrappers for fortran MPI calls
C

C
C MPI_Send
C
      subroutine my_mpi_send(buf, count, dattype, dest, tag, comm, ierr)            

      implicit none
      include 'mpif.h'
      double precision buf(*)
      integer count, dattype, dest, tag, comm, ierr
      integer rank

      call mpi_comm_rank(comm, rank, ierr)      

      print *, rank, 'MPI_send: count =', count,
     1 ', dest =', dest, ', tag =', tag

      call mpi_send(buf, count, dattype, dest, tag, comm, ierr)

      return
      end

C
C MPI_Isend
C
      subroutine my_mpi_isend(buf, count, dattype, dest, tag, comm, req, 
     +                        ierr)            

      implicit none
      include 'mpif.h'
      double precision buf(*)
      integer count, dattype, dest, tag, comm, ierr, req
      integer rank

      call mpi_comm_rank(comm, rank, ierr)      

      print *, rank, 'MPI_Isend: count =', count,
     1 ', dest =', dest, ', tag =', tag

      call mpi_isend(buf, count, dattype, dest, tag, comm, req, ierr)

      return
      end

C
C MPI_Get_count
C
      subroutine my_mpi_get_count(stat,dattype,nbytes,ierr)

      implicit none
      include 'mpif.h'
      integer stat(MPI_STATUS_SIZE),dattype,nbytes,ierr
      integer rank

      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

C      print *, rank, 'MPI_Get_count: source =', stat(MPI_SOURCE), 
C     1 ', tag =', stat(MPI_TAG)

      call mpi_get_count(stat,dattype,nbytes,ierr)

      print *, rank, 'MPI_Get_count: source =', stat(MPI_SOURCE), 
     1', tag =', stat(MPI_TAG), ', count=', nbytes

      return
      end
