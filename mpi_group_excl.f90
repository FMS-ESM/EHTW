include "mpif.h"
  implicit none
  integer group_world, odd_group, even_group
  integer i, p, Neven, Nodd, nonmembers(0:7), ierr

  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
  call MPI_Comm_group(MPI_COMM_WORLD, group_world, ierr)

  Neven = (p + 1)/2   ! processes of MPI_COMM_WORLD are divided
  Nodd = p - Neven    ! into odd- and even-numbered groups
  do i=0,Neven - 1    ! "nonmembers" are even-numbered procs
    nonmembers(i) = 2*i
  enddo

  call MPI_Group_excl(group_world, Neven, nonmembers, odd_group, ierr)

