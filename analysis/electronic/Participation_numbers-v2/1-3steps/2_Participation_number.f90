program main
implicit none
integer :: i,ierror,in
real*8  :: x,sum,sum2
integer :: count
integer :: nnstep
integer :: ids

!participation number using all atoms

open(1,file="normalized_contribution_atom.dat")
open(10,file="Participation_number.dat")
do
  read(1,*,iostat=ierror)nnstep,count,ids
  if(ierror /= 0)exit
  sum  = 0.d0
  sum2 = 0.d0
  do i=1,count
     read(1,*)in,x
     sum  = sum + x
     sum2 = sum2 + x**2
  end do
  write(10,'(I6,1x,2(1x,ES10.3),I5)')nnstep,sum,1.d0/sum2,ids
end do
close(1)
close(10)

stop
end program





