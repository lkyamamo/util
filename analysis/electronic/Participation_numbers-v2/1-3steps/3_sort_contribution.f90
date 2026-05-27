program main
implicit none
integer :: i,j,ierror
real*8  :: x,sum
integer :: nnstep,count
real*8  :: val(1000),val2,val3(1000)
integer :: nval(1000),nval2
integer :: id(1000),id2
integer :: ids

open(1,file="normalized_contribution_atom.dat")
open(10,file="normalized_contribution_atom_sort.dat")
do
  read(1,*,iostat=ierror)nnstep,count,ids
  if(ierror /= 0)exit
  do i=1,count
     read(1,*)nval(i),val(i),val3(i),id(i)
  end do
  do i=1,count
     do j=i+1,count
      if(val(i) < val(j))then
         val2 = val(i)
         val(i) = val(j)
         val(j) = val2
         nval2 = nval(i)
         nval(i) = nval(j)
         nval(j) = nval2
         id2 = id(i)
         id(i) = id(j)
         id(j) = id2
      end if
     end do
  end do
  write(10,*)nnstep,count,ids
  do i=1,count
     write(10,'(I6,1x,2(1x,ES10.3),1x,I3)')nval(i),val(i),val(i)*100.d0,id(i)
  end do
end do
close(1)
close(10)

stop
end program




