program main
implicit none
integer :: i,ierror,ierror2,in,il,in2,id
real*8  :: x,sum,y(3,1000)
integer :: id2(1000)
integer :: flag = 1  !1: electron, 2: hole 
integer :: nnstep,count,count2
integer :: nnnstep
integer :: nn
integer :: ids

!if(flag == 1)then
   open(1,file="../elec.dat")
!else if(flag == 2)then
!   open(1,file="../hole.dat")
!end if 
open(10,file="contribution_sum.dat")
do
  read(1,*,iostat=ierror)nnstep,count,ids
  if(ierror /= 0)exit
  sum = 0.d0
  do i=1,count
     read(1,*)in,il,nn
     sum = sum + dble(nn)
  end do
  write(10,'(I6,1x,ES10.3,I4,I5)')nnstep,sum,in,ids
  write(*,*)"sum=",sum,ids
end do
close(1)
close(10)


!if(flag == 1)then
   open(1,file="../elec.dat")
!else if(flag == 2)then
!   open(1,file="../hole.dat")
!end if
open(10,file="contribution_sum.dat")
open(11,file="normalized_contribution_l.dat")
do
  read(1,*,iostat=ierror)nnstep,count,ids
  if(ierror /= 0)exit
  read(10,*)nnnstep,sum,count2
  !write(11,*)nnstep,count,count2
  write(11,*)nnstep,count,count2,ids
  do i=1,count
     read(1,*)in,il,nn,id
     write(11,'(I4,I4,1x,ES10.3,1x,I3)')in,il,dble(nn)/sum &!*100.d0
                                 &,id
  end do
end do
close(1)
close(10)
close(11)

open(11,file="normalized_contribution_l.dat")
open(12,file="normalized_contribution_atom.dat")
do
  read(11,*,iostat=ierror)nnstep,count,count2,ids
  if(ierror /= 0)exit
  write(12,*)nnstep,count2,ids
  !write(4649,*)nnstep,count2
  do i=1,count
     read(11,*)y(1:3,i),id2(i)
  end do
  sum = 0.d0
  do i=1,count
     if(i == count)then
       sum = sum + y(3,i)
       write(12,'(I4,1x,2(1x,ES10.3),1x,I3)')nint(y(1,i)),sum,sum*100.d0,id2(i)
     else
       if(y(1,i+1) /= y(1,i))then
          sum = sum + y(3,i)     
          write(12,'(I4,1x,2(1x,ES10.3),1x,I3)')nint(y(1,i)),sum,sum*100.d0,id2(i)
          sum = 0.d0
       else if(y(1,i+1) == y(1,i))then
          sum = sum + y(3,i)     
       end if 
     end if
  end do

  !in2 = 1
  !sum = 0.d0
  !do i=1,100000
  !   read(11,'(I4,I4,1x,PE10.3)')in,il,x
  !   !write(4649,*)in,il,x
  !   if(in /= in2)then
  !      write(12,'(I4,1x,PE10.3)')in2,sum
  !      !write(4649,'(I4,1x,PE10.3)')in2,sum
  !      sum = 0.d0
  !      in2 = in
  !      if(il > 5)then
  !        write(4649,*)"il=",il,nnstep,i
  !        BACKSPACE 11
  !        exit
  !      end if
  !   else if(in == in2)then
  !      sum = sum + x
  !      !write(4649,*)x,sum
  !      !write(4649,*)i,count
  !   end if
  !end do
  !write(12,'(I4,1x,PE10.3)')in,sum
end do
close(11)
close(12)

!open(1,file="fort.3333")
!sum = 0.d0
!do
!  read(1,*,iostat=ierror)in,x
!  if(ierror /= 0)exit 
!    sum = sum + x
!end do

!write(*,*)"sum =",sum


stop
end program




