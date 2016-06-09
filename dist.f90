subroutine get_dist(molx)
!subroutine get_dist(nat,npair,echo,xyz,dist)
use parm
implicit none
real(8) dx,dy,dz
logical echo
integer istat

type(molecule) molx

!if(npair.gt.thrP) then
! simple version
! loop over all unique pairs
  k=1
  do i=1,nat-1
   do j=i+1,nat
   dx=molx%xyz(1,i)-molx%xyz(1,j)
   dy=molx%xyz(2,i)-molx%xyz(2,j)
   dz=molx%xyz(3,i)-molx%xyz(3,j)
   molx%dist(k)=sqrt(dx*dx+dy*dy+dz*dz)
!    print*, molx%dist(k)
   k=k+1
   enddo
  enddo
  k=k-1

!else

!k=0
!do i=1,nat
! do j=1,nat
! if(i.gt.j.or.i.eq.j) cycle
! k=k+1
! dx=xyz(1,i)-xyz(1,j)
! dy=xyz(2,i)-xyz(2,j)
! dz=xyz(3,i)-xyz(3,j)
! dist(k)=dx*dx+dy*dy+dz*dz
! enddo
!enddo

! slower and more memory ... 
!if(istat.ne.0) stop 'allocation error @ index vectors'

! make index vectors
!k=0
!ki=0
!kj=0
!do i=1,nat-1
! do j=i+1,nat
! k=k+1
! kj(k)=j
! ki(k)=i
! enddo
!enddo

!omp parallel do private(i,j,k,dx,dy,dz)
!do k=1,npair
!i=ki(k)
!j=kj(k)
! dx=molx%xyz(1,i)-molx%xyz(1,j)
! dy=molx%xyz(2,i)-molx%xyz(2,j)
! dz=molx%xyz(3,i)-molx%xyz(3,j)
! molx%dist(k)=dx*dx+dy*dy+dz*dz
!enddo
!omp end parallel do

!endif


!if(echo) then
!write(*,'(2x,a,F6.2,a)') 'Keeping: ',(dble(k)/dble(npair)*100d0),' % of atom pairs after distance screening'
!write(*,'(2x,a,I)') 'npair before screening : ' ,npair
!npair=k
!write(*,'(2x,a,I)') 'npair after  screening : ' ,npair
!endif
end subroutine
