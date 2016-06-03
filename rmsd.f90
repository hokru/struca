!*********************************
!* A Tool for Structure Analysis *
!*********************************

subroutine run_drmsd(xyz1,diA,xyz2,diB)
use omp_lib
use constant, only: au2ang
use parm, only: npair,nat
implicit none
integer i,j,ntr,k
real(8) xyz1(3,nat),dia(npair)
real(8) xyz2(3,nat),dib(npair)
real(8) e,w1,w0,t0,t1
real(8) mem1,mem2,mem3

call cpu_time(t0)
w0=omp_get_wtime ()


print*, '******************'
print*, '* distance RMSDs *'
print*, '******************'



!$omp parallel
call omp_set_dynamic( .true. )
ntr=omp_get_num_threads()
!$omp end parallel

if(ntr.gt.int(npair)) then
  call omp_set_num_threads(int(npair))
  ntr=omp_get_num_threads()
endif

write(*,*)
write(*,'(2x,a,I12)')    'Threads               : ', ntr
write(*,*)
mem1=2d0*(dble(npair)*8d0)/(1024d0**3) ! diA,diB
mem2=5d0*(dble(nat)*8d0)/(1024d0**3)  ! xyz,iat,ifrez
mem3=2d0*(dble(npair)*8d0)/(1024d0**3) ! ki,kj
write(*,'(2x,a,I12)')    'nat                   : ' ,nat
write(*,'(2x,a,I12)')    'npair                 : ' ,npair
write(*,'(2x,a,F9.4)') 'least Memory (Gb)     : ' ,mem1+mem2+mem3
write(*,*)

!print*,xyz1,xyz2,dia,dib
!xyz1=xyz1*au2ang
!xyz2=xyz2*au2ang

!$omp parallel do private(k)
do k=1,npair
 diA(k)=sqrt(diA(k)*au2ang)
 diB(k)=sqrt(diB(k)*au2ang)
enddo
!$omp end parallel do

e=0
!$OMP parallel do reduction( + : e ) &
!$OMP private(k)
do k=1,npair
  e=e+((diA(k)-diB(k))*(diA(k)-diB(k)))
enddo
!$OMP end parallel do 
e=(e/dble(npair)) 

write(*,'(''dRMSD: [ang] '',F8.4)') e

call cpu_time(t1)
w1=omp_get_wtime()

!call prtim(6,t1-t0,'t','gCP ')
call prtim(6,w1-w0,'w',' ')
!print*,'Speed-up % :',(w1-w0)/(t1-t0)
end





subroutine get_rmsd(nat,xyz1,xyz2)
implicit none
real(8) xyz1(3,nat),xyz2(3,nat)
real(8) rab,e,dx,dy,dz
integer i, nat

e=0
do i=1,nat
e=e+  ( (xyz1(1,i)-xyz2(1,i))**2+(xyz1(2,i)-xyz2(2,i))**2+(xyz1(3,i)-xyz2(3,i))**2 )
enddo

e=sqrt(e/nat)

write(*,'(''RMSD: [ang] '',F8.4)') e

end subroutine
