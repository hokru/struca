! ASSUMES ANGSTROM!
subroutine bondmatrix(xyz,iat,bond)
use parm, only: nat,i,j
use atomdata, only: rcov
use constant, only: au2ang
implicit none
integer bond(nat,nat)
real(8) r1,r2,rab
real(8), intent(in) :: xyz(3,nat)
integer, intent(in) :: iat(nat)


! full symmetric bond matrix
bond=0
do i=1,nat
 do j=1,nat
   if(i==j) cycle
   r1=rcov(iat(i))+rcov(iat(j))*0.5
   r2=rab(xyz(1,i),xyz(1,j))
   if(abs(r1-r2) < 0.5) bond(i,j)=1
 enddo
enddo

! call prmati(6,bond,nat,nat,'bond matrix')
end