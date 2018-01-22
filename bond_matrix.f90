! ASSUMES ANGSTROM!
subroutine bondmatrix(xyz,iat,bond)
use parm, only: nat,i,j
use logic, only: bondf
use atomdata, only: rcov,rcov3
use constant, only: au2ang
implicit none
integer bond(nat,nat)
real(8) r1,r2,rab,fac_bond
real(8), intent(in) :: xyz(3,nat)
integer, intent(in) :: iat(nat)


! full symmetric bond matrix
!bond=0
!do i=1,nat
! do j=1,nat
!   if(i==j) cycle
!   r1=rcov(iat(i))+rcov(iat(j))*0.5
!   r2=rab(xyz(1,i),xyz(1,j))
!   if(abs(r1-r2) < 0.5) bond(i,j)=1
! enddo
!enddo


fac_bond=bondf
! full symmetric bond matrix
bond=0
do i=1,nat
 do j=1,nat
   if(i==j) cycle
!   r1=(rcov(iat(i))+rcov(iat(j))) ! in A
   r1=(rcov3(iat(i))+rcov3(iat(j))) ! in A
   r2=rab(xyz(1,i),xyz(1,j))
   if(r2 <= r1*fac_bond) bond(i,j)=1
 enddo
enddo


! call prmati(6,bond,nat,nat,'bond matrix')
end




subroutine count_fragments(xyz,nat,bond)
implicit none
real(8) xyz(3,nat)
integer bond(nat,nat)
integer i,j,nfrag,nat

nfrag=1
do i=1,nat
 do j=1,nat
  if(bond(i,j)==1) then
  endif
 enddo
enddo


end subroutine
