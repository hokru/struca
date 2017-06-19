!*********************************
!* A Tool for Structure Analysis *
!*********************************

subroutine run_drmsd(xyz1,diA,xyz2,diB)
use constant, only: au2ang
use parm, only: npair,nat
implicit none
integer i,j,ntr,k
real(8) xyz1(3,nat),dia(npair)
real(8) xyz2(3,nat),dib(npair)
real(8) e,w1,w0,t0,t1
real(8) mem1,mem2,mem3

do k=1,npair
 diA(k)=sqrt(diA(k)*au2ang)
 diB(k)=sqrt(diB(k)*au2ang)
enddo

e=0
do k=1,npair
  e=e+((diA(k)-diB(k))*(diA(k)-diB(k)))
enddo
e=(e/dble(npair)) 

write(*,'(''distance RMSD: [ang] '',F8.4)') e
end




! standard RMSD
! if this does not fit to the qfit RMSD, it means the alignment is not working correctly!
subroutine get_rmsd(nat,xyz1,xyz2,iat)
implicit none
real(8) xyz1(3,nat),xyz2(3,nat)
real(8) rab,e,dx,dy,dz,eheavy
integer i, nat,iat(nat)

print*,' conventional RMSD between 2 sets of coordinates'

e=0
eheavy=0
do i=1,nat
e=e+  ( (xyz1(1,i)-xyz2(1,i))**2+(xyz1(2,i)-xyz2(2,i))**2+(xyz1(3,i)-xyz2(3,i))**2 )
if (iat(i)==1) cycle
eheavy=eheavy+  ( (xyz1(1,i)-xyz2(1,i))**2+(xyz1(2,i)-xyz2(2,i))**2+(xyz1(3,i)-xyz2(3,i))**2 )
enddo

e=sqrt(e/nat)
eheavy=sqrt(eheavy/nat)

write(*,'(2x,''RMSD: [ang] '',F8.4)') e
write(*,'(2x,''RMSD (heavy): [ang] '',F8.4)') eheavy

end subroutine
