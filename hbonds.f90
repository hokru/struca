! hydrogen bonding analysis for a molecule

subroutine hbonds(molx)
use parm
use atomdata, only: el
use logic, only:thresh_hbr,thresh_hba
implicit none
real(8) ang,anggrad
integer acc(6),bond(nat,nat),n_hb
type(molecule) molx
real(8) rthr,athr,a1,a2
real(8) rab,r
real(8) dummy(max_dummy)

print*,''
print*,'********************'
print*,'* H-BOND ANALYSIS  *'
print*,'********************'
print*,''
print*,'type: Y-H---X'
print*,'  X/Y=N,O,F,P,S,Cl'
print*,''

! bonding matrix
call bondmatrix(molx%xyz,molx%iat,bond)

! distance and angle criteria (simple!)
rthr=thresh_hbr ! Ang
athr=thresh_hba ! deg D-H--A
a1=180d0-athr
a2=180d0+athr
print*,'thresholds:'
print'(a,F8.1)',     ' max  Hbond length [A]   :',rthr
print'(a,2(F8.1,x))',' Hbond angle range [deg] :', a1,a2
print*,''

n_hb=0
! list of active sites: N,O,F,P,S,Cl
acc(1)=7
acc(2)=8
acc(3)=9
acc(4)=15
acc(5)=16
acc(6)=17
i=0;j=0;k=0
write(*,'(a)'),'         atom      atom    H--X [A]  Y-H--X [deg]'
do i=1,nat ! Y
 if(.not.any(acc==molx%iat(i))) cycle
 do j=1,nat ! H
   if(i==j) cycle 
   if(bond(i,j)/=1) cycle
   if(molx%iat(j)/=1) cycle 
   do k=1,nat 
     if(k==i.or.k==j) cycle
     if(.not.any(acc==molx%iat(k))) cycle
       r=rab(molx%xyz(:,j),molx%xyz(:,k))
       call  angle(molx%xyz,i,j,k,ang,anggrad)
       if(r<=rthr.and.anggrad>=a1.and.anggrad<=a2) then
       n_hb=n_hb+1
         write(*,'(a,x,I5,''['',a2,'']'',I5,''['',a2,'']'',x,F8.4,x,F8.1)'),'H-bond',j,el(molx%iat(i)),k,el(molx%iat(j)),r,anggrad
!         print*,i,j,k,r,anggrad
       endif
   enddo
 enddo
enddo

print*,''
print'(a,I4)',' # H-bonds: ',n_hb
print*,''

end subroutine


! H-bond difference

subroutine delta_hbonds(mol1,mol2)
use parm
use atomdata, only: el
use logic, only:thresh_hbr,thresh_hba
implicit none
real(8) ang1,anggrad1
real(8) ang2,anggrad2
integer acc(6),bond(nat,nat),n_hb
type(molecule) mol1,mol2
real(8) rthr,athr,a1,a2
real(8) rab,r1,r2,da,dr
real(8) dummy(max_dummy)
real(8), allocatable :: tvec(:)


print*,''
print*,'********************'
print*,'* H-BOND ANALYSIS  *'
print*,'********************'
print*,''
print*,'type: Y-H---X'
print*,'  X/Y=N,O,F,P,S,Cl'
print*,''

! bonding matrix
call bondmatrix(mol1%xyz,mol1%iat,bond)

! distance and angle criteria (simple!)
rthr=3.0 ! Ang
athr=120 ! deg D-H--A
rthr=thresh_hbr
athr=thresh_hba
a1=180d0-athr
a2=180d0+athr
print*,'thresholds:'
print'(a,F8.1)',     ' max  Hbond length [A]   :',rthr
print'(a,2(F8.1,x))',' Hbond angle range [deg] :', a1,a2
print*,''

n_hb=0
! list of active sites: N,O,F,P,S,Cl
acc(1)=7
acc(2)=8
acc(3)=9
acc(4)=15
acc(5)=16
acc(6)=17
i=0;j=0;k=0
write(*,'(a)'),'                atom      atom    d(H--X) [A]  d(Y-H--X) [deg]'
do i=1,nat ! Y
 if(.not.any(acc==mol1%iat(i))) cycle
 do j=1,nat ! H
   if(i==j) cycle
   if(bond(i,j)/=1) cycle
   if(mol1%iat(j)/=1) cycle
   do k=1,nat
     if(k==i.or.k==j) cycle
     if(.not.any(acc==mol1%iat(k))) cycle
       r1=rab(mol1%xyz(:,j),mol1%xyz(:,k))
       r2=rab(mol2%xyz(:,j),mol2%xyz(:,k))
       call  angle(mol1%xyz,i,j,k,ang1,anggrad1)
       call  angle(mol2%xyz,i,j,k,ang2,anggrad2)
       da=anggrad1-anggrad2
       dr=r1-r2
       if(r1<=rthr.and.anggrad1>=a1.and.anggrad1<=a2) then
       n_hb=n_hb+1
       dummy(n_hb)=dr
         write(*,'(a,x,I5,''['',a2,'']'',I5,''['',a2,'']'',x,F8.4,x,F8.1)'),'delta_H-bond',j,el(mol1%iat(i)),k,el(mol1%iat(j)),dr,da
!         print*,i,j,k,r,anggrad
       endif
   enddo
 enddo
enddo

print*,''
print'(a,I4)',' # H-bonds: ',n_hb
print*,''
if(n_hb>0) then
 allocate(tvec(n_hb))
 tvec(1:n_hb)=dummy(1:n_hb)
 print'(a,F8.1)',' mean abs Hbond deviation [A] : ', sum(abs(tvec))/dble(n_hb)
 print'(a,F8.1)',' mean Hbond deviation [A]     : ', sum(tvec)/dble(n_hb)
 print'(a,F8.1)',' max (+) deviation [A]        : ', maxval(tvec)
 print'(a,F8.1)',' max (-) deviation [A]        : ', minval(tvec)
endif
end subroutine

