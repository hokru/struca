subroutine composition(nat,iat)
implicit none
integer nat, iat(nat)
integer nel,comp(107),c,k,i,ic(107)
character(2) esym,to_upper
character(5) num

comp=0
ic=0

do k=1,107
c=0
 do i=1,nat
  if(iat(i)==k) then
    c=c+1
    comp(k)=c
    ic(k)=i
  endif
 enddo
enddo

write(*,'(a,$)') ' element composition : '
do k=1,107
if(ic(k)/=0) then
write(num,'(I5)') comp(k)
write(*,'(a,$)') trim(adjustl(to_upper(esym(k))))//trim(adjustl(num))//' '
endif
enddo
write(*,'(a)')''

end subroutine

subroutine getCOM(com,nat,xyz,iat,orient)
use atomdata, only: ams
use constant, only: au2ang
implicit none
integer nat, iat(nat)
real(8) xyz(3,nat)
real(8) com(3), mmass
integer i,j
logical orient

xyz=xyz/au2ang
mmass=0d0
com=0d0
do i=1,nat
mmass=mmass+ams(iat(i))
 do j=1,3
  com(j)=com(j)+xyz(j,i)*ams(iat(i))
 enddo
enddo

com=com/mmass


!write(*,'(3x,a,3F12.5)') ' molecular mass      : ',mmass
write(*,'(3x,a,3F12.5)') ' center of mass (A)  : ',com(1:3)*au2ang
!write(*,'(3x,a,3F12.5)') 'center of mass (au)',com(1:3)

if(orient) then
! move molecule to COM
do i=1,nat
 do j=1,3
    xyz(j,i)=xyz(j,i)-com(j) 
 enddo 
enddo
endif
!call wrxyz('cema_xopt.xyz')
!endif


end


! z**2+y**2    -xy        -xz     
! -xy         z**2+x**2   -yz
! -xz           -yz        x*2+y**2
! calculate moment of inertia and principle axis
subroutine getIntertia(nat,iat,xyz,orient,rot)
use atomdata, only: ams
use constant, only: Planck, pc_c, pi, bohr2m, amu2kg,au2ang
implicit none
integer i,j,k,nat,iat(nat)
real(8) xyz(3,nat)
real(8) mom(3,3),paxis(3),e(3),rot(3)
real(8) x,y,z,m
real(8) conv,s,ddot
real(8) coord(3,nat)
logical orient

! Conversion factor from moments to rotational constants.
  conv= Planck / (8.0d0 * pi *pi * pc_c)
! Add factor to put moments into SI units - give result in wavenumbers.
  conv=conv / (bohr2m * bohr2m * amu2kg * 100.d0)
!  conv=conv * au2ang*au2ang




mom=0d0
!xyz=xyz/au2ang
do i=1,nat
x=xyz(1,i)
y=xyz(2,i)
z=xyz(3,i)
m=ams(iat(i))
mom(1,1)=mom(1,1)+(z**2+y**2)*m
mom(1,2)=mom(1,2)-x*y*m
mom(1,3)=mom(1,3)-x*z*m
mom(2,2)=mom(2,2)+(z**2+x**2)*m
mom(2,3)=mom(2,3)-y*z*m
mom(3,3)=mom(3,3)+(x**2+y**2)*m
enddo

!print*,xyz(1:3,1)
mom(2,1)=mom(1,2)
mom(3,1)=mom(1,3)
mom(3,2)=mom(2,3)

!do i=1,3
!write(*,'(3F14.6)') mom(i,1:3)
!enddo

! diag, mom contains now eigenvectors
call DiagSM(3,mom,e)
!print*, mom(1:3,1)
!print*, mom(1:3,2)
!print*, mom(1:3,3)

!mom(1:3,2)=-mom(1:3,2)
!mom(1:3,3)=-mom(1:3,3)

do i=1,3
if(e(i)<1e-5) then
  rot=0d0
else
  rot(i)=conv/e(i)
endif
enddo

if(e(1)==0) then
 write(*,'(3x,a)') '  linear molecule!'
 write(*,'(3x,a,2(F12.5,a))') ' Rotational constants: A= ***  B= ',rot(2) ,' C= ',rot(3), ' [cm^-1]'
 rot=rot*pc_c/10000.d0
 write(*,'(3x,a,2(F12.5,a))') ' Rotational constants: A= ***  B= ',rot(2) ,' C= ',rot(3), ' [Mhz]'
else
 write(*,'(3x,a,3(F12.5,a))') ' Rotational constants:  A= ', rot(1),' B= ',rot(2) ,' C= ',rot(3), ' [cm^-1]'
 rot=rot*pc_c/10000.d0
 write(*,'(3x,a,3(F12.5,a))') ' Rotational constants:  A= ', rot(1),' B= ',rot(2) ,' C= ',rot(3), ' [Mhz]'
endif





if(orient) then
print*, '   ** rotating molecule to principle axis frame  **'
! handedness
s=mom(1,1)*(mom(2,2)*mom(3,3)-mom(3,2)*mom(2,3)) +  &               
  mom(1,2)*(mom(2,3)*mom(3,1)-mom(2,1)*mom(3,3)) +  &               
  mom(1,3)*(mom(2,1)*mom(3,2)-mom(2,2)*mom(3,1))  
!print*,'s',s                 
! invert if left-handed
if(s<0) then
  do i=1,3
  mom(i,1)=-mom(i,1)
  enddo
endif



coord=xyz
! rotate to principle axis frame
do i=1,nat
 do j=1,3
 xyz(j,i)=ddot(3,coord(1,i),1,mom(1,j),1)
 enddo
enddo

!xyz=xyz/au2ang
!print*, '   ** writing cema_coord.xyz  **'
!call wrxyz(iat,nat,xyz,'cema_coord.xyz')

endif


end subroutine
